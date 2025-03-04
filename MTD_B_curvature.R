
# Install the necessary packages
install.packages("dplyr")
install.packages("ggplot2")
install.packages("ggforce")

#Load the necessary libraries
library(dplyr)
library(ggplot2)
library(ggforce)

# The WT, FAP20_PACRG and invivo data in separated folders.
# Adapt the directions to the folder containing the data.
dir_wt = # your folder directory (eg. "N:/Commun/Umut/Proteins/")
dir_prot = # your folder directory (eg. "N:/Commun/Umut/Proteins/")
dir_invivo = # your folder directory (eg. "N:/Commun/Umut/Proteins/")


list_dir = c(dir_wt, dir_prot, dir_invivo )
condition_umut = c("WT", "FAP20_PACRG", "dir_invivo" )


# All files---------------------------
first = TRUE

for (i in seq_along(list_dir)){
  
  dir = list_dir[i]
  print(dir)
  condition_seq = condition_umut[i]
  print(condition_seq)

  setwd(dir)
  files_umut = list.files(pattern = 'Position', recursive = TRUE)
  print(files_umut)

  
  # Loop over the files in the directory for intensity ----------------------
  for (f in files_umut) {
    
    # Read the file using the helper function
    dat = read.table(f, header = TRUE, sep = "", dec = ".")
    
    # Extract cell type from filename
    t_number = sub(".* of ", "", f)
    experiment_umut = sub("/.*$", "",t_number)
    Position_beam = sub(".*_(Position[0-9]+)_.*", "\\1",f)
    micro_doublet = sub(".*_(MTD[0-9]+)_.*", "\\1",f)
    section = sub(".*_(sec[0-9]+)\\.txt$", "\\1",f)
    Doublet_pos = sub(".*_(right|left)_.*", "\\1",f)
    
    # Add columns for 'Cell_Type' and 'experiment'
    dat = dat %>%
      mutate(tomo_number = t_number, 
             experiment = experiment_umut,
             position = Position_beam,
             MTD = micro_doublet,
             reslice = section,
             doublet_position = Doublet_pos,
             condition = condition_seq
             )
    
    dat$tomo_number = gsub(".txt", "", dat$tomo_number)
   
    
    print(f)
    print(dat$condition)
    
    
    # If it's the first file, initialize the final dataframe
    if (first) {
      dat_coordinates = dat
      first = FALSE
    } else {
      # Ensure columns are the same
      common_columns = intersect(names(dat_coordinates), names(dat))
      dat = dat[, common_columns, drop = FALSE]
      dat_coordinates = dat_coordinates[, common_columns, drop = FALSE]
      
      # Combine data
      dat_coordinates = rbind(dat_coordinates, dat)
      
    }
  }
}
dat_coordinates$doublet_position = ifelse(dat_coordinates$doublet_position %in% c("left", "right"), dat_coordinates$doublet_position, NA)
dat_coordinates$condition =as.factor(dat_coordinates$condition)
levels(dat_coordinates$condition)

dim(dat_coordinates)
dat_coord2 = dat_coordinates # This is the full raw database
dim(dat_coord2)

dat_coord2$key = paste(dat_coord2$position, dat_coord2$MTD, dat_coord2$doublet_position, sep = "_")
dat_coord2$key = as.factor(dat_coord2$key)







# Normalize the full database------------------------------------------------------------

# Colors used for ploting
color_circle = c("#009966")
inner_color = c("#FFFFFF")
color_lines = c("#3399CC")
wt_protein_color = c( "black", "tomato1",  "slategray2" )



# Function to normalize a single dataset
normalize_dataset <- function(data, target_angle, target_distance) {
  # Translate: Set Point 0 at (0, 0)
  x0 <- data$x[1]
  y0 <- data$y[1]
  data$x <- data$x - x0
  data$y <- data$y - y0
  # Calculate rotation angle for Point 1
  x1 <- data$x[2]
  y1 <- data$y[2]
  rotation_angle <- atan2(y1, x1) - target_angle
  # Rotate all points
  rotated_coords <- rotate_points(data$x, data$y, -rotation_angle)
  data$x <- rotated_coords$x
  data$y <- rotated_coords$y
  # Scale all points so Point 1 is at the target distance
  current_distance <- sqrt(data$x[2]^2 + data$y[2]^2)
  scale_factor <- target_distance / current_distance
  data$x <- data$x * scale_factor
  data$y <- data$y * scale_factor
  return(data)
}

# Function to rotate points
rotate_points <- function(x, y, rotation_angle) {
  x_rot <- x * cos(rotation_angle) - y * sin(rotation_angle)
  y_rot <- x * sin(rotation_angle) + y * cos(rotation_angle)
  data.frame(x = x_rot, y = y_rot)
}

# Parameters
target_angle <- pi / 4  # 45° in radians
target_distance <- 1    # Fixed distance for Point 1 from Point 0



# Process all dataset by grouping and applying normalization
aligned_data <- dat_coord2 %>%
  group_by(condition, tomo_number, ) %>%
  group_modify(~ normalize_dataset(.x, target_angle, target_distance)) %>%
  ungroup()

aligned_data = as.data.frame(aligned_data)
aligned_data$key = paste(aligned_data$position, aligned_data$MTD, aligned_data$doublet_position, sep = "_")
aligned_data$key = as.factor(aligned_data$key)


# When the point 7 of one group is x<0 and y>0, write negative in a new column
aligned_data <- aligned_data %>%
  group_by(tomo_number)%>%
  rowwise()%>%
  mutate(
    # Identify counter-clockwise groups
    oriantation = ifelse(point == 7 & x < 0 & y > 0, "negative", "positive" ))
  

# When "negative" is found within a group, write TRUE in a new column 
aligned_data <- aligned_data %>%
  group_by(tomo_number)%>%
  mutate(
    # Identify counter-clockwise groups
    is_counter_clockwise = any(oriantation == "negative"))%>%
  rowwise()


# Reflect all points in the aligned_data dataset
df_reflected <- aligned_data %>%
  group_by(tomo_number)%>%
  mutate(
    # Identify counter-clockwise groups
    is_counter_clockwise = any(oriantation == "negative"))%>%
  rowwise() %>%
  mutate(temp_x = ifelse(is_counter_clockwise == TRUE, x, y),# Temporarily store x
         temp_y = ifelse(is_counter_clockwise == TRUE, y, x),# Temporarily store y
         x = temp_y,       # Assign y to x
         y = temp_x) %>%  # Assign stored x to y
  select(-temp_x,-is_counter_clockwise, -temp_y)  # Drop the temporary column

df_reflected = as.data.frame(df_reflected) # This is the data after geometric reflection of the counter clockwise microtubules



# Plot all the microtubules from WT ant Protein conditions (Not currated data)
ggplot(df_reflected, aes(x = x, y = y)) +
  geom_path(aes(group = tomo_number), linewidth = 0.5, color = color_lines, alpha = 0.4) +
  geom_circle(aes(x0=0, y0=0, r=1 ), fill = inner_color, color = color_circle, linewidth = 1, inherit.aes=FALSE) +
  theme_bw()+
  theme(panel.grid.major =element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        strip.text.x = element_text(size = 4,face="bold"),
        legend.position = "none",
  )+
  facet_wrap(~condition)




# Data averaging----------------------------------------------------------------

# Compute the average trajectory by averaging across all datasets at each index for non currated data
average_trajectory <- df_reflected %>%
  group_by(condition, tomo_number ) %>%
  mutate(index = row_number(),
         condition = condition) %>%  # Create an index to match the data points in each group
  ungroup() %>%
  group_by(condition, index ) %>%
  summarize(mean_x = mean(x), mean_y = mean(y)) %>%
  ungroup()
average_trajectory = average_trajectory[1:24,]



#Plot the WT, Protein, and invivo condition with the line average together

ggplot(df_reflected, aes(x = x, y = y)) +
  geom_path(aes(group = tomo_number), linewidth = 0.8, color = color_lines, alpha = 0.4) +
  geom_path(data = average_trajectory, aes(x = mean_x, y = mean_y), 
            color = "black", linetype = "solid", linewidth = 4)+ 
  geom_circle(aes(x0=0, y0=0, r=1 ), fill = inner_color, color = color_circle, linewidth = 6, inherit.aes=FALSE) +
  theme_bw()+
  theme(panel.grid.major =element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        strip.text.x = element_text(size = 12,face="bold"),
        legend.position = "none",
  )+
  facet_wrap(~ condition)





# Data curration------------------------------------------

# Make two keys to isolate the microtubules and the individual sections per microtubules 
df_reflected$reslice = as.factor(df_reflected$reslice)
df_reflected$key = paste(df_reflected$position, df_reflected$MTD, df_reflected$doublet_position, sep = "_")
df_reflected$key_reslice = paste(df_reflected$position, df_reflected$MTD, df_reflected$doublet_position, df_reflected$reslice , sep = "_")

df_reflected$key = as.factor(df_reflected$key)
df_reflected$key_reslice = as.factor(df_reflected$key_reslice)


# Plot the condition to observe each microtubule individually and found the once with weird sections
# Adapt the condition for "FAP20_PACRG" or "dir_invivo"
ggplot(subset(df_reflected, condition == "WT" ), aes(x = x, y = y)) +
  geom_path(aes(group = tomo_number), linewidth = 0.5, color = color_lines, alpha = 0.4) +
  geom_circle(aes(x0=0, y0=0, r=1 ), fill = inner_color, color = color_circle, linewidth = 1, inherit.aes=FALSE) +
  theme_bw()+
  theme(panel.grid.major =element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        strip.text.x = element_text(size = 4,face="bold"),
        legend.position = "none",
  )+
  facet_wrap(~key)



# plot the separated sections for the selected microtubule to identify the weird section using the key_reslice
# Adapt the condition for "FAP20_PACRG" or "dir_invivo" and the key for the desired position
ggplot(subset(df_reflected, condition == "WT" & key == "Position11_MTD1_right"), aes(x = x, y = y)) +
  geom_path(aes(group = key_reslice), linewidth = 0.5, color = color_lines, alpha = 0.4) +
  geom_circle(aes(x0=0, y0=0, r=1 ), fill = inner_color, color = color_circle, linewidth = 1, inherit.aes=FALSE) +
  theme_bw()+
  theme(panel.grid.major =element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        strip.text.x = element_text(size = 4,face="bold"),
        legend.position = "none",
  )+
  facet_wrap(~key_reslice)


# Make a vector with the sections to remove from the analysis
to_remove = c("Position20_MTD1_NA_sec2", "Position29_MTD1_NA_sec5", "Position29_MTD1_NA_sec7", "Position30_MTD1_right_sec26", "Position30_MTD1_right_sec28", 
              "Position11_MTD1_right_sec3")





# Make CURRATED database from the reflected data -------------------------------
# This dataset will be cleaned-up for unwanted section identified above

# Make a new dataframe where the weird data were removed
df_reflected_currated = subset(df_reflected, !key_reslice %in% to_remove)



# Summarise the data by counting the number of microtubule sections per condition
summary_table <- df_reflected_currated %>%
  group_by(condition) %>%
  summarise(
    Number_of_microtubules = n_distinct(key), # Count unique microtubules
    unique_sections = n_distinct(key_reslice), # Count unique tomo_number
    total_coordinates_B_microtubule = n()- unique_sections# Count total number of points measured on the B microtubule
  )

print(summary_table)




# Compute the average trajectory by averaging across all datasets at each index using the currated data
average_trajectory_currated <- df_reflected_currated %>%
  group_by(condition, tomo_number ) %>%
  mutate(index = row_number(),
         condition = condition) %>%  # Create an index to match the data points in each group
  ungroup() %>%
  group_by(condition, index ) %>%
  summarize(mean_x = mean(x), mean_y = mean(y)) %>%
  ungroup()

average_trajectory_currated = average_trajectory_currated[1:24,] # This removes the only point 9 that was measure
average_trajectory_currated = as.data.frame(average_trajectory_currated)



# Compute the average trajectory by averaging across one chosen WT microtubule 
subset_MTD_WT = subset(df_reflected_currated, condition == "WT" & key == "Position27_MTD1_NA")

average_WT <- subset_MTD_WT %>%
  group_by(condition, tomo_number ) %>%
  mutate(index = row_number(),
         condition = condition) %>%  # Create an index to match the data points in each group
  ungroup() %>%
  group_by(condition, index ) %>%
  summarize(mean_x = mean(x), mean_y = mean(y)) %>%
  ungroup()
average_WT = as.data.frame(average_WT)

# Plot the choosen WT microtubule doublet with the avaerage of all its sections
ggplot(subset(df_reflected_currated, condition == "WT" & key == "Position27_MTD1_NA"), aes(x = x, y = y)) +
  geom_path(aes(group = key_reslice), linewidth = 1, color = color_lines, alpha = 0.4) +
  geom_path(data = average_WT, aes(x = mean_x, y = mean_y), 
            color = "black", linetype = "longdash", linewidth = 2)+ 
  geom_circle(aes(x0=0, y0=0, r=1 ), fill = inner_color, color = color_circle, linewidth = 1, inherit.aes=FALSE) +
  theme_bw()+
  theme(panel.grid.major =element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        strip.text.x = element_text(size = 4,face="bold"),
        legend.position = "none",
  )





# Compute the average trajectory by averaging across one choosen FAP20_PARG microtubule 
subset_MTD_FAP20_PACRG = subset(df_reflected_currated, condition == "FAP20_PACRG" & key == "Position23_MTD1_NA")

average_FAP20_PACRG <- subset_MTD_FAP20_PACRG %>%
  group_by(condition, tomo_number ) %>%
  mutate(index = row_number(),
         condition = condition) %>%  # Create an index to match the data points in each group
  ungroup() %>%
  group_by(condition, index ) %>%
  summarize(mean_x = mean(x), mean_y = mean(y)) %>%
  ungroup()
average_FAP20_PACRG = as.data.frame(average_FAP20_PACRG)

# Plot the choosen WT microtubule doublet with the avaerage of all its sections
# This plot was used for publication
ggplot(subset_MTD_FAP20_PACRG, aes(x = x, y = y)) +
  geom_path(aes(group = key_reslice), linewidth = 1, color = "tomato1", alpha = 0.4) +
  geom_path(data = average_FAP20_PACRG, aes(x = mean_x, y = mean_y), 
            color = "black", linetype = "longdash", linewidth = 2)+ 
  geom_circle(aes(x0=0, y0=0, r=1 ), fill = inner_color, color = color_circle, linewidth = 1, inherit.aes=FALSE) +
  theme_bw()+
  theme(panel.grid.major =element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        strip.text.x = element_text(size = 4,face="bold"),
        legend.position = "none",
  )




#Plot the WT and FAP20_PACRG condition with the line average together
# This plot was used for publication
ggplot(subset(df_reflected_currated, condition != "dir_invivo") , aes(x = x, y = y, color = condition)) +
  geom_path(aes(group = tomo_number), linewidth = 0.8, alpha = 0.02) +
  geom_path(data = subset(average_trajectory_currated, condition != "dir_invivo"), aes(x = mean_x, y = mean_y), 
            color = "black", linetype = "longdash", linewidth = 2)+ 
  geom_circle(aes(x0=0, y0=0, r=1 ), fill = inner_color, color = color_circle, linewidth = 2, inherit.aes=FALSE) +
  scale_color_manual(values = c("FAP20_PACRG" = "tomato1", "WT" = "#3399CC" ))+
  theme_bw()+
  theme(panel.grid.major =element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        strip.text.x = element_text(size = 12,face="bold"),
        legend.position = "none",
  )+
  facet_wrap(~ condition)




# Plot the averaged conditions (WT, invivo, FAP20_PACRG) together without all the measurments
# This plot was used for publication
ggplot(average_trajectory_currated, aes(x = mean_x, y = mean_y)) +
  geom_path(aes(group = condition, color = condition), linetype = "solid", linewidth = 4, alpha = 1) +
  geom_circle(aes(x0=0, y0=0, r=1 ), fill = inner_color, color = color_circle, linewidth = 4, inherit.aes=FALSE) +
  scale_color_manual(values = wt_protein_color) +
  theme_bw()+
  theme(panel.grid.major =element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        strip.text.x = element_text(size = 4,face="bold"),
        legend.position = "none",
  )







