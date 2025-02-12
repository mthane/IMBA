##### binning data #####
# the binning data functions allow binning data among a certain variable

#' Title
#'
#' @param df 
#' @param variable 
#' @param width 
#' @param frame_interval 
#' @param distance_to_odor_interval 
#' @param Abs_HC_Angle_interval 
#' @param Abs_bearing_angle 
#' @param radius 
#' @param frame_rate 
#' @param direction 
#'
#' @return
#' 
#' @import dplyr
#' @import data.table
#' @export
#'c
 
binning_data <- function(df,
                         variable,
                         width = 10,
                         frame_interval=c(0,2880),
                         distance_to_odor_interval=c(0,Inf),
                         Abs_HC_Angle_interval=c(0,360),
                         Abs_bearing_angle=c(0,360),
                         spinepoint_y_6_interval =c(-50,50),
                         spinepoint_x_6_interval =c(-50,50),
                         abs_y_angle_interval =c(0,360),
                         radius=42.5,
                         frame_rate= 16,
                         direction = NULL,
                         group_condition="group_condition"
) {
  tic()
  message(paste("binning data ... ", variable , sep = " "))
  
  if (variable == "frame") {
    width = width*frame_rate
  }
  
  if (variable == "time_to_trans_to") {#new after publication
    width = width*frame_rate
  }
  
  if (variable == "time_to_trans_from") {#new after publication
    width = width*frame_rate
  }
  

  radius = as.numeric(radius)
  df <- df %>% filter(
    frame %in% seq(from = frame_interval[1] * frame_rate,
                   to = frame_interval[2] * frame_rate) &
      distance_to_odor >= distance_to_odor_interval[1] &
      distance_to_odor <= distance_to_odor_interval[2] &
      abs(bearing_angle) >= Abs_bearing_angle[1] &
      abs(bearing_angle) <= Abs_bearing_angle[2] &
      spinepoint_y_6_conv >= spinepoint_y_6_interval[1] &
      spinepoint_y_6_conv <= spinepoint_y_6_interval[2] &
      spinepoint_x_6_conv >= spinepoint_x_6_interval[1] &
      spinepoint_x_6_conv <= spinepoint_x_6_interval[2] &
      abs_y_angle >= abs_y_angle_interval[1] &
      abs_y_angle <= abs_y_angle_interval[2] 
  )
  
  #direction filter
  if(direction=="forwards"){
    df <- as.data.frame(df)%>%
      filter(run_direction==1)
  }
  if(direction=="backwards"){
    df <- as.data.frame(df)%>%
      filter(run_direction==-1)
  }
  if(direction=="both"){
    df <- as.data.frame(df)%>%
      filter(run_direction==1|run_direction==-1)
  }
  
  
  data <-binning_data_ts(df,variable,width)%>%
    left_join(binning_data_hc(df, variable, width, Abs_HC_Angle_interval,frame_rate=frame_rate,group_condition)) %>%
    left_join(binning_data_is(df, variable,  width,group_condition)) %>%
    left_join(binning_data_run_speed(df, variable, width,group_condition)) %>%
    left_join(binning_data_preference(df, variable, radius, width,group_condition))
  
  
  if (variable == "frame") {
    data$midpoint = data$midpoint / frame_rate
  }
  
  if (variable == "time_to_trans_to") {#new after publication
    data$midpoint = data$midpoint / frame_rate
  }
  
  if (variable == "time_to_trans_from") {#new after publication
    data$midpoint = data$midpoint / frame_rate
  }
  

  it = toc()
  exectime <- it$toc - it$tic
  message(paste("data binning... ", exectime , " seconds elapsed", sep =
                  " "))
  data
  
}



#' Title
#'
#' @param df 
#' @param variable 
#' @param width 
#' @param frame_interval 
#' @param distance_to_odor_interval 
#' @param Abs_HC_Angle_interval 
#' @param Abs_bearing_angle 
#' @param radius 
#' @param frame_rate 
#' @param grouping 
#' @param direction 
#'
#' @return
#' @export
#' @import dplyr
#' @import data.table
#' 
 
binning_data_grouped <- function(df,
                                 variable,
                                 width = 10,
                                 frame_interval=c(0,2880),
                                 distance_to_odor_interval=c(0,Inf),
                                 Abs_HC_Angle_interval=c(0,360),
                                 Abs_bearing_angle=c(0,360),
                                 spinepoint_y_6_interval =c(-50,50),
                                 spinepoint_x_6_interval =c(-50,50),
                                 abs_y_angle_interval =c(0,360),
                                 radius=42.5,
                                 frame_rate= 16,
                                 grouping="id",
                                 direction = NULL,group_condition="group_condition"
) {
  tic()
  message(paste("binning data ... ", variable , sep = " "))
  radius = as.numeric(radius)
  
  df <- df %>% filter(
    frame %in% seq(from = frame_interval[1] * frame_rate,
                   to = frame_interval[2] * frame_rate) &
      distance_to_odor >= distance_to_odor_interval[1] &
      distance_to_odor <= distance_to_odor_interval[2] &
      abs(bearing_angle) >= Abs_bearing_angle[1]&
      abs(bearing_angle) <= Abs_bearing_angle[2]&
      spinepoint_y_6_conv >= spinepoint_y_6_interval[1] &
      spinepoint_y_6_conv <= spinepoint_y_6_interval[2] &
      spinepoint_x_6_conv >= spinepoint_x_6_interval[1] &
      spinepoint_x_6_conv <= spinepoint_x_6_interval[2] &
      abs_y_angle >= abs_y_angle_interval[1] &
      abs_y_angle <= abs_y_angle_interval[2] 
  )
  
  #direction filter
  if(direction=="forwards"){
    df <- as.data.frame(df)%>%
      filter(run_direction==1)
  }
  if(direction=="backwards"){
    df <- as.data.frame(df)%>%
      filter(run_direction==-1)
  }
  if(direction=="both"){
    df <- as.data.frame(df)%>%
      filter(run_direction==1|run_direction==-1)
  }
  
  
  data <-binning_data_ts_grouped(df,variable,width,grouping=grouping,group_condition)%>%
    left_join(binning_data_hc_grouped(df, variable, width, grouping=grouping,Abs_HC_Angle_interval,frame_rate=frame_rate,group_condition)) %>%
    left_join(binning_data_is_grouped(df, variable, width,grouping=grouping,group_condition)) %>%
    left_join(binning_data_run_speed_grouped(df, variable, width,grouping=grouping,group_condition)) %>%
    left_join(binning_data_preference_grouped(df, variable, radius, width,grouping=grouping,group_condition))
  
  if (variable == "frame") {
    data$midpoint = data$midpoint / frame_rate
  }
  
  if (variable == "time_to_trans_to") {#new after publication
    data$midpoint = data$midpoint / frame_rate
  }
  
  if (variable == "time_to_trans_from") {#new after publication
    data$midpoint = data$midpoint / frame_rate
  }
  

  it = toc()
  exectime <- it$toc - it$tic
  message(paste("data binning... ", exectime , " seconds elapsed", sep =
                  " "))
  data
  
}

#' Title
#'
#' @param df 
#' @param variable 
#' @param width 
#'
#' @return
#' 
#' @import dplyr
#' @import data.table
#' @import gmodels
#' 
#' @export
 
binning_data_ts <- function(df,
                            variable,
                            width = 10,
                            group_condition = "group_condition") {
  message("binning TS...")
  # Get breaks and midpoints
  var <- df[[variable]]
  var_min <- min(var, na.rm = TRUE)
  var_max <- max(var, na.rm = TRUE)
  breaks <- seq(var_min, var_max, width)
  midpoints <- (breaks[2:length(breaks)] + breaks[1:(length(breaks) - 1)]) / 2
  
  # Filter data
  filtered_data <- df
  
  # Do binning
  bin <- cut(filtered_data[[variable]], breaks, labels = FALSE)
  bindf1 <- data.frame(bin = bin)
  
  if (length(midpoints) != length(sort(unique(bin)))) {
    stop("Binning failed. Try a higher width.")
  }
  
  bindf2 <- data.frame(bin = sort(unique(bin)), midpoints = midpoints)
  bindf <- left_join(bindf1, bindf2, by = "bin")
  binned_data <- filtered_data %>%
    mutate(bin = bindf$bin,
           midpoint = bindf$midpoints)
  
  # Aggregation
  data <-binned_data %>%
    group_by(across(all_of(group_condition)), bin, midpoint) %>%
    summarise(
      bearing_angle.lowCI = ci(bearing_angle, na.rm = TRUE)[2],
      bearing_angle.hiCI = ci(bearing_angle, na.rm = TRUE)[3],
      bearing_angle.mean = mean(bearing_angle, na.rm = TRUE),
      abs_bearing_angle.lowCI = ci(abs(bearing_angle), na.rm = TRUE)[2],
      abs_bearing_angle.hiCI = ci(abs(bearing_angle), na.rm = TRUE)[3],
      abs_bearing_angle.mean = mean(abs(bearing_angle), na.rm = TRUE),
      
      heading_angle.lowCI = ci(heading_angle, na.rm = TRUE)[2],
      heading_angle.hiCI = ci(heading_angle, na.rm = TRUE)[3],
      heading_angle.mean = mean(heading_angle, na.rm = TRUE),
      abs_heading_angle.lowCI = ci(abs(heading_angle), na.rm = TRUE)[2],
      abs_heading_angle.hiCI = ci(abs(heading_angle), na.rm = TRUE)[3],
      abs_heading_angle.mean = mean(abs(heading_angle), na.rm = TRUE),
      
      y_angle.lowCI = ci(y_angle, na.rm = TRUE)[2],
      y_angle.hiCI = ci(y_angle, na.rm = TRUE)[3],
      y_angle.mean = mean(y_angle, na.rm = TRUE),
      abs_y_angle.lowCI = ci(abs(y_angle), na.rm = TRUE)[2],
      abs_y_angle.hiCI = ci(abs(y_angle), na.rm = TRUE)[3],
      abs_y_angle.mean = mean(abs(y_angle), na.rm = TRUE),
      
      bending_angle.lowCI = ci(bending_angle, na.rm = TRUE)[2],
      bending_angle.hiCI = ci(bending_angle, na.rm = TRUE)[3],
      bending_angle.mean = mean(bending_angle, na.rm = TRUE),
      
      abs_bending_angle.lowCI = ci(abs(bending_angle), na.rm = TRUE)[2],
      abs_bending_angle.hiCI = ci(abs(bending_angle), na.rm = TRUE)[3],
      abs_bending_angle.mean = mean(abs(bending_angle), na.rm = TRUE),
      
      head_vector_y.lowCI = ci(head_vector_extended_y, na.rm = TRUE)[2],
      head_vector_y.hiCI = ci(head_vector_extended_y, na.rm = TRUE)[3],
      head_vector_y.mean = mean(head_vector_extended_y, na.rm = TRUE),
      
      head_vector_x.lowCI = ci(head_vector_extended_x, na.rm = TRUE)[2],
      head_vector_x.hiCI = ci(head_vector_extended_x, na.rm = TRUE)[3],
      head_vector_x.mean = mean(head_vector_extended_x, na.rm = TRUE),
      
      head_vector_angular_speed.lowCI = ci(head_vector_angular_speed, na.rm = TRUE)[2],
      head_vector_angular_speed.hiCI = ci(head_vector_angular_speed, na.rm = TRUE)[3],
      head_vector_angular_speed.mean = mean(head_vector_angular_speed, na.rm = TRUE),
      
      # abs_head_vector_angular_speed.lowCI = ci(abs_head_vector_angular_speed, na.rm = TRUE)[2],
      # abs_head_vector_angular_speed.hiCI = ci(abs_head_vector_angular_speed, na.rm = TRUE)[3],
      # abs_head_vector_angular_speed.mean = mean(abs_head_vector_angular_speed, na.rm = TRUE),
      
      tail_vector_angular_speed.lowCI = ci(tail_vector_angular_speed, na.rm = TRUE)[2],
      tail_vector_angular_speed.hiCI = ci(tail_vector_angular_speed, na.rm = TRUE)[3],
      tail_vector_angular_speed.mean = mean(tail_vector_angular_speed, na.rm = TRUE),
      
      # abs_tail_vector_angular_speed.lowCI = ci(abs_tail_vector_angular_speed, na.rm = TRUE)[2],
      # abs_tail_vector_angular_speed.hiCI = ci(abs_tail_vector_angular_speed, na.rm = TRUE)[3],
      # abs_tail_vector_angular_speed.mean = mean(abs_tail_vector_angular_speed, na.rm = TRUE),
      
      head_vel_forward.lowCI = ci(head_vel_forward, na.rm = TRUE)[2],
      head_vel_forward.hiCI = ci(head_vel_forward, na.rm = TRUE)[3],
      head_vel_forward.mean = mean(head_vel_forward, na.rm = TRUE),
      
      tail_vel_forward.lowCI = ci(tail_vel_forward, na.rm = TRUE)[2],
      tail_vel_forward.hiCI = ci(tail_vel_forward, na.rm = TRUE)[3],
      tail_vel_forward.mean = mean(tail_vel_forward, na.rm = TRUE),
      
      midpoint_speed.lowCI = ci(midpoint_speed, na.rm = TRUE)[2],
      midpoint_speed.hiCI = ci(midpoint_speed, na.rm = TRUE)[3],
      midpoint_speed.mean = mean(midpoint_speed, na.rm = TRUE),
      
      midpoint_distance.lowCI = ci(midpoint_distance, na.rm = TRUE)[2],
      midpoint_distance.hiCI = ci(midpoint_distance, na.rm = TRUE)[3],
      midpoint_distance.mean = mean(midpoint_distance, na.rm = TRUE),
      
      distance_to_odor.lowCI = ci(distance_to_odor, na.rm = TRUE)[2],
      distance_to_odor.hiCI = ci(distance_to_odor, na.rm = TRUE)[3],
      distance_to_odor.mean = mean(distance_to_odor, na.rm = TRUE),
      
      spine_length.lowCI = ci(spine_length, na.rm = TRUE)[2],
      spine_length.hiCI = ci(spine_length, na.rm = TRUE)[3],
      spine_length.mean = mean(spine_length, na.rm = TRUE),
      
      head_vel_forward_bl.lowCI = ci(head_vel_forward / spine_length.mean, na.rm = TRUE)[2],
      head_vel_forward_bl.hiCI = ci(head_vel_forward / spine_length.mean, na.rm = TRUE)[3],
      head_vel_forward_bl.mean = mean(head_vel_forward / spine_length.mean, na.rm = TRUE),
      
      tail_vel_forward_bl.lowCI = ci(tail_vel_forward / spine_length.mean, na.rm = TRUE)[2],
      tail_vel_forward_bl.hiCI = ci(tail_vel_forward / spine_length.mean, na.rm = TRUE)[3],
      tail_vel_forward_bl.mean = mean(tail_vel_forward / spine_length.mean, na.rm = TRUE),
      
      midpoint_speed_bl.lowCI = ci(midpoint_speed / spine_length.mean, na.rm = TRUE)[2],
      midpoint_speed_bl.hiCI = ci(midpoint_speed / spine_length.mean, na.rm = TRUE)[3],
      midpoint_speed_bl.mean = mean(midpoint_speed / spine_length.mean, na.rm = TRUE),
      
      midpoint_distance_bl.lowCI = ci(midpoint_distance / spine_length.mean, na.rm = TRUE)[2],
      midpoint_distance_bl.hiCI = ci(midpoint_distance / spine_length.mean, na.rm = TRUE)[3],
      midpoint_distance_bl.mean = mean(midpoint_distance / spine_length.mean, na.rm = TRUE)
    )
  print(data)
  data
}

#' Title
#'
#' @param df 
#' @param variable 
#' @param width 
#' @param grouping 
#'
#' @return
#' @export
#' @import dplyr
#' @import data.table
#' @import gmodels
 
binning_data_ts_grouped <- function(df,
                                    variable,
                                    width = 10,
                                    grouping="id",
                                    group_condition="group_condition"
){
  message("binning TS...")
  
  # Get breaks and midpoints
  var <- df[[variable]]
  var_min <- min(var, na.rm = TRUE)
  var_max <- max(var, na.rm = TRUE)
  breaks <- seq(var_min, var_max, width)
  midpoints <- (breaks[2:length(breaks)] + breaks[1:(length(breaks) - 1)]) / 2
  
  # Filter data
  filtered_data <- df
  
  # Do binning
  bin <- cut(filtered_data[[variable]], breaks, labels = FALSE)
  bindf1 <- data.frame(bin = bin)
  
  if(length(midpoints) != length(sort(unique(bin)))) {
    stop("Binning failed. Try a higher width.")
  }
  
  bindf2 <- data.frame(bin = sort(unique(bin)), midpoints = midpoints)
  bindf <- left_join(bindf1, bindf2, by = "bin")
  
  binned_data <- filtered_data %>%
    mutate(bin = factor(bindf$bin, levels = sort(unique(bindf$bin))),  # Ensure bin is a factor
           midpoint = bindf$midpoints)
  
  binned_data %>%
    group_by_at(c(group_condition, grouping, "bin", "midpoint")) %>%
    summarise(
      bearing_angle.mean = mean(bearing_angle, na.rm = TRUE),
      abs_bearing_angle.mean = mean(abs(bearing_angle), na.rm = TRUE),
      heading_angle.mean = mean(heading_angle, na.rm = TRUE),
      abs_heading_angle.mean = mean(abs(heading_angle), na.rm = TRUE),
      bending_angle.mean = mean(bending_angle, na.rm = TRUE),
      abs_bending_angle.mean = mean(abs(bending_angle), na.rm = TRUE),
      head_vector_angular_speed.mean = mean(head_vector_angular_speed, na.rm = TRUE),
      tail_vector_angular_speed.mean = mean(tail_vector_angular_speed, na.rm = TRUE),
      head_vel_forward.mean = mean(head_vel_forward, na.rm = TRUE),
      tail_vel_forward.mean = mean(tail_vel_forward, na.rm = TRUE),
      # abs_head_vel_forward.mean = mean(abs_head_vel_forward, na.rm = TRUE),
      # abs_tail_vel_forward.mean = mean(abs_tail_vel_forward, na.rm = TRUE),
      midpoint_speed.mean = mean(midpoint_speed, na.rm = TRUE),
      midpoint_distance.mean = mean(midpoint_distance, na.rm = TRUE),
      distance_to_odor.mean = mean(distance_to_odor, na.rm = TRUE),
      spine_length.mean = mean(spine_length, na.rm = TRUE),
      head_vel_forward_bl.mean = mean(head_vel_forward / spine_length, na.rm = TRUE),
      tail_vel_forward_bl.mean = mean(tail_vel_forward / spine_length, na.rm = TRUE),
      midpoint_speed_bl.mean = mean(midpoint_speed / spine_length, na.rm = TRUE),
      midpoint_distance_bl.mean = mean(midpoint_distance / spine_length, na.rm = TRUE)
    ) %>%
    group_by_at(c(group_condition, "bin", "midpoint")) %>%
    summarise(
      N = n(),
      bearing_angle.lowCI = ci(bearing_angle.mean, na.rm = TRUE)[2],
      bearing_angle.hiCI = ci(bearing_angle.mean, na.rm = TRUE)[3],
      bearing_angle.mean = mean(bearing_angle.mean, na.rm = TRUE),
      abs_bearing_angle.lowCI = ci(abs_bearing_angle.mean, na.rm = TRUE)[2],
      abs_bearing_angle.hiCI = ci(abs_bearing_angle.mean, na.rm = TRUE)[3],
      abs_bearing_angle.mean = mean(abs_bearing_angle.mean, na.rm = TRUE),
      heading_angle.lowCI = ci(heading_angle.mean, na.rm = TRUE)[2],
      heading_angle.hiCI = ci(heading_angle.mean, na.rm = TRUE)[3],
      heading_angle.mean = mean(heading_angle.mean, na.rm = TRUE),
      abs_heading_angle.lowCI = ci(abs_heading_angle.mean, na.rm = TRUE)[2],
      abs_heading_angle.hiCI = ci(abs_heading_angle.mean, na.rm = TRUE)[3],
      abs_heading_angle.mean = mean(abs_heading_angle.mean, na.rm = TRUE),
      bending_angle.lowCI = ci(bending_angle.mean, na.rm = TRUE)[2],
      bending_angle.hiCI = ci(bending_angle.mean, na.rm = TRUE)[3],
      bending_angle.mean = mean(bending_angle.mean, na.rm = TRUE),
      abs_bending_angle.lowCI = ci(abs_bending_angle.mean, na.rm = TRUE)[2],
      abs_bending_angle.hiCI = ci(abs_bending_angle.mean, na.rm = TRUE)[3],
      abs_bending_angle.mean = mean(abs_bending_angle.mean, na.rm = TRUE),
      head_vector_angular_speed.lowCI = ci(head_vector_angular_speed.mean, na.rm = TRUE)[2],
      head_vector_angular_speed.hiCI = ci(head_vector_angular_speed.mean, na.rm = TRUE)[3],
      head_vector_angular_speed.mean = mean(head_vector_angular_speed.mean, na.rm = TRUE),
      # abs_head_vector_angular_speed.lowCI = ci(abs_head_vector_angular_speed.mean, na.rm = TRUE)[2],
      # abs_head_vector_angular_speed.hiCI = ci(abs_head_vector_angular_speed.mean, na.rm = TRUE)[3],
      # abs_head_vector_angular_speed.mean = mean(abs_head_vector_angular_speed.mean, na.rm = TRUE),
      tail_vector_angular_speed.lowCI = ci(tail_vector_angular_speed.mean, na.rm = TRUE)[2],
      tail_vector_angular_speed.hiCI = ci(tail_vector_angular_speed.mean, na.rm = TRUE)[3],
      tail_vector_angular_speed.mean = mean(tail_vector_angular_speed.mean, na.rm = TRUE),
      # abs_tail_vector_angular_speed.lowCI = ci(abs_tail_vector_angular_speed.mean, na.rm = TRUE)[2],
      # abs_tail_vector_angular_speed.hiCI = ci(abs_tail_vector_angular_speed.mean, na.rm = TRUE)[3],
      # abs_tail_vector_angular_speed.mean = mean(abs_tail_vector_angular_speed.mean, na.rm = TRUE),
      head_vel_forward.lowCI = ci(head_vel_forward.mean, na.rm = TRUE)[2],
      head_vel_forward.hiCI = ci(head_vel_forward.mean, na.rm = TRUE)[3],
      head_vel_forward.mean = mean(head_vel_forward.mean, na.rm = TRUE),
      tail_vel_forward.lowCI = ci(tail_vel_forward.mean, na.rm = TRUE)[2],
      tail_vel_forward.hiCI = ci(tail_vel_forward.mean, na.rm = TRUE)[3],
      tail_vel_forward.mean = mean(tail_vel_forward.mean, na.rm = TRUE),
      midpoint_speed.lowCI = ci(midpoint_speed.mean, na.rm = TRUE)[2],
      midpoint_speed.hiCI = ci(midpoint_speed.mean, na.rm = TRUE)[3],
      midpoint_speed.mean = mean(midpoint_speed.mean, na.rm = TRUE),
      midpoint_distance.lowCI = ci(midpoint_distance.mean, na.rm = TRUE)[2],
      midpoint_distance.hiCI = ci(midpoint_distance.mean, na.rm = TRUE)[3],
      midpoint_distance.mean = mean(midpoint_distance.mean, na.rm = TRUE),
      distance_to_odor.lowCI = ci(distance_to_odor.mean, na.rm = TRUE)[2],
      distance_to_odor.hiCI = ci(distance_to_odor.mean, na.rm = TRUE)[3],
      distance_to_odor.mean = mean(distance_to_odor.mean, na.rm = TRUE),
      spine_length.lowCI = ci(spine_length.mean, na.rm = TRUE)[2],
      spine_length.hiCI = ci(spine_length.mean, na.rm = TRUE)[3],
      spine_length.mean = mean(spine_length.mean, na.rm = TRUE),
      head_vel_forward_bl.lowCI = ci(head_vel_forward_bl.mean, na.rm = TRUE)[2],
      head_vel_forward_bl.hiCI = ci(head_vel_forward_bl.mean, na.rm = TRUE)[3],
      head_vel_forward_bl.mean = mean(head_vel_forward_bl.mean, na.rm = TRUE),
      tail_vel_forward_bl.lowCI = ci(tail_vel_forward_bl.mean, na.rm = TRUE)[2],
      tail_vel_forward_bl.hiCI = ci(tail_vel_forward_bl.mean, na.rm = TRUE)[3],
      tail_vel_forward_bl.mean = mean(tail_vel_forward_bl.mean, na.rm = TRUE),
      midpoint_speed_bl.lowCI = ci(midpoint_speed_bl.mean, na.rm = TRUE)[2],
      midpoint_speed_bl.hiCI = ci(midpoint_speed_bl.mean, na.rm = TRUE)[3],
      midpoint_speed_bl.mean = mean(midpoint_speed_bl.mean, na.rm = TRUE),
      midpoint_distance_bl.lowCI = ci(midpoint_distance_bl.mean, na.rm = TRUE)[2],
      midpoint_distance_bl.hiCI = ci(midpoint_distance_bl.mean, na.rm = TRUE)[3],
      midpoint_distance_bl.mean = mean(midpoint_distance_bl.mean, na.rm = TRUE)
    )
}

#' Title
#'
#' @param df 
#' @param variable 
#' @param width 
#' @param Abs_HC_Angle_interval 
#' @param frame_rate 
#'
#' @return
#' @export
#' @import dplyr
#' @import data.table
#' @import gmodels
 
binning_data_hc <- function(df,
                            variable,
                            width = 10,
                            Abs_HC_Angle_interval = c(20, 360),
                            frame_rate = 16,
                            group_condition = "group_condition") {
  
  message("binning HC data")
  
  # Get breaks and midpoints
  var <- df[[variable]]
  var_min <- min(var, na.rm = TRUE)
  var_max <- max(var, na.rm = TRUE)
  breaks <- seq(var_min, var_max, width)
  midpoints <- (breaks[2:length(breaks)] + breaks[1:(length(breaks) - 1)]) / 2
  
  # Filter data
  filtered_data <- df %>%
    filter(
      is.na(Abs_HC_angle) |
        (Abs_HC_angle >= Abs_HC_Angle_interval[1] &
           Abs_HC_angle <= Abs_HC_Angle_interval[2])
    )
  
  # Perform binning
  bin <- cut(filtered_data[[variable]], breaks, labels = FALSE)
  bindf1 <- data.frame(bin = bin)
  
  if (length(midpoints) != length(sort(unique(bin)))) {
    stop("Binning failed. Try a higher width.")
  }
  
  bindf2 <- data.frame(bin = sort(unique(bin)), midpoints = midpoints)
  bindf <- left_join(bindf1, bindf2, by = "bin")
  
  binned_data <- filtered_data %>%
    mutate(
      bin = bindf$bin,
      midpoint = bindf$midpoints
    )
  
  # Aggregation
  binned_data %>%
    group_by(.data[[group_condition]], bin, midpoint) %>%
    summarise(
      N = n(),
      HC_rate = (sum(HCs, na.rm = TRUE) / n()) * frame_rate,
      
      HC_angle.lowCI = ci(HC_angle, na.rm = TRUE)[2],
      HC_angle.hiCI = ci(HC_angle, na.rm = TRUE)[3],
      HC_angle.mean = mean(HC_angle, na.rm = TRUE),
      
      Abs_HC_angle.lowCI = ci(Abs_HC_angle, na.rm = TRUE)[2],
      Abs_HC_angle.hiCI = ci(Abs_HC_angle, na.rm = TRUE)[3],
      Abs_HC_angle.mean = mean(Abs_HC_angle, na.rm = TRUE),
      
      HC_reorientation.lowCI = ci(HC_reorientation, na.rm = TRUE)[2],
      HC_reorientation.hiCI = ci(HC_reorientation, na.rm = TRUE)[3],
      HC_reorientation.mean = mean(HC_reorientation, na.rm = TRUE),
      
      HC_rate_modulation = (sum(HCs[odor_orientation == "towards"], na.rm = TRUE) -
                              sum(HCs[odor_orientation == "away"], na.rm = TRUE)) /
        (sum(HCs[odor_orientation == "towards"], na.rm = TRUE) +
           sum(HCs[odor_orientation == "away"], na.rm = TRUE))
    )
}


#' Title
#'
#' @param df 
#' @param variable 
#' @param width 
#' @param grouping 
#' @param Abs_HC_Angle_interval 
#' @param frame_rate 
#'
#' @return
#' @export
#' @import dplyr
#' @import data.table
#' @import gmodels
binning_data_hc_grouped <- function(df,
                                    variable,
                                    width = 10,
                                    grouping = "id",
                                    Abs_HC_Angle_interval = c(20, 360),
                                    frame_rate = 16,
                                    group_condition = "group_condition") {
  
  message("binning HC data")
  
  # Get breaks and midpoints
  var <- df[[variable]]
  var_min <- min(var, na.rm = TRUE)
  var_max <- max(var, na.rm = TRUE)
  breaks <- seq(var_min, var_max, width)
  midpoints <- (breaks[2:length(breaks)] + breaks[1:(length(breaks) - 1)]) / 2
  
  # Filter data
  filtered_data <- df %>%
    filter(
      is.na(Abs_HC_angle) |
        (Abs_HC_angle >= Abs_HC_Angle_interval[1] &
           Abs_HC_angle <= Abs_HC_Angle_interval[2])
    )
  
  # Perform binning
  bin <- cut(filtered_data[[variable]], breaks, labels = FALSE)
  bindf1 <- data.frame(bin = bin)
  
  if (length(midpoints) != length(sort(unique(bin)))) {
    stop("Binning failed. Try a higher width.")
  }
  
  bindf2 <- data.frame(bin = sort(unique(bin)), midpoints = midpoints)
  bindf <- left_join(bindf1, bindf2, by = "bin")
  
  binned_data <- filtered_data %>%
    mutate(
      bin = bindf$bin,
      midpoint = bindf$midpoints
    )
  
  # Aggregation with grouping
  binned_data %>%
    group_by(.data[[group_condition]], .data[[grouping]], bin, midpoint) %>%
    summarise(
      N = n(),
      HC_rate = (sum(HCs, na.rm = TRUE) / n()) * frame_rate,
      HC_angle.mean = mean(HC_angle, na.rm = TRUE),
      Abs_HC_angle.mean = mean(Abs_HC_angle, na.rm = TRUE),
      HC_reorientation.mean = mean(HC_reorientation, na.rm = TRUE),
      
      HC_rate_modulation = (sum(HCs[odor_orientation == "towards"], na.rm = TRUE) -
                              sum(HCs[odor_orientation == "away"], na.rm = TRUE)) /
        (sum(HCs[odor_orientation == "towards"], na.rm = TRUE) +
           sum(HCs[odor_orientation == "away"], na.rm = TRUE))
    ) %>%
    group_by(.data[[group_condition]], bin, midpoint) %>%
    summarise(
      HC_rate.lowCI = ci(HC_rate, na.rm = TRUE)[2],
      HC_rate.hiCI = ci(HC_rate, na.rm = TRUE)[3],
      HC_rate.mean = mean(HC_rate, na.rm = TRUE),
      
      HC_angle.lowCI = ci(HC_angle.mean, na.rm = TRUE)[2],
      HC_angle.hiCI = ci(HC_angle.mean, na.rm = TRUE)[3],
      HC_angle.mean = mean(HC_angle.mean, na.rm = TRUE),
      
      Abs_HC_angle.lowCI = ci(Abs_HC_angle.mean, na.rm = TRUE)[2],
      Abs_HC_angle.hiCI = ci(Abs_HC_angle.mean, na.rm = TRUE)[3],
      Abs_HC_angle.mean = mean(Abs_HC_angle.mean, na.rm = TRUE),
      
      HC_reorientation.lowCI = ci(HC_reorientation.mean, na.rm = TRUE)[2],
      HC_reorientation.hiCI = ci(HC_reorientation.mean, na.rm = TRUE)[3],
      HC_reorientation.mean = mean(HC_reorientation.mean, na.rm = TRUE),
      
      HC_rate_modulation.lowCI = ci(HC_rate_modulation, na.rm = TRUE)[2],
      HC_rate_modulation.hiCI = ci(HC_rate_modulation, na.rm = TRUE)[3],
      HC_rate_modulation.mean = mean(HC_rate_modulation, na.rm = TRUE)
    )
}
#' Title
#'
#' @param df 
#' @param variable 
#' @param width 
#'
#' @return
#' @export
#' @import dplyr
#' @import data.table
#' @import gmodels
 
binning_data_is <- function(df,
                            variable,
                            width = 10,
                            group_condition = "group_condition") {
  
  message("binning IS data")
  
  # Get breaks and midpoints
  var <- df[[variable]]
  var_min <- min(var, na.rm = TRUE)
  var_max <- max(var, na.rm = TRUE)
  breaks <- seq(var_min, var_max, width)
  midpoints <- (breaks[2:length(breaks)] + breaks[1:length(breaks) - 1]) / 2
  
  # Perform binning
  filtered_data <- df  # No filter on step_boolean is applied here
  bin <- cut(filtered_data[[variable]], breaks, labels = FALSE, right = TRUE)
  bindf1 <- data.frame(bin = bin)
  
  # Check if binning was successful
  if (length(midpoints) != length(sort(unique(bin)))) {
    stop("Binning failed. Try a higher width.")
  }
  
  # Create data frame with bin midpoints
  bindf2 <- data.frame(bin = sort(unique(bin)), midpoints = midpoints)
  bindf <- left_join(bindf1, bindf2, by = "bin")
  
  # Add bin and midpoint columns to the filtered data
  binned_data <- filtered_data %>%
    mutate(bin = bindf$bin,
           midpoint = bindf$midpoints)
  
  # Aggregation with summarisation of relevant variables
  binned_data %>%
    group_by(across(all_of(group_condition)), bin, midpoint) %>%
    summarise(
      IS_angle.lowCI = ci(IS_angle, na.rm = TRUE)[2],
      IS_angle.hiCI = ci(IS_angle, na.rm = TRUE)[3],
      IS_angle.mean = mean(IS_angle, na.rm = TRUE),
      
      Abs_IS_angle.lowCI = ci(Abs_IS_angle, na.rm = TRUE)[2],
      Abs_IS_angle.hiCI = ci(Abs_IS_angle, na.rm = TRUE)[3],
      Abs_IS_angle.mean = mean(Abs_IS_angle, na.rm = TRUE),
      
      IS_reorientation.lowCI = ci(IS_reorientation, na.rm = TRUE)[2],
      IS_reorientation.hiCI = ci(IS_reorientation, na.rm = TRUE)[3],
      IS_reorientation.mean = mean(IS_reorientation, na.rm = TRUE),
      
      IS_distance.lowCI = ci(IS_distance, na.rm = TRUE)[2],
      IS_distance.hiCI = ci(IS_distance, na.rm = TRUE)[3],
      IS_distance.mean = mean(IS_distance, na.rm = TRUE),
      
      IS_distance_bl.lowCI = ci(IS_distance_bl, na.rm = TRUE)[2],
      IS_distance_bl.hiCI = ci(IS_distance_bl, na.rm = TRUE)[3],
      IS_distance_bl.mean = mean(IS_distance_bl, na.rm = TRUE),
      
      IS_interval.lowCI = ci(IS_interval, na.rm = TRUE)[2],
      IS_interval.hiCI = ci(IS_interval, na.rm = TRUE)[3],
      IS_interval.mean = mean(IS_interval, na.rm = TRUE)
    )
}


#' Title
#'
#' @param df 
#' @param variable 
#' @param width 
#' @param grouping 
#'
#' @return
#' @export
#' @import dplyr
#' @import data.table
#' @import gmodels
 
binning_data_is_grouped <- function(df,
                                    variable,
                                    width = 10,
                                    grouping = "id",
                                    group_condition = "group_condition") {
  
  message("binning IS data")
  
  # Get breaks and midpoints
  var <- df[[variable]]
  var_min <- min(var, na.rm = TRUE)
  var_max <- max(var, na.rm = TRUE)
  breaks <- seq(var_min, var_max, width)
  midpoints <- (breaks[2:length(breaks)] + breaks[1:length(breaks) - 1]) / 2
  
  # Perform binning
  filtered_data <- df  # No filter on step_boolean applied here
  bin <- cut(filtered_data[[variable]], breaks, labels = FALSE)
  bindf1 <- data.frame(bin = bin)
  
  # Check if binning was successful
  if (length(midpoints) != length(sort(unique(bin)))) {
    stop("Binning failed. Try a higher width.")
  }
  
  # Create data frame with bin midpoints
  bindf2 <- data.frame(bin = sort(unique(bin)), midpoints = midpoints)
  bindf <- left_join(bindf1, bindf2, by = "bin")
  
  # Add bin and midpoint columns to the filtered data
  binned_data <- filtered_data %>%
    mutate(bin = bindf$bin,
           midpoint = bindf$midpoints)
  
  # First summarization: group by `group_condition`, `grouping`, `bin`, and `midpoint`
  first_summary <- binned_data %>%
    group_by(across(all_of(group_condition)), across(all_of(grouping)), bin, midpoint) %>%
    summarise(
      IS_angle.mean = mean(IS_angle, na.rm = TRUE),
      Abs_IS_angle.mean = mean(Abs_IS_angle, na.rm = TRUE),
      IS_reorientation.mean = mean(IS_reorientation, na.rm = TRUE),
      IS_distance.mean = mean(IS_distance, na.rm = TRUE),
      IS_distance_bl.mean = mean(IS_distance_bl, na.rm = TRUE),
      IS_interval.mean = mean(IS_interval, na.rm = TRUE)
    )
  
  # Second summarization: group by `group_condition`, `bin`, and `midpoint`
  final_summary <- first_summary %>%
    group_by(across(all_of(group_condition)), bin, midpoint) %>%
    summarise(
      IS_angle.lowCI = ci(IS_angle.mean, na.rm = TRUE)[2],
      IS_angle.hiCI = ci(IS_angle.mean, na.rm = TRUE)[3],
      IS_angle.mean = mean(IS_angle.mean, na.rm = TRUE),
      
      Abs_IS_angle.lowCI = ci(Abs_IS_angle.mean, na.rm = TRUE)[2],
      Abs_IS_angle.hiCI = ci(Abs_IS_angle.mean, na.rm = TRUE)[3],
      Abs_IS_angle.mean = mean(Abs_IS_angle.mean, na.rm = TRUE),
      
      IS_reorientation.lowCI = ci(IS_reorientation.mean, na.rm = TRUE)[2],
      IS_reorientation.hiCI = ci(IS_reorientation.mean, na.rm = TRUE)[3],
      IS_reorientation.mean = mean(IS_reorientation.mean, na.rm = TRUE),
      
      IS_distance.lowCI = ci(IS_distance.mean, na.rm = TRUE)[2],
      IS_distance.hiCI = ci(IS_distance.mean, na.rm = TRUE)[3],
      IS_distance.mean = mean(IS_distance.mean, na.rm = TRUE),
      
      IS_distance_bl.lowCI = ci(IS_distance_bl.mean, na.rm = TRUE)[2],
      IS_distance_bl.hiCI = ci(IS_distance_bl.mean, na.rm = TRUE)[3],
      IS_distance_bl.mean = mean(IS_distance_bl.mean, na.rm = TRUE),
      
      IS_interval.lowCI = ci(IS_interval.mean, na.rm = TRUE)[2],
      IS_interval.hiCI = ci(IS_interval.mean, na.rm = TRUE)[3],
      IS_interval.mean = mean(IS_interval.mean, na.rm = TRUE)
    )
  
  return(final_summary)
}



#' Title
#'
#' @param df 
#' @param variable 
#' @param width 
#'
#' @return
#' @export
#' @import dplyr
#' @import data.table
#' @import gmodels
 
binning_data_run_speed <- function(df,
                                   variable,
                                   width = 10,
                                   group_condition = "group_condition") {
  
  message("binning run_speed data")
  
  # Get breaks and midpoints
  var <- df[[variable]]
  var_min <- min(var, na.rm = TRUE)
  var_max <- max(var, na.rm = TRUE)
  breaks <- seq(var_min, var_max, width)
  midpoints <- (breaks[2:length(breaks)] + breaks[1:length(breaks) - 1]) / 2
  
  # Perform binning
  filtered_data <- df  # No filter on run_flag applied here
  bin <- cut(filtered_data[[variable]], breaks, labels = FALSE)
  bindf1 <- data.frame(bin = bin)
  
  # Check if binning was successful
  if (length(midpoints) != length(sort(unique(bin)))) {
    stop("Binning failed. Try a higher width.")
  }
  
  # Create data frame with bin midpoints
  bindf2 <- data.frame(bin = sort(unique(bin)), midpoints = midpoints)
  bindf <- left_join(bindf1, bindf2, by = "bin")
  
  # Add bin and midpoint columns to the filtered data
  binned_data <- filtered_data %>%
    mutate(bin = bindf$bin,
           midpoint = bindf$midpoints)
  
  # Aggregation step
  binned_data %>%
    group_by(across(all_of(group_condition)), bin, midpoint) %>%
    summarise(
      run_speed.lowCI = ci(run_speed, na.rm = TRUE)[2],
      run_speed.hiCI = ci(run_speed, na.rm = TRUE)[3],
      run_speed.mean = mean(run_speed, na.rm = TRUE),
      
      run_speed_bl.lowCI = ci(run_speed_bl, na.rm = TRUE)[2],
      run_speed_bl.hiCI = ci(run_speed_bl, na.rm = TRUE)[3],
      run_speed_bl.mean = mean(run_speed_bl, na.rm = TRUE),
      
      run_speed_modulation = (mean(run_speed[odor_orientation == "towards"], na.rm = TRUE) -
                                mean(run_speed[odor_orientation == "away"], na.rm = TRUE)) /
        (mean(run_speed[odor_orientation == "towards"], na.rm = TRUE) +
           mean(run_speed[odor_orientation == "away"], na.rm = TRUE))
    )
}

#' Title
#'
#' @param df 
#' @param variable 
#' @param width 
#' @param grouping 
#'
#' @return
#' @export
#' @import dplyr
#' @import data.table
#' @import gmodels
 
#' 
binning_data_run_speed_grouped <- function(df,
                                           variable,
                                           width = 10,
                                           grouping = "id",
                                           group_condition = "group_condition") {
  
  message("binning run_speed data")
  
  # Get breaks and midpoints
  var <- df[[variable]]
  var_min <- min(var, na.rm = TRUE)
  var_max <- max(var, na.rm = TRUE)
  breaks <- seq(var_min, var_max, width)
  midpoints <- (breaks[2:length(breaks)] + breaks[1:length(breaks) - 1]) / 2
  
  # Perform binning
  filtered_data <- df  # No filter on run_flag applied here
  bin <- cut(filtered_data[[variable]], breaks, labels = FALSE)
  bindf1 <- data.frame(bin = bin)
  
  # Check if binning was successful
  if (length(midpoints) != length(sort(unique(bin)))) {
    stop("Binning failed. Try a higher width.")
  }
  
  # Create data frame with bin midpoints
  bindf2 <- data.frame(bin = sort(unique(bin)), midpoints = midpoints)
  bindf <- left_join(bindf1, bindf2, by = "bin")
  
  # Add bin and midpoint columns to the filtered data
  binned_data <- filtered_data %>%
    mutate(bin = bindf$bin,
           midpoint = bindf$midpoints)
  
  # Aggregation step
  binned_data %>%
    group_by(across(all_of(group_condition)), across(all_of(grouping)), bin, midpoint) %>%
    summarise(
      run_speed = mean(run_speed, na.rm = TRUE),
      run_speed_bl = mean(run_speed_bl, na.rm = TRUE),
      run_speed_modulation = (mean(run_speed[odor_orientation == "towards"], na.rm = TRUE) -
                                mean(run_speed[odor_orientation == "away"], na.rm = TRUE)) /
        (mean(run_speed[odor_orientation == "towards"], na.rm = TRUE) +
           mean(run_speed[odor_orientation == "away"], na.rm = TRUE))
    ) %>%
    group_by(across(all_of(group_condition)), bin, midpoint) %>%
    summarise(
      run_speed.lowCI = ci(run_speed, na.rm = TRUE)[2],
      run_speed.hiCI = ci(run_speed, na.rm = TRUE)[3],
      run_speed.mean = mean(run_speed, na.rm = TRUE),
      
      run_speed_bl.lowCI = ci(run_speed_bl, na.rm = TRUE)[2],
      run_speed_bl.hiCI = ci(run_speed_bl, na.rm = TRUE)[3],
      run_speed_bl.mean = mean(run_speed_bl, na.rm = TRUE),
      
      run_speed_modulation.lowCI = ci(run_speed_modulation, na.rm = TRUE)[2],
      run_speed_modulation.hiCI = ci(run_speed_modulation, na.rm = TRUE)[3],
      run_speed_modulation.mean = mean(run_speed_modulation, na.rm = TRUE)
    )
}

#' Title
#'
#' @param df 
#' @param variable 
#' @param radius 
#' @param width 
#'
#' @return
#' @export
#' @import dplyr
#' @import data.table
#' @import gmodels
 
#' 
binning_data_preference <- function(df, variable, radius, width = 10, group_condition = "group_condition") {
  
  message("binning preference data")
  
  # Get breaks and midpoints for the binning process
  var <- df[[variable]]
  var_min <- min(var, na.rm = TRUE)
  var_max <- max(var, na.rm = TRUE)
  breaks <- seq(var_min, var_max, width)
  midpoints <- (breaks[2:length(breaks)] + breaks[1:length(breaks) - 1]) / 2
  
  # Perform binning
  filtered_data <- df  # No additional filtering applied here
  bin <- cut(filtered_data[[variable]], breaks, labels = FALSE)
  bindf1 <- data.frame(bin = bin)
  
  # Check if binning was successful
  if (length(midpoints) != length(sort(unique(bin)))) {
    stop("Binning failed. Try a higher width.")
  }
  
  # Create data frame with bin midpoints
  bindf2 <- data.frame(bin = sort(unique(bin)), midpoints = midpoints)
  bindf <- left_join(bindf1, bindf2, by = "bin")
  
  # Add bin and midpoint columns to the filtered data
  binned_data <- filtered_data %>%
    mutate(bin = bindf$bin,
           midpoint = bindf$midpoints)
  
  # Aggregation step
  binned_data %>%
    group_by(across(all_of(group_condition)), bin, midpoint) %>%
    summarise(
      preference = mean(pref, na.rm = TRUE),
      
      preference_dist.lowCI = -1 * rescale(
        ci(distance_to_odor, na.rm = TRUE)[2],
        to = c(-1, 1),
        from = c(0, radius * 2)
      ),
      
      preference_dist.hiCI = -1 * rescale(
        ci(distance_to_odor, na.rm = TRUE)[3],
        to = c(-1, 1),
        from = c(0, radius * 2)
      ),
      
      preference_dist.mean = -1 * rescale(
        mean(distance_to_odor, na.rm = TRUE),
        to = c(-1, 1),
        from = c(0, radius * 2)
      ),
      
      odor_speed.lowCI = ci(odor_speed, na.rm = TRUE)[2],
      odor_speed.hiCI = ci(odor_speed, na.rm = TRUE)[3],
      odor_speed.mean = mean(odor_speed, na.rm = TRUE),
      
      ratio_towards = sum(odor_orientation == "towards") / (
        sum(odor_orientation == "towards") + sum(odor_orientation == "away")
      )
    )
}

#' Title
#'
#' @param df 
#' @param variable 
#' @param radius 
#' @param width 
#' @param grouping 
#'
#' @return
#' @export
#' @import dplyr
#' @import data.table
#' @import gmodels
 
#' 
binning_data_preference_grouped <- function(df, variable, radius, width = 10, grouping = "id") {
  
  message("Binning preference data")
  
  # Get breaks and midpoints for the binning process
  var <- df[[variable]]
  var_min <- min(var, na.rm = TRUE)
  var_max <- max(var, na.rm = TRUE)
  breaks <- seq(var_min, var_max, width)
  midpoints <- (breaks[2:length(breaks)] + breaks[1:length(breaks) - 1]) / 2
  
  # Perform binning
  filtered_data <- df
  bin <- cut(filtered_data[[variable]], breaks, labels = FALSE)
  bindf1 <- data.frame(bin = bin)
  
  # Check if binning was successful
  if (length(midpoints) != length(sort(unique(bin)))) {
    stop("Binning failed. Try a higher width.")
  }
  
  # Create data frame with bin midpoints
  bindf2 <- data.frame(bin = sort(unique(bin)), midpoints = midpoints)
  bindf <- left_join(bindf1, bindf2, by = "bin")
  
  # Add bin and midpoint columns to the filtered data
  binned_data <- filtered_data %>%
    mutate(bin = bindf$bin,
           midpoint = bindf$midpoints)
  
  # Aggregation - first group by the specified 'grouping' variable
  aggregated_data <- binned_data %>%
    group_by(across(all_of(grouping)), bin, midpoint) %>%
    summarise(
      preference = mean(pref, na.rm = TRUE),
      
      preference_dist = -1 * rescale(
        mean(distance_to_odor, na.rm = TRUE),
        to = c(-1, 1),
        from = c(0, radius * 2)
      ),
      
      odor_speed = mean(odor_speed, na.rm = TRUE),
      
      ratio_towards = sum(odor_orientation == "towards", na.rm = TRUE) / (
        sum(odor_orientation == "towards", na.rm = TRUE) +
          sum(odor_orientation == "away", na.rm = TRUE)
      ),
      .groups = 'drop'
    )
  
  # Further aggregation for mean and CI of other variables
  aggregated_data %>%
    group_by(across(all_of(grouping)), bin, midpoint) %>%
    summarise(
      preference.mean = mean(preference, na.rm = TRUE),
      preference.lowCI = ci(preference, na.rm = TRUE)[2],
      preference.hiCI = ci(preference, na.rm = TRUE)[3],
      
      preference_dist.mean = mean(preference_dist, na.rm = TRUE),
      preference_dist.lowCI = ci(preference_dist, na.rm = TRUE)[2],
      preference_dist.hiCI = ci(preference_dist, na.rm = TRUE)[3],
      
      odor_speed.mean = mean(odor_speed, na.rm = TRUE),
      odor_speed.lowCI = ci(odor_speed, na.rm = TRUE)[2],
      odor_speed.hiCI = ci(odor_speed, na.rm = TRUE)[3],
      
      ratio_towards.mean = mean(ratio_towards, na.rm = TRUE),
      .groups = 'drop'
    )
}
