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
#'
 
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
                         direction = NULL
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
    left_join(binning_data_hc(df, variable, width, Abs_HC_Angle_interval,frame_rate=frame_rate)) %>%
    left_join(binning_data_is(df, variable,  width)) %>%
    left_join(binning_data_run_speed(df, variable, width)) %>%
    left_join(binning_data_preference(df, variable, radius, width))
  
  
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
                                 direction = NULL
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
  
  
  data <-binning_data_ts_grouped(df,variable,width,grouping=grouping)%>%
    left_join(binning_data_hc_grouped(df, variable, width, grouping=grouping,Abs_HC_Angle_interval,frame_rate=frame_rate)) %>%
    left_join(binning_data_is_grouped(df, variable, width,grouping=grouping)) %>%
    left_join(binning_data_run_speed_grouped(df, variable, width,grouping=grouping)) %>%
    left_join(binning_data_preference_grouped(df, variable, radius, width,grouping=grouping))
  
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
 
binning_data_ts <-   function(df,
                              variable,
                              width = 10){
  message("binning TS...")
  #get breaks and midpoints
  var <- df[[variable]]
  var_min <- min(var,na.rm=T)
  var_max <- max(var,na.rm=T)
  breaks <- seq(var_min,var_max,width)
  midpoints <- (breaks[2:length(breaks)]+breaks[1:length(breaks)-1])/2
  
  #filter data 
  filtered_data <- df
  
  #do binning
  bin <- cut(filtered_data[[variable]],breaks,labels=F)
  bindf1 <- data.frame(list(bin=bin))
  
  if(length(midpoints)!=length(sort(unique(bin)))){
    stop("Binning failed.Try a higher width.")
  }
  bindf2 <- data.frame(list(bin = sort(unique(bin)),midpoints = midpoints))
  bindf <- left_join(bindf1,bindf2,by="bin")
  binned_data <-  filtered_data%>%
    mutate(bin = bindf$bin,
           midpoint = bindf$midpoints
    )
  
  #aggregation
  binned_data %>%
    group_by(group_condition,bin,midpoint) %>%
    summarise(
      
      bearing_angle.lowCI = ci(bearing_angle, na.rm = T)[2],
      bearing_angle.hiCI = ci(bearing_angle, na.rm = T)[3],
      bearing_angle.mean = mean(bearing_angle,na.rm=T),
      abs_bearing_angle.lowCI = ci(abs(bearing_angle), na.rm = T)[2],
      abs_bearing_angle.hiCI = ci(abs(bearing_angle), na.rm = T)[3],
      abs_bearing_angle.mean = mean(abs(bearing_angle),na.rm=T),
      
      heading_angle.lowCI = ci(heading_angle, na.rm = T)[2],
      heading_angle.hiCI = ci(heading_angle, na.rm = T)[3],
      heading_angle.mean = mean(heading_angle,na.rm=T),
      abs_heading_angle.lowCI = ci(abs(heading_angle), na.rm = T)[2],
      abs_heading_angle.hiCI = ci(abs(heading_angle), na.rm = T)[3],
      abs_heading_angle.mean = mean(abs(heading_angle),na.rm=T),

      y_angle.lowCI = ci(y_angle, na.rm = T)[2],
      y_angle.hiCI = ci(y_angle, na.rm = T)[3],
      y_angle.mean = mean(y_angle,na.rm=T),
      abs_y_angle.lowCI = ci(abs(y_angle), na.rm = T)[2],
      abs_y_angle.hiCI = ci(abs(y_angle), na.rm = T)[3],
      abs_y_angle.mean = mean(abs(y_angle),na.rm=T),
      
      bending_angle.lowCI = ci(bending_angle, na.rm = T)[2],
      bending_angle.hiCI = ci(bending_angle, na.rm = T)[3],
      bending_angle.mean = mean(bending_angle,na.rm=T),
      
      abs_bending_angle.lowCI = ci(abs(bending_angle), na.rm = T)[2],
      abs_bending_angle.hiCI = ci(abs(bending_angle), na.rm = T)[3],
      abs_bending_angle.mean = mean(abs(bending_angle),na.rm=T),
      
      head_vector_y.lowCI = ci(head_vector_extended_y, na.rm = T)[2],
      head_vector_y.hiCI = ci(head_vector_extended_y, na.rm = T)[3],
      head_vector_y.mean = mean(head_vector_extended_y,na.rm=T),
      
      head_vector_x.lowCI = ci(head_vector_extended_x, na.rm = T)[2],
      head_vector_x.hiCI = ci(head_vector_extended_x, na.rm = T)[3],
      head_vector_x.mean = mean(head_vector_extended_x,na.rm=T),
      
      head_vector_angular_speed.lowCI = ci(head_vector_angular_speed, na.rm = T)[2],
      head_vector_angular_speed.hiCI = ci(head_vector_angular_speed, na.rm = T)[3],
      head_vector_angular_speed.mean = mean(head_vector_angular_speed,na.rm=T),
      
      abs_head_vector_angular_speed.lowCI = ci(abs_head_vector_angular_speed, na.rm = T)[2],
      abs_head_vector_angular_speed.hiCI = ci(abs_head_vector_angular_speed, na.rm = T)[3],
      abs_head_vector_angular_speed.mean = mean(abs_head_vector_angular_speed,na.rm=T),
      
      tail_vector_angular_speed.lowCI = ci(tail_vector_angular_speed, na.rm = T)[2],
      tail_vector_angular_speed.hiCI = ci(tail_vector_angular_speed, na.rm = T)[3],
      tail_vector_angular_speed.mean = mean(tail_vector_angular_speed,na.rm=T),
      
      abs_tail_vector_angular_speed.lowCI = ci(abs_tail_vector_angular_speed, na.rm = T)[2],
      abs_tail_vector_angular_speed.hiCI = ci(abs_tail_vector_angular_speed, na.rm = T)[3],
      abs_tail_vector_angular_speed.mean = mean(abs_tail_vector_angular_speed,na.rm=T),
      
      
      head_vel_forward.lowCI = ci(head_vel_forward, na.rm = T)[2],
      head_vel_forward.hiCI = ci(head_vel_forward, na.rm = T)[3],
      head_vel_forward.mean = mean(head_vel_forward,na.rm=T),
      
      tail_vel_forward.lowCI = ci(tail_vel_forward, na.rm = T)[2],
      tail_vel_forward.hiCI = ci(tail_vel_forward, na.rm = T)[3],
      tail_vel_forward.mean = mean(tail_vel_forward,na.rm=T),
      
      midpoint_speed.lowCI = ci(midpoint_speed, na.rm = T)[2],
      midpoint_speed.hiCI = ci(midpoint_speed, na.rm = T)[3],
      midpoint_speed.mean = mean(midpoint_speed,na.rm=T),
      
      midpoint_distance.lowCI = ci(midpoint_distance, na.rm = T)[2],
      midpoint_distance.hiCI = ci(midpoint_distance, na.rm = T)[3],
      midpoint_distance.mean = mean(midpoint_distance,na.rm=T),
      
      distance_to_odor.lowCI = ci(distance_to_odor, na.rm = T)[2],
      distance_to_odor.hiCI = ci(distance_to_odor, na.rm = T)[3],
      distance_to_odor.mean = mean(distance_to_odor,na.rm=T),
      
      spine_length.lowCI = ci(spine_length,na.rm=T)[2],
      spine_length.hiCI = ci(spine_length,na.rm=T)[3],
      spine_length.mean = mean(spine_length,na.rm=T),
      
      head_vel_forward_bl.lowCI = ci(head_vel_forward/spine_length.mean, na.rm = T)[2],
      head_vel_forward_bl.hiCI = ci(head_vel_forward/spine_length.mean, na.rm = T)[3],
      head_vel_forward_bl.mean = mean(head_vel_forward/spine_length.mean,na.rm=T),
      
      tail_vel_forward_bl.lowCI = ci(tail_vel_forward/spine_length.mean, na.rm = T)[2],
      tail_vel_forward_bl.hiCI = ci(tail_vel_forward/spine_length.mean, na.rm = T)[3],
      tail_vel_forward_bl.mean = mean(tail_vel_forward/spine_length.mean,na.rm=T),
      
      midpoint_speed_bl.lowCI = ci(midpoint_speed/spine_length.mean, na.rm = T)[2],
      midpoint_speed_bl.hiCI = ci(midpoint_speed/spine_length.mean, na.rm = T)[3],
      midpoint_speed_bl.mean = mean(midpoint_speed/spine_length.mean,na.rm=T),
      
      midpoint_distance_bl.lowCI = ci(midpoint_distance/spine_length.mean, na.rm = T)[2],
      midpoint_distance_bl.hiCI = ci(midpoint_distance/spine_length.mean, na.rm = T)[3],
      midpoint_distance_bl.mean = mean(midpoint_distance/spine_length.mean,na.rm=T)
      
      
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
 
binning_data_ts_grouped <-   function(df,
                                      variable,
                                      width = 10,
                                      grouping="id"
){
  message("binning TS...")
  
  #get breaks and midpoints
  var <- df[[variable]]
  var_min <- min(var,na.rm=T)
  var_max <- max(var,na.rm=T)
  breaks <- seq(var_min,var_max,width)
  midpoints <- (breaks[2:length(breaks)]+breaks[1:length(breaks)-1])/2
  
  #filter data 
  filtered_data <- df
  
  #do binning
  bin <- cut(filtered_data[[variable]],breaks,labels=F)
  bindf1 <- data.frame(list(bin=bin))
  if(length(midpoints)!=length(sort(unique(bin)))){
    stop("Binning failed.Try a higher width.")
  }
  bindf2 <- data.frame(list(bin = sort(unique(bin)),midpoints = midpoints))
  bindf <- left_join(bindf1,bindf2,by="bin")
  binned_data <-  filtered_data%>%
    mutate(bin = bindf$bin,
           midpoint = bindf$midpoints
    )
  
  
  binned_data %>%
    group_by_("group_condition",grouping,"bin","midpoint") %>%
    summarise(
      bearing_angle.mean = mean(bearing_angle,na.rm=T),
      abs_bearing_angle.mean = mean(abs(bearing_angle),na.rm=T),
      heading_angle.mean = mean(heading_angle,na.rm=T),
      abs_heading_angle.mean = mean(abs(heading_angle),na.rm=T),
      bending_angle.mean = mean(bending_angle,na.rm=T),
      abs_bending_angle.mean = mean(abs(bending_angle),na.rm=T),
      head_vector_angular_speed.mean = mean(head_vector_angular_speed,na.rm=T),
      tail_vector_angular_speed.mean = mean(tail_vector_angular_speed,na.rm=T),
      head_vel_forward.mean = mean(head_vel_forward,na.rm=T),
      tail_vel_forward.mean = mean(tail_vel_forward,na.rm=T),
      abs_head_vel_forward.mean = mean(abs_head_vel_forward,na.rm=T),
      abs_tail_vel_forward.mean = mean(abs_tail_vel_forward,na.rm=T),
      
      
      midpoint_speed.mean = mean(midpoint_speed,na.rm=T),
      midpoint_distance.mean = mean(midpoint_distance,na.rm=T),
      distance_to_odor.mean = mean(distance_to_odor,na.rm=T),
      spine_length.mean = mean(spine_length,na.rm=T),
      head_vel_forward_bl.mean = mean(head_vel_forward/spine_length,na.rm=T),
      tail_vel_forward_bl.mean = mean(tail_vel_forward/spine_length.mean,na.rm=T),
      midpoint_speed_bl.mean = mean(midpoint_speed/spine_length.mean,na.rm=T),
      midpoint_distance_bl.mean = mean(midpoint_distance/spine_length.mean,na.rm=T)
    )%>%
    group_by(group_condition,bin,midpoint)%>%
    summarise(
      N= n(),
      bearing_angle.lowCI = ci(bearing_angle.mean, na.rm = T)[2],
      bearing_angle.hiCI = ci(bearing_angle.mean, na.rm = T)[3],
      bearing_angle.mean = mean(bearing_angle.mean,na.rm=T),
      
      abs_bearing_angle.lowCI = ci(abs_bearing_angle.mean, na.rm = T)[2],
      abs_bearing_angle.hiCI = ci(abs_bearing_angle.mean, na.rm = T)[3],
      abs_bearing_angle.mean = mean(abs_bearing_angle.mean,na.rm=T),
      
      heading_angle.lowCI = ci(heading_angle.mean, na.rm = T)[2],
      heading_angle.hiCI = ci(heading_angle.mean, na.rm = T)[3],
      heading_angle.mean = mean(heading_angle.mean,na.rm=T),
      
      
      abs_heading_angle.lowCI = ci(abs_heading_angle.mean, na.rm = T)[2],
      abs_heading_angle.hiCI = ci(abs_heading_angle.mean, na.rm = T)[3],
      abs_heading_angle.mean = mean(abs_heading_angle.mean,na.rm=T),
      
      bending_angle.lowCI = ci(bending_angle.mean, na.rm = T)[2],
      bending_angle.hiCI = ci(bending_angle.mean, na.rm = T)[3],
      bending_angle.mean = mean(bending_angle.mean,na.rm=T),
      
      
      abs_bending_angle.lowCI = ci(abs_bending_angle.mean, na.rm = T)[2],
      abs_bending_angle.hiCI = ci(abs_bending_angle.mean, na.rm = T)[3],
      abs_bending_angle.mean = mean(abs_bending_angle.mean,na.rm=T),
      
      head_vector_angular_speed.lowCI = ci(head_vector_angular_speed.mean, na.rm = T)[2],
      head_vector_angular_speed.hiCI = ci(head_vector_angular_speed.mean, na.rm = T)[3],
      head_vector_angular_speed.mean = mean(head_vector_angular_speed.mean,na.rm=T),
      
      abs_head_vector_angular_speed.lowCI = ci(abs_head_vector_angular_speed.mean, na.rm = T)[2],
      abs_head_vector_angular_speed.hiCI = ci(abs_head_vector_angular_speed.mean, na.rm = T)[3],
      abs_head_vector_angular_speed.mean = mean(abs_head_vector_angular_speed.mean,na.rm=T),
      
      tail_vector_angular_speed.lowCI = ci(tail_vector_angular_speed.mean, na.rm = T)[2],
      tail_vector_angular_speed.hiCI = ci(tail_vector_angular_speed.mean, na.rm = T)[3],
      tail_vector_angular_speed.mean = mean(tail_vector_angular_speed.mean,na.rm=T),
      
      abs_tail_vector_angular_speed.lowCI = ci(abs_tail_vector_angular_speed.mean, na.rm = T)[2],
      abs_tail_vector_angular_speed.hiCI = ci(abs_tail_vector_angular_speed.mean, na.rm = T)[3],
      abs_tail_vector_angular_speed.mean = mean(abs_tail_vector_angular_speed.mean,na.rm=T),
      
      
      head_vel_forward.lowCI = ci(head_vel_forward.mean, na.rm = T)[2],
      head_vel_forward.hiCI = ci(head_vel_forward.mean, na.rm = T)[3],
      head_vel_forward.mean = mean(head_vel_forward.mean,na.rm=T),
      
      tail_vel_forward.lowCI = ci(tail_vel_forward.mean, na.rm = T)[2],
      tail_vel_forward.hiCI = ci(tail_vel_forward.mean, na.rm = T)[3],
      tail_vel_forward.mean = mean(tail_vel_forward.mean,na.rm=T),
      
      midpoint_speed.lowCI = ci(midpoint_speed.mean, na.rm = T)[2],
      midpoint_speed.hiCI = ci(midpoint_speed.mean, na.rm = T)[3],
      midpoint_speed.mean = mean(midpoint_speed.mean,na.rm=T),
      
      midpoint_distance.lowCI = ci(midpoint_distance.mean, na.rm = T)[2],
      midpoint_distance.hiCI = ci(midpoint_distance.mean, na.rm = T)[3],
      midpoint_distance.mean = mean(midpoint_distance.mean,na.rm=T),
      
      distance_to_odor.lowCI = ci(distance_to_odor.mean, na.rm = T)[2],
      distance_to_odor.hiCI = ci(distance_to_odor.mean, na.rm = T)[3],
      distance_to_odor.mean = mean(distance_to_odor.mean,na.rm=T),
      
      spine_length.lowCI = ci(spine_length.mean,na.rm=T)[2],
      spine_length.hiCI = ci(spine_length.mean,na.rm=T)[3],
      spine_length.mean = mean(spine_length.mean,na.rm=T),
      
      head_vel_forward_bl.lowCI = ci(head_vel_forward_bl.mean, na.rm = T)[2],
      head_vel_forward_bl.hiCI = ci(head_vel_forward_bl.mean, na.rm = T)[3],
      head_vel_forward_bl.mean = mean(head_vel_forward_bl.mean,na.rm=T),
      
      tail_vel_forward_bl.lowCI = ci(tail_vel_forward_bl.mean, na.rm = T)[2],
      tail_vel_forward_bl.hiCI = ci(tail_vel_forward_bl.mean, na.rm = T)[3],
      tail_vel_forward_bl.mean = mean(tail_vel_forward_bl.mean,na.rm=T),
      
      midpoint_speed_bl.lowCI = ci(midpoint_speed_bl.mean, na.rm = T)[2],
      midpoint_speed_bl.hiCI = ci(midpoint_speed_bl.mean, na.rm = T)[3],
      midpoint_speed_bl.mean = mean(midpoint_speed_bl.mean,na.rm=T),
      
      midpoint_distance_bl.lowCI = ci(midpoint_distance_bl.mean, na.rm = T)[2],
      midpoint_distance_bl.hiCI = ci(midpoint_distance_bl.mean, na.rm = T)[3],
      midpoint_distance_bl.mean = mean(midpoint_distance_bl.mean,na.rm=T)
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
 
binning_data_hc <-
  function(df,
           variable,
           width = 10,
           Abs_HC_Angle_interval = c(20, 360),
           frame_rate=16
  ) {
    
    message("binning HC data")
    #get breaks and midpoints
    var <- df[[variable]]
    var_min <- min(var,na.rm=T)
    var_max <- max(var,na.rm=T)
    breaks <- seq(var_min,var_max,width)
    midpoints <- (breaks[2:length(breaks)]+breaks[1:length(breaks)-1])/2
    
    #filter data 
    filtered_data <- df %>%
      filter(
        is.na(Abs_HC_angle)  |
          Abs_HC_angle >= Abs_HC_Angle_interval[1] &
          Abs_HC_angle <= Abs_HC_Angle_interval[2]
      )
    
    #do binning
    bin <- cut(filtered_data[[variable]],breaks,labels=F)
    bindf1 <- data.frame(list(bin=bin))
    if(length(midpoints)!=length(sort(unique(bin)))){
      stop("Binning failed.Try a higher width.")
    }
    bindf2 <- data.frame(list(bin = sort(unique(bin)),midpoints = midpoints))
    bindf <- left_join(bindf1,bindf2,by="bin")
    binned_data <-  filtered_data%>%
      mutate(bin = bindf$bin,
             midpoint = bindf$midpoints
      )
    
    #aggregation
    binned_data %>%
      group_by(group_condition,bin,midpoint) %>%
      summarise(
        N=n(),
        HC_rate = (sum(HCs) / n())*frame_rate,
        
        HC_angle.lowCI = ci(HC_angle, na.rm = T)[2],
        HC_angle.hiCI = ci(HC_angle, na.rm = T)[3],
        HC_angle.mean = mean(HC_angle, na.rm = T),
        
        Abs_HC_angle.lowCI = ci(Abs_HC_angle, na.rm = T)[2],
        Abs_HC_angle.hiCI = ci(Abs_HC_angle, na.rm = T)[3],
        Abs_HC_angle.mean = mean(Abs_HC_angle, na.rm = T),
        
        HC_reorientation.lowCI = ci(HC_reorientation, na.rm = T)[2],
        HC_reorientation.hiCI = ci(HC_reorientation, na.rm = T)[3],
        HC_reorientation.mean = mean(HC_reorientation, na.rm = T),
        
        HC_rate_modulation=
          (sum(HCs[odor_orientation=="towards"],na.rm=T)-
             sum(HCs[odor_orientation=="away"],na.rm=T))/
          (sum(HCs[odor_orientation=="towards"],na.rm=T)+
             sum(HCs[odor_orientation=="away"],na.rm=T))
        
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
 
binning_data_hc_grouped <-
  function(df,
           variable,
           width = 10,
           grouping="id",
           Abs_HC_Angle_interval = c(20, 360),
           frame_rate=16
  ) {
    
    message("binning HC data")
    
    #get breaks and midpoints
    var <- df[[variable]]
    var_min <- min(var,na.rm=T)
    var_max <- max(var,na.rm=T)
    breaks <- seq(var_min,var_max,width)
    midpoints <- (breaks[2:length(breaks)]+breaks[1:length(breaks)-1])/2
    
    #filter data 
    filtered_data <- df %>%
      filter(
        is.na(Abs_HC_angle)  |
          Abs_HC_angle >= Abs_HC_Angle_interval[1] &
          Abs_HC_angle <= Abs_HC_Angle_interval[2]
      )
    
    #do binning
    bin <- cut(filtered_data[[variable]],breaks,labels=F)
    bindf1 <- data.frame(list(bin=bin))
    if(length(midpoints)!=length(sort(unique(bin)))){
      stop("Binning failed.Try a higher width.")
    }
    bindf2 <- data.frame(list(bin = sort(unique(bin)),midpoints = midpoints))
    bindf <- left_join(bindf1,bindf2,by="bin")
    binned_data <-  filtered_data%>%
      mutate(bin = bindf$bin,
             midpoint = bindf$midpoints
      )
    
    #aggregation
    binned_data %>%
      group_by_("group_condition",grouping,"bin","midpoint") %>%
      summarise(
        N=n(),
        HC_rate = (sum(HCs) / n())*frame_rate,
        HC_angle.mean = mean(HC_angle,na.rm=T),
        Abs_HC_angle.mean = mean(Abs_HC_angle,na.rm=T),
        HC_reorientation.mean = mean(HC_reorientation,na.rm=T),
        
        HC_rate_modulation=
          (sum(HCs[odor_orientation=="towards"],na.rm=T)-
             sum(HCs[odor_orientation=="away"],na.rm=T))/
          (sum(HCs[odor_orientation=="towards"],na.rm=T)+
             sum(HCs[odor_orientation=="away"],na.rm=T))
        , by = c(grouping,"bin"))%>%
      group_by(group_condition,bin,midpoint)%>%
      summarise(
        HC_rate.lowCI = ci(HC_rate, na.rm = T)[2],
        HC_rate.hiCI = ci(HC_rate, na.rm = T)[3],
        HC_rate.mean = mean(HC_rate,na.rm=T),
        
        HC_angle.lowCI = ci(HC_angle.mean, na.rm = T)[2],
        HC_angle.hiCI = ci(HC_angle.mean, na.rm = T)[3],
        HC_angle.mean = mean(HC_angle.mean,na.rm=T),
        
        Abs_HC_angle.lowCI = ci(Abs_HC_angle.mean, na.rm = T)[2],
        Abs_HC_angle.hiCI = ci(Abs_HC_angle.mean, na.rm = T)[3],
        Abs_HC_angle.mean = mean(Abs_HC_angle.mean,na.rm=T),
        
        HC_reorientation.lowCI = ci(HC_reorientation.mean, na.rm = T)[2],
        HC_reorientation.hiCI = ci(HC_reorientation.mean, na.rm = T)[3],
        HC_reorientation.mean = mean(HC_reorientation.mean,na.rm=T),
        
        HC_rate_modulation.lowCI = ci(HC_rate_modulation, na.rm = T)[2],
        HC_rate_modulation.hiCI = ci(HC_rate_modulation, na.rm = T)[3],
        HC_rate_modulation.mean = mean(HC_rate_modulation,na.rm=T),
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
                            width = 10) {
  
  message("binning IS data")
  
  #get breaks and midpoints
  var <- df[[variable]]
  var_min <- min(var,na.rm=T)
  var_max <- max(var,na.rm=T)
  breaks <- seq(var_min,var_max,width)
  midpoints <- (breaks[2:length(breaks)]+breaks[1:length(breaks)-1])/2
  
  #do binning
  filtered_data <-  df #%>%filter(step_boolean == 1) 
  bin <- cut(filtered_data[[variable]],breaks,labels=F,right=T)
  bindf1 <- data.frame(list(bin=bin))
  if(length(midpoints)!=length(sort(unique(bin)))){
    stop("Binning failed.Try a higher width.")
  }
  bindf2 <- data.frame(list(bin = sort(unique(bin)),midpoints = midpoints))
  bindf <- left_join(bindf1,bindf2,by="bin")
  
  binned_data <-  filtered_data%>%
    mutate(bin = bindf$bin,
           midpoint = bindf$midpoints
    )
  #aggregation
  binned_data%>%
    group_by(group_condition,bin,midpoint) %>%
    summarise(
      IS_angle.lowCI = ci(IS_angle, na.rm = T)[2],
      IS_angle.hiCI = ci(IS_angle, na.rm = T)[3],
      IS_angle.mean = mean(IS_angle,na.rm=T),
      
      Abs_IS_angle.lowCI = ci(Abs_IS_angle, na.rm = T)[2],
      Abs_IS_angle.hiCI = ci(Abs_IS_angle, na.rm = T)[3],
      Abs_IS_angle.mean = mean(Abs_IS_angle,na.rm=T),
      
      
      IS_reorientation.lowCI = ci(IS_reorientation, na.rm = T)[2],
      IS_reorientation.hiCI = ci(IS_reorientation, na.rm = T)[3],
      IS_reorientation.mean = mean(IS_reorientation,na.rm=T),
      
      IS_distance.lowCI = ci(IS_distance, na.rm = T)[2],
      IS_distance.hiCI = ci(IS_distance, na.rm = T)[3],
      IS_distance.mean = mean(IS_distance,na.rm=T),
      
      IS_distance_bl.lowCI = ci(IS_distance_bl, na.rm = T)[2],
      IS_distance_bl.hiCI = ci(IS_distance_bl, na.rm = T)[3],
      IS_distance_bl.mean = mean(IS_distance_bl,na.rm=T),
      
      IS_interval.lowCI = ci(IS_interval, na.rm = T)[2],
      IS_interval.hiCI = ci(IS_interval, na.rm = T)[3],
      IS_interval.mean = mean(IS_interval,na.rm=T)
      
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
                                    grouping="id"
) {
  
  message("binning IS data")
  
  #get breaks and midpoints
  var <- df[[variable]]
  var_min <- min(var,na.rm=T)
  var_max <- max(var,na.rm=T)
  breaks <- seq(var_min,var_max,width)
  midpoints <- (breaks[2:length(breaks)]+breaks[1:length(breaks)-1])/2
  
  #do binning
  filtered_data <-  df #%>%filter(step_boolean == 1) 
  bin <- cut(filtered_data[[variable]],breaks,labels=F)
  bindf1 <- data.frame(list(bin=bin))
  if(length(midpoints)!=length(sort(unique(bin)))){
    stop("Binning failed.Try a higher width.")
  }
  bindf2 <- data.frame(list(bin = sort(unique(bin)),midpoints = midpoints))
  bindf <- left_join(bindf1,bindf2,by="bin")
  binned_data <-  filtered_data%>%
    mutate(bin = bindf$bin,
           midpoint = bindf$midpoints
    )
  
  binned_data %>%
    group_by_("group_condition",grouping,"bin","midpoint")%>%
    summarise(
      IS_angle.mean = mean(IS_angle,na.rm=T),
      
      Abs_IS_angle.mean = mean(Abs_IS_angle,na.rm=T),
      IS_reorientation.mean = mean(IS_reorientation,na.rm=T),
      IS_distance.mean = mean(IS_distance,na.rm=T),
      IS_distance_bl.mean = mean(IS_distance_bl,na.rm=T),
      IS_interval.mean = mean(IS_interval,na.rm=T)
    )%>%
    group_by(group_condition,bin,midpoint)%>%
    summarise(
      IS_angle.lowCI = ci(IS_angle.mean, na.rm = T)[2],
      IS_angle.hiCI = ci(IS_angle.mean, na.rm = T)[3],
      IS_angle.mean = mean(IS_angle.mean,na.rm=T),
      
      Abs_IS_angle.lowCI = ci(Abs_IS_angle.mean, na.rm = T)[2],
      Abs_IS_angle.hiCI = ci(Abs_IS_angle.mean, na.rm = T)[3],
      Abs_IS_angle.mean = mean(Abs_IS_angle.mean,na.rm=T),
      
      IS_reorientation.lowCI = ci(IS_reorientation.mean, na.rm = T)[2],
      IS_reorientation.hiCI = ci(IS_reorientation.mean, na.rm = T)[3],
      IS_reorientation.mean = mean(IS_reorientation.mean,na.rm=T),
      
      IS_distance.lowCI = ci(IS_distance.mean, na.rm = T)[2],
      IS_distance.hiCI = ci(IS_distance.mean, na.rm = T)[3],
      IS_distance.mean = mean(IS_distance.mean,na.rm=T),
      
      IS_distance_bl.lowCI = ci(IS_distance_bl.mean, na.rm = T)[2],
      IS_distance_bl.hiCI = ci(IS_distance_bl.mean, na.rm = T)[3],
      IS_distance_bl.mean = mean(IS_distance_bl.mean,na.rm=T),
      
      IS_interval.lowCI = ci(IS_interval.mean, na.rm = T)[2],
      IS_interval.hiCI = ci(IS_interval.mean, na.rm = T)[3],
      IS_interval.mean = mean(IS_interval.mean,na.rm=T)
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
 
binning_data_run_speed <-
  function(df,
           variable,
           width = 10) {
    message("binning run_speed data")
    
    #get breaks and midpoints
    var <- df[[variable]]
    var_min <- min(var,na.rm=T)
    var_max <- max(var,na.rm=T)
    breaks <- seq(var_min,var_max,width)
    midpoints <- (breaks[2:length(breaks)]+breaks[1:length(breaks)-1])/2
    
    #do binning
    filtered_data <-  df #%>%filter(run_flag == T)
    bin <- cut(filtered_data[[variable]],breaks,labels=F)
    bindf1 <- data.frame(list(bin=bin))
    if(length(midpoints)!=length(sort(unique(bin)))){
      stop("Binning failed.Try a higher width.")
    }
    bindf2 <- data.frame(list(bin = sort(unique(bin)),midpoints = midpoints))
    bindf <- left_join(bindf1,bindf2,by="bin")
    binned_data <-  filtered_data%>%
      mutate(bin = bindf$bin,
             midpoint = bindf$midpoints
      )
    
    # aggregation
    binned_data %>%
      group_by(group_condition,bin,midpoint) %>%
      summarise(
        run_speed.lowCI = ci(run_speed, na.rm = T)[2],
        run_speed.hiCI = ci(run_speed, na.rm = T)[3],
        run_speed.mean = mean(run_speed,na.rm=T),
        
        run_speed_bl.lowCI = ci(run_speed_bl, na.rm = T)[2],
        run_speed_bl.hiCI = ci(run_speed_bl, na.rm = T)[3],
        run_speed_bl.mean = mean(run_speed_bl,na.rm=T),
        
        runs_peed_modulation = 
          (mean(run_speed[odor_orientation=="towards"],na.rm=T)-
             mean(run_speed[odor_orientation=="away"],na.rm=T))/(
               mean(run_speed[odor_orientation=="towards"],na.rm=T)+
                 mean(run_speed[odor_orientation=="away"],na.rm=T)
               
             )
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
binning_data_run_speed_grouped <-
  function(df,
           variable,
           width = 10,
           grouping="id"
  ) {
    message("binning run_speed data")
    #get breaks and midpoints
    var <- df[[variable]]
    var_min <- min(var,na.rm=T)
    var_max <- max(var,na.rm=T)
    breaks <- seq(var_min,var_max,width)
    midpoints <- (breaks[2:length(breaks)]+breaks[1:length(breaks)-1])/2
    
    #do binning
    filtered_data <-  df #%>%filter(run_flag == T)
    bin <- cut(filtered_data[[variable]],breaks,labels=F)
    bindf1 <- data.frame(list(bin=bin))
    if(length(midpoints)!=length(sort(unique(bin)))){
      stop("Binning failed.Try a higher width.")
    }
    bindf2 <- data.frame(list(bin = sort(unique(bin)),midpoints = midpoints))
    bindf <- left_join(bindf1,bindf2,by="bin")
    binned_data <-  filtered_data%>%
      mutate(bin = bindf$bin,
             midpoint = bindf$midpoints
      )
    
    # aggregation
    binned_data %>%
      group_by_("group_condition",grouping,"bin","midpoint") %>%
      summarise(
        run_speed = mean(run_speed,na.rm=T),
        run_speed_bl = mean(run_speed_bl,na.rm=T),
        run_speed_modulation =
          (mean(run_speed[odor_orientation=="towards"],na.rm=T)-
             mean(run_speed[odor_orientation=="away"],na.rm=T))/
          (mean(run_speed[odor_orientation=="towards"],na.rm=T)+
             mean(run_speed[odor_orientation=="away"],na.rm=T))
      )%>%
      group_by(group_condition,bin,midpoint)%>%
      summarise(
        run_speed.lowCI = ci(run_speed, na.rm = T)[2],
        run_speed.hiCI = ci(run_speed, na.rm = T)[3],
        run_speed.mean = mean(run_speed,na.rm=T),
        
        run_speed_bl.lowCI = ci(run_speed_bl, na.rm = T)[2],
        run_speed_bl.hiCI = ci(run_speed_bl, na.rm = T)[3],
        run_speed_bl.mean = mean(run_speed_bl,na.rm=T),
        
        run_speed_modulation.lowCI = ci(run_speed_modulation, na.rm = T)[2],
        run_speed_modulation.hiCI = ci(run_speed_modulation, na.rm = T)[3],
        run_speed_modulation.mean = mean(run_speed_modulation,na.rm=T)
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
binning_data_preference <-
  function(df, variable, radius, width = 10) {
    
    message("binning preference data")
    #get breaks and midpoints
    var <- df[[variable]]
    var_min <- min(var,na.rm=T)
    var_max <- max(var,na.rm=T)
    breaks <- seq(var_min,var_max,width)
    midpoints <- (breaks[2:length(breaks)]+breaks[1:length(breaks)-1])/2
    
    #do binning
    filtered_data <-  df
    bin <- cut(filtered_data[[variable]],breaks,labels=F)
    bindf1 <- data.frame(list(bin=bin))
    if(length(midpoints)!=length(sort(unique(bin)))){
      stop("Binning failed.Try a higher width.")
    }
    bindf2 <- data.frame(list(bin = sort(unique(bin)),midpoints = midpoints))
    bindf <- left_join(bindf1,bindf2,by="bin")
    binned_data <-  filtered_data%>%
      mutate(bin = bindf$bin,
             midpoint = bindf$midpoints
      )
    
    # aggregation
    binned_data %>%
      group_by(group_condition,bin,midpoint)%>%
      summarise(
        
        preference = mean(pref,na.rm=T),
        preference_dist.lowCI = -1 * rescale(
          ci(distance_to_odor, na.rm = T)[2],
          to = c(-1, 1),
          from = c(0, radius * 2)
        ),
        preference_dist.hiCI = -1 * rescale(
          ci(distance_to_odor, na.rm = T)[3],
          to = c(-1, 1),
          from = c(0, radius * 2)
        ),
        preference_dist.mean = -1 * rescale(
          mean(distance_to_odor, na.rm = TRUE),
          to = c(-1, 1),
          from = c(0, radius * 2)
        ),
        odor_speed.lowCI = ci(odor_speed, na.rm = T)[2],
        odor_speed.hiCI = ci(odor_speed, na.rm = T)[3],
        
        odor_speed.mean = mean(odor_speed, na.rm = T),
        ratio_towards = sum(odor_orientation=="towards")/(
          sum(odor_orientation == "towards")+
            sum(odor_orientation=="away"))
        
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
binning_data_preference_grouped <-
  function(df, variable, radius, width = 10, grouping="id") {
    
    message("binning preference data")
    #get breaks and midpoints
    var <- df[[variable]]
    var_min <- min(var,na.rm=T)
    var_max <- max(var,na.rm=T)
    breaks <- seq(var_min,var_max,width)
    midpoints <- (breaks[2:length(breaks)]+breaks[1:length(breaks)-1])/2
    
    #do binning
    filtered_data <-  df
    bin <- cut(filtered_data[[variable]],breaks,labels=F)
    bindf1 <- data.frame(list(bin=bin))
    if(length(midpoints)!=length(sort(unique(bin)))){
      stop("Binning failed.Try a higher width.")
    }
    bindf2 <- data.frame(list(bin = sort(unique(bin)),midpoints = midpoints))
    bindf <- left_join(bindf1,bindf2,by="bin")
    binned_data <-  filtered_data%>%
      mutate(bin = bindf$bin,
             midpoint = bindf$midpoints
      )
    
    # aggregation
    binned_data %>%
      group_by_("group_condition",grouping,"bin","midpoint") %>%
      summarise(
        preference = mean(pref,na.rm=T),
        
        preference_dist = -1 * rescale(
          mean(distance_to_odor, na.rm = TRUE),
          to = c(-1, 1),
          from = c(0, radius * 2)
        ),
        odor_speed = mean(odor_speed, na.rm = T),
        ratio_towards = sum(odor_orientation=="towards",na.rm=T)/(
          sum(odor_orientation == "towards",na.rm=T)+
            sum(odor_orientation=="away",na.rm=T))
        
      )%>%
      group_by(group_condition,bin,midpoint)%>%
      summarise(
        preference.mean = mean(preference,na.rm=T),
        preference.lowCI = ci(preference,na.rm=T)[2],
        preference.hiCI = ci(preference,na.rm=T)[3],
        
        preference_dist.mean = mean(preference_dist,na.rm=T),
        preference_dist.lowCI = ci(preference_dist,na.rm=T)[2],
        preference_dist.hiCI = ci(preference_dist,na.rm=T)[3],
        
        odor_speed.mean = mean(odor_speed,na.rm=T),
        odor_speed.lowCI = ci(odor_speed,na.rm=T)[2],
        odor_speed.hiCI = ci(odor_speed,na.rm=T)[3],
        
        ratio_towards.mean = mean(ratio_towards,na.rm=T),
        ratio_towards.lowCI = ci(ratio_towards,na.rm=T)[2],
        ratio_towards.hiCI = ci(ratio_towards,na.rm=T)[3]
        
        
      )
  }

