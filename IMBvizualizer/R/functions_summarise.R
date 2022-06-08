##### SUMMARY: data.table #####

#' @param data 
#'
#' @param grouped_by 
#' @param frame_interval 
#' @param distance_to_odor_interval 
#' @param Abs_HC_Angle_interval 
#' @param Abs_bearing_angle_interval 
#' @param radius 
#' @param video_length 
#' @param frame_rate 
#' @param direction 
#' @param threshold 
#' @import dplyr
#' @import data.table
#' @export
create_summarized_analysis_dt <- function(data,
                                          grouped_by="trial",
                                          frame_interval=c(0,180),
                                          distance_to_odor_interval=c(0,2000),
                                          Abs_HC_Angle_interval = c(20,360),
                                          Abs_bearing_angle_interval =c(0,180),
                                          radius=42.5,
                                          video_length=180,
                                          frame_rate=16,
                                          direction = NULL,
                                          threshold=90){
  gc()
  data$group_condition <- paste(data$group,data$condition,sep="-")
  
  if (grouped_by == "trial") {
    gb =  c("trial", "group", "condition","group_condition")
  }
  
  if (grouped_by == "id") {
    gb =  c("id", "group", "condition","group_condition")
    idx <- data[,.(
      valid =.N>frame_rate*threshold
    ),id]
    idx <- idx[valid==TRUE]
    data <- data[id %in% idx$id ]
  }
  message("filter data...")

  if(!is.null(direction)){
    #direction filter
    if(direction=="forwards"){
      data <- as.data.frame(data)%>%
        filter(run_direction==1)
    }
    if(direction=="backwards"){
      data <- as.data.frame(data)%>%
        filter(run_direction==-1)
    }
    if(direction=="both"){
      data <- as.data.frame(data)%>%
        filter(run_direction==1|run_direction==-1)
    }
  }

  filtered_data <- as.data.frame(data) %>%
    dplyr::filter(na_flag_flip == F) %>%
    dplyr::filter(
      frame %in% seq(from = frame_interval[1] * frame_rate,
                     to = frame_interval[2] * frame_rate) &
        distance_to_odor >= distance_to_odor_interval[1] &
        distance_to_odor <= distance_to_odor_interval[2] &
        abs(bearing_angle) >= Abs_bearing_angle_interval[1] &
        abs(bearing_angle) <= Abs_bearing_angle_interval[2]
    )%>%as.data.table()
  
  if (length(distance_to_odor_interval) < 2) {
    distance_to_odor_interval = c(0, as.numeric(radius) * 2)
  }
  
 
  test = as.data.table(data.frame(list(id=unique(data$id))))

  DT_list <- list(
    
    summarise_TS_per_group_dt(filtered_data, gb),
    summarise_HC_per_group_dt(filtered_data, gb,Abs_HC_Angle_interval,frame_rate),
    summarise_RUN_per_group_dt(filtered_data, gb),
    summarise_PREF_per_group_dt(filtered_data, gb, radius),
    summarise_TRACK_per_group_dt(filtered_data,gb,frame_rate),
    aggregate_groups(data,gb)
  )
  result = Reduce(function(...) merge(...,by=gb, all = TRUE),DT_list)
  df <- as.data.frame(result)
  df <- df %>%
    mutate(
      HC_rate =coalesce(HC_rate,0),
      HCs = coalesce(HCs,0)
    )%>%
    select(!test)
  as.data.table(df)
  
}

#' Title
#'
#' @param data 
#' @param grouped_by 
#'
#' @return
#' @export
#' @import data.table
#' @examples
aggregate_groups<- function(data,grouped_by){
  message("summarising TRACK-data ...")
  #data <- as.disk.frame(data)
  data[,.(
    test = sum(spine_length,na.rm=T)
  )
  ,by=grouped_by]
  
}

#' @param data 
#'
#' @param grouped_by 
#' @import dplyr
#' @import data.table
#' @export
summarise_TS_per_group_dt <- function(data,
                                      grouped_by){
  message("summarising TS-data ...")
  #data <- as.disk.frame(data)
  data[,.(
    bearing_angle_mean = mean(bearing_angle,na.rm=TRUE),
    bearing_angle_var = var(bearing_angle,na.rm=TRUE),
    heading_angle_mean = mean(heading_angle,na.rm=TRUE),
    heading_angle_var = var(heading_angle,na.rm=TRUE),
    bending_angle_mean = mean(bending_angle,na.rm=TRUE),
    bending_angle_var = var(bending_angle,na.rm=TRUE),
    
    abs_bending_angle_mean = mean(abs(bending_angle),na.rm=T),
    abs_bending_angle_var = var(abs(bending_angle),na.rm=T),
    abs_bearing_angle_mean = mean(abs(bearing_angle),na.rm=T),
    abs_bearing_angle_var = var(abs(bearing_angle),na.rm=T),
    abs_heading_angle_mean = mean(abs(heading_angle),na.rm=T),
    abs_heading_angle_var = var(abs(heading_angle),na.rm=T),
    
    head_vector_angular_speed_mean = mean(head_vector_angular_speed,na.rm=TRUE),
    head_vector_angular_speed_var = var(head_vector_angular_speed,na.rm=TRUE),
    abs_head_vector_angular_speed_mean = mean(abs(head_vector_angular_speed),na.rm=TRUE),
    abs_head_vector_angular_speed_var = var(abs(head_vector_angular_speed),na.rm=TRUE),
    head_vector_angular_acc_mean = mean(diff(head_vector_angular_speed),na.rm=TRUE),
    head_vector_angular_acc_var = var(diff(head_vector_angular_speed),na.rm=TRUE),
    
    
    head_vel_forward_mean = mean(head_vel_forward,na.rm=TRUE),
    head_vel_forward_var = var(head_vel_forward,na.rm=TRUE),
    head_vel_forward_bl_mean = mean(head_vel_forward/mean(spine_length,na.rm=T), na.rm = TRUE),
    head_vel_forward_bl_var = var(head_vel_forward/mean(spine_length,na.rm=T), na.rm = TRUE),
    head_acc_forward_mean = mean(diff(head_vel_forward),na.rm=TRUE),
    head_acc_forward_var = var(diff(head_vel_forward),na.rm=TRUE),
    
    tail_vector_angular_speed_mean = mean(tail_vector_angular_speed,na.rm=TRUE),
    tail_vector_angular_speed_var = var(tail_vector_angular_speed,na.rm=TRUE),
    abs_tail_vector_angular_speed_mean = mean(abs(tail_vector_angular_speed),na.rm=TRUE),
    abs_tail_vector_angular_speed_var = var(abs(tail_vector_angular_speed),na.rm=TRUE),
    tail_vector_angular_acc_mean = mean(diff(tail_vector_angular_speed),na.rm=TRUE),
    tail_vector_angular_acc_var = var(diff(tail_vector_angular_speed),na.rm=TRUE),
    
    tail_vel_forward_mean = mean(tail_vel_forward,na.rm=TRUE),
    tail_vel_forward_var = var(tail_vel_forward,na.rm=TRUE),
    tail_vel_forward_bl_mean = mean(tail_vel_forward/mean(spine_length,na.rm=T), na.rm = TRUE),
    tail_vel_forward_bl_var = var(tail_vel_forward/mean(spine_length,na.rm=T), na.rm = TRUE),
    tail_acc_forward_mean = mean(diff(tail_vel_forward),na.rm=TRUE),
    tail_acc_forward_var = var(diff(tail_vel_forward),na.rm=TRUE),
    
    midpoint_speed_mean = mean(midpoint_speed,na.rm=TRUE),
    midpoint_speed_var = var(midpoint_speed,na.rm=TRUE),
    larvae_speed_mean =  mean(midpoint_speed/mean(spine_length,na.rm=TRUE),na.rm=TRUE),
    larvae_speed_var =  var(midpoint_speed/mean(spine_length,na.rm=TRUE),na.rm=TRUE),
    
    
    midpoint_acc_mean = mean(diff(midpoint_speed),na.rm=TRUE),
    midpoint_acc_var = var(diff(midpoint_speed),na.rm=TRUE),
    midpoint_distance_mean = mean(max(midpoint_distance),na.rm=TRUE),
    larvae_distance_mean =  mean(max(midpoint_distance)/mean(spine_length,na.rm=TRUE),na.rm=TRUE),
    
    distance_to_odor_mean = mean(distance_to_odor,na.rm=TRUE),
    distance_to_sp_mean = mean(distance_to_sp,na.rm=TRUE),
    distance_to_sp_max =max(distance_to_sp,na.rm=TRUE)
    
    
  )
  ,by=grouped_by]
  
}

#' @param data 
#'
#' @param grouped_by 
#' @param frame_rate 
#' @import data.table
#' @export
summarise_TRACK_per_group_dt <- function(data,
                                         grouped_by,
                                         frame_rate
                                         ){
  message("summarising TRACK-data ...")
  #data <- as.disk.frame(data)
  data[,.(
    area_mean = mean(area,na.rm=TRUE),
    area_var = var(area,na.rm=TRUE),
    grey_mean = mean(grey,na.rm=TRUE),
    grey_var = var(grey,na.rm=TRUE),
    length_mean = mean(spine_length,na.rm=TRUE),
    length_var = var(spine_length,na.rm=TRUE),
    width_mean = mean(width,na.rm=TRUE),
    width_var = var(width,na.rm=TRUE),
    perimeter_mean = mean(perimeter,na.rm=TRUE),
    perimeter_var = var(perimeter,na.rm=TRUE),
    ncollisions = sum(as.logical(collision_flag),na.rm=TRUE),
    neg_flank = frame[tail(which(flanks==-1),1)]/frame_rate,
    pos_flank = frame[tail(which(flanks==1),1)]/frame_rate
  )
  ,by=grouped_by]
  
}


#' @param data 
#'
#' @param grouped_by 
#' @param Abs_HC_Angle_interval 
#' @param frame_rate 
#' @import dplyr
#' @import data.table
#' @export
summarise_HC_per_group_dt <- function(data,
                                      grouped_by,
                                      Abs_HC_Angle_interval=c(20,360),
                                      frame_rate = 16
                                      ) {
  
  message("summarising HC-data ...")
  
  df1 <- data[,
              .(
                frames_towards = sum(odor_orientation=="towards",na.rm=T),
                frames_away = sum(odor_orientation=="away",na.rm=T),
                frames = .N
              ),grouped_by]
  
  df2 <- data[
    Abs_HC_angle >= Abs_HC_Angle_interval[1] &
      Abs_HC_angle <= Abs_HC_Angle_interval[2],
    .(
        HC_angle_mean = mean(HC_angle,na.rm=TRUE),
        HC_angle_var = var(HC_angle,na.rm=TRUE),
        
        Abs_HC_angle_mean = mean(Abs_HC_angle,na.rm=TRUE),
        Abs_HC_angle_var = var(Abs_HC_angle,na.rm=TRUE),
        
        HC_reorientation_mean = mean(HC_reorientation,na.rm=TRUE),
        HC_reorientation_var = var(HC_reorientation,na.rm=TRUE),
        HCs_towards = sum(HCs[odor_orientation=="towards"],na.rm=TRUE),
        HCs_away= sum(HCs[odor_orientation=="away"],na.rm=TRUE),
        
        HCs = sum(HCs,na.rm=T)
      ),grouped_by]

  df <- merge(df1,df2,all=T)
  df$HCs[is.na(df$HCs)] <- 0
  df$HCs_towards[is.na(df$HCs_towards)] <- 0
  df$HCs_away[is.na(df$HCs_away)] <- 0
  
  
  df <- df[,':='(
    HC_rate_modulation = ((HCs_away/frames_away)-(HCs_towards/frames_towards))/
      ((HCs_away/frames_away)+(HCs_towards/frames_towards)),
    HC_rate = (HCs/frames)*frame_rate
  )
  ]
  
  df
  #result
}

#' @param data 
#'
#' @param grouped_by 
#' @import dplyr
#' @import data.table
#' @export
summarise_RUN_per_group_dt <- function(data,
                                       grouped_by) {
  
  message("summarising RUN-data ...")
  data[step_boolean == TRUE ,.(
    run_speed_mean = mean(run_speed, na.rm = TRUE),
    run_speed_var = var(run_speed, na.rm = TRUE),
    
    run_speed_bl_mean = mean(run_speed/mean(spine_length,na.rm=T), na.rm = TRUE),
    run_speed_bl_var = var(run_speed/mean(spine_length,na.rm=T), na.rm = TRUE),
    

    run_speed_modulation = 
      (mean(run_speed[odor_orientation=="towards"],na.rm=TRUE)-
         mean(run_speed[odor_orientation=="away"],na.rm=TRUE))/
      (mean(run_speed[odor_orientation=="towards"],na.rm=TRUE)+
         mean(run_speed[odor_orientation=="away"],na.rm=TRUE)),
    
    
    IS_angle_mean = mean(IS_angle, na.rm = TRUE),
    IS_angle_var = var(IS_angle, na.rm = TRUE),
    
    Abs_IS_angle_mean = mean(Abs_IS_angle, na.rm = TRUE),
    Abs_IS_angle_var = var(Abs_IS_angle, na.rm = TRUE),
    
    IS_reorientation_mean = mean(IS_reorientation, na.rm = TRUE),
    IS_reorientation_var = var(IS_reorientation, na.rm = TRUE),
    
    IS_interval_mean = mean(IS_interval, na.rm = TRUE),
    IS_interval_var = var(IS_interval, na.rm = TRUE),
    
    IS_distance_mean = mean(IS_distance, na.rm = TRUE),
    IS_distance_var = var(IS_distance, na.rm =TRUE),
    
    IS_distance_bl_mean = mean(IS_distance_bl,na.rm=T),
    IS_distance_bl_var = var(IS_distance_bl,na.rm=T),
    
    step_extr_mean = mean(step_extr,na.rm=T),
    
    ratio_minima = 1-sum(gt_1min,na.rm=T)/.N
    
  ),by=grouped_by]
  
}


#' @param data 
#'
#' @param grouped_by 
#' @param radius 
#'
#' @export
summarise_PREF_per_group_dt <-function(data, grouped_by, radius) {
    
    message("summarising PREF-data ...")
    preference_dist_se <- function(df, radius, grouped_by = "id") {
      df <- na.omit(df, cols="distance_to_odor")
      data <- df[,.(
        preference_dist_se =
          distance_to_odor[1] - tail(distance_to_odor, 1)
      ),by=c("id","trial")]
      
      data %>% select(!trial)
      if (grouped_by[1] == "trial") {
        data <- data[,.(
          preference_dist_se = mean(preference_dist_se,na.rm=TRUE)
        ),trial]
      }
      data
    }
    
    df1 <- data[,.(
      preference = (sum(spinepoint_y_6_conv>0,na.rm=TRUE)-sum(spinepoint_y_6_conv<0,na.rm=TRUE))/.N,
      pref_dist = -1 * rescale(
        mean(distance_to_odor, na.rm = TRUE),
        to = c(-1, 1),
        from = c(0, as.numeric(radius) * 2)),
      ratio_towards = sum(odor_orientation=="towards")/.N,
      odor_speed = mean(odor_speed,na.rm=T)
      
    ),by = grouped_by]
    
    df2=preference_dist_se(data, radius, grouped_by)
    merge(df1,df2,all=T)
  }
