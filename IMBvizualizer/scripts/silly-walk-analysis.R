library(data.table)
library(dplyr)
library(illav)

summarise_RUN_per_group_dt_ <- function(data,
                                        grouped_by) {
  
  message("summarising RUN-data ...")
  data[step_boolean == TRUE|step_boolean_bw ==TRUE ,.(
    run_speed_mean = mean(run_speed, na.rm = TRUE),
    run_speed_var = var(run_speed, na.rm = TRUE),
    
    run_speed_bw_mean = mean(run_speed_bw, na.rm = TRUE),
    run_speed_bw_var = var(run_speed_bw, na.rm = TRUE),
    
    run_speed_bl_mean = mean(run_speed/mean(spine_length,na.rm=T), na.rm = TRUE),
    run_speed_bl_var = var(run_speed/mean(spine_length,na.rm=T), na.rm = TRUE),
    
    run_speed_bl_bw_mean = mean(run_speed_bw/mean(spine_length,na.rm=T), na.rm = TRUE),
    run_speed_bl_bw_var = var(run_speed_bw/mean(spine_length,na.rm=T), na.rm = TRUE),
    
    
    run_speed_modulation = 
      (mean(run_speed[odor_orientation=="towards"],na.rm=TRUE)-
         mean(run_speed[odor_orientation=="away"],na.rm=TRUE))/
      (mean(run_speed[odor_orientation=="towards"],na.rm=TRUE)+
         mean(run_speed[odor_orientation=="away"],na.rm=TRUE)),
    
    
    IS_angle_mean = mean(IS_angle, na.rm = TRUE),
    IS_angle_var = var(IS_angle, na.rm = TRUE),
    
    IS_angle_bw_mean = mean(IS_angle_bw, na.rm = TRUE),
    IS_angle_bw_var = var(IS_angle_bw, na.rm = TRUE),
    
    Abs_IS_angle_mean = mean(Abs_IS_angle, na.rm = TRUE),
    Abs_IS_angle_var = var(Abs_IS_angle, na.rm = TRUE),
    
    Abs_IS_angle_bw_mean = mean(Abs_IS_angle_bw, na.rm = TRUE),
    Abs_IS_angle_bw_var = var(Abs_IS_angle_bw, na.rm = TRUE),
    
    IS_reorientation_mean = mean(IS_reorientation, na.rm = TRUE),
    IS_reorientation_var = var(IS_reorientation, na.rm = TRUE),
    
    IS_reorientation_bw_mean = mean(IS_reorientation_bw, na.rm = TRUE),
    IS_reorientation_bw_var = var(IS_reorientation_bw, na.rm = TRUE),
    
    IS_interval_mean = mean(IS_interval, na.rm = TRUE),
    IS_interval_var = var(IS_interval, na.rm = TRUE),
    
    IS_interval_bw_mean = mean(IS_interval_bw, na.rm = TRUE),
    IS_interval_bw_var = var(IS_interval_bw, na.rm = TRUE),
    
    IS_distance_mean = mean(IS_distance, na.rm = TRUE),
    IS_distance_var = var(IS_distance, na.rm =TRUE),
    
    IS_distance_bw_mean = mean(IS_distance_bw, na.rm = TRUE),
    IS_distance_bw_var = var(IS_distance_bw, na.rm =TRUE),
    
    IS_distance_bl_mean = mean(IS_distance_bl,na.rm=T),
    IS_distance_bl_var = var(IS_distance_bl,na.rm=T),
    
    IS_distance_bl_bw_mean = mean(IS_distance_bl_bw,na.rm=T),
    IS_distance_bl_bw_var = var(IS_distance_bl_bw,na.rm=T),
    
    step_extr_mean = mean(step_extr,na.rm=T),
    
    is_sillyStep = sum(step_extr,na.rm=T)>1
    
  ),by=grouped_by]
  
}


data <- readRDS("C:/MThane/TestData//CIRL-blind-2021-11-17-analysis/CIRL-blind-2021-11-17_analyzed.rds")
data <- data$data %>%
  group_by(id)%>%
  mutate(step_count = cumsum(step_boolean))

ids <- data.frame(id=sort(rep(unique(data$id),3)))

TS <- illav::summarise_TS_per_group_dt(as.data.table(data), c("group","condition","trial","id","step_count"))

TS <- as.data.frame(TS)

RUN <-  summarise_RUN_per_group_dt_(as.data.table(data), c("id","step_count"))%>%
  as.data.frame(TS)


result <- left_join(TS,RUN,by=c("id","step_count"))%>% 
  mutate (group_condition = as.factor(paste(group,condition,is_sillyStep,sep="-")))


#levels(result$group_condition) <- c("no-step","more-then-two-minima","zero-minima","one-minimum", "two-minima")


fwrite(result, "sillyWalk_summarySteps.csv")




data <- data %>%
  group_by(id)%>%
  mutate(step_count = cumsum(step_boolean))



ids <- data.frame(id=sort(rep(unique(data$id),3)))

TS <- summarise_TS_per_group_dt(as.data.table(data), c("group","condition","trial","id","minima"))

TS <- as.data.frame(TS)%>% 
  complete(id, nesting(minima))

RUN <-  summarise_RUN_per_group_dt(as.data.table(data), c("group","condition","trial","id","minima"))%>%
  as.data.frame(TS)%>% 
  complete(id, nesting(minima))


result <- left_join(TS,RUN,by=c("group","condition","trial","id","minima"))%>% 
  mutate (group_condition = as.factor(minima))


levels(result$group_condition) <- c("no-step","more-then-two-minima","zero-minima","one-minimum", "two-minima")


fwrite(result, "sillyWalk_summaryIds.csv")



