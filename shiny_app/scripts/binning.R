library(zoo)
library(gmodels)

binning_data_run_speed_grouped <-
  function(df,
           variable,
           width = 10,
           grouping="id"
  ) {
    message("binning run_speed data")
    df %>%
      filter(run_flag == 1) %>%
      filter(!is.na(run_speed)) %>%
      mutate(bin = cut_width(.[[variable]], width=length,labels=F)) %>%
      mutate(bin = sapply(bin, get_midpoint))%>%
      group_by_(grouping,"bin")%>%
      summarise(
        run_speed = mean(run_speed,na.rm=T),
        run_speed_bw = mean(run_speed,na.rm=T),
        run_speed_bl = mean(run_speed_bl,na.rm=T),
        run_speed_bl_bw = mean(run_speed_bl_bw,na.rm=T),
        run_speed_modulation =
          (mean(midpoint_speed[odor_orientation=="towards"],na.rm=T)-
             mean(midpoint_speed[odor_orientation=="away"],na.rm=T))/
          (mean(midpoint_speed[odor_orientation=="towards"],na.rm=T)+
             mean(midpoint_speed[odor_orientation=="away"],na.rm=T))
      )%>%
      group_by(bin)%>%
      summarise(
        run_speed.lowCI = ci(run_speed, na.rm = T)[2],
        run_speed.hiCI = ci(run_speed, na.rm = T)[3],
        run_speed.mean = mean(run_speed,na.rm=T),
        
        run_speed_bw.lowCI = ci(run_speed_bw, na.rm = T)[2],
        run_speed_bw.hiCI = ci(run_speed_bw, na.rm = T)[3],
        run_speed_bw.mean = mean(run_speed_bw,na.rm=T),
        
        run_speed_bl.lowCI = ci(run_speed_bl, na.rm = T)[2],
        run_speed_bl.hiCI = ci(run_speed_bl, na.rm = T)[3],
        run_speed_bl.mean = mean(run_speed_bl,na.rm=T),
        
        run_speed_bl_bw.lowCI = ci(run_speed_bl_bw, na.rm = T)[2],
        run_speed_bl_bw.hiCI = ci(run_speed_bl_bw, na.rm = T)[3],
        run_speed_bl_bw.mean = mean(run_speed_bl_bw,na.rm=T),
        
        run_speed_modulation.lowCI = ci(run_speed_modulation, na.rm = T)[2],
        run_speed_modulation.hiCI = ci(run_speed_modulation, na.rm = T)[3],
        run_speed_modulation.mean = mean(run_speed_modulation,na.rm=T)
      )%>%mutate(bin_id=seq(1,nrow(.)))
    
  }



data <- readRDS("C:/MThane/TestData//CIRL-blind-2021-11-17-analysis/CIRL-blind-2021-11-17_analyzed.rds")
binning_data_run_speed_grouped(data$data,variable = "frame")


bn <-data$data %>%
  filter(!is.na(run_speed)) %>%
  mutate(bin = cut(frame,seq(0,180,10)))
View(bn)
