
##### SINGLE MODE: PLOTS #####

#' Title
#'
#' @param track 
#' @param variable 
#' @param interval 
#' @param bins 
#' @param use_density 
#'
#' @return
#' @import ggplot2
#' @import dplyr
#' @import plotly
#' @export
#'
plot_single_histogram <-
  function(track,
           variable,
           interval,
           bins,
           use_density = F
           ) {
    interval = interval*16
    data <- track %>%
      filter(frame %in% seq(from = interval[1],
                            to = interval[2]))%>%
      mutate(id = factor(id))
    
    p <- ggplot(data, aes_string(variable,fill="id")) +

      ylab("Number of observations")+
      theme_classic()
    if(use_density){
      p <- p + geom_density(alpha=0.5)
    }else{
      p <- p + geom_histogram(
        bins = bins,
        alpha = 0.5,
        na.rm = T
      ) 
    }
    
    ggplotly(p)
    
  }


#####  line plots #####
#' Title
#'
#' @param track 
#' @param variables 
#' @param SINGLE_VARTABLE 
#' @param interval 
#' @param scale 
#' @param use_limits 
#' @param limits 
#' @param frame_rate 
#' @import ggplot2
#' @import dplyr
#' @import plotly
#' @return
#' @export
#'
plot_timeseries <- function(track,
                            variables,
                            SINGLE_VARTABLE,
                            interval,
                            scale=T,
                            use_limits=F,
                            limits=c(-5,5),
                            frame_rate = 16
){
  
  
  #prepare data
  
  interval = interval*frame_rate
  if(scale){
    data <- track%>%
      mutate_at(variables, funs(c(scale(.))))
  }else{
    data <- track
  }
  
  ### changing labels to variable names
  varNames <- SINGLE_VARTABLE %>%
    filter(variable %in% variables)%>%
    select(name,variable)%>%
    na.omit()
  varNames$idx <- match(varNames$variable,colnames(data))
  varNames <- na.omit(varNames)
  colnames(data)[varNames$idx] <- unlist(varNames$name)
  
  
  data <- data%>%
    mutate(time = frame/frame_rate,
           id_gc = paste(group_condition,id,sep="-"))%>%
    filter(frame >= interval[1] &
             frame <= interval[2])
  
  data_pivot <- data %>%
    pivot_longer(cols=varNames$name,
                 names_to="variable")%>%
    select(id_gc,
           time,
           value,
           variable)%>%
    na.omit()
  min_time = min(data_pivot$time,na.rm=T)
  max_time = max(data_pivot$time,na.rm=T)
  #build plot
  manual_colors = c("#0000FF", "#000000",RColorBrewer::brewer.pal(8,"Set1"))
  
  p <- ggplot(data_pivot, aes(x = time, y = value,color=variable)) +
    geom_line(alpha=0.8)+
    
    scale_x_continuous(breaks=seq(min_time,max_time,20))+
    xlab("Time in seconds")+
    theme_classic()+ 
    scale_color_manual(values=manual_colors)+
    #scale_color_viridis(discrete = TRUE, option = "D")+
    facet_grid(rows = vars(id_gc))+
    scale_x_continuous(breaks = scales::pretty_breaks(n = 10))
  
  
  #adding vlines  
  for (id_ in unique(data$id_gc)){
    data_ <- data %>%filter(id_gc==id_ &
                              frame >= interval[1] &
                              frame <= interval[2])
    
    start <- na.omit(data$time)[1]
    # left HCs
    df_hc_left <- data.frame(idx = data_$frame[which(as.logical(data_$HCs_left))]/frame_rate)
    df_hc_left$id_gc <- rep(id_,nrow(df_hc_left))
    df_hc_left$direction <- rep("red",nrow(df_hc_left))
    # right HCs
    df_hc_right <- data.frame(idx = data_$frame[which(as.logical(data_$HCs_right))]/frame_rate)
    df_hc_right$id_gc <- rep(id_,nrow(df_hc_right))
    df_hc_right$direction <- rep("green",nrow(df_hc_right))
    
    switching_point_negative <-data.frame(
      list(list(idx =data_$frame[which(data_$flanks==-1)]/frame_rate,
                id_gc = rep(id_,length(which(data_$flanks==-1)))
                )))
    
    
    switching_point_positive <-data.frame(
      list(list(idx =data_$frame[which(data_$flanks==1)]/frame_rate,
                id_gc = rep(id_,length(which(data_$flanks==1)))
      )))
    
    df_hc <- rbind(df_hc_left,df_hc_right)
    if(nrow(df_hc_right)>0){
      p<- p+ geom_vline(
        data = df_hc_right,
        aes(group=id_gc,xintercept=idx),
        color="green",
        alpha = 0.1)
    }
    if(nrow(df_hc_left)>0){
      p<- p+ geom_vline(
        data = df_hc_left,
        aes(group=id_gc,xintercept=idx),
        color="red",
        alpha = 0.1)
    }
    df_steps <- data%>%
      filter(step_boolean==T)%>%
      mutate(idx = (1+frame)/frame_rate)%>% 
      select(idx,id_gc)
# new after Thane et al 2023: added a +1 to correct for a shifted display
    
    if(nrow(df_steps)>0){
      p<- p+ geom_vline(
        data = df_steps,
        aes(group=id_gc,xintercept=idx),
        alpha = 0.1)
    }
    
    #add vlines for switching point
    
    # if(nrow(switching_point_positive)>0){
    #   p <-  p+ geom_vline(data=switching_point_positive,
    #                       aes(group= id_gc,
    #                           xintercept = idx))
    # }
    # if(nrow(switching_point_negative)>0){
    #   p <-  p+ geom_vline(data=switching_point_negative,
    #                       aes(group= id_gc,
    #                           xintercept = idx),
    #                       color="blue",
    #                       linetype="dashed")
    # }
    
  }
  
  #adding hlines
  if ("head_vector_angular_speed" %in% variables) {
    p <- p + geom_hline(yintercept = 35,
                        linetype = "dotted",
                        color = "red") +
      geom_hline(yintercept = -35,
                 linetype = "dotted",
                 color = "green")
  }
  
  if ("tail_vector_angular_speed" %in% variables) {
    p <- p + geom_hline(
      yintercept = 45,
      linetype = "dotted",
      color = "blue",
      alpha = 0.3) + 
      geom_hline(
        yintercept =-45,
        linetype = "dotted",
        color = "blue",
        alpha = 0.3)
  }
  
  if(use_limits){
    p<- p+ ylim(limits)
    
  }
  pl <-ggplotly(p,height=100+length(unique(track$id))*250)%>% 
    layout(legend = list(orientation = 'h',y=5),margin=list(t=50))
  pl
}



##### BOXPLOTS

#' Title
#'
#' @param data 
#' @param variable 
#' @param groupColorValues 
#' @param groupNameValues 
#' @param remove_outliers 
#' @param ylim 
#' @param show_zero 
#' @param height 
#' @param width 
#' @param margin 
#' @param titlefontsize 
#' @param axisfontsize 
#' @param violin 
#'
#' @return
#' @import ggplot2
#' @import dplyr
#' @import plotly
#' @export
#'
boxplot_per_group <-
  function(data,
           variable,
           groupColorValues,
           groupNameValues,
           remove_outliers = F,
           ylim = NULL,
           show_zero = FALSE,
           height=200,
           width=200,
           margin= list(l=50,r=50,b=50,t=50,pad=4),
           titlefontsize = 16,
           axisfontsize = 12,
           #sel_row=NULL,
           violin=FALSE
  ){
    # if(!is.null(sel_row)){
    #   tickColors <- ifelse(unique(data$group_condition) %in%  
    #                          c(sel_row$X1, sel_row$X2), "red", "black")
    #   
    # }else{
    #   tickColors="black"
    # }
    v <- variable$variable
    #names(groupColorValues) <- groupNameValues
    if (v %in% c(
      "HC_angle",
      "HC_rate_modulation",
      "HC_reorientation",
      "preference",
      "pref_dist",
      "odor_speed"
    )) {
      show_zero = TRUE
    }
    if (v %in% c("preference", "preference_dist")) {
      ylim = c(-1, 1)
    }
    
    names(groupColorValues) <- groupNameValues
    
    #show number of observations
    plot <-
      ggplot(data,
             aes_string(x = "group_condition",
                        y = v,
                        fill="group_condition")) +
      
      labs(title=paste(variable$name," [",variable$unit,"]",sep=""),
           x = "Experimental condition",
           y = paste(variable$name," [",variable$unit,"]",sep="")
      )+
      stat_n_text() 
    
    if(remove_outliers){
      outlier_shape = NA
      plot <- plot + geom_boxplot()
      
    }else{
      outlier_shape=1
    }
    
    if(violin){
      plot <- plot +    
        geom_violin(width=3.0,alpha=0.5,outlier.shape = outlier_shape)  +
        geom_boxplot(width=0.2,alpha=0.8, outlier.shape = outlier_shape)
      # stat_summary(
      #   fun.data = "mean_sdl",  fun.args = list(mult = 1), 
      #   geom = "pointrange", color = "black"
      # )
    }else{
      plot <- plot + geom_boxplot( outlier.shape = outlier_shape)
    }
    
    
    
    # showing zero line
    if (show_zero) {
      plot <- plot + geom_hline(yintercept = 0,
                                alpha = 0.5,
                                linetype = "dashed")
    }
    if(!is.null(unlist(groupColorValues))){
      plot <- plot + scale_color_manual(values= groupColorValues)
    }
    
    if (!is.null(ylim)) {
      plot <- plot +
        ylim(ylim[1], ylim[2])
    }
    
    if(variable$unit==""){
      plot <- plot+labs(
        
        title = variable$name,
        
        y = variable$name
      )
    }
    plot <- plot+
      theme_classic()+
      scale_x_discrete(labels = function(x) str_wrap(x, width = 10))+
      theme(axis.text.x = element_text(angle = 45),
            
            legend.position ="none", 
            legend.spacing.x = unit(3.0, 'cm')
      )
    
    
    # Need to modify the plotly object and make outlier points have opacity equal to 0
    plot <- plotly_build(plot)
    if(remove_outliers){
      for(i in 1:length(plot$x$data)) {
        plot$x$data[[i]]$marker$opacity = 0
      }
    }
    
    plot%>%
      config(
        toImageButtonOptions = list(
          format = "svg",
          filename = paste0(variable$variable,"_boxplot"),
          width = 600,
          height = 700
        )
      )
  }



#' Title
#'
#' @param data 
#' @param x 
#' @param y 
#' @param radius 
#' @param binning_vartable 
#' @param npoints 
#' @param circ 
#' @param binning_mode 
#' @param use_polar 
#'
#' @return
#' @import ggplot2
#' @import dplyr
#' @import plotly
#' @import stats
#' @export
#'
binning_data_plot <- function(data,
                              x,
                              y,
                              radius,
                              binning_vartable,
                              npoints = 11,
                              circ=FALSE,
                              binning_mode = "dish",
                              use_polar = T
){
  message("binning data ...")
  
  if(x=="bearing_angle"){
    circ=TRUE
  }
  # if confidence interval can be calculated
  ci_vars <- binning_vartable%>%
    filter(has_ci==TRUE)%>%
    select(variable)%>%
    unlist()
  
  xvar <- binning_vartable %>%
    filter(variable==x)%>%
    select(name,unit)
  xvariable_name <-  paste(xvar$name," [",xvar$unit,"] ",sep="")
  
  yvar <- binning_vartable %>%
    filter(variable==y)%>%
    select(name,unit,description)
  yvariable_name <-  paste(yvar$name," [",yvar$unit,"] ",sep="")
  
  
  
  if(x=="frame"){
    xvariable_name <- "Time [s]"
  }
  
  if(x=="time_to_trans_to"){
    xvariable_name <- "Time until next transition [s]"
  }
  
  if(x=="time_to_trans_from"){
    xvariable_name <- "Time until next transition [s]"
  }
  
  if(x=="spinepoint_y_6_conv"){
    xvariable_name <- "Y coordinate"
  }
  
  if(x=="spinepoint_x_6_conv"){
    xvariable_name <- "X coordinate"
  }
  
  if(y %in% ci_vars | binning_mode!="all"){
    data <- data %>%
      group_by(group_condition) %>%
      summarise(
        N = N,
        midpoint = midpoint,
        variable = stats::filter(!!as.symbol(paste0(y,".mean")),
                                 rep(1, npoints) / npoints,
                                 circular = circ),
        lower = as.numeric(
          stats::filter(
            !!as.symbol(paste0(y, ".lowCI")),
            rep(1, npoints) / npoints,
            circular = circ
          )
        ),
        upper = as.numeric(
          stats::filter(
            !!as.symbol(paste0(y, ".hiCI")),
            rep(1, npoints) / npoints,
            circular = circ
          )
        )
        
        
      ) %>%
      ungroup() %>%
      na.omit()
    
  }else{
    data <- data %>%
      group_by(group_condition) %>%
      summarise(
        N = N,
        midpoint = midpoint,
        variable = stats::filter(!!as.symbol(y),
                                 rep(1, npoints) / npoints, circular = circ))
    
    
  }
  
  message("binning finished!")
  
  if(x=="bearing_angle"& use_polar){
    
    fig <- plot_ly(
      data, color = ~group_condition,
      theta= ~midpoint, r = ~variable,
      type = 'scatterpolar',
      mode = 'lines',
      subplot = 'polar'
    )
    
    
    # if confidence interval can be calculated
    if(y %in% ci_vars| binning_mode!="all"){
      fig <- fig %>% add_trace(
        theta = ~midpoint, r = ~lower,
        line = list(
          dash = 'dot'
        )
      )
      fig <- fig %>% add_trace( 
        theta = ~midpoint, r = ~upper,
        line = list(
          dash = 'dot'
        )
      )
      
      
    }
  }else{
    fig <- plot_ly(data, x = ~midpoint, y = ~variable, color = ~group_condition,text=~N)
    fig %>% add_lines( line = list(
      x = ~midpoint, y = ~variable,
      dash = "dot"
    ))
    
    # if confidence interval can be calculated
    if(y %in% ci_vars| binning_mode!="all"){
      g.uni <- unique(data$group_condition)
      for(i in g.uni){
        fig <- fig  %>%
          add_fun(function(plot) {
            plot %>% 
              filter(group_condition == i) %>%
              add_ribbons(x = ~midpoint,
                          ymin = ~lower,
                          ymax = ~upper
              ) %>%
              add_lines(
                x = ~midpoint,
                y = ~variable,
                line = list(
                  dash = "dot"
                )
              )
          }
          )
      }
    }

  }
  
  m <- list(
    l = 50,
    r = 50,
    b = 100,
    t = 100,
    pad = 4
  )
  
  fig <- fig%>%
    layout(title= paste(yvariable_name,xvariable_name, sep = " with respect to "),
           xaxis = list(title = xvariable_name, showline = TRUE,showgrid = F),
           yaxis = list(title = yvariable_name, showline = TRUE,showgrid = F),
           template = "ggplot2",
           margin=m
    )%>%
    config(
      toImageButtonOptions = list(
        format = "svg",
        filename = paste(x,y,"lineplot",sep="_"),
        width = 600,
        height = 700
      )
    )
  
  if(x=="bearing_angle"){
    fig <- fig%>%layout( 
      
      legend = list(orientation = 'h'),
      polar = list(
        radialaxis = list(
          angle = 90
        ),
        angularaxis = list(
          direction = 'clockwise',
          period = 6
        )
      )
      
    )
  }
  fig
}


#' Title
#'
#' @param data 
#' @param x 
#' @param y 
#' @param z 
#' @param nbins 
#' @param filter_radius 
#' @param xlabel 
#' @param ylabel 
#' @param ratio_fixed 
#' @param limits 
#' @param radius 
#' 
#' @return
#' 
#' @import ggplot2
#' @import dplyr
#' @import plotly
#' @import imagine
#' @import reshape2
#' @export
#'
create2dHeatmap <- function(data,x,y,z,
                            nbins=50,
                            filter_radius = 5,
                            xlabel="",ylabel="",
                            ratio_fixed = F,
                            limits=c(0,100),
                            radius=NA
){
  #create normal 2d binning plot and extract data
  if(z=="HC_rate"){
    z="HCs"
  }
  
  if(y=="time"| x=="time"){
    data <- data %>% mutate(time=frame/16)
  }
  
  if(z=="abs_heading_angle"){
    data = data%>%mutate(abs_heading_angle = abs(heading_angle))
  }
  if(z=="abs_bearing_angle"){
    data = data%>%mutate(abs_bearing_angle = abs(bearing_angle))
  }
  if(z=="abs_bending_angle"){
    data = data%>%mutate(abs_bending_angle = abs(bending_angle))
  }
  if(z=="Abs_HC_angle"){
    data = data%>%mutate(Abs_HC_angle = abs(Abs_HC_angle))
  }
  
  if(z=="Abs_IS_angle"){
    data = data%>%mutate(Abs_HC_angle = abs(Abs_IS_angle))
  }
  
  p <- ggplot(data, aes_string(x, y, z = z))+
    stat_summary_2d(bins = nbins)
  df <- ggplot_build(p)$data[[1]]%>%
    dplyr::select("x","y","value")
  colnames(df) <-c("x","y","z")
  # convert plot into matrix
  out <- dcast(setDT(df), x ~ y, value.var ="z", sep = "",drop=F)%>%
    as.data.frame()
  M <- out %>% select(!x)%>%as.matrix()
  rownames(M)<- out$x
  
  #filter matrix using imagine
  filtered     <- medianFilter(X = M,radius=filter_radius,times=10)
  colnames(filtered)<- colnames(M)
  rownames(filtered)<- rownames(M)
  df <-   melt(filtered)
  p<-ggplot(df)+
    geom_raster(aes(Var1,Var2,fill=value))+
    scale_fill_viridis_c(limits=limits,
                         #breaks= seq(limits[1],limits[2],limits[2]/5),
                         oob=squish)+
    theme_classic()+
    xlab(xlabel)+
    ylab(ylabel)
  
  if(!is.na(radius)){
    p <- p+annotate("path",
                    x = 0 + radius * cos(seq(0, 2 * pi, length.out = 100)),
                    y = 0 + radius * sin(seq(0, 2 * pi, length.out = 100))) 
  }
  
  if(ratio_fixed){
    p <-p+
      coord_fixed()
  }
  
  ggplotly(p)%>%
    config(
      toImageButtonOptions = list(
        format = "svg",
        filename = paste0("heatmap"),
        width = 600,
        height = 700
      )
    )
  
}


empty_plot <- function(message){
  
  p <- ggplot(data.frame(c(0,0),c(0,0))) +
    theme_classic() +
    geom_text(aes(0,0,label=message)) +
    ggtitle("Empty Plot") #optional, but safer in case another theme is applied later
    
  ggplotly(p)
  
}

