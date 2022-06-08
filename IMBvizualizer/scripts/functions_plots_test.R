
##### SINGLE MODE: PLOTS #####

plot_single_histogram <-
  function(track,
           variable,
           interval,
           bins) {
    interval = interval*16
    data <- track %>%
      filter(frame %in% seq(from = interval[1],
                            to = interval[2]))
    p <- ggplot(data, aes_string(variable)) +
      geom_histogram(
        bins = bins,
        color = 'grey',
        alpha = 0.7,
        na.rm = T
      ) +
      ylab("Number of observations")+
      theme_classic()
    ggplotly(p)
    
  }


#####  line plots #####
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
  
  #build plot
  p <- ggplot(data_pivot, aes(x = time, y = value,color=variable)) +
    geom_line()+
    scale_x_continuous(breaks=seq(0,180,30))+
    xlab("Time in seconds")+
    theme_classic()+ 
    scale_color_brewer(palette="Set2")+
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
      filter(step_boolean==T | step_boolean_bw==T)%>%
      mutate(idx = frame/frame_rate)%>%
      select(idx,id_gc)
  
    if(nrow(df_steps)>0){
      p<- p+ geom_vline(
        data = df_steps,
        aes(group=id_gc,xintercept=idx),
        alpha = 0.1)
    }
    
  }
  p <- p
  
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
           sel_row=NULL,
           violin=FALSE
  ){
    if(!is.null(sel_row)){
      tickColors <- ifelse(unique(data$group_condition) %in%  
                             c(sel_row$X1, sel_row$X2), "red", "black")
      
    }else{
      tickColors="black"
    }
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



binning_data_plot <- function(data,
                              x,
                              y,
                              radius,
                              binning_vartable,
                              npoints = 11,
                              circ=FALSE,
                              binning_mode = "dish"
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
  
  if(y %in% ci_vars | binning_mode!="all"){
    data <- data %>%
      group_by(group_condition) %>%
      summarise(
        N = N,
        bin = bin,
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
        bin = bin,
        variable = stats::filter(!!as.symbol(y),
                                 rep(1, npoints) / npoints, circular = circ))
    
    
  }
  
  message("binning finished!")
  
  if(x!="bearing_angle"){
    fig <- plot_ly(data, x = ~bin, y = ~variable, color = ~group_condition,text=~N)
    fig %>% add_lines( line = list(
      x = ~bin, y = ~variable,
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
              add_ribbons(x = ~bin,
                          ymin = ~lower,
                          ymax = ~upper
              ) %>%
              add_lines(
                x = ~bin,
                y = ~variable,
                line = list(
                  dash = "dot"
                )
              )
          }
          )
      }
    }
  }else{
    
    
    fig <- plot_ly(
      data, color = ~group_condition,
      theta= ~bin, r = ~variable,
      type = 'scatterpolar',
      mode = 'lines',
      subplot = 'polar'
    )
    
    
    # if confidence interval can be calculated
    if(y %in% ci_vars| binning_mode!="all"){
      fig <- fig %>% add_trace(
        theta = ~bin, r = ~lower,
        line = list(
          dash = 'dot'
        )
      )
      fig <- fig %>% add_trace( 
        theta = ~bin, r = ~upper,
        line = list(
          dash = 'dot'
        )
      )
      
      
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

