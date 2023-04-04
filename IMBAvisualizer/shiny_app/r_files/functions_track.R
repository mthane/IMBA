

#' Title
#'
#' @param track 
#' @param in_columns 
#' @param odor_pos 
#'
#' @return
#' @export
#' @import stringr
#' @import dplyr
#' 
#' 
calculate_rotate_columns <- function(track, in_columns, odor_pos) {
  track <- track
  odor_vector_len = distance2(c(0, 0), odor_pos)
  angle <- clockwise_angle(c(0, odor_vector_len), odor_pos)
  odor_a_rot  = data.frame(rep(0, nrow(track)), rep(odor_vector_len, nrow(track)))
  colnames(odor_a_rot) <- c("odor_a_x_rot", "odor_a_y_rot")
  columns = list()
  i = 1
  for (col in in_columns) {
    col_x <- unlist(stringr::str_split(col, " "))[1]
    col_y <- unlist(stringr::str_split(col, " "))[2]
    x = track[col_x]
    y = track[col_y]
    x_center = 0
    y_center = 0
    M <- cbind(x, y)
    alpha <- deg2rad(angle)
    rotm <-
      matrix(c(cos(alpha), sin(alpha),-sin(alpha), cos(alpha)), ncol = 2)
    M2 <- t(rotm %*% t(M))
    data <- M2#as.data.frame(M2)
    colnames(data) <-
      c(paste0(col_x, "_rot"), paste0(col_y, "_rot"))
    columns[[i]] <- data
    i = i + 1
  }
  data <- bind_cols(c(columns, odor_a_rot))
  
  return(data)
  
}

#' Title
#'
#' @param track
#' @param frame_rate the frame rate of the tracked videos
#'
#' @return
#' @export
#' @import RcppRoll


calculate_convolve_columns <-
  function(track,
           in_columns,
           column_names,
           frame_rate = 16) {
    df <-
      roll_mean(as.matrix(track[in_columns]), round(frame_rate * 0.3), fill = NA)
    colnames(df) <- column_names
    return(df)
  }






#' Title
#'
#' @param track 
#' @param head_start_pos 
#' @param head_end_pos 
#' @param tail_start_pos 
#' @param tail_end_pos 
#'
#' @return
#' @export
#'
calculate_head_tail_vectors <- function(track,
                                        head_start_pos = c("spinepoint_x_9_conv",
                                                           "spinepoint_y_9_conv"),
                                        head_end_pos =
                                          c("spinepoint_x_11_conv",
                                            "spinepoint_y_11_conv"),
                                        tail_start_pos =  c("spinepoint_x_2_conv",
                                                            "spinepoint_y_2_conv"),
                                        tail_end_pos = c("spinepoint_x_6_conv",
                                                         "spinepoint_y_6_conv")) {
  ### calculating head vector and tail vector ###
  ### from Emmanouil Paisios et. al 2017 ###
  tail_vector <-
    cbind(track[tail_end_pos[1]] - track[tail_start_pos[1]],
          track[tail_end_pos[2]] - track[tail_start_pos[2]])
  tail_vector <-
    tail_vector / colSums(t((as.matrix(tail_vector))) ** 2)
  head_vector <-
    cbind(track[head_end_pos[1]] - track[head_start_pos[1]],
          track[head_end_pos[2]] - track[head_start_pos[2]])
  
  df <- cbind(tail_vector, head_vector)
  colnames(df) <-
    c("tail_vector_x",
      "tail_vector_y",
      "head_vector_x",
      "head_vector_y")
  
  df
}



#' Title
#'
#' @param track 
#'
#' @return
#' @export
#'
calculate_spine_length_sum <- function(track) {
  track <- track
  distances <- rep(0, nrow(track))
  for (i in 1:11) {
    A = cbind(track[paste0("spinepoint_x_", i)] - track[paste0("spinepoint_x_", i +
                                                                 1)])
    B = cbind(track[paste0("spinepoint_y_", i)] - track[paste0("spinepoint_y_", i +
                                                                 1)])
    distance <- sqrt(A ** 2 + B ** 2)
    distances <- distances + distance
    #dist_sum <- dist_sum + distance
  }
  unlist(distances)
}


#' Title
#'
#' @param df 
#' @param frame_rate 
#'
#' @return
#' @export
#'
calculate_speed_direction_odor <- function(df, frame_rate = 16) {
  odor_vector_x = df$odor_a_x_rot - df$spinepoint_x_6_conv
  odor_vector_y = df$odor_a_y_rot - df$spinepoint_y_6_conv
  
  odor_vector <- cbind(odor_vector_x, odor_vector_y)
  
  odor_vector <-
    odor_vector / colSums(t((as.matrix(odor_vector))) ** 2)
  
  diff_x = c(diff(as.numeric(unlist(df["spinepoint_x_6_conv"]))), NA)
  diff_y = c(diff(as.numeric(unlist(df["spinepoint_y_6_conv"]))), NA)
  
  odor_speed <-
    frame_rate * (odor_vector[, 1] * diff_x + odor_vector[, 2] * diff_y)
  odor_speed
}


#' Title
#'
#' @param df 
#' @param frame_rate 
#'
#' @return
#' @export
#'
#' @examples
calculate_head_vector_angular_speed <-
  function(df, frame_rate = 16) {
    ### from Emmanouil Paisios et. al 2017 ###
    #   # head_vector_angular_speed
    # calculate differnce of angles
    len = length(df$head_vector_x)
    dt = 1 / frame_rate
    M = cbind(df$head_vector_x[1:len - 1],
              df$head_vector_y[1:len - 1],
              df$head_vector_x[2:len],
              df$head_vector_y[2:len])
    
    head_vector_angular_speed <-
      frame_rate * clockwise_angle_v(M[, 1], M[, 2], M[, 3], M[, 4])
    
    
    npoints = round(0.45 * frame_rate)
    x = head_vector_angular_speed
    y = rep(1, npoints) / npoints
    if (length(x) > length(y) & frame_rate > 8) {
      head_vector_angular_speed <- stats::filter(x, y)
    }
    
    head_vector_angular_speed = c(head_vector_angular_speed, NA)
    head_vector_angular_speed
  }


#' Title
#'
#' @param df 
#' @param frame_rate 
#'
#' @return
#' @export
#'
#' @examples
calculate_tail_vector_angular_speed <-
  function(df, frame_rate = 16) {
    ### from Emmanouil Paisios et. al 2017 ###
    #   # head_vector_angular_speed
    # calculate differnce of angles
    len = length(df$tail_vector_x)
    M = cbind(df$tail_vector_x[1:len - 1],
              df$tail_vector_y[1:len - 1],
              df$tail_vector_x[2:len],
              df$tail_vector_y[2:len])
    tail_vector_angular_speed <-
      frame_rate * clockwise_angle_v(M[, 1], M[, 2], M[, 3], M[, 4])
    c(tail_vector_angular_speed, NA)
  }


#' Title
#'
#' @param df 
#' @param frame_rate 
#'
#' @return
#' @export
#'
#' @examples
calculate_head_vel_forward <- function(df, frame_rate = 16) {
  ### from Emmanouil Paisios et. al 2017 ###
  diff_x = c(diff(as.numeric(unlist(df["spinepoint_x_12_conv"]))), NA)
  diff_y = c(diff(as.numeric(unlist(df["spinepoint_y_12_conv"]))), NA)
  frame_rate * (df['tail_vector_x'] * diff_x + df['tail_vector_x'] *
                  diff_y) %>% unlist()
  
}

#' Title
#'
#' @param df 
#' @param frame_rate 
#'
#' @return
#' @export
#'
#' @examples
calculate_tail_vel_forward <- function(df, frame_rate = 16) {
  ### from Emmanouil Paisios et. al 2017 ###
  diff_x = c(diff(as.numeric(unlist(df["spinepoint_x_1_conv"]))), NA)
  diff_y = c(diff(as.numeric(unlist(df["spinepoint_y_1_conv"]))), NA)
  tail_vel_forward <-
    frame_rate * (df['tail_vector_x'] * diff_x + df['tail_vector_y'] * diff_y)
  npoints = round(0.3 * frame_rate)
  if (frame_rate > 8) {
    stats::filter(tail_vel_forward,
                  rep(1 / npoints, npoints))
  } else{
    tail_vel_forward
  }
  
  
}


#' Title
#'
#' @param track 
#'
#' @return
#' @export
#' @import stats
#'
#' @examples
calculate_na_flag_ht_flip <- function(track) {
  is_invalid <- abs(track$head_vector_angular_speed) > 180 |
    track$tail_vector_angular_speed > 180
  as.logical(as.integer(stats::filter(as.integer(is_invalid), rep(1, 5))))
}



#' Title
#'
#' @param df 
#'
#' @return
#' @export
#'
calculate_bearing_angle  <-
  function(df) {
    odor_vector_x <- df$odor_a_x_rot - df$spinepoint_x_6_conv
    odor_vector_y <- df$odor_a_y_rot - df$spinepoint_y_6_conv
    
    M <-
      data.frame(df$tail_vector_x,
                 df$tail_vector_y,
                 odor_vector_x,
                 odor_vector_y)
    clockwise_angle_v(M[, 1], M[, 2], M[, 3], M[, 4])
  }


#' Title
#'
#' @param df 
#'
#' @return
#' @export
#'
calculate_heading_angle <-
  function(df) {
    M <-
      df[c(
        'head_vector_x',
        'head_vector_y',
        'spinepoint_x_9_conv',
        'spinepoint_y_9_conv',
        'odor_a_x_rot',
        'odor_a_y_rot'
      )]
    odor_vector_x = M$odor_a_x_rot - M$spinepoint_x_9_conv
    odor_vector_y = M$odor_a_y_rot - M$spinepoint_y_9_conv
    
    
    
    clockwise_angle_v(M[, 1], M[, 2], odor_vector_x, odor_vector_y)
  }


#' Title
#'
#' @param tail_vector_x 
#' @param tail_vector_y 
#' @param head_vector_extended_x 
#' @param head_vector_extended_y 
#' @param frame_rate 
#'
#' @return
#' @export
#'
calculate_bending_angle <- function(tail_vector_x,
                                    tail_vector_y,
                                    head_vector_extended_x,
                                    head_vector_extended_y,
                                    frame_rate) {
  if (frame_rate > 8) {
    stats::filter(
      -clockwise_angle_v(
        tail_vector_x,
        tail_vector_y,
        head_vector_extended_x,
        head_vector_extended_y
      ),
      rep(1,  round(0.45 * frame_rate)) /  round(0.45 * frame_rate)
    )
  } else{
    -clockwise_angle_v(tail_vector_x,
                       tail_vector_y,
                       head_vector_extended_x,
                       head_vector_extended_y)
  }
  
}

#' Title
#'
#' @param df 
#' @param frame_rate 
#' @param threshold_head_vector_angular_speed 
#' @param threshold_tail_vector_angular_speed 
#' @param min_hc_size 
#' @param max_hc_size 
#'
#' @return
#' @export
#' @import dplyr
calculate_HC_variables <-  function(df,
                                    frame_rate = 16,
                                    threshold_head_vector_angular_speed = 35,
                                    threshold_tail_vector_angular_speed = 45,
                                    min_hc_size = 0,
                                    max_hc_size = Inf
                                    ) {
  # ## from Emmanouil Paisios et. al 2017 ###
  HCs_left <-
    df$head_vector_angular_speed > threshold_head_vector_angular_speed
  HCs_right <-
    df$head_vector_angular_speed < (-threshold_head_vector_angular_speed)
  
  HC_index <- (HCs_right | HCs_left) #& HC_condition_tail
  
  HC_index[is.na(HC_index)] <- FALSE
  
  if (all(!HC_index)) {
    return(data.frame(
      list(
        HCs = rep(0, nrow(df)),
        HC_flag = rep(0, nrow(df)),
        HCs_right = rep(0, nrow(df)),
        HCs_left = rep(0, nrow(df)),
        HC_angle = rep(NA, nrow(df)),
        Abs_HC_angle = rep(NA, nrow(df)),
        HC_reorientation = rep(NA, nrow(df))
      )
    ))
  }
  
  if (HC_index[1] == TRUE) {
    HC_index[1] = FALSE
  }
  if (tail(HC_index, 1) == TRUE) {
    HC_index[tail(HC_index, 1)] = FALSE
  }
  
  HC_index_change <- c(diff(HC_index), NA)
  
  HC_start_idx <- which(HC_index_change == 1)
  HC_end_idx <- which(HC_index_change == -1)
  
  HC_intervals <-
    data.frame(start = HC_start_idx , end = HC_end_idx)
  
  
  true_HCs <- function(x) {
    values <- abs(df$head_vector_angular_speed[x[1]:x[2]])
    
    return(x[1] + which(values == max(values, na.rm = T))[1])
  }
  
  tail_condition <- function(x) {
    max(abs(df$tail_vector_angular_speed[x[1]:x[2]]), na.rm = T) <
      threshold_tail_vector_angular_speed &
      min(abs(df$tail_vector_angular_speed[x[1]:x[2]]), na.rm = T) > -5
  }
  
  size_condition <- function(x) {
    hc_size <- abs(df$bending_angle[x[1]]-df$bending_angle[x[2]])
    hc_size>min_hc_size & hc_size < max_hc_size
  }
  
  
  HC_intervals$tail_condition <-
    apply(HC_intervals, 1, function(x)
      tail_condition(x))
  
  HC_intervals$size_condition <-
    apply(HC_intervals, 1, function(x)
      size_condition(x))
  
  
  HC_intervals$HC_true <-
    apply(HC_intervals, 1, function(x)
      true_HCs(x))
  
  HC_intervals <- HC_intervals %>%
    filter(tail_condition == TRUE)%>%
    filter(size_condition == TRUE)
  
  HCs <- rep(0, nrow(df))
  if (nrow(HC_intervals) > 0) {
    HCs[unlist(HC_intervals$HC_true)] = 1
    
  }
  
  if (nrow(HC_intervals) > 0) {
    HC_index <- c()
    for (i in 1:nrow(HC_intervals)) {
      start <- HC_intervals[i, "start"]
      end <- HC_intervals[i, "end"]
      
      #if(!is.na(start) & !is.na(end)) {
      HC_index <- c(HC_index, seq(start, end))
      #}
    }
    
  }
  HC_flag <- rep(0, nrow(df))
  HC_flag[HC_index] <- 1
  
  
  
  HC_angle <-
    interval_difference(df,
                        HC_intervals$start,
                        HC_intervals$end,
                        "bending_angle",
                        HC_intervals$HC_true,
                        function(x, y)
                          y - x
                        )
  Abs_HC_angle <-
    interval_difference(df,
                        HC_intervals$start,
                        HC_intervals$end,
                        "bending_angle",
                        HC_intervals$HC_true,
                        function(x, y)
                          abs(y - x)
                        )
  
  HC_reorientation <-
    interval_difference(df,
                        HC_intervals$start,
                        HC_intervals$end,
                        "heading_angle",
                        HC_intervals$HC_true,
                        function(x, y)
                          abs(x) - abs(y)
                        
                        )
  
  
  HC_left = which(HC_flag == 1 & df$head_vector_angular_speed > 0)
  HC_right = which(HC_flag == 1 & df$head_vector_angular_speed < 0)
  
  HCs_right <- rep(0, nrow(df))
  HCs_right[HC_right] <- 1
  
  HCs_left <- rep(0, nrow(df))
  HCs_left[HC_left] <- 1
  
  df <-
    data.frame(HCs,
               HC_flag,
               HCs_right,
               HCs_left,
               HC_angle,
               Abs_HC_angle,
               HC_reorientation)
  return(df)
}



#' Title
#'
#' @param df 
#' @param frame_rate 
#'
#' @return
#' @export
#' @import stats
#' @examples
calculate_run_flag <- function(df, frame_rate = 16) {

  HC_flag <- df$HCs_right | df$HCs_left
  
  if (length(HC_flag) > 2 * 1.5 * frame_rate) {
    HC_flag <-
      stats::filter(HC_flag, rep(1, 2 * 1.5 * frame_rate))#df$HCs_left & df$HCs_right
  }
  HC_flag <- as.logical(HC_flag)
  ! HC_flag
  
}



calculate_step_variables <- function(df,
                                     frame_rate = 16,
                                     threshold_tail_speed_foward = 0.6,
                                     
                                     threshold_tail_speed_backward = -0.4) {
  # ## from Emmanouil Paisios et. al 2017 ###
  
  npoints = round(0.3 * frame_rate)
  
#' Title
#'
#' @param df 
#' @param threshold_tail_speed_foward 
#' @param threshold_tail_speed_backward 
#' @param kind 
#'
#' @return
#' @export
#'
  calculate_steps <- function(df,
                              threshold_tail_speed_foward = 0.6 ,
                              threshold_tail_speed_backward = -0.4,
                              kind = "forward") {
    
    step_forward <- diff(sign(diff(df$tail_vel_forward))) == -2
    step_forward[which(df$tail_vel_forward < threshold_tail_speed_foward)] <- FALSE
    step_forward[is.na(step_forward)] <- FALSE
    
    step_backward <- rollapply(df$tail_vel_forward, 5, function(x) which.min(x) == 2)
    step_backward[which(df$tail_vel_forward > threshold_tail_speed_backward)] <-FALSE
    step_backward[is.na(step_backward)] <- FALSE
    
    
    step_boolean = rep(F,nrow(df))
    step_boolean[which(step_forward)] <- T
    step_boolean[which(step_backward)] <- T
    
    step_start_idx = list()
    step_end_idx = list()
    run_direction <- as.numeric(df$run_flag)
    
    if (length(which(step_boolean)) > 1) {
      run_idx = which(df$run_flag)
      step_idx_all <- which(step_boolean)
      for (i in seq(1, length(step_idx_all) - 1)) {
        start <- step_idx_all[i]
        end <- step_idx_all[i + 1]
        
        if(step_forward[start]){
          run_direction[start:end]=1
        }
        
        if(step_backward[start]){
          run_direction[start:end]=-1
        }
        
        if ((all(c(start,end) %in%  run_idx))) {
          step_start_idx[[i]] <- start
          step_end_idx[[i]] <- end
        }
      }
    }
    step_boolean = rep(F, nrow(df))
    step_boolean[unlist(step_end_idx)] = T
    step_boolean <- as.logical(step_boolean)
    
    run_direction[which(!df$run_flag)] <-0
    run_direction[is.na(run_direction)] <- 0
    
    return(
      list(
        step_boolean = step_boolean,
        step_forward = step_forward,
        step_backward = step_backward,
        step_start_idx = unlist(step_start_idx),
        step_end_idx = unlist(step_end_idx),
        run_direction = run_direction
      )
    )
  }
  
  steps <- calculate_steps(df,
                           threshold_tail_speed_foward,
                           threshold_tail_speed_backward)
  step_boolean <- steps$step_boolean
  step_start_idx <- steps$step_start_idx
  step_end_idx <- steps$step_end_idx
  run_direction <- steps$run_direction
  
  if (length(step_start_idx) > 2) {
    IS_angle <- interval_angle(
      df,
      step_start_idx,
      step_end_idx,
      "tail_vector_x",
      "tail_vector_y",
      step_start_idx
    )
    Abs_IS_angle <- abs(IS_angle)
    
    
    IS_reorientation <-
      interval_difference(df,
                          step_start_idx,
                          step_end_idx,
                          "heading_angle",
                          step_start_idx,
                          function(x, y)
                            abs(x) - abs(y)
                          
      )
    
    IS_interval <-
      interval_difference(df,
                          step_start_idx,
                          step_end_idx,
                          "frame",
                          step_start_idx,
                          function(x, y)
                            (y - x) / frame_rate
      )
    
    IS_distance <-
      interval_difference(df,
                          step_start_idx,
                          step_end_idx,
                          "midpoint_distance",
                          step_start_idx,
                          function(x, y)
                            y - x)
    
    run_speed <-
      interval_mean(df,
                    step_start_idx,
                    step_end_idx,
                    "midpoint_speed",
                    step_start_idx)
    
    extr <- interval_extr(df,
                          step_start_idx,
                          step_end_idx,
                          step_start_idx)
    
    step_idx <- extr$step_index
    step_extr <- extr$step_extr
    gt_1min <- extr$gt_1min
    minima <- extr$minima
    
  } else{
    step_boolean = rep(F, nrow(df))
    IS_angle = rep(NA, nrow(df))
    Abs_IS_angle = rep(NA, nrow(df))
    IS_reorientation = rep(NA, nrow(df))
    IS_interval = rep(NA, nrow(df))
    IS_distance = rep(NA, nrow(df))
    run_speed = rep(NA, nrow(df))
    step_idx = rep(NA, nrow(df))
    step_extr = rep(NA, nrow(df))
    gt_1min = rep(NA, nrow(df))
    minima = rep(NA, nrow(df))
    run_direction= rep(NA,nrow(df))
  }
  return(
    data.frame(
      run_direction,
      step_boolean,
      IS_angle,
      Abs_IS_angle,
      IS_reorientation,
      IS_interval,
      IS_distance,
      run_speed,
      step_idx,
      step_extr,
      gt_1min,
      minima
    )
  )
  
}




#' Title
#'
#' @param track 
#' @param frame_rate 
#' @param seconds 
#'
#' @return
#' @export
#' @import RcppRolll


calculate_switching_point <- function(track,frame_rate=16,seconds=10){
  print(track)
  if(nrow(track)<frame_rate*seconds*2){
    rep(NA,nrow(track))
  }else{
    c(NA,diff(sign(roll_mean(track$tail_vel_forward,frame_rate*seconds,fill=NA))))==2
  }
}


#' Title
#'
#' @param track 
#' @param frame_rate 
#' @param seconds_filter 
#' @param seconds_interval 
#'
#' @return
#' @import RcppRoll
#' @export
#'
calculate_flanks <- function(track,
                             frame_rate=16,
                             seconds_filter=5,
                             seconds_interval=10){
  
  
  tail_vel <- track$tail_vel_forward_bl
  tail_vel <- roll_median(tail_vel,frame_rate*seconds_filter,fill=NA)
  idx_zero <- tail_vel<0.03 &tail_vel>-0.03
  tail_vel[idx_zero] <- 0
  tail_vel <- roll_mean(tail_vel,frame_rate*seconds_filter,fill=NA)
  sig <- sign(tail_vel)
  d <-  c(NA,diff(sign(sig)))
  print(d)
  edges = which(d!=0 & !is.na(d))
  to_neg <- c()
  to_pos <- c()
  if(length(edges>0)){
    for(i in seq(1,length(edges))){
      edge_idx <- edges[i]
      before_idx = edge_idx-frame_rate*seconds_interval
      after_idx = edge_idx+frame_rate*seconds_interval
      if(before_idx <0){
        before_idx = 0
      }
      before <- sig[seq(before_idx,edge_idx-1)]
      after <- sig[seq(edge_idx+1,after_idx)]
      before[is.na(before)]<- 0
      after[is.na(after)]<- 0
      offset = 0#(frame_rate*seconds_filter)/2
      
      #detect positive to negative or zero to negative
      if((d[edge_idx]==-2 |d[edge_idx]==-1)& 
         all(before==1|before==0)&all(after==-1)){
        to_neg <- c(to_neg,edge_idx-offset)
      }
      
      #detect negative or zero to positive
      if((d[edge_idx]==2 |d[edge_idx]==1)&
        all(before==-1|before==0)&all(after==1)){
        to_pos <- c(to_pos,edge_idx-offset)
      }
      
    }
  }

  flanks <- rep(0,nrow(track))
  if(!is.null(to_pos)){
    flanks[to_pos] <- 1
  }
  if(!is.null(to_neg)){
    flanks[to_neg] <- -1
  }
  
  flanks
}



#' Title
#'
#' @param track 
#' @param odor_pos 
#' @param frame_rate 
#' @param threshold_head_vector_angular_speed 
#' @param threshold_tail_vector_angular_speed 
#' @param min_hc_size 
#' @param max_hc_size 
#' @param spines_pairs 
#' @param contours_pairs 
#' @param drop_contours 
#' @param drop_collisions 
#' @param drop_curls 
#'
#' @return
#' @import dplyr
#' @import data.table
#' @export
#'
process_track <- function(track,
                          odor_pos,
                          frame_rate = 16,
                          threshold_head_vector_angular_speed = 35,
                          threshold_tail_vector_angular_speed = 45,
                          min_hc_size = 0,
                          max_hc_size = Inf,
                          spines_pairs = paste(paste("spinepoint_x", seq(1, 12), sep = "_"),
                                               paste("spinepoint_y", seq(1, 12), sep = "_")),
                          contours_pairs = paste(paste("contourpoint_x", seq(1, 22), sep = "_"),
                                                 paste("contourpoint_y", seq(1, 22), sep = "_")),
                          drop_contours = FALSE,
                          drop_collisions=TRUE,
                          drop_curls = TRUE
){
  spines <-
    unlist(str_split(spines_pairs, " "))
  
  contours <-
    unlist(str_split(contours_pairs, " "))


  
  
  df <- track %>%
    mutate(spine_length = calculate_spine_length_sum(track))
  if (!drop_contours) {
    df <-  df %>%
      cbind(calculate_rotate_columns(., spines_pairs, odor_pos)) %>%
      select(!c("odor_a_x_rot", "odor_a_y_rot")) %>%
      cbind(
        calculate_convolve_columns(
          .,
          in_columns = paste(spines, "rot", sep = "_"),
          column_names = paste(spines, "conv", sep = "_"),
          frame_rate = frame_rate
        )
      ) %>%
      cbind(calculate_rotate_columns(., contours_pairs, odor_pos)) %>%
      cbind(
        calculate_convolve_columns(
          .,
          in_columns = paste(contours, "rot", sep = "_"),
          column_names = paste(contours, "conv", sep = "_"),
          frame_rate = frame_rate
        )
      )%>%
      #remove original spines and rotated columns
      select(!spines)%>%
      select(!paste(spines, "rot", sep = "_"))%>%
      select(!contours)%>%
      select(!paste(contours, "rot", sep = "_"))
  } else{
    df <- df %>%
      cbind(calculate_rotate_columns(., spines_pairs, odor_pos)) %>%
      cbind(
        calculate_convolve_columns(
          .,
          in_columns = paste(spines, "rot", sep = "_"),
          column_names = paste(spines, "conv", sep = "_"),
          frame_rate = frame_rate
        )
      ) %>%     #remove original spines and rotated columns
      select(!spines)%>%
      select(!paste(spines, "rot", sep = "_"))
  }
  
  df <- df %>%
    cbind(calculate_head_tail_vectors(.)) %>%
    mutate(midpoint_speed = difference(spinepoint_x_6_conv, spinepoint_y_6_conv) *
             frame_rate) %>%
    mutate(midpoint_speed_bl = midpoint_speed / mean(spine_length, na.rm =
                                                       T)) %>%
    mutate(midpoint_distance = cumsum(coalesce((
      midpoint_speed / frame_rate
    ), 0)) + midpoint_speed / frame_rate * 0) %>%
    mutate(midpoint_distance_bl = midpoint_distance / mean(spine_length, na.rm =
                                                             T)) %>%
    mutate(distance_to_sp =
             sqrt((spinepoint_x_6_conv - spinepoint_x_6_conv[5]) ^ 2 +
                    (spinepoint_y_6_conv - spinepoint_y_6_conv[5]) ^ 2
             )) %>%
    mutate(
      head_vector_angular_speed = calculate_head_vector_angular_speed(., frame_rate),
      tail_vector_angular_speed = calculate_tail_vector_angular_speed(., frame_rate)
    ) %>%
    mutate(
      abs_head_vector_angular_speed =abs(head_vector_angular_speed),
      abs_tail_vector_angular_speed =abs(tail_vector_angular_speed)
    )%>%
    mutate(
      head_vel_forward = calculate_head_vel_forward(., frame_rate),
      tail_vel_forward = calculate_tail_vel_forward(., frame_rate)
    ) %>%
    mutate(
      head_vel_forward_bl = head_vel_forward / mean(spine_length, na.rm = T),
      tail_vel_forward_bl = tail_vel_forward / mean(spine_length, na.rm = T)
      
    ) %>%
    mutate(na_flag_flip = calculate_na_flag_ht_flip(.)) %>%
    mutate(
      head_vector_extended_x = spinepoint_x_11_conv - spinepoint_x_5_conv,
      head_vector_extended_y = spinepoint_y_11_conv - spinepoint_y_5_conv
    ) %>%
    mutate(
      bending_angle = calculate_bending_angle(
        tail_vector_x,
        tail_vector_y,
        head_vector_extended_x,
        head_vector_extended_y,
        frame_rate
      )
    ) %>%
    mutate(abs_bending_angle = abs(bending_angle))%>%
    mutate(bearing_angle = calculate_bearing_angle(.),
           heading_angle = calculate_heading_angle(.)) %>%
    mutate(abs_bearing_angle = abs(bearing_angle),
           abs_heading_angle = abs(heading_angle)) %>%
    mutate(
      distance_to_odor =
        sqrt((spinepoint_x_6_conv - odor_a_x_rot) ^ 2 +
               (spinepoint_y_6_conv - odor_a_y_rot) ^ 2
        ),
      odor_orientation = ifelse(bearing_angle <= 90, "towards", "away"),
      odor_speed = calculate_speed_direction_odor(., frame_rate),
      pref = ifelse(spinepoint_y_5_conv > 0, 1,-1)
    ) 
  
  # filter out NA frames
  
  if(drop_collisions){
    df[which(df$collision_flag==1),c(paste("spinepoint_x", seq(1, 12), sep = "_"),
                                           paste("spinepoint_y", seq(1, 12), sep = "_"),
                                           paste("contourpoint_x", seq(1, 22), sep = "_"),
                                           paste("contourpoint_y", seq(1, 22), sep = "_"),
                                     
                                     "head_vector_x","head_vector_y", "tail_vector_x", "tail_vector_y",
                                     "midpoint_speed","midpoint_speed_bl","midpoint_distance","midpoint_distance_bl",
                                     "distance_from_sp","head_vector_angular_speed","abs_head_vector_angular_speed",
                                     "tail_vector_angular_speed","abs_tail_vector_angular_speed",
                                     "head_vel_forward", "head_vel_forward_bl","tail_vel_forward","tail_vel_forward_bl",
                                     "head_vector_extended_x","head_vector_extended_y",
                                     "bending_angle","bearing_angle","heading_angle","distance_to_odor","odor_orientation","odor_speed",
                                     "pref"
                                     )] <- NA
  }
  if(drop_curls){
    df[which(df$collision_flag==2),c(paste("spinepoint_x", seq(1, 12), sep = "_"),
                                     paste("spinepoint_y", seq(1, 12), sep = "_"),
                                     paste("contourpoint_x", seq(1, 22), sep = "_"),
                                     paste("contourpoint_y", seq(1, 22), sep = "_"),
                                     
                                     "head_vector_x","head_vector_y", "tail_vector_x", "tail_vector_y",
                                     "midpoint_speed","midpoint_speed_bl","midpoint_distance","midpoint_distance_bl",
                                     "distance_from_sp","head_vector_angular_speed","abs_head_vector_angular_speed",
                                     "tail_vector_angular_speed","abs_tail_vector_angular_speed",
                                     "head_vel_forward", "head_vel_forward_bl","tail_vel_forward","tail_vel_forward_bl",
                                     "head_vector_extended_x","head_vector_extended_y",
                                     "bending_angle","bearing_angle","heading_angle","distance_to_odor","odor_orientation","odor_speed",
                                     "pref")] <- NA
  }
    df <- df %>%
    cbind(calculate_HC_variables(., 
                                 frame_rate = frame_rate,
                                 threshold_head_vector_angular_speed = threshold_head_vector_angular_speed,
                                 threshold_tail_vector_angular_speed = threshold_tail_vector_angular_speed,
                                 min_hc_size = min_hc_size,
                                 max_hc_size = max_hc_size)) %>%
    mutate(run_flag = calculate_run_flag(., frame_rate)) %>%
    cbind(calculate_step_variables(., frame_rate))%>%
    mutate(
      IS_distance_bl = IS_distance / mean(spine_length, na.rm =T),
      run_speed_bl = run_speed / mean(spine_length, na.rm =T)
    )%>%
    #mutate(back_to_for = calculate_switching_point(.,frame_rate,30))
    mutate(flanks = calculate_flanks(.,frame_rate,5,5))
}
