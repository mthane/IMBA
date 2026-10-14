
##### HELP FUNCTIONS #####

distance2 <- function(x1, x2)
  sqrt(sum((x1 - x2) ^ 2))

rad2deg <- function(rad) {
  (rad * 180) / (pi)
}


deg2rad <- function(deg) {
  (deg / 180) * (pi)
}


#' Title
#'
#' @param n 
#'
#' @return
#' @export
#'
#' @import shiny
#' @examples
linebreaks <- function(n) {
  HTML(strrep(br(), n))
}
#' Title
#'
#' @param a vector
#' @param b vector
#'
#' @description Calculates clockwise angle between two vectors.
#' Positive means clockwise change and negative counter clockwise.
#' @return clockwise angle between vector a and b
#'
#'
clockwise_angle <- function(a, b) {
  angle <- atan2(b[2], b[1]) - atan2(a[2], a[1])
  return(-rad2deg((angle - (abs(
    angle
  ) > pi) * sign(angle) * 2 * pi)))
}


#' Title
#'
#' @param a1 
#' @param a2 
#' @param b1 
#' @param b2 
#'
#' @return
#' @export
#' 
clockwise_angle_v <- function(a1, a2, b1, b2) {
  angle <- atan2(b2, b1) - atan2(a2, a1)
  return(-rad2deg((angle - (abs(
    angle
  ) > pi) * sign(angle) * 2 * pi)))
}


#' Title
#'
#' @param vector 
#' @param angle 
#'
#' @return
#' @export
#'
rotate_vector_clockwise <- function(vector, angle) {
  c(cos(angle) * vector[1] + sin(angle) * vector[2],-sin(angle) * vector[1] + cos(angle) * vector[2])
}


#' Title
#'
#' @param column_x 
#' @param column_y 
#'
#' @return
#' @export
#'
#' @examples
difference <- function(column_x, column_y) {
  column_x <- as.numeric(unlist(column_x))
  column_y <- as.numeric(unlist(column_y))
  vx <- as.numeric(c(NA, diff(column_x)))
  vy <- as.numeric(c(NA, diff(column_y)))
  return(sqrt(vx ^ 2 + vy ^ 2))
}

#' Title
#'
#' @param x 
#' @param n 
#'
#' @return
#' @export
#' @import stats
#' @examples
ma <- function(x, n = 3) {
  if (length(x) < n) {
    x
  } else{
    stats::filter(x, rep(1 / n, n), sides = 2)
  }
}

#' Title
#'
#' @param df 
#' @param start_idx 
#' @param end_idx 
#' @param measurement 
#' @param centers 
#' @param FUN 
#'
#' @return
#' @export
#'
#' @examples
interval_difference <-
  function(df,
           start_idx,
           end_idx,
           measurement,
           centers,
           FUN
           ) {
    values <-
      FUN(df[start_idx, measurement], df[end_idx, measurement])
    #if HC in between
    result <- rep(NA, nrow(df))
    if (length(values) > 0){
      result[unlist(centers)] <- values
    }
    result
  }


#' Title
#'
#' @param df 
#' @param start_idx 
#' @param end_idx 
#' @param vector_x 
#' @param vector_y 
#' @param centers 
#'
#' @return
#' @export
#'
interval_angle <-
  function(df,
           start_idx,
           end_idx,
           vector_x,
           vector_y,
           centers) {
    values <-
      clockwise_angle_v(
                        df[end_idx, vector_x], df[end_idx,vector_y],
                        df[start_idx, vector_x], df[start_idx,vector_y]
                        )
    result <- rep(NA, nrow(df))
    if (length(values) > 0) {
      result[unlist(centers)] <- values
    }
    result
  }

#' Title
#'
#' @param df 
#' @param start_idx 
#' @param end_idx 
#' @param measurement 
#' @param centers 
#'
#' @return
#' @export
#'
interval_mean <-
  function(df,
           start_idx,
           end_idx,
           measurement,
           centers) {
    values <- list()
    steps <- data.frame(start_idx, end_idx)
    if (is.null(start_idx) | is.null(end_idx)) {
      values <- NA
    } else{
      #iterate over the steps
      for (idx in  1:nrow(steps)) {
        step_measure <-
          df[steps[idx,]$start_idx:steps[idx,]$end_idx, measurement]
        values[[idx]] <- mean(step_measure, na.rm = T)
      }
    }
    values <- unlist(values)
    result <- rep(NA, nrow(df))
    if (length(values) > 0) {
      result[unlist(centers)] <- values
    }
    result
  }

interval_max <- #new after Thane et al. 2023
  function(df,
           start_idx,
           end_idx,
           measurement,
           centers) {
    values <- list()
    steps <- data.frame(start_idx, end_idx)
    if (is.null(start_idx) | is.null(end_idx)) {
      values <- NA
    } else{
      #iterate over the steps
      for (idx in  1:nrow(steps)) {
        step_measure <-
          df[steps[idx,]$start_idx:steps[idx,]$end_idx, measurement]
        values[[idx]] <- max(step_measure, na.rm = T)
      }
    }
    values <- unlist(values)
    result <- rep(NA, nrow(df))
    if (length(values) > 0) {
      result[unlist(centers)] <- values
    }
    result
  }

interval_maxindex <- #new after Thane et al. 2023
  function(df,
           start_idx,
           end_idx,
           measurement,
           centers) {
    values <- list()
    steps <- data.frame(start_idx, end_idx)
    if (is.null(start_idx) | is.null(end_idx)) {
      values <- NA
    } else{
      #iterate over the steps
      for (idx in  1:nrow(steps)) {
        step_measure <-
          df[steps[idx,]$start_idx:steps[idx,]$end_idx, measurement]
        values[[idx]] <- which.max(step_measure)
      }
    }
    values <- unlist(values)
    result <- rep(NA, nrow(df))
    if (length(values) > 0) {
      result[unlist(centers)] <- values
    }
    result
  }

#' Title
#'
#' @param df 
#' @param start_idx 
#' @param end_idx 
#' @param centers 
#' @param window 
#'
#' @return
#' @export
#' @import zoo
#' 
interval_extr <-
  function(df,
           start_idx,
           end_idx,
           centers,
           window = 5) {
    values <- list()
    
    steps <- data.frame(start_idx, end_idx)
    
    ## boolean for the "minima" > 1 minimum
    gt_1min <- rep(NA, nrow(df))
    minima <- rep(NA, nrow(df))
    step_extr <- rep(NA, nrow(df))
    step_index <- rep(NA, nrow(df))
    total = 0
    if (is.null(start_idx) | is.null(end_idx)) {
      values <- NA
    } else{
      #iterate over the steps
      for (idx in  1:nrow(steps)) {
        step_index[steps[idx,]$start_idx:steps[idx,]$end_idx] <- idx
        step_measure <-
          df[steps[idx,]$start_idx:steps[idx,]$end_idx, "tail_vel_forward"]
        
        nmin = NA
        if (length(step_measure) > 5) {
          x = ma(diff(step_measure)) #ma: moving average
          xz <- as.zoo(x) #zoo: ordered objects
          extr <-
            c(NA, rollapply(xz, window, function(x)
              which.min(x) == 2), NA)
          nmin <- sum(extr, na.rm = T)
          values[[idx]] <- nmin
        } else{
          xz <- as.zoo(step_measure)
          if(length(step_measure)>2){
            extr <-
              c(NA, rollapply(xz, 3, function(x)
                which.min(x) == 2), NA)
            
            nmin <- sum(extr, na.rm = T)
            values[[idx]] <- nmin
          }else{
            values[[idx]] <- NA
          }
        }
        total = total + length(step_measure)
        if (!is.na(nmin)) {
          if (nmin == 0) {
            minima[steps[idx,]$start_idx:steps[idx,]$end_idx] <- "0"
          }
          if (nmin == 1) {
            minima[steps[idx,]$start_idx:steps[idx,]$end_idx] <- "1"
            gt_1min[steps[idx,]$start_idx:steps[idx,]$end_idx] <-
              F
            
          } else{
            gt_1min[steps[idx,]$start_idx:steps[idx,]$end_idx] <- T
          }
          if (nmin == 2) {
            minima[steps[idx,]$start_idx:steps[idx,]$end_idx] <-  "2"
          }
          if (nmin > 2) {
            minima[steps[idx,]$start_idx:steps[idx,]$end_idx] <-  ">2"
          }
        } else{
          minima[steps[idx,]$start_idx:steps[idx,]$end_idx] <- NA
          
          gt_1min[steps[idx,]$start_idx:steps[idx,]$end_idx] <- NA
        }
        
        
      }
    }
    
    values <- unlist(values)
    gt_1min <- unlist(gt_1min)
    minima <- unlist(minima)
    if (length(values) > 0) {
      step_extr[unlist(centers)] <- values
    }
    
    data.frame(step_index, step_extr, gt_1min, minima)
  }

head_extr <- #new after Thane et al. 2023
  function(df,
           start_idx,
           end_idx,
           centers,
           window = 5) {
    values <- list()
    
    steps <- data.frame(start_idx, end_idx)
    
    ## boolean for the "minima" > 1 minimum
    gh_3max <- rep(NA, nrow(df))
    hmaxima <- rep(NA, nrow(df))
    step_hmax <- rep(NA, nrow(df))
    step_index <- rep(NA, nrow(df))
    total = 0
    if (is.null(start_idx) | is.null(end_idx)) {
      values <- NA
    } else{
      #iterate over the steps
      for (idx in  1:nrow(steps)) {
        step_index[steps[idx,]$start_idx:steps[idx,]$end_idx] <- idx
        step_measure <-
          df[steps[idx,]$start_idx:steps[idx,]$end_idx, "head_vel_forward"]
        
        nmax = NA
        if (length(step_measure) > 5) {
          x = ma(diff(step_measure)) #ma: moving average
          xz <- as.zoo(x) #zoo: ordered objects
          extr <-
            c(NA, rollapply(xz, window, function(x)
              which.max(x) == 2), NA)  #list with all frames of IS with local maxima = T
          nmax <- sum(extr, na.rm = T) #number of local maxima in IS
          values[[idx]] <- nmax
        } else{
          xz <- as.zoo(step_measure)
          if(length(step_measure)>2){
            extr <-
              c(NA, rollapply(xz, 3, function(x)
                which.max(x) == 2), NA)
            
            nmax <- sum(extr, na.rm = T)
            values[[idx]] <- nmax
          }else{
            values[[idx]] <- NA
          }
        }
        total = total + length(step_measure)
        if (!is.na(nmax)) {
            if (nmax == 0 | nmax == 1 | nmax == 2 | nmax == 3) {
            hmaxima[steps[idx,]$start_idx:steps[idx,]$end_idx] <- nmax
            gh_3max[steps[idx,]$start_idx:steps[idx,]$end_idx] <-
              F
            
          } else{
            hmaxima[steps[idx,]$start_idx:steps[idx,]$end_idx] <-  ">3"
            gh_3max[steps[idx,]$start_idx:steps[idx,]$end_idx] <- T
          }
        } else{
          hmaxima[steps[idx,]$start_idx:steps[idx,]$end_idx] <- NA
          gh_3max[steps[idx,]$start_idx:steps[idx,]$end_idx] <- NA
        }
        
        
      }
    }
    
    values <- unlist(values)
    gh_3max <- unlist(gh_3max)
    hmaxima <- unlist(hmaxima)
    if (length(values) > 0) {
      step_hmax[unlist(centers)] <- values
    }
    
    data.frame(step_hmax, gh_3max, hmaxima)
  }

#' Title
#'
#' @param df 
#' @param variables 
#'
#' @return
#' @export
#' @import data.table
#' @import stats
#' 
uTests <- function(df,variables) {

  combs <-
    data.frame(combinations(length(unique(
      df$group_condition
    )), 2, unique(df$group_condition)))
  comb_dfs <-list()
  i=1
  
  print(combs)
  for (row in 1:nrow(combs)){
    name <- paste(combs$X1[row],combs$X2[row],sep="-vs-")
    dft <- df%>%
      filter(group_condition==combs$X1[row] | group_condition==combs$X2[row])
    pvalues = list()
    j=1
    for (v in variables) {
      wilcox <-wilcox.test(as.formula(paste0(v, "~ group_condition")), data = dft)
      pvalues[j] <- wilcox$p.value
      j= j +1
    }
    df_i <- data.frame(rep(name,length(variables)),variables,as.numeric(unlist(pvalues)))
    
    colnames(df_i) <- c("names","variable","pvalue")
    comb_dfs[[i]] <- df_i%>% arrange(pvalue)
    i=i+1
  }
  rbindlist(comb_dfs)
}





#' Title
#'
#' @param df 
#' @param variable 
#'
#' @return
#' @export
#' @import gtools
#' 
pvalues <- function(df, variable) {
  #df$group_condition <- paste(df$group, df$condition, sep = "-")
  combs <-
    data.frame(combinations(length(unique(
      df$group_condition
    )), 2, unique(df$group_condition), set = FALSE))
  
  pvalues = list()
  for (row in 1:nrow(combs)) {
    data <- df %>% filter(group_condition == combs[row, "X1"] |
                            group_condition == combs[row, "X2"])
    
    #  grouping factor must have exactly 2 levels
    if (length(unique(data$group_condition)) == 2) {
      wilcox <-
        wilcox.test(as.formula(paste0(variable, "~ group_condition")), data = data)
      pvalues[[row]] <- wilcox$p.value
    } else{
      pvalues[[row]] <- NA
    }
    
  }
  cbind(combs, data.frame(pvalues = unlist(pvalues)))
}


#' Title
#'
#' @param containerId 
#' @param expr 
#'
#' @return
#' @export
#'
withConsoleRedirect <- function(containerId, expr) {
  ##### https://gist.github.com/jcheng5/3830244757f8ca25d4b00ce389ea41b3 #####
  # Change type="output" to type="message" to catch stderr
  # (messages, warnings, and errors) instead of stdout.
  txt <- capture.output(results <- expr, type = "message")
  if (length(txt) > 0) {
    insertUI(
      paste0("#", containerId),
      where = "beforeEnd",
      ui = paste0(txt, "\n", collapse = "")
    )
  }
  results
}
