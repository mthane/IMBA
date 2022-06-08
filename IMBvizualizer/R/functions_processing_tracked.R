#' Title
#'
#' @param input_dir the directory containing the tracked data
#' @param output_dir where the analysis should be saved
#' @param save_as the folder name under which the analysis should be saved
#' @param radius the radius of the petri dish in mm
#' @param drop_contours if contours should be dropped (contours will be useful for the animation of larvae)
#' @param frame_rate the frame rate of the tracked videos
#'
#' @return data.frame with the analyzed data set
#' @export
#' @import tictoc
process_data_dt <- function(input_dir,
                            output_dir,
                            save_as,
                            radius,
                            drop_contours,
                            frame_rate,
                            video_length,
                            threshold_head_vector_angular_speed = 35,
                            threshold_tail_vector_angular_speed = 45,
                            min_hc_size = 0,
                            max_hc_size = Inf
) {
  tic()
  message("begin...")
  start <- Sys.time() # get start time
  
  radius <- as.numeric(radius)
  
  data <- create_dataset_from_dir(input_dir,
                                  threshold = 30,
                                  drop_contours = drop_contours)
  
  message("finished preprocessing")
  data <- analyze_data(data,
                       output_dir,
                       save_as,
                       radius,
                       drop_contours,
                       frame_rate,
                       video_length,
                       threshold_head_vector_angular_speed,
                       threshold_tail_vector_angular_speed,
                       min_hc_size,
                       max_hc_size
  )
  
  it <- toc()
  exectime <- it$toc - it$tic
  message(paste0(exectime , "seconds elapsed"))
  return(data)
}

##### CREATE DATASETS #####

#' Title
#'
#' @param input_dir 
#' @param threshold 
#' @param drop_contours 
#'
#' @return
#' @export
#' @import dplyr
#' @import data.table
#' @import tictoc
#'
create_dataset_from_dir <-
  function(input_dir,
           threshold = 30,
           drop_contours = F) {
    #has to be the experiment folder
    track_dirs <- dir(input_dir, "[[:digit:]].csv$", recursive = T)

    spines <-
      unlist(str_split(paste(
        paste("spinepoint_x", seq(1, 12), sep = "_"),
        paste("spinepoint_y", seq(1, 12), sep = "_")
      ), " "))
    contours <-
      unlist(str_split(paste(
        paste("contourpoint_x", seq(1, 22), sep = "_"),
        paste("contourpoint_y", seq(1, 22), sep = "_")
      ), " "))
    
    column_names <- c(
      "frame",
      spines,
      "centoid_x",
      "centroid_y",
      "blob_orient",
      "area",
      "grey",
      "spine_length",
      "width",
      "perimeter",
      "collision_flag"
    )
    column_names_contours <- c(
      "frame",
      spines,
      contours,
      "centoid_x",
      "centroid_y",
      "blob_orient",
      "area",
      "grey",
      "spine_length",
      "width",
      "perimeter",
      "collision_flag"
    )
    
    
    tracks <- list()
    tcount <- 1
    for (track_dir in track_dirs) {
      percent = round(tcount / length(track_dirs) * 100)
      message(
        paste(
          "PREPROCESSING   ",
          "track: ",
          track_dir,
          "tracks:",
          length(track_dirs),
          "progress:",
          percent,
          "%",
          "time: "
        )
      )
      
      
      spl <-
        unlist(str_split(paste(input_dir, track_dir, sep = "/"), "/"))
      experiment <- tail(spl, 5)
      group <- tail(spl, 4)[1]
      condition <- tail(spl, 3)[1]
      trial <- tail(spl, 2)[1]
      track <- tail(spl, 1)[1]
      
      base_dir <-
        paste(spl[seq(1, length(spl) - 4)], collapse = "/")
      
      meta_dir <-
        paste(c(
          base_dir,
          group,
          condition,
          trial,
          "vidAndLogs/metadata.txt"
        ),
        collapse = "/")
      
      if ((file.info(meta_dir)$size == 0)) {
        odor_location = NA
      } else{
        metadata <- read.delim(meta_dir)
        odor_location = unlist(str_split(metadata[6, ], "="))[2]
      }
      
      if (is.na(odor_location)) {
        odor_xy = c(0, 0)
      }
      else{
        odor_xy = as.numeric(unlist(str_split(odor_location, ",")))
      }
      
      if (drop_contours) {
        track <-
          as.data.frame(fread(
            paste(input_dir, track_dir, sep = "/"),
            drop = 26:69,
            header = F
          ))
        
        colnames(track) <- column_names
      } else{
        track <-
          as.data.frame(fread(paste(input_dir, track_dir, sep = "/"), header = F))
        colnames(track) <- column_names_contours
        track[contours] <- sapply(track[contours], as.numeric)
      }
      id <- rep(tcount, nrow(track))
      name <- rep(tail(spl, 1), nrow(track))
      trial <- rep(tail(spl, 2)[1], nrow(track))
      condition <- rep(tail(spl, 3)[1], nrow(track))
      group <- rep(tail(spl, 4)[1], nrow(track))
      odor_a_x <- rep(odor_xy[1], nrow(track))
      odor_a_y <- rep(odor_xy[2], nrow(track))
      track <-
        cbind(track, id, name, trial, condition, group, odor_a_x, odor_a_y)
      track[spines] <- sapply(track[spines], as.numeric)
      
      
      
      if (nrow(track > threshold)) {
        tracks[[tcount]] = track
        
      }
      tcount = tcount + 1
    }
    data <- bind_rows(tracks)
    
    #fwrite(data, paste(dir, paste0(save_as, "_raw_data.csv"), sep ="/"))
    return(data)
  }



#' Title
#'
#' @param data the data.frame which is the output of the function create_dataset_from_dir
#' @param output_dir the directory where the analysis should be saved
#' @param save_as the name under which the analysis csv-data should be saved
#' @param radius the radius of the petri dish
#' @param drop_contours if contours should be dropped (contours will be useful for the animation of larvae)
#' @param frame_rate the frame rate of the tracked videos
#'
#' @return
#' @import dplyr
#' @import data.table
#' @export
#'

analyze_data <-
  function(data,
           output_dir = ".",
           save_as = "analyzed",
           radius = 42.5,
           drop_contours = F,
           frame_rate = 16,
           video_length=180,
           threshold_head_vector_angular_speed = 35,
           threshold_tail_vector_angular_speed = 45,
           min_hc_size = 0,
           max_hc_size = Inf
           ) {
    tic()
    data <-  as.data.table(data)
    track_ids = unique(data$id)
    
    
    spines <-
      unlist(str_split(paste(
        paste("spinepoint_x", seq(1, 12), sep = "_"),
        paste("spinepoint_y", seq(1, 12), sep = "_")
      ), " "))
    
    spines_c <-
      unlist(str_split(paste(
        paste("spinepoint_x", seq(1, 12), "conv", sep = "_"),
        paste("spinepoint_y", seq(1, 12), "conv", sep = "_")
      ), " "))
    contours <-
      unlist(str_split(paste(
        paste("contourpoint_x", seq(1, 22), sep = "_"),
        paste("contourpoint_y", seq(1, 22), sep = "_")
      ), " "))
    
    tcount = 1
    tracks = list()
    
    ntracks = length(track_ids)
    time_per_group = 0
    for (id_ in track_ids) {
      start_interval <- Sys.time()
      
      track <- as.data.frame(data[id == id_])
      if (nrow(track) > round(frame_rate / 0.5)) {
        percent = round(tcount / length(track_ids) * 100)
        message(
          paste(
            "ANALYSIS ",
            "track: ",
            id_,
            "name: ",
            unique(track$name),
            "trial: ",
            unique(track$trial),
            "condition: ",
            unique(track$condition),
            "group: ",
            unique(track$group),
            "tracks:",
            ntracks,
            "progress:",
            percent,
            "%",
            sep = " "
          )
        )
        odor_pos = c(unique(track$odor_a_x), unique(track$odor_a_y))
        track <- process_track(track, odor_pos,
                               drop_contours = drop_contours,
                               frame_rate= frame_rate,
                               threshold_head_vector_angular_speed = threshold_head_vector_angular_speed,
                               threshold_tail_vector_angular_speed = threshold_tail_vector_angular_speed,
                               min_hc_size = min_hc_size,
                               max_hc_size = max_hc_size
                               )
        tracks[[tcount]] <- as.data.table(track)
        end_interval <- Sys.time() - start_interval
        if (tcount %% 200 == 0) {
          gc()
        }
        
      }
      tcount = tcount + 1
      
    }
    
    message("bind rows... this can take some minutes")
    data = rbindlist(tracks, fill = TRUE)
    message(paste0("write data... this can take some minutes"))
    fwrite(data, paste(output_dir,
                       paste0(save_as, "_analyzed.csv"),
                       sep = "/"))
    
    upload = list(
      timestamp = Sys.time(),
      output_dir = output_dir,
      save_as = save_as,
      radius = radius,
      drop_contours = drop_contours,
      frame_rate = frame_rate,
      data = data,
      video_length=video_length,
      threshold_head_vector_angular_speed = threshold_head_vector_angular_speed,
      threshold_tail_vector_angular_speed = threshold_tail_vector_angular_speed,
      min_hc_size = min_hc_size,
      max_hc_size = max_hc_size
    )
    
    saveRDS(upload,paste(output_dir,
                         paste0(save_as, "_analyzed.rds"),
                         sep = "/"))
    it <- toc()
    exectime <- it$toc - it$tic
    message(paste0(exectime , "seconds elapsed"))
    return(data)
  }

