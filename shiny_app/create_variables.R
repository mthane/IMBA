library(dplyr)
library(tidyr)
library(stringr)

list_names <- function(df,cat){
  df <- df %>% 
    filter(category==cat)%>%
    select(variable,name)
  variables <- df$variable
  names(variables) <- df$name
  variables
}


##### SINGLE
SINGLE_VARTABLE = read.csv("variables/singleModeVariables.csv", encoding="UTF-8")


SPINES <-
  unlist(str_split(paste(
    paste("spinepoint_x", seq(1, 12), sep = "_"),
    paste("spinepoint_y", seq(1, 12), sep = "_")
  ), " "))
SPINES_C <-
  unlist(str_split(paste(
    paste("spinepoint_x", seq(1, 12),"conv", sep = "_"),
    paste("spinepoint_y", seq(1, 12), "conv",sep = "_")
  ), " "))


TS_VARIABLES_SINGLE <-
  SINGLE_VARTABLE %>%
  filter(category=="TS")%>%
  select("variable","name")%>%
  spread(name, variable)%>%
  as.list()



HC_VARIABLES_SINGLE <-
  SINGLE_VARTABLE %>%
  filter(category=="HC")%>%
  select("variable","name")%>%
  spread(name, variable)%>%
  as.list()

RUN_VARIABLES_SINGLE =
  SINGLE_VARTABLE %>%
  filter(category=="IS")%>%
  select("variable","name")%>%
  spread(name, variable)%>%
  as.list()

PREF_VARIABLES_SINGLE <-
  SINGLE_VARTABLE %>%
  filter(category=="PREF")%>%
  select("variable","name")%>%
  spread(name, variable)%>%
  as.list()

SINGLE_VARIABLES =
  list(`TS-Variables` = TS_VARIABLES_SINGLE,
       `HC-Variables` = HC_VARIABLES_SINGLE,
       `Run-Variables` = RUN_VARIABLES_SINGLE,
       `Pref-Variables` = PREF_VARIABLES_SINGLE
  )

##### BINNING #####

BINNING_VARTABLE = read.csv("variables/binningVariables.csv")
BINNING_VARIABLES =
  list(`TS-Variables` = list_names(BINNING_VARTABLE,"TS"),
       `HC-Variables` = list_names(BINNING_VARTABLE,"HC"),
       `Run-Variables` =list_names(BINNING_VARTABLE,"RUN"),
       `Pref-Variables` =list_names(BINNING_VARTABLE,"PREF")
  )

HEATMAP_VARTABLE = read.csv("variables/heatmapVariables.csv", encoding="UTF-8")
HEATMAP_VARIABLES =
  list(`TS-Variables` = list_names(HEATMAP_VARTABLE,"TS"),
       `HC-Variables` = list_names(HEATMAP_VARTABLE,"HC"),
       `Run-Variables` =list_names(HEATMAP_VARTABLE,"RUN"),
       `Pref-Variables` =list_names(HEATMAP_VARTABLE,"PREF")
  )




#####  SUMMARY #####
SUMMARY_VARTABLE = read.csv("variables/summaryVariables.csv", encoding="UTF-8")%>%
  filter(variable!="")
View(SUMMARY_VARTABLE)
TS_VARIABLES <-  list_names(SUMMARY_VARTABLE,"TS")
HC_VARIABLES <-  list_names(SUMMARY_VARTABLE,"HC")
RUN_VARIABLES <-  list_names(SUMMARY_VARTABLE,"RUN")
print(RUN_VARIABLES)
PREF_VARIABLES <-  list_names(SUMMARY_VARTABLE,"PREF")
TRACK_VARIABLES <-  list_names(SUMMARY_VARTABLE,"TRACK")





SUMMARY_VARIABLES =
  list(`TS-Variables` = list_names(SUMMARY_VARTABLE,"TS"),
       `HC-Variables` = list_names(SUMMARY_VARTABLE,"HC"),
       `Run-Variables` =list_names(SUMMARY_VARTABLE,"RUN"),
       `Pref-Variables` =list_names(SUMMARY_VARTABLE,"PREF"),
       `Tracking-Variables` =list_names(SUMMARY_VARTABLE,"TRACK")
  )



SUMMARY_VARIABLES_2 <- unlist(unname(c(list_names(SUMMARY_VARTABLE,"TS"),
                                       list_names(SUMMARY_VARTABLE,"HC"),
                                       list_names(SUMMARY_VARTABLE,"RUN"),
                                       list_names(SUMMARY_VARTABLE,"PREF"),
                                       list_names(SUMMARY_VARTABLE,"TRACK")
)))


VARIABLES <- list(
  SPINES = SPINES,
  SPINES_C = SPINES_C,
  TS_VARIABLES_SINGLE = TS_VARIABLES_SINGLE,
  HC_VARIABLES_SINGLE = HC_VARIABLES_SINGLE,
  RUN_VARIABLES_SINGLE = RUN_VARIABLES_SINGLE,
  SINGLE_VARTABLE = SINGLE_VARTABLE,
  SINGLE_VARIABLES = SINGLE_VARIABLES,
  BINNING_VARTABLE = BINNING_VARTABLE,
  BINNING_VARIABLES = BINNING_VARIABLES,
  HEATMAP_VARTABLE = HEATMAP_VARTABLE,
  HEATMAP_VARIABLES = HEATMAP_VARIABLES,
  TS_VARIABLES = TS_VARIABLES,
  HC_VARIABLES = HC_VARIABLES,
  RUN_VARIABLES = RUN_VARIABLES,
  PREF_VARIABLES = PREF_VARIABLES,
  TRACK_VARIABLES = TRACK_VARIABLES,
  SUMMARY_VARTABLE = SUMMARY_VARTABLE,
  SUMMARY_VARIABLES = SUMMARY_VARIABLES,
  SUMMARY_VARIABLES_2 = SUMMARY_VARIABLES_2
)
saveRDS(VARIABLES,"global_variables.Rds")
