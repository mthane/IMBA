source("functions_track_test.R")
source("functions_help_test.R")
source("functions_processing_tracked_test.R")
source("functions_summarise_test.R")
source("functions_binning_test.R")
library(tictoc)
library(data.table)
library(RcppRoll)
library(dplyr)
library(stats)
library(stringr)
library(zoo)
library(ggplot2)
library(gmodels)
library(scales)

library(mclust)
data = fread("C:/MThane/TestData/experiment-2021-12-17-analysis/experiment-2021-12-17_analyzed.csv")
colnames(data)
ids = unique(data$id)
ids
aggr <- summarise_TRACK_per_group_dt_test(data,"id")
View(aggr)

dataM <- data%>%
  select(midpoint_speed_bl,
         tail_vel_forward_bl,
         head_vector_angular_speed,
         tail_vector_angular_speed,
         head_vector_x,
         head_vector_y,
         tail_vector_x,
         tail_vector_y,
         bending_angle,
         bending_angle,
         spine_length,
         perimeter,
         width
         )%>%
  na.omit()

model <- Mclust(dataM)

model2 <- Mclust(dataM,G = 5)
df_class <- data%>%
  select(group,
         frame,
         id,
        midpoint_speed_bl,
         tail_vel_forward_bl,
         head_vector_angular_speed,
         tail_vector_angular_speed,
         head_vector_x,
         head_vector_y,
         tail_vector_x,
         tail_vector_y,
         bending_angle,
         bending_angle,
         spine_length,
         perimeter,
         width
  )%>%
  na.omit()%>%
  mutate(cluster = factor(model2$classification))%>%
  mutate(group = factor(group))


# 
for (i in c(12)){#unique(df_class$id)){
  track <- df_class %>%
    filter(id==i)
  intercept1 = track$frame[which(track$cluster==1)]
  intercept2 = track$frame[which(track$cluster==2)]
  intercept3 =  track$frame[which(track$cluster==3)]
  
  track %>%
    ggplot(aes(frame,tail_vel_forward_bl))+
    geom_line()+
    geom_vline(xintercept = intercept1,color="green",alpha=0.1)+
    geom_vline(xintercept = intercept2,color="blue",alpha=0.1)+
    geom_vline(xintercept = intercept3,color="red",alpha=0.1)
  
  
  ggsave(paste0(i,"gmm3.png"))
}


track <- df_class %>%
  filter(id==4)
intercept1 = track$frame[which(track$cluster==1)]
intercept2 = track$frame[which(track$cluster==2)]
intercept3 =  track$frame[which(track$cluster==3)]
intercept4 =  track$frame[which(track$cluster==4)]
intercept5 =  track$frame[which(track$cluster==5)]

track %>%
  ggplot(aes(frame,tail_vel_forward_bl))+
  geom_line()+
  geom_line(aes(frame,head_vector_angular_speed),data=track)+
  geom_vline(xintercept = intercept1,color="green",alpha=0.1)+
  geom_vline(xintercept = intercept2,color="blue",alpha=0.1)+
  geom_vline(xintercept = intercept3,color="red",alpha=0.1)+
  geom_vline(xintercept = intercept4,color="orange",alpha=0.1)+
  geom_vline(xintercept = intercept5,color="yellow",alpha=0.1)


df_tsne <- df_class%>%
  mutate(id = factor(id))%>%
  select(
         group,
         id,
         tail_vel_forward_bl,
         midpoint_speed_bl,
         head_vector_angular_speed,
         head_vector_x,
         head_vector_y,
         tail_vector_x,
         tail_vector_y)

df_tsne <- df_tsne%>%
  filter(id %in% unique(df_tsne$id)[1:50])

library(stats)

M <- df_tsne%>%select(!c("cluster","group","id"))%>%as.matrix()
components <- prcomp(t(M))
dft <- as.data.frame(components$rotation)
df_test <- cbind(df_tsne,dft)

df_test%>%
  ggplot(aes(PC1,PC2,color=cluster))+
  geom_line(alpha=0.3)


df_test
library(randomForest)
df_tsne <- droplevels(df_tsne)

rf1 <-randomForest(id~.,data=df_tsne%>%select(!group))
rf1

rf1$importance


rf1$confusion

library(dplyr)
library(ggplot2)
data%>%
  mutate(time=frame/16)%>%
  select(time,flanks)%>%
  filter(flanks!=0)%>%
  mutate(flanks= factor(flanks))%>%
  ggplot(aes(x=time,fill=flanks))+
  geom_density(alpha=0.3)
ggsave("switch_point_density.png")


track <- fread("C:/MThane/TestData/Tracked_2021-10-22/1_2/EM_Rewarded/lightbox2-2021-10-22_15_16_04/2.csv")
track%>%
  ggplot(aes(V1,V2))+
  geom_line(
    
    
  )

View(track[which(track$V78!=0)])



# Dummy data
x <- LETTERS[1:20]
y <- paste0("var", seq(1,20))
data <- expand.grid(X=x, Y=y)
data$Z <- runif(400, 0, 5)

