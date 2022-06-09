#!/usr/bin/env Rscript

#Sebastian Glaess
#Masterstudiengang Statistik
#Otto-von-Guericke-Universitaet Magdeburg


FORCING = 0				#1: Force an assignment for an illogically sound combination of assignments
					#0: Stick to the high quality


### TIME ###
#start <- Sys.time ()


### INITIATE NUMBERS FOR DECISIONS ###
#gesamt=0				#number of all collisions where 2-bundle -> 2 seperate larvae
#erkannt=0				#number of all solved collisions
#nichterkannt=0			#number of non-solved collisions
#erkannt.normal=0			#number of all solved collisions where before and after collisions >=20 frames for each larva
#nichterkannt.normal=0		#number of non-solved collisions where before and after collisions >=20 frames for each larva


##################################
### LOAD LIBRARIES, SET PATH #####
##################################
#library(xlsx)
#library(FactoMineR)
#C:\Users\sebas_000\Desktop\DATEN 27.01.2016\AM_Rewarded
#C:\Users\Sebastian\Desktop\TRACKED_ON_13_01_2016M\Train_FRU_Test_T1\AM_Rewarded
#C:\Users\sebas_000\Desktop\Backup AM_Rewarded\AM_Rewarded
#folder.path <- file.path("C:", "Users", "sebas_000", "Desktop", "Backup AM_Rewarded", "AM_Rewarded")
#setwd(folder.path)
#foldernames <- dir(pattern="box")
#foldernames


##################################
### INITIATE VARIABLES ###########
##################################
number.all.results <- NULL
number.detected <- NULL
#micro.coll <- NULL
#number.micro.coll <- NULL		#Amount of multiple collisions within the video (without NA's)
#coll.micro.coll <- NULL		#Amount of collisions within all multiple collisions per video (without NA's)
#micro.coll.fulldata <- NULL		#Full data frame with NA's
						#NA in out2/frames.out2 and frames.out1==0 -> new larva enters collisions
						#NA in frames.out1 -> larvae split, but <20 frames til video ends or larva leaves tracking area
						#NA in out1/out2/frames.out1/frames.out2 -> larvae together while end of video or together outside
result.table <- NULL

number.3coll <- NULL
number.3coll.1solved <- NULL
number.3coll.3solved <- NULL
result3.table <- NULL
halfvideo <- 0				#Count of larvae with at least 1440 frames (over all videos in the folder)
guess <- NULL				#Number of guessed amound of larvae on a petridish
allguess <- NULL				#Guessed amount of larvae on the petridishes (over all videos in the folder)
model2x2 <- NULL				#Number of solved 2x2 collisions by Manos model
single <- NULL				#Number of single larvae tracks on a petridish
nframes1 <- 0				#Amount of all frames (over all videos in the folder)

eigenvalues <- NULL
expl.var <- NULL

number2x2 <- 0				#Number of all 2x2 colls, all frames accepted

n_2x2_assign <- 0 			#Amount of assignments for 2x2 colls
n_3x3_assign <- 0				#Amount of assignments for 3x3 colls
n_nx1_assign <- 0				#Amount of assignments for nx1 colls
n_1x1_assign <- 0				#Amount of 1x1 assignments (not real statistics...)
n_2x1_assign <- 0
n_3x1_assign <- 0
n_4x1_assign <- 0
n_5x1_assign <- 0
n_assignments <- 0 			#Amount of assignments in total

assignments <- data.frame(Incoming=numeric(),Outgoing=numeric(),TypeOfCollision=character(),
                          Duration=numeric(),Forced=numeric())
#Forced==1: combination of assignments is not logically sound, 
#Forced==0: logically sound combination of assignments
#Notice: Nx1 is always logically sound!

##################################
### FUNCTIONS ####################
##################################

#####Correlations between the variables#####
#mosthighlycorrelated and calcpc
#http://little-book-of-r-for-multivariate-analysis.readthedocs.org/en/latest/src/multivariateanalysis.html
mosthighlycorrelated <- function(mydataframe,numtoreport)
  {
     cormatrix <- cor(mydataframe, method="spearman")				#Spearman f�r nonparamatrical data 
     diag(cormatrix) <- 0
     cormatrix[lower.tri(cormatrix)] <- 0
     fm <- as.data.frame(as.table(cormatrix))
     names(fm) <- c("First.Variable", "Second.Variable","Correlation")
     head(fm[order(abs(fm$Correlation),decreasing=T),],n=numtoreport)
  }

#####New component from PCA##################
calcpc <- function(variables,loadings)
  {
     as.data.frame(variables)
     numsamples <- nrow(variables)		#number of rows
     pc <- numeric(numsamples)		#vector
     numvariables <- length(variables)	#number of variables
     for (i in 1:numsamples)			#calculate value for the prcomp for everything
     {
        valuei <- 0
        for (j in 1:numvariables)
        {
           valueij <- variables[i,j]
           loadingj <- loadings[j]
           valuei <- valuei + (valueij * loadingj)
        }
        pc[i] <- valuei
     }
     return(pc)
  }

print("FOO")
sayHello <- function() {
		print('hello')
}


box <- commandArgs(trailingOnly = TRUE)
print(box)
#box.path <- file.path(box.path, i)
box.path <- box
log.path <- file.path(box.path, "vidAndLogs")
setwd(box.path)

#####Load Tracking Data##############
names <-dir(pattern="*.csv")
alldata <- NULL
for(i in names)
  {
     tab = read.csv(file=i, header=F, sep=",")
     tab["csv_name"] <- basename(i)
     tab <- tab[,-c(2:68,72)]
     csvnumber <- sapply(strsplit(tab[,11], split='.', fixed=TRUE), function(x) (x[1]))
     csvnumber <- strtoi(csvnumber)
     tab["csv_number"] <- csvnumber
     colnames(tab) <- c("Frame","Centroid_Y","Centroid_X","Orientation","Area",
				"Grey_Sum","Length","Width","Perimeter","Indicator","csv_Name","csv_Number")
     alldata <- rbind(alldata, tab)
  }
alldata["Grey_Average"] <- alldata$Grey_Sum/alldata$Area
alldata <- alldata[order(alldata$csv_Number, alldata$Frame),]
alldata <- subset(alldata, alldata$Indicator==0)			#Model approach data (2), tracking error (1) excluded

needed <- alldata[,-c(1:4,6,10:11)]						#Leave only Area(5), Length(7), Width(8), 
											#Perimeter(9), csv_Number(12), Grey_Average(13)

nframes1 <- nrow(needed)+nframes1

#####Load Collisions Data############	
setwd(log.path)
collisions <- read.csv("collisions.csv", header=F, sep=",")			
colnames(collisions) <- c("csv_Number", "New1", "New2", "New3", "Origin", "Bundle")
collisions$n_frames <- NULL	
collisions <- collisions[order(collisions$csv_Number, collisions$New1, collisions$New2),]

for(i in nrow(collisions):2)
  {
     if(collisions$csv_Number[i]==collisions$csv_Number[i-1])
       {
          collisions <- collisions[-c(i-1),]
       }
  }

for(i in 1:max(collisions$csv_Number))
  {
     if(is.na(match(i,collisions$csv_Number)))
       {
          collisions <- rbind(collisions, c(i,NA,NA,NA,NA,0,NA))
       }
  }
collisions <- collisions[order(collisions$csv_Number),]
for(i in collisions$csv_Number)
  {
     collisions$n_frames[i] <- nrow(subset(needed, needed$csv_Number==i))
  }	
#if collisions$csv_Number[i] does not exist -> collisions$n_frames[i]=0


collisions$Start <- rep(NA,nrow(collisions))		#Starting frame
collisions$End <- rep(NA,nrow(collisions))		#Ending frame --- will be used for border and odor cup decisions
for(i in 1:nrow(collisions))
  {
     if(!is.na(subset(alldata, alldata$csv_Number==i)$Frame[1]) 
     && !is.na(subset(alldata, alldata$csv_Number==i)$Frame[nrow(subset(alldata, alldata$csv_Number==i))])) 
       {
          collisions$Start[i]=min(subset(alldata,alldata$csv_Number==i)$Frame)
          collisions$End[i]=max(subset(alldata,alldata$csv_Number==i)$Frame)     
       }
  }


######################################################
### TESTING CORRELATIONS AND PCA #####################
######################################################
setwd(box.path)

#####Correlations####################
spearman.corr <- mosthighlycorrelated(needed[,1:4],6)
#write.xlsx(spearman.corr, row.names = FALSE, "Spearman_Correlations.xlsx")

#####PCA#############################
#pca1 <- PCA(needed[,1:4])
#PrComp <- pca1$eig 					#(Eigenvalue = sdev^2) Kaisers Criterion: Use all components with Eigenvalue > 1
#PrComp.corr <- dimdesc(pca1, axes=c(1,2)) 	#Describes Correlations between PrComp and Variables
#write.xlsx(PrComp, "PrComp.xlsx")
#write.xlsx(PrComp.corr, "PrCompCorr.xlsx")

standardise.needed <- as.data.frame(scale(needed[,1:4]))  		#Leave out the Grey_Average
needed.pca <- prcomp(standardise.needed)                
summary(needed.pca)
#screeplot(needed.pca, type="lines")
#needed.pca$rotation[,1] 			
#sum((needed.pca$rotation[,1])^2)
first.comp <- calcpc(standardise.needed, needed.pca$rotation[,1])
needed$PrComp1 <- first.comp
#second.comp <- calcpc(standardise.needed, needed.pca$rotation[,2])
#needed$PrComp2 <- second.comp
eigenvalues <- rbind(eigenvalues,(needed.pca$sdev)^2)
expl.var <- rbind(expl.var,(needed.pca$sdev)^2/length(needed.pca$sdev))
#PCA.results <- cbind(c(1:4), needed.pca$sdev, eigenvalues, expl.var)
#colnames(PCA.results) <- c("PrComp","sdev", "eigenvalue", "expl.var")
#write.xlsx(PCA.results, row.names = FALSE, "PCA_sdev_eig_explvar.xlsx")
#spearman.corr <- mosthighlycorrelated(cbind(needed[,1:4],first.comp),10)


#Plotting stuff... 
if(T){
jpeg(file="saving_plot1.jpeg")
boxplot(cbind(expl.var,rep(1,nrow(eigenvalues))), xaxt="n", yaxt="n", ylim=c(0,1), horizontal=T, border=c("black","blue","darkgreen","brown","white"),at=c(1,1,1,1,.8), names=F)
mtext("eigenvalues",side=1, line=2.25) 
axis(1, at=seq(0,1,1))
mtext("explained variance",side=3, line=2.25) 
axis(3, at=seq(0,1,.1))
legend((-.28+4*0.003),.5,"PC4",bty="n",cex=.8, text.col="brown")
legend((-.15+4*0.015),.5,"PC3",bty="n",cex=.8,text.col="darkgreen")
legend((-.25+4*0.117),.5,"PC2",bty="n",cex=.8,text.col="blue")
legend((-.25+4*0.89),.5,"PC1",bty="n",cex=.8,text.col="black")
dev.off()
}

######################################################
### TEST OF NORMALITY ################################	#no need for this, we stick with non parametric tests
######################################################
#START_IF_FALSE
if(FALSE){

setwd(box.path)

#####Test by view####################
#L�nge der tracks bestimmen
tracklength <- data.frame(unique(needed$csv_Number))
colnames(tracklength) = c("csv_Number")
tracklength$length <- 0
for(i in 1:length(unique(needed$csv_Number))){
	tracklength$length[i] <- nrow(subset(needed, needed$csv_Number == unique(needed$csv_Number)[i]))
	}

par(mfrow=c(1,2))
qqnorm(subset(needed$PrComp1, needed$csv_Number==159), xlim=c(-5,5), ylim=c(-5,5))
abline(0,1)
qqnorm(subset(needed$PrComp1, needed$csv_Number==169), xlim=c(-5,5), ylim=c(-5,5))
abline(0,1)


#####Shapiro-Wilk####################				#needs between 3 and 5000 values (max frames per track: 2880)
Shapiro_P <- as.data.frame(rep(NA,max(needed$csv_Number)))	#Test on normal dist. with first comp. of PCA
for(i in unique(needed$csv_Number))
  {
     Shapiro_P[i,1] <- i     
     Shapiro_P[i,2] <- shapiro.test(subset(needed$PrComp1, needed$csv_Number==i))$p.value
     Shapiro_P[i,3] <- nrow(subset(needed, needed$csv_Number==i))
  }
names(Shapiro_P) <- c("csv-Number", "p-value", "samplesize")
Shapiro_P <- na.omit(Shapiro_P)
#write.xlsx(Shapiro_P, row.names =FALSE, "Shapiro.xlsx")

#####Komogorov-Smirnov###############				#Test of Normality (fails at n>1000)
#ks.test(subset(needed$Length, needed$csv_Number==1), "pnorm")

}
#END_IF_FALSE


######################################################
### 2x2 collisions ###################################
######################################################
#no cross testing
setwd(box.path)

#####Wilcoxon Mann-Whitney Rank Sum Test#####	
Wilcox_P <- array(NA, c(0,8))				#only 2-bundle collisions with >=10 frames per larva
for(i in collisions$csv_Number)
  {
     laufdf <- rbind(subset(collisions,collisions$New1==i),subset(collisions,collisions$csv_Number==i),
			   subset(collisions,collisions$csv_Number==collisions$New1[i]),
			   subset(collisions,collisions$csv_Number==collisions$New2[i]))

     #all 2x2 colls:
     if(nrow(laufdf)==5 && laufdf$Bundle[1]==1 && laufdf$Bundle[2]==1 && laufdf$Bundle[3]==2 && laufdf$Bundle[4]==1 && laufdf$Bundle[5]==1)
	 {number2x2 = number2x2 + 1}

     if(nrow(laufdf)==5 && laufdf$Bundle[1]==1 && laufdf$Bundle[2]==1 && laufdf$Bundle[3]==2 && laufdf$Bundle[4]==1 && laufdf$Bundle[5]==1
		&& nrow(subset(needed,needed$csv_Number==laufdf$csv_Number[1]))>=10
		&& nrow(subset(needed,needed$csv_Number==laufdf$csv_Number[2]))>=10
		&& nrow(subset(needed,needed$csv_Number==laufdf$New1[3]))>=10
		&& nrow(subset(needed,needed$csv_Number==laufdf$New2[3]))>=10)
	 {  
          in1.frames <- nrow(subset(needed,needed$csv_Number==laufdf$csv_Number[1]))
          in2.frames <- nrow(subset(needed,needed$csv_Number==laufdf$csv_Number[2]))
          out1.frames <- nrow(subset(needed,needed$csv_Number==laufdf$New1[3]))
          out2.frames <- nrow(subset(needed,needed$csv_Number==laufdf$New2[3]))
          max.frames <- min(in1.frames, in2.frames, out1.frames, out2.frames, 45)

	    datas = rbind(subset(needed, needed$csv_Number==laufdf$csv_Number[1])[(in1.frames-max.frames+1):(in1.frames),],
				subset(needed, needed$csv_Number==laufdf$csv_Number[2])[(in2.frames-max.frames+1):(in2.frames),],
		  		subset(needed, needed$csv_Number==laufdf$New1[3])[1:max.frames,],
		  		subset(needed, needed$csv_Number==laufdf$New2[3])[1:max.frames,])

	    if(kruskal.test(PrComp1~csv_Number, data=datas)$p.value <= 0.05)
          {
	    p.11 <- wilcox.test(PrComp1~csv_Number, 
		data = rbind(subset(needed, needed$csv_Number==laufdf$csv_Number[1])[(in1.frames-max.frames+1):(in1.frames),],
                         subset(needed, needed$csv_Number==laufdf$New1[3])[1:max.frames,]))
	    p.12 <- wilcox.test(PrComp1~csv_Number, 
		data = rbind(subset(needed, needed$csv_Number==laufdf$csv_Number[1])[(in1.frames-max.frames+1):(in1.frames),],
                         subset(needed, needed$csv_Number==laufdf$New2[3])[1:max.frames,]))
	    p.21 <- wilcox.test(PrComp1~csv_Number, 
		data = rbind(subset(needed, needed$csv_Number==laufdf$csv_Number[2])[(in2.frames-max.frames+1):(in2.frames),],
                         subset(needed, needed$csv_Number==laufdf$New1[3])[1:max.frames,]))
	    p.22 <- wilcox.test(PrComp1~csv_Number, 
		data = rbind(subset(needed, needed$csv_Number==laufdf$csv_Number[2])[(in2.frames-max.frames+1):(in2.frames),],
                         subset(needed, needed$csv_Number==laufdf$New2[3])[1:max.frames,]))
	    Wilcox_P <- rbind(Wilcox_P, c(laufdf$csv_Number[1], laufdf$csv_Number[2], laufdf$New1[3], laufdf$New2[3], 
		p.11$p.value, p.12$p.value, p.21$p.value, p.22$p.value))
	    }
	    else{Wilcox_P=rbind(Wilcox_P, c(laufdf$csv_Number[1], laufdf$csv_Number[2], laufdf$New1[3], laufdf$New2[3],
		  		0,0,0,0))
	    }
	 }
  }
Wilcox_P <- as.data.frame(Wilcox_P)

wilcox.p <- Wilcox_P	

number.all.results <- rbind(number.all.results, nrow(wilcox.p))
wilcox.p$detected <- rep(0, nrow(wilcox.p))
wilcox.p$assignm1 <- rep(NA, nrow(wilcox.p))
wilcox.p$assignm2 <- rep(NA, nrow(wilcox.p))

colnames(wilcox.p) <- c("In1", "In2", "Out1", "Out2", "p.In1.Out1", "p.In1.Out2", "p.In2.Out1", "p.In2.Out2", 
				"detected", "assignm1", "assignm2")

if(nrow(wilcox.p)!=0)
{
#detection:	1 -> pattern ok (2 highest values are I1O1/I2O2 or I1O2/I2O1)
#		0 -> pattern not correct or values are to close to each other
for(i in 1:nrow(wilcox.p))
  {
     if(((max(wilcox.p[i,5:8])!=0) &&
        (wilcox.p$p.In1.Out1[i]==max(wilcox.p[i,5:8]) 
	    && wilcox.p$p.In2.Out2[i]==max(wilcox.p[i,6:8]))||
	  (wilcox.p$p.In2.Out2[i]==max(wilcox.p[i,5:8])
          && wilcox.p$p.In1.Out1[i]==max(wilcox.p[i,5:7]))||
	  (wilcox.p$p.In1.Out2[i]==max(wilcox.p[i,5:8])
	    && wilcox.p$p.In2.Out1[i]==max(wilcox.p[i,c(5,7:8)]))||
	  (wilcox.p$p.In2.Out1[i]==max(wilcox.p[i,5:8])
	    && wilcox.p$p.In1.Out2[i]==max(wilcox.p[i,c(5:6,8)])))
        && !is.na(min(sort(wilcox.p[i,5:8])[3]/sort(wilcox.p[i,5:8])[2], sort(wilcox.p[i,5:8])[2]/sort(wilcox.p[i,5:8])[3])<0.001)
	  && if(!is.na(min(sort(wilcox.p[i,5:8])[3]/sort(wilcox.p[i,5:8])[2], sort(wilcox.p[i,5:8])[2]/sort(wilcox.p[i,5:8])[3])<0.001))
               {
                 min(sort(wilcox.p[i,5:8])[3]/sort(wilcox.p[i,5:8])[2], sort(wilcox.p[i,5:8])[2]/sort(wilcox.p[i,5:8])[3])<0.001
               })

	 {
	    wilcox.p$detected[i] <- 1
	 }
     #Duration of the collision
     Duration = collisions$Start[wilcox.p$Out1[i]] - collisions$End[wilcox.p$In1[i]] - 1

     #Forced assignments:
     if(FORCING==1)
      {
	 if(wilcox.p$detected[i]==0 && (wilcox.p$p.In1.Out1[i]==max(wilcox.p[i,5:8]) || wilcox.p$p.In2.Out2[i]==max(wilcox.p[i,5:8])))
       {
          wilcox.p$assignm1[i] <- paste(wilcox.p$In1[i], "->", wilcox.p$Out1[i])						#ASSIGNMENTS
	    wilcox.p$assignm2[i] <- paste(wilcox.p$In2[i], "->", wilcox.p$Out2[i])   
	    assignments <- rbind(assignments, data.frame('Incoming'=wilcox.p$In1[i], 'Outgoing'=wilcox.p$Out1[i], 
					'TypeOfCollision'="2x2", 'Duration'=Duration, 'Forced'=1))
	    assignments <- rbind(assignments, data.frame('Incoming'=wilcox.p$In2[i], 'Outgoing'=wilcox.p$Out2[i], 
					'TypeOfCollision'="2x2", 'Duration'=Duration, 'Forced'=1))         
	    needed$csv_Number[needed$csv_Number==wilcox.p$In1[i]] <- wilcox.p$Out1[i]						#COMBINE THE TRACKS, InTrack gets OutNumber
          needed$csv_Number[needed$csv_Number==wilcox.p$In2[i]] <- wilcox.p$Out2[i]
          collisions$n_frames[wilcox.p$In1[i]] <- nrow(subset(needed,needed$csv_Number==wilcox.p$In1[i]))		#ADJUST FRAMENUMBERS
          collisions$n_frames[wilcox.p$Out1[i]] <- nrow(subset(needed,needed$csv_Number==wilcox.p$Out1[i]))
          collisions$n_frames[wilcox.p$In2[i]] <- nrow(subset(needed,needed$csv_Number==wilcox.p$In2[i]))
          collisions$n_frames[wilcox.p$Out2[i]] <- nrow(subset(needed,needed$csv_Number==wilcox.p$Out2[i]))
          collisions$End[wilcox.p$In1[i]] <- subset(collisions, collisions$csv_Number==wilcox.p$Out1[i])$End		#ADJUST STARTING AND ENDING FRAME
          collisions$End[wilcox.p$In2[i]] <- subset(collisions, collisions$csv_Number==wilcox.p$Out2[i])$End
          collisions$Start[wilcox.p$Out1[i]] <- subset(collisions, collisions$csv_Number==wilcox.p$In1[i])$Start	
          collisions$Start[wilcox.p$Out2[i]] <- subset(collisions, collisions$csv_Number==wilcox.p$In2[i])$Start
       }
	 if(wilcox.p$detected[i]==0 && (wilcox.p$p.In1.Out2[i]==max(wilcox.p[i,5:8]) || wilcox.p$p.In2.Out1[i]==max(wilcox.p[i,5:8])))
       {
          wilcox.p$assignm1[i] <- paste(wilcox.p$In1[i], "->", wilcox.p$Out2[i])
	    wilcox.p$assignm2[i] <- paste(wilcox.p$In2[i], "->", wilcox.p$Out1[i])
	    assignments <- rbind(assignments, data.frame('Incoming'=wilcox.p$In1[i], 'Outgoing'=wilcox.p$Out2[i], 
					'TypeOfCollision'="2x2", 'Duration'=Duration, 'Forced'=1))
	    assignments <- rbind(assignments, data.frame('Incoming'=wilcox.p$In2[i], 'Outgoing'=wilcox.p$Out1[i], 
					'TypeOfCollision'="2x2", 'Duration'=Duration, 'Forced'=1))
          needed$csv_Number[needed$csv_Number==wilcox.p$In1[i]] <- wilcox.p$Out2[i]
          needed$csv_Number[needed$csv_Number==wilcox.p$In2[i]] <- wilcox.p$Out1[i]
          collisions$n_frames[wilcox.p$In1[i]] <- nrow(subset(needed,needed$csv_Number==wilcox.p$In1[i]))
          collisions$n_frames[wilcox.p$Out1[i]] <- nrow(subset(needed,needed$csv_Number==wilcox.p$Out1[i]))
          collisions$n_frames[wilcox.p$In2[i]] <- nrow(subset(needed,needed$csv_Number==wilcox.p$In2[i]))
          collisions$n_frames[wilcox.p$Out2[i]] <- nrow(subset(needed,needed$csv_Number==wilcox.p$Out2[i]))
          collisions$End[wilcox.p$In1[i]] <- subset(collisions, collisions$csv_Number==wilcox.p$Out2[i])$End	
          collisions$End[wilcox.p$In2[i]] <- subset(collisions, collisions$csv_Number==wilcox.p$Out1[i])$End
          collisions$Start[wilcox.p$Out1[i]] <- subset(collisions, collisions$csv_Number==wilcox.p$In2[i])$Start	
          collisions$Start[wilcox.p$Out2[i]] <- subset(collisions, collisions$csv_Number==wilcox.p$In1[i])$Start
       }
      }

     #Logically sound assignments:
     if(wilcox.p$detected[i]==1 && (wilcox.p$p.In1.Out1[i]==max(wilcox.p[i,5:8]) || wilcox.p$p.In2.Out2[i]==max(wilcox.p[i,5:8])))
       {
          wilcox.p$assignm1[i] <- paste(wilcox.p$In1[i], "->", wilcox.p$Out1[i])						#ASSIGNMENTS
	    wilcox.p$assignm2[i] <- paste(wilcox.p$In2[i], "->", wilcox.p$Out2[i])   
	    assignments <- rbind(assignments, data.frame('Incoming'=wilcox.p$In1[i], 'Outgoing'=wilcox.p$Out1[i], 
					'TypeOfCollision'="2x2", 'Duration'=Duration, 'Forced'=0))
	    assignments <- rbind(assignments, data.frame('Incoming'=wilcox.p$In2[i], 'Outgoing'=wilcox.p$Out2[i], 
					'TypeOfCollision'="2x2", 'Duration'=Duration, 'Forced'=0))
          needed$csv_Number[needed$csv_Number==wilcox.p$In1[i]] <- wilcox.p$Out1[i]						#COMBINE THE TRACKS, InTrack gets OutNumber
          needed$csv_Number[needed$csv_Number==wilcox.p$In2[i]] <- wilcox.p$Out2[i]
          collisions$n_frames[wilcox.p$In1[i]] <- nrow(subset(needed,needed$csv_Number==wilcox.p$In1[i]))		#ADJUST FRAMENUMBERS
          collisions$n_frames[wilcox.p$Out1[i]] <- nrow(subset(needed,needed$csv_Number==wilcox.p$Out1[i]))
          collisions$n_frames[wilcox.p$In2[i]] <- nrow(subset(needed,needed$csv_Number==wilcox.p$In2[i]))
          collisions$n_frames[wilcox.p$Out2[i]] <- nrow(subset(needed,needed$csv_Number==wilcox.p$Out2[i]))
          collisions$End[wilcox.p$In1[i]] <- subset(collisions, collisions$csv_Number==wilcox.p$Out1[i])$End		#ADJUST STARTING AND ENDING FRAME
          collisions$End[wilcox.p$In2[i]] <- subset(collisions, collisions$csv_Number==wilcox.p$Out2[i])$End
          collisions$Start[wilcox.p$Out1[i]] <- subset(collisions, collisions$csv_Number==wilcox.p$In1[i])$Start	
          collisions$Start[wilcox.p$Out2[i]] <- subset(collisions, collisions$csv_Number==wilcox.p$In2[i])$Start
       }
     if(wilcox.p$detected[i]==1 && (wilcox.p$p.In1.Out2[i]==max(wilcox.p[i,5:8]) || wilcox.p$p.In2.Out1[i]==max(wilcox.p[i,5:8])))
       {
          wilcox.p$assignm1[i] <- paste(wilcox.p$In1[i], "->", wilcox.p$Out2[i])
	    wilcox.p$assignm2[i] <- paste(wilcox.p$In2[i], "->", wilcox.p$Out1[i])
	    assignments <- rbind(assignments, data.frame('Incoming'=wilcox.p$In1[i], 'Outgoing'=wilcox.p$Out2[i], 
					'TypeOfCollision'="2x2", 'Duration'=Duration, 'Forced'=0))
	    assignments <- rbind(assignments, data.frame('Incoming'=wilcox.p$In2[i], 'Outgoing'=wilcox.p$Out1[i], 
					'TypeOfCollision'="2x2", 'Duration'=Duration, 'Forced'=0))
          needed$csv_Number[needed$csv_Number==wilcox.p$In1[i]] <- wilcox.p$Out2[i]
          needed$csv_Number[needed$csv_Number==wilcox.p$In2[i]] <- wilcox.p$Out1[i]
          collisions$n_frames[wilcox.p$In1[i]] <- nrow(subset(needed,needed$csv_Number==wilcox.p$In1[i]))
          collisions$n_frames[wilcox.p$Out1[i]] <- nrow(subset(needed,needed$csv_Number==wilcox.p$Out1[i]))
          collisions$n_frames[wilcox.p$In2[i]] <- nrow(subset(needed,needed$csv_Number==wilcox.p$In2[i]))
          collisions$n_frames[wilcox.p$Out2[i]] <- nrow(subset(needed,needed$csv_Number==wilcox.p$Out2[i]))
          collisions$End[wilcox.p$In1[i]] <- subset(collisions, collisions$csv_Number==wilcox.p$Out2[i])$End	
          collisions$End[wilcox.p$In2[i]] <- subset(collisions, collisions$csv_Number==wilcox.p$Out1[i])$End
          collisions$Start[wilcox.p$Out1[i]] <- subset(collisions, collisions$csv_Number==wilcox.p$In2[i])$Start	
          collisions$Start[wilcox.p$Out2[i]] <- subset(collisions, collisions$csv_Number==wilcox.p$In1[i])$Start
       }
  }
}
number.detected <- rbind(number.detected, sum(wilcox.p$detected))

#setwd(log.path)
for(i in 1:nrow(wilcox.p)){
	if(wilcox.p$detected[i]==0){wilcox.p$detected[i]="no"}
	else{wilcox.p$detected[i]="yes"}}
#write.xlsx(wilcox.p, row.names =FALSE, "results.xlsx")
wilcox.p

n_2x2_assign <- n_2x2_assign + 2*nrow(subset(wilcox.p, wilcox.p$detected=="yes"))
n_assignments <- n_assignments + 2*nrow(subset(wilcox.p, wilcox.p$detected=="yes"))
if(FORCING==1){
  n_2x2_assign <- n_2x2_assign + 2*nrow(subset(wilcox.p, wilcox.p$detected=="no"))
  n_assignments <- n_assignments + 2*nrow(subset(wilcox.p, wilcox.p$detected=="no"))
		 }


#################################################
### 3x3 collisions ##############################
#################################################
#Results have to fit in one out of 6 patterns 
#If they don't fit, test for 1 assignment instead of 3 (maybe one larva is strongly diffenrent to the others)
#cross testing added

setwd(box.path)
									#only 3-bundle collisions with >=10 frames per larva
bundle3 <- subset(collisions, collisions$Bundle==3)
wilcox3 <- array(NA, c(0,15))
retest.p <- array(NA, c(0,6))

for(i in bundle3$csv_Number)					#generate p-values for 9 possible assignments
  {
     s0 <- subset(bundle3,bundle3$csv_Number==i)									
     s1 <- subset(collisions,collisions$New1==i)									
     s2 <- subset(collisions,collisions$New1==subset(s1,s1$Bundle==2)$csv_Number)				
     s3 <- subset(collisions,collisions$csv_Number==subset(bundle3,bundle3$csv_Number==i)$New1)	
     s4 <- subset(collisions,collisions$csv_Number==subset(bundle3,bundle3$csv_Number==i)$New2)
     s5 <- subset(collisions,collisions$csv_Number==subset(rbind(s3,s4),rbind(s3,s4)$Bundle==2)$New1)
     s6 <- subset(collisions,collisions$csv_Number==subset(rbind(s3,s4),rbind(s3,s4)$Bundle==2)$New2)
     laufdf3 <- rbind(s2,s1,s0,s3,s4,s5,s6)
     inandout <- subset(laufdf3,laufdf3$Bundle==1) 

     if(nrow(inandout)==6 
		&& isTRUE(inandout$n_frames[1]>=10)
		&& isTRUE(inandout$n_frames[2]>=10)
		&& isTRUE(inandout$n_frames[3]>=10)
		&& isTRUE(inandout$n_frames[4]>=10)
		&& isTRUE(inandout$n_frames[5]>=10)
		&& isTRUE(inandout$n_frames[6]>=10))
       {
          in1.frames <- inandout$n_frames[1]
          in2.frames <- inandout$n_frames[2]
          in3.frames <- inandout$n_frames[3]
          out1.frames <- inandout$n_frames[4]
          out2.frames <- inandout$n_frames[5]
          out3.frames <- inandout$n_frames[6]
          max.frames <- min(in1.frames,in2.frames,in3.frames,out1.frames,out2.frames,out3.frames,45)

          in1.data <- subset(needed, needed$csv_Number==inandout$csv_Number[1])[(in1.frames-max.frames+1):(in1.frames),]
          in2.data <- subset(needed, needed$csv_Number==inandout$csv_Number[2])[(in2.frames-max.frames+1):(in2.frames),]
          in3.data <- subset(needed, needed$csv_Number==inandout$csv_Number[3])[(in3.frames-max.frames+1):(in3.frames),]
          out1.data <- subset(needed, needed$csv_Number==inandout$csv_Number[4])[1:max.frames,]
          out2.data <- subset(needed, needed$csv_Number==inandout$csv_Number[5])[1:max.frames,]
          out3.data <- subset(needed, needed$csv_Number==inandout$csv_Number[6])[1:max.frames,]

	    Duration1A <- inandout[4,]$Start - inandout[1,]$End - 1
	    Duration2A <- inandout[4,]$Start - inandout[2,]$End - 1
	    Duration3A <- inandout[4,]$Start - inandout[3,]$End - 1
          Duration1B <- inandout[5,]$Start - inandout[1,]$End - 1
          Duration2B <- inandout[5,]$Start - inandout[2,]$End - 1	
          Duration3B <- inandout[5,]$Start - inandout[3,]$End - 1
	    Duration1C <- inandout[6,]$Start - inandout[1,]$End - 1
          Duration2C <- inandout[6,]$Start - inandout[2,]$End - 1
          Duration3C <- inandout[6,]$Start - inandout[3,]$End - 1

	    datas <- rbind(in1.data, in2.data, in3.data, out1.data, out2.data, out3.data)
          
          if(kruskal.test(PrComp1~csv_Number, data=datas)$p.value <= 0.05)
            {
               p.1A <- wilcox.test(PrComp1~csv_Number, data=rbind(in1.data,out1.data))$p.value
               p.1B <- wilcox.test(PrComp1~csv_Number, data=rbind(in1.data,out2.data))$p.value
               p.1C <- wilcox.test(PrComp1~csv_Number, data=rbind(in1.data,out3.data))$p.value
               p.2A <- wilcox.test(PrComp1~csv_Number, data=rbind(in2.data,out1.data))$p.value
               p.2B <- wilcox.test(PrComp1~csv_Number, data=rbind(in2.data,out2.data))$p.value
               p.2C <- wilcox.test(PrComp1~csv_Number, data=rbind(in2.data,out3.data))$p.value
               p.3A <- wilcox.test(PrComp1~csv_Number, data=rbind(in3.data,out1.data))$p.value
               p.3B <- wilcox.test(PrComp1~csv_Number, data=rbind(in3.data,out2.data))$p.value
               p.3C <- wilcox.test(PrComp1~csv_Number, data=rbind(in3.data,out3.data))$p.value  
               wilcox3 <- rbind(wilcox3, c(inandout$csv_Number[1:6],p.1A,p.1B,p.1C,p.2A,p.2B,p.2C,p.3A,p.3B,p.3C))
               p.12 <- wilcox.test(PrComp1~csv_Number, data=rbind(in1.data,in2.data))$p.value
               p.13 <- wilcox.test(PrComp1~csv_Number, data=rbind(in1.data,in3.data))$p.value
               p.23 <- wilcox.test(PrComp1~csv_Number, data=rbind(in2.data,in3.data))$p.value
               p.AB <- wilcox.test(PrComp1~csv_Number, data=rbind(out1.data,out2.data))$p.value
               p.AC <- wilcox.test(PrComp1~csv_Number, data=rbind(out1.data,out3.data))$p.value
               p.BC <- wilcox.test(PrComp1~csv_Number, data=rbind(out2.data,out3.data))$p.value
		   retest.p <- rbind(retest.p, c(p.12,p.13,p.23,p.AB,p.AC,p.BC))
            }   
	    else
            {
		   wilcox3 = rbind(wilcox3, c(inandout$csv_Number[1:6],0,0,0,0,0,0,0,0,0))
		   retest.p = rbind(retest.p, c(0,0,0,0,0,0))
	      }    
 
       }
  }
wilcox3 <- as.data.frame(wilcox3)
colnames(wilcox3) <- c("In1","In2","In3","Out1","Out2","Out3","p.1-A","p.1-B","p.1-C","p.2-A","p.2-B","p.2-C","p.3-A","p.3-B","p.3-C")
wilcox3$decided <- rep(0, nrow(wilcox3))
wilcox3$assignm1 <- rep(NA, nrow(wilcox3))
wilcox3$assignm2 <- rep(NA, nrow(wilcox3))
wilcox3$assignm3 <- rep(NA, nrow(wilcox3))

if(nrow(wilcox3)!=0)
{
for(i in 1:nrow(wilcox3))				
  {												#Assign matches for all 3 larvae
     if(max(wilcox3[i,7:15])!=0
		&& ((p.1A == max(wilcox3[i,7:9]) && p.2B == max(wilcox3[i,10:12]) && p.3C == max(wilcox3[i,13:15])) ||
		    (p.1A == max(wilcox3[i,7:9]) && p.2C == max(wilcox3[i,10:12]) && p.3B == max(wilcox3[i,13:15])) ||
		    (p.1B == max(wilcox3[i,7:9]) && p.2A == max(wilcox3[i,10:12]) && p.3C == max(wilcox3[i,13:15])) ||
		    (p.1B == max(wilcox3[i,7:9]) && p.2C == max(wilcox3[i,10:12]) && p.3A == max(wilcox3[i,13:15])) ||
		    (p.1C == max(wilcox3[i,7:9]) && p.2A == max(wilcox3[i,10:12]) && p.3B == max(wilcox3[i,13:15])) ||
		    (p.1C == max(wilcox3[i,7:9]) && p.2B == max(wilcox3[i,10:12]) && p.3A == max(wilcox3[i,13:15])))
		&& !is.na(min(max(wilcox3[i,7:9])/sort(wilcox3[i,7:9])[2], sort(wilcox3[i,7:9])[2]/max(wilcox3[i,7:9]))<0.001)
		&& !is.na(min(max(wilcox3[i,10:12])/sort(wilcox3[i,10:12])[2], sort(wilcox3[i,10:12])[2]/max(wilcox3[i,10:12]))<0.001)
		&& !is.na(min(max(wilcox3[i,13:15])/sort(wilcox3[i,13:15])[2], sort(wilcox3[i,13:15])[2]/max(wilcox3[i,13:15]))<0.001)
		&& 
                (  ((p.1A==max(wilcox3[7:9]) && p.1A > max(p.12,p.13,p.AB,p.AC) &&
                     if(!is.na(min(min(p.1A,max(p.12,p.13,p.AB,p.AC))/max(p.1A,max(p.12,p.13,p.AB,p.AC)), max(p.1A,max(p.12,p.13,p.AB,p.AC))/min(p.1A,max(p.12,p.13,p.AB,p.AC))) < 0.001))
                       {min(min(p.1A,max(p.12,p.13,p.AB,p.AC))/max(p.1A,max(p.12,p.13,p.AB,p.AC)), max(p.1A,max(p.12,p.13,p.AB,p.AC))/min(p.1A,max(p.12,p.13,p.AB,p.AC))) < 0.001})
                 || (p.1B==max(wilcox3[7:9]) && p.1B > max(p.12,p.13,p.AB,p.BC) &&
                     if(!is.na(min(min(p.1B,max(p.12,p.13,p.AB,p.BC))/max(p.1B,max(p.12,p.13,p.AB,p.BC)), max(p.1B,max(p.12,p.13,p.AB,p.BC))/min(p.1B,max(p.12,p.13,p.AB,p.BC))) < 0.001))
                       {min(min(p.1B,max(p.12,p.13,p.AB,p.BC))/max(p.1B,max(p.12,p.13,p.AB,p.BC)), max(p.1B,max(p.12,p.13,p.AB,p.BC))/min(p.1B,max(p.12,p.13,p.AB,p.BC))) < 0.001})
                 || (p.1C==max(wilcox3[7:9]) && p.1C > max(p.12,p.13,p.AC,p.BC) &&
                     if(!is.na(min(min(p.1C,max(p.12,p.13,p.AC,p.BC))/max(p.1C,max(p.12,p.13,p.AC,p.BC)), max(p.1C,max(p.12,p.13,p.AC,p.BC))/min(p.1C,max(p.12,p.13,p.AC,p.BC))) < 0.001))
                       {min(min(p.1C,max(p.12,p.13,p.AC,p.BC))/max(p.1C,max(p.12,p.13,p.AC,p.BC)), max(p.1C,max(p.12,p.13,p.AC,p.BC))/min(p.1C,max(p.12,p.13,p.AC,p.BC))) < 0.001})
                   )
                 &&((p.2A==max(wilcox3[10:12]) && p.2A > max(p.12,p.23,p.AB,p.AC) &&
                     if(!is.na(min(min(p.2A,max(p.12,p.23,p.AB,p.AC))/max(p.2A,max(p.12,p.23,p.AB,p.AC)), max(p.2A,max(p.12,p.23,p.AB,p.AC))/min(p.2A,max(p.12,p.23,p.AB,p.AC))) < 0.001))
                       {min(min(p.2A,max(p.12,p.23,p.AB,p.AC))/max(p.2A,max(p.12,p.23,p.AB,p.AC)), max(p.2A,max(p.12,p.23,p.AB,p.AC))/min(p.2A,max(p.12,p.23,p.AB,p.AC))) < 0.001})
                 || (p.2B==max(wilcox3[10:12]) && p.2B > max(p.12,p.23,p.AB,p.BC) &&
                     if(!is.na(min(min(p.2B,max(p.12,p.23,p.AB,p.BC))/max(p.2B,max(p.12,p.23,p.AB,p.BC)), max(p.2B,max(p.12,p.23,p.AB,p.BC))/min(p.2B,max(p.12,p.23,p.AB,p.BC))) < 0.001))
                       {min(min(p.2B,max(p.12,p.23,p.AB,p.BC))/max(p.2B,max(p.12,p.23,p.AB,p.BC)), max(p.2B,max(p.12,p.23,p.AB,p.BC))/min(p.2B,max(p.12,p.23,p.AB,p.BC))) < 0.001})
                 || (p.2C==max(wilcox3[10:12]) && p.2C > max(p.12,p.23,p.AC,p.BC) &&
                     if(!is.na(min(min(p.2C,max(p.12,p.23,p.AC,p.BC))/max(p.2C,max(p.12,p.23,p.AC,p.BC)), max(p.2C,max(p.12,p.23,p.AC,p.BC))/min(p.2C,max(p.12,p.23,p.AC,p.BC))) < 0.001))
                       {min(min(p.2C,max(p.12,p.23,p.AC,p.BC))/max(p.2C,max(p.12,p.23,p.AC,p.BC)), max(p.2C,max(p.12,p.23,p.AC,p.BC))/min(p.2C,max(p.12,p.23,p.AC,p.BC))) < 0.001})
                   )             
                 &&((p.3A==max(wilcox3[13:15]) && p.3A > max(p.13,p.23,p.AB,p.AC) &&
                     if(!is.na(min(min(p.3A,max(p.13,p.23,p.AB,p.AC))/max(p.3A,max(p.13,p.23,p.AB,p.AC)), max(p.3A,max(p.13,p.23,p.AB,p.AC))/min(p.3A,max(p.13,p.23,p.AB,p.AC))) < 0.001))
                       {min(min(p.3A,max(p.13,p.23,p.AB,p.AC))/max(p.3A,max(p.13,p.23,p.AB,p.AC)), max(p.3A,max(p.13,p.23,p.AB,p.AC))/min(p.3A,max(p.13,p.23,p.AB,p.AC))) < 0.001})
                 || (p.3B==max(wilcox3[13:15]) && p.3B > max(p.13,p.23,p.AB,p.BC) &&
                     if(!is.na(min(min(p.3B,max(p.13,p.23,p.AB,p.BC))/max(p.3B,max(p.13,p.23,p.AB,p.BC)), max(p.3B,max(p.13,p.23,p.AB,p.BC))/min(p.3B,max(p.13,p.23,p.AB,p.BC))) < 0.001))
                       {min(min(p.3B,max(p.13,p.23,p.AB,p.BC))/max(p.3B,max(p.13,p.23,p.AB,p.BC)), max(p.3B,max(p.13,p.23,p.AB,p.BC))/min(p.3B,max(p.13,p.23,p.AB,p.BC))) < 0.001})
                 || (p.3C==max(wilcox3[13:15]) && p.3C > max(p.13,p.23,p.AC,p.BC) &&
                     if(!is.na(min(min(p.3C,max(p.13,p.23,p.AC,p.BC))/max(p.3C,max(p.13,p.23,p.AC,p.BC)), max(p.3C,max(p.13,p.23,p.AC,p.BC))/min(p.3C,max(p.13,p.23,p.AC,p.BC))) < 0.001))
                       {min(min(p.3C,max(p.13,p.23,p.AC,p.BC))/max(p.3C,max(p.13,p.23,p.AC,p.BC)), max(p.3C,max(p.13,p.23,p.AC,p.BC))/min(p.3C,max(p.13,p.23,p.AC,p.BC))) < 0.001})
                   )
		   )
	 )
       {
          wilcox3$decided[i] <- 3
          if(p.1A == max(wilcox3[i,7:9]))
            {wilcox3$assignm1[i] <- paste(wilcox3$In1[i], "->", wilcox3$Out1[i])						#Assignment
             assignments <- rbind(assignments, data.frame('Incoming'=wilcox3$In1[i], 'Outgoing'=wilcox3$Out1[i], 
					'TypeOfCollision'="3x3", 'Duration'=Duration1A, 'Forced'=0))
             needed$csv_Number[needed$csv_Number==wilcox3$In1[i]] <- wilcox3$Out1[i]					#Combine tracks (use outgoing nr)
		 collisions$n_frames[wilcox3$In1[i]] <- nrow(subset(needed,needed$csv_Number==wilcox3$In1[i]))		#Adjust frame numbers
             collisions$n_frames[wilcox3$Out1[i]] <- nrow(subset(needed,needed$csv_Number==wilcox3$Out1[i]))
		 collisions$End[wilcox3$In1[i]] <- collisions$End[wilcox3$Out1[i]]						#Adjust ending and starting frame
		 collisions$Start[wilcox3$Out1[i]] <- collisions$Start[wilcox3$In1[i]]}
          if(p.1B == max(wilcox3[i,7:9]))
            {wilcox3$assignm1[i] <- paste(wilcox3$In1[i], "->", wilcox3$Out2[i])
             assignments <- rbind(assignments, data.frame('Incoming'=wilcox3$In1[i], 'Outgoing'=wilcox3$Out2[i], 
					'TypeOfCollision'="3x3", 'Duration'=Duration1B, 'Forced'=0))
             needed$csv_Number[needed$csv_Number==wilcox3$In1[i]] <- wilcox3$Out2[i]
             collisions$n_frames[wilcox3$In1[i]] <- nrow(subset(needed,needed$csv_Number==wilcox3$In1[i]))
             collisions$n_frames[wilcox3$Out2[i]] <- nrow(subset(needed,needed$csv_Number==wilcox3$Out2[i]))
		 collisions$End[wilcox3$In1[i]] <- collisions$End[wilcox3$Out2[i]]						
		 collisions$Start[wilcox3$Out2[i]] <- collisions$Start[wilcox3$In1[i]]}
          if(p.1C == max(wilcox3[i,7:9]))
            {wilcox3$assignm1[i] <- paste(wilcox3$In1[i], "->", wilcox3$Out3[i])
             assignments <- rbind(assignments, data.frame('Incoming'=wilcox3$In1[i], 'Outgoing'=wilcox3$Out3[i], 
					'TypeOfCollision'="3x3", 'Duration'=Duration1C, 'Forced'=0))
             needed$csv_Number[needed$csv_Number==wilcox3$In1[i]] <- wilcox3$Out3[i]
             collisions$n_frames[wilcox3$In1[i]] <- nrow(subset(needed,needed$csv_Number==wilcox3$In1[i]))
             collisions$n_frames[wilcox3$Out3[i]] <- nrow(subset(needed,needed$csv_Number==wilcox3$Out3[i]))
		 collisions$End[wilcox3$In1[i]] <- collisions$End[wilcox3$Out3[i]]						
		 collisions$Start[wilcox3$Out3[i]] <- collisions$Start[wilcox3$In1[i]]}
          if(p.2A == max(wilcox3[i,10:12]))
            {wilcox3$assignm2[i] <- paste(wilcox3$In2[i], "->", wilcox3$Out1[i])
             assignments <- rbind(assignments, data.frame('Incoming'=wilcox3$In2[i], 'Outgoing'=wilcox3$Out1[i], 
					'TypeOfCollision'="3x3", 'Duration'=Duration2A, 'Forced'=0))
             needed$csv_Number[needed$csv_Number==wilcox3$In2[i]] <- wilcox3$Out1[i]
             collisions$n_frames[wilcox3$In2[i]] <- nrow(subset(needed,needed$csv_Number==wilcox3$In2[i]))
             collisions$n_frames[wilcox3$Out1[i]] <- nrow(subset(needed,needed$csv_Number==wilcox3$Out1[i]))
		 collisions$End[wilcox3$In2[i]] <- collisions$End[wilcox3$Out1[i]]						
		 collisions$Start[wilcox3$Out1[i]] <- collisions$Start[wilcox3$In2[i]]}
          if(p.2B == max(wilcox3[i,10:12]))
            {wilcox3$assignm2[i] <- paste(wilcox3$In2[i], "->", wilcox3$Out2[i])
             assignments <- rbind(assignments, data.frame('Incoming'=wilcox3$In2[i], 'Outgoing'=wilcox3$Out2[i], 
					'TypeOfCollision'="3x3", 'Duration'=Duration2B, 'Forced'=0))
             needed$csv_Number[needed$csv_Number==wilcox3$In2[i]] <- wilcox3$Out2[i]
             collisions$n_frames[wilcox3$In2[i]] <- nrow(subset(needed,needed$csv_Number==wilcox3$In2[i]))
             collisions$n_frames[wilcox3$Out2[i]] <- nrow(subset(needed,needed$csv_Number==wilcox3$Out2[i]))
		 collisions$End[wilcox3$In2[i]] <- collisions$End[wilcox3$Out2[i]]						
		 collisions$Start[wilcox3$Out2[i]] <- collisions$Start[wilcox3$In2[i]]}
          if(p.2C == max(wilcox3[i,10:12]))
            {wilcox3$assignm2[i] <- paste(wilcox3$In2[i], "->", wilcox3$Out3[i])
             assignments <- rbind(assignments, data.frame('Incoming'=wilcox3$In2[i], 'Outgoing'=wilcox3$Out3[i], 
					'TypeOfCollision'="3x3", 'Duration'=Duration2C, 'Forced'=0))
             needed$csv_Number[needed$csv_Number==wilcox3$In2[i]] <- wilcox3$Out3[i]
             collisions$n_frames[wilcox3$In2[i]] <- nrow(subset(needed,needed$csv_Number==wilcox3$In2[i]))
             collisions$n_frames[wilcox3$Out3[i]] <- nrow(subset(needed,needed$csv_Number==wilcox3$Out3[i]))
		 collisions$End[wilcox3$In2[i]] <- collisions$End[wilcox3$Out3[i]]						
		 collisions$Start[wilcox3$Out3[i]] <- collisions$Start[wilcox3$In2[i]]}
          if(p.3A == max(wilcox3[i,13:15]))
            {wilcox3$assignm3[i] <- paste(wilcox3$In3[i], "->", wilcox3$Out1[i])
             assignments <- rbind(assignments, data.frame('Incoming'=wilcox3$In3[i], 'Outgoing'=wilcox3$Out1[i], 
					'TypeOfCollision'="3x3", 'Duration'=Duration3A, 'Forced'=0))
             needed$csv_Number[needed$csv_Number==wilcox3$In3[i]] <- wilcox3$Out1[i]
             collisions$n_frames[wilcox3$In3[i]] <- nrow(subset(needed,needed$csv_Number==wilcox3$In3[i]))
             collisions$n_frames[wilcox3$Out1[i]] <- nrow(subset(needed,needed$csv_Number==wilcox3$Out1[i]))
		 collisions$End[wilcox3$In3[i]] <- collisions$End[wilcox3$Out1[i]]						
		 collisions$Start[wilcox3$Out1[i]] <- collisions$Start[wilcox3$In3[i]]}
          if(p.3B == max(wilcox3[i,13:15]))
            {wilcox3$assignm3[i] <- paste(wilcox3$In3[i], "->", wilcox3$Out2[i])
             assignments <- rbind(assignments, data.frame('Incoming'=wilcox3$In3[i], 'Outgoing'=wilcox3$Out2[i], 
					'TypeOfCollision'="3x3", 'Duration'=Duration3B, 'Forced'=0))
             needed$csv_Number[needed$csv_Number==wilcox3$In3[i]] <- wilcox3$Out2[i]
             collisions$n_frames[wilcox3$In3[i]] <- nrow(subset(needed,needed$csv_Number==wilcox3$In3[i]))
             collisions$n_frames[wilcox3$Out2[i]] <- nrow(subset(needed,needed$csv_Number==wilcox3$Out2[i]))
		 collisions$End[wilcox3$In3[i]] <- collisions$End[wilcox3$Out2[i]]						
		 collisions$Start[wilcox3$Out2[i]] <- collisions$Start[wilcox3$In3[i]]}
          if(p.3C == max(wilcox3[i,13:15]))
            {wilcox3$assignm3[i] <- paste(wilcox3$In3[i], "->", wilcox3$Out3[i])
             assignments <- rbind(assignments, data.frame('Incoming'=wilcox3$In3[i], 'Outgoing'=wilcox3$Out3[i], 
					'TypeOfCollision'="3x3", 'Duration'=Duration3C, 'Forced'=0))
             needed$csv_Number[needed$csv_Number==wilcox3$In3[i]] <- wilcox3$Out3[i]
             collisions$n_frames[wilcox3$In3[i]] <- nrow(subset(needed,needed$csv_Number==wilcox3$In3[i]))
             collisions$n_frames[wilcox3$Out3[i]] <- nrow(subset(needed,needed$csv_Number==wilcox3$Out3[i]))
		 collisions$End[wilcox3$In3[i]] <- collisions$End[wilcox3$Out3[i]]						
		 collisions$Start[wilcox3$Out3[i]] <- collisions$Start[wilcox3$In3[i]]}
       }

     array.right <- array(0,c(3,3))							#Assign matches for one larva
     dimnames(array.right) <- list(c(1,2,3),c("A","B","C"))
     if(p.1A == max(wilcox3[i,7:9])){array.right[1,1] <- 1}
     if(p.1B == max(wilcox3[i,7:9])){array.right[1,2] <- 1}
     if(p.1C == max(wilcox3[i,7:9])){array.right[1,3] <- 1}
     if(p.2A == max(wilcox3[i,10:12])){array.right[2,1] <- 1}
     if(p.2B == max(wilcox3[i,10:12])){array.right[2,2] <- 1}
     if(p.2C == max(wilcox3[i,10:12])){array.right[2,3] <- 1}
     if(p.3A == max(wilcox3[i,13:15])){array.right[3,1] <- 1}
     if(p.3B == max(wilcox3[i,13:15])){array.right[3,2] <- 1}
     if(p.3C == max(wilcox3[i,13:15])){array.right[3,3] <- 1}
     
     in.nr <- NA
     out.nr <- NA
     in.column <- NA
     in.row <- NA
     for(k in 1:3)
       {
          if(sum(array.right[,k])==1)
            {
               out.nr <- inandout$csv_Number[k+3]
               in.column <- k
            }
       }
     for(k in 1:3)
       {
          if(!is.na(in.column) && array.right[k,in.column]==1)
            {
               in.nr <- inandout$csv_Number[k]
               in.row <- k
            }
       }

     if(!is.na(in.nr) && !is.na(out.nr)
		&& !is.na(min(max(wilcox3[i,(7+(in.row-1)*3):(9+(in.row-1)*3)])/sort(wilcox3[i,(7+(in.row-1)*3):(9+(in.row-1)*3)])[2],
                              sort(wilcox3[i,(7+(in.row-1)*3):(9+(in.row-1)*3)])[2]/max(wilcox3[i,(7+(in.row-1)*3):(9+(in.row-1)*3)]))<=0.001)
       	&&
                       ((p.1A==max(wilcox3[1:3]) && p.1A > max(p.12,p.13,p.AB,p.AC) &&
                         if(!is.na(min(min(p.1A,max(p.12,p.13,p.AB,p.AC))/max(p.1A,max(p.12,p.13,p.AB,p.AC)), max(p.1A,max(p.12,p.13,p.AB,p.AC))/min(p.1A,max(p.12,p.13,p.AB,p.AC))) < 0.001))
                           {min(min(p.1A,max(p.12,p.13,p.AB,p.AC))/max(p.1A,max(p.12,p.13,p.AB,p.AC)), max(p.1A,max(p.12,p.13,p.AB,p.AC))/min(p.1A,max(p.12,p.13,p.AB,p.AC))) < 0.001})
                     || (p.2B==max(wilcox3[4:6]) && p.2B > max(p.12,p.23,p.AB,p.BC) &&
                         if(!is.na(min(min(p.2B,max(p.12,p.23,p.AB,p.BC))/max(p.2B,max(p.12,p.23,p.AB,p.BC)), max(p.2B,max(p.12,p.23,p.AB,p.BC))/min(p.2B,max(p.12,p.23,p.AB,p.BC))) < 0.001))
                           {min(min(p.2B,max(p.12,p.23,p.AB,p.BC))/max(p.2B,max(p.12,p.23,p.AB,p.BC)), max(p.2B,max(p.12,p.23,p.AB,p.BC))/min(p.2B,max(p.12,p.23,p.AB,p.BC))) < 0.001})
                     || (p.3C==max(wilcox3[7:9]) && p.3C > max(p.13,p.23,p.AC,p.BC) &&
                         if(!is.na(min(min(p.3C,max(p.13,p.23,p.AC,p.BC))/max(p.3C,max(p.13,p.23,p.AC,p.BC)), max(p.3C,max(p.13,p.23,p.AC,p.BC))/min(p.3C,max(p.13,p.23,p.AC,p.BC))) < 0.001))
                           {min(min(p.3C,max(p.13,p.23,p.AC,p.BC))/max(p.3C,max(p.13,p.23,p.AC,p.BC)), max(p.3C,max(p.13,p.23,p.AC,p.BC))/min(p.3C,max(p.13,p.23,p.AC,p.BC))) < 0.001})
                       ))
	 {
          wilcox3$decided[i] <- 1
          wilcox3$assignm1[i] <- paste(in.nr, "->", out.nr)

	    if(in.nr==inandout[1,]$csv_Number && out.nr==inandout[4,]$csv_Number){Duration <- Duration1A}
	    if(in.nr==inandout[1,]$csv_Number && out.nr==inandout[5,]$csv_Number){Duration <- Duration1B}
	    if(in.nr==inandout[1,]$csv_Number && out.nr==inandout[6,]$csv_Number){Duration <- Duration1C}
	    if(in.nr==inandout[2,]$csv_Number && out.nr==inandout[4,]$csv_Number){Duration <- Duration2A}
	    if(in.nr==inandout[2,]$csv_Number && out.nr==inandout[5,]$csv_Number){Duration <- Duration2B}
	    if(in.nr==inandout[2,]$csv_Number && out.nr==inandout[6,]$csv_Number){Duration <- Duration2C}
	    if(in.nr==inandout[3,]$csv_Number && out.nr==inandout[4,]$csv_Number){Duration <- Duration3A}
	    if(in.nr==inandout[3,]$csv_Number && out.nr==inandout[5,]$csv_Number){Duration <- Duration3B}
	    if(in.nr==inandout[3,]$csv_Number && out.nr==inandout[6,]$csv_Number){Duration <- Duration3C}

          assignments <- rbind(assignments, data.frame('Incoming'=in.nr, 'Outgoing'=out.nr, 
					'TypeOfCollision'="3x3", 'Duration'=Duration, 'Forced'=0))
          needed$csv_Number[needed$csv_Number==in.nr] <- out.nr
          collisions$n_frames[in.nr] <- nrow(subset(needed,needed$csv_Number==in.nr))
          collisions$n_frames[out.nr] <- nrow(subset(needed,needed$csv_Number==out.nr))
          collisions$End[in.nr] <- collisions$End[out.nr]
          collisions$Start[out.nr] <- collisions$Start[in.nr]
       }
  }
}
wilcox3

n_3x3_assign <- n_3x3_assign + sum(wilcox3$decided)
n_assignments <- n_assignments + sum(wilcox3$decided)


#number.3coll <- rbind(number.3coll, nrow(wilcox3))
#number.3coll.1solved <- rbind(number.3coll.1solved, sum(wilcox3$decided==1))
#number.3coll.3solved <- rbind(number.3coll.3solved, sum(wilcox3$decided==3))
#if(nrow(wilcox3)!=0)
#  {
#    write.xlsx(wilcox3, row.names =FALSE, "Wilcox.3collisions.xlsx")
#  }


#################################################
### Nx1 collisions ##############################
#################################################
#Das ist noch bei weitem nicht perfekt. Sobald ein 2-bundle in den Rand geht, 
#kann ich nichts mehr l�sen, sobald zu Beginn != 0 im Rand sind, kann ich auch 
#nichts machen... 

#cross testing added
setwd(box.path)
assigntable <- array(NA,c(0,2))
colnames(assigntable) <- c("assignment","incoming")
n_dependend = 0

for(origin in 0:-2)
{

new=1
while(new==1){
new=0

inborder <- subset(collisions, collisions$New1==origin)
inborder <- inborder[order(inborder$End),]
inborder$Way <- rep(1, nrow(inborder))
inborder$Frame <- inborder$End
outborder <- subset(collisions, collisions$Origin==origin)
outborder <- outborder[order(outborder$Start),]
outborder$Way <- rep(2, nrow(outborder))
outborder$Frame <- outborder$Start

inout <- rbind(inborder,outborder)
inout <- inout[,-c(2:9)]
inout <- inout[order(inout$Way, decreasing=TRUE),]
inout <- inout[order(inout$Frame),]

inside <- NULL

if(sum(rbind(inborder,outborder)$Bundle)==nrow(inout)){

if(nrow(inout)>0)
{
for(i in 1:nrow(inout))
  {
     if(inout$Way[i]==1)
       {inside <- rbind(inside,inout[i,])}
     if(inout$Way[i]==2 && is.null(inside)==FALSE)
	 {
          if(nrow(inside)==1)
            {
		assigntable <- rbind(assigntable, c(paste(inside$csv_Number[1], "->", inout$csv_Number[i]),1))
		Duration1x1 <- inout$Frame[i] - inside$Frame[1] - 1
            assignments <- rbind(assignments, data.frame('Incoming'=inside$csv_Number[1], 'Outgoing'=inout$csv_Number[i], 
					'TypeOfCollision'="1x1", 'Duration'=Duration1x1, 'Forced'=0))
		needed$csv_Number[needed$csv_Number==inside$csv_Number[1]] <- inout$csv_Number[i]
      	collisions$n_frames[inside$csv_Number[1]] <- nrow(subset(needed,needed$csv_Number==inside$csv_Number[1]))
        	collisions$n_frames[inout$csv_Number[i]] <- nrow(subset(needed,needed$csv_Number==inout$csv_Number[i]))
		collisions$Origin[inout$csv_Number[i]] <- collisions$Origin[inside$csv_Number[1]]
		collisions$New1[inside$csv_Number[1]] <- collisions$New1[inout$csv_Number[i]]
        	collisions$End[inside$csv_Number[1]] <- subset(collisions, collisions$csv_Number==inout$csv_Number[i])$End	
        	collisions$Start[inout$csv_Number[i]] <- subset(collisions, collisions$csv_Number==inside$csv_Number[1])$Start
		collisions[inside$csv_Number[1],] <- c(inside$csv_Number[1],NA,NA,NA,NA,1,0,NA,NA)
		if(inout$Way[i-1]==1){n_1x1_assign = n_1x1_assign+1}
		else{n_dependend = n_dependend+1}
		inside <- inside[-c(1),]
		new=1
		}
	    if(nrow(inside)>=2 && nrow(inside)<=5)
            {
            p.vals = rep(0,5)
            retest.p = rep(0,10)

		in1.frames = (subset(inborder, inborder$csv_Number==inside$csv_Number[1])$n_frames)
		Duration1 = inout$Frame[i] - inside$Frame[1] - 1
		in2.frames = (subset(inborder, inborder$csv_Number==inside$csv_Number[2])$n_frames)
		Duration2 = inout$Frame[i] - inside$Frame[2] - 1
		out.frames = (subset(outborder, outborder$csv_Number==inout$csv_Number[i])$n_frames)
		max.frames = min(in1.frames,in2.frames,out.frames,45)

		in3.data = NULL
		if(nrow(inside)>=3){
                 in3.frames = (subset(inborder, inborder$csv_Number==inside$csv_Number[3])$n_frames)
		     Duration3 = inout$Frame[i] - inside$Frame[3] - 1
                 max.frames = min(max.frames,in3.frames)
		  }
		in4.data = NULL
		if(nrow(inside)>=4){
                 in4.frames = (subset(inborder, inborder$csv_Number==inside$csv_Number[4])$n_frames)
		     Duration4 = inout$Frame[i] - inside$Frame[4] - 1
                 max.frames = min(max.frames,in4.frames)
		  }
		in5.data = NULL
		if(nrow(inside)>=5){
                 in5.frames = (subset(inborder, inborder$csv_Number==inside$csv_Number[5])$n_frames)
		     Duration5 = inout$Frame[i] - inside$Frame[5] - 1
                 max.frames = min(max.frames,in5.frames)
		  }

	    if(max.frames >=10)
	    {
	      in1.data <- subset(needed, needed$csv_Number==inside$csv_Number[1])[(in1.frames-max.frames+1):(in1.frames),]
            in2.data <- subset(needed, needed$csv_Number==inside$csv_Number[2])[(in2.frames-max.frames+1):(in2.frames),]
            out.data <- subset(needed, needed$csv_Number==inout$csv_Number[i])[1:max.frames,]
		if(nrow(inside)>=3){
            in3.data = subset(needed, needed$csv_Number==inside$csv_Number[3])[(in3.frames-max.frames+1):(in3.frames),]
           				 }
		if(nrow(inside)>=4){
		in4.data = subset(needed, needed$csv_Number==inside$csv_Number[4])[(in4.frames-max.frames+1):(in4.frames),]
					 }
		if(nrow(inside)>=5){
 		in5.data = subset(needed, needed$csv_Number==inside$csv_Number[5])[(in5.frames-max.frames+1):(in5.frames),]
					 }

		datas <- rbind(in1.data,in2.data,in3.data,in4.data,in5.data,out.data)

		if(kruskal.test(PrComp1~csv_Number, data=datas)$p.value <=0.05)
		  {
		      p.1A <- wilcox.test(PrComp1~csv_Number, data=rbind(in1.data,out.data))$p.value
                  p.2A <- wilcox.test(PrComp1~csv_Number, data=rbind(in2.data,out.data))$p.value
			p.3A <- 0
			if(nrow(inside)>=3){
			p.3A <- wilcox.test(PrComp1~csv_Number, data=rbind(in3.data,out.data))$p.value
						 }
			p.4A <- 0
			if(nrow(inside)>=4){
			p.4A <- wilcox.test(PrComp1~csv_Number, data=rbind(in4.data,out.data))$p.value
						 }
			p.5A <- 0
			if(nrow(inside)>=5){
			p.5A <- wilcox.test(PrComp1~csv_Number, data=rbind(in5.data,out.data))$p.value
						 }
			wilcox.nx1 <- c(p.1A, p.2A, p.3A, p.4A, p.5A)

			retest.p[1] <- wilcox.test(PrComp1~csv_Number, data=rbind(in1.data,in2.data))$p.value
			if(nrow(inside)>=3){
			retest.p[2] <- wilcox.test(PrComp1~csv_Number, data=rbind(in1.data,in3.data))$p.value
			retest.p[3] <- wilcox.test(PrComp1~csv_Number, data=rbind(in2.data,in3.data))$p.value
						 }
			if(nrow(inside)>=4){
			retest.p[4] <- wilcox.test(PrComp1~csv_Number, data=rbind(in1.data,in4.data))$p.value
			retest.p[5] <- wilcox.test(PrComp1~csv_Number, data=rbind(in2.data,in4.data))$p.value
			retest.p[6] <- wilcox.test(PrComp1~csv_Number, data=rbind(in3.data,in4.data))$p.value
						 }
			if(nrow(inside)>=5){
			retest.p[7] <- wilcox.test(PrComp1~csv_Number, data=rbind(in1.data,in5.data))$p.value
			retest.p[8] <- wilcox.test(PrComp1~csv_Number, data=rbind(in2.data,in5.data))$p.value
			retest.p[9] <- wilcox.test(PrComp1~csv_Number, data=rbind(in3.data,in5.data))$p.value
			retest.p[10] <- wilcox.test(PrComp1~csv_Number, data=rbind(in4.data,in5.data))$p.value
						 }
			if(wilcox.nx1[1]==max(wilcox.nx1)){retest.max=max(retest.p[1],retest.p[2],retest.p[4],retest.p[7])}
			if(wilcox.nx1[2]==max(wilcox.nx1)){retest.max=max(retest.p[1],retest.p[3],retest.p[5],retest.p[8])}
			if(wilcox.nx1[3]==max(wilcox.nx1)){retest.max=max(retest.p[2],retest.p[3],retest.p[6],retest.p[9])}
			if(wilcox.nx1[4]==max(wilcox.nx1)){retest.max=max(retest.p[4],retest.p[5],retest.p[6],retest.p[10])}
			if(wilcox.nx1[5]==max(wilcox.nx1)){retest.max=max(retest.p[7],retest.p[8],retest.p[9],retest.p[10])}

			if(wilcox.nx1[1]==max(wilcox.nx1)){Duration = Duration1}
			if(wilcox.nx1[2]==max(wilcox.nx1)){Duration = Duration2}
			if(wilcox.nx1[3]==max(wilcox.nx1)){Duration = Duration3}
			if(wilcox.nx1[4]==max(wilcox.nx1)){Duration = Duration4}
			if(wilcox.nx1[5]==max(wilcox.nx1)){Duration = Duration5}
		  }
		else
		  {
		  wilcox.nx1 = c(0,0,0,0,0)
		  retest.p = c(0,0,0,0,0,0,0,0,0,0)
		  }
		
		if(nrow(inside)==2)
		  {
		    if(((max(wilcox.nx1[1:2])!=0)          
			&& (wilcox.nx1[1]!=wilcox.nx1[2])
                  && !is.na(min(min(wilcox.nx1[1:2])/max(wilcox.nx1[1:2]), max(wilcox.nx1[1:2])/min(wilcox.nx1[1:2])) < 0.001)
                  && if(!is.na(min(min(wilcox.nx1[1:2])/max(wilcox.nx1[1:2]), max(wilcox.nx1[1:2])/min(wilcox.nx1[1:2])) < 0.001))
                       {min(min(wilcox.nx1[1:2])/max(wilcox.nx1[1:2]), max(wilcox.nx1[1:2])/min(wilcox.nx1[1:2])) < 0.001})
                  && 
                    ( 
                      (max(wilcox.nx1[1:2])>retest.p[1])
                      && if(!is.na(min(min(max(wilcox.nx1[1:2]),retest.p[1])/max(max(wilcox.nx1[1:2]),retest.p[1]), max(max(wilcox.nx1[1:2]),retest.p[1])/min(max(wilcox.nx1[1:2]),retest.p[1])) < 0.001))
                       {min(min(max(wilcox.nx1[1:2]),retest.p[1])/max(max(wilcox.nx1[1:2]),retest.p[1]), max(max(wilcox.nx1[1:2]),retest.p[1])/min(max(wilcox.nx1[1:2]),retest.p[1])) < 0.001}
                    )
                   )
		     {
			 for(j in 1:2){
			 if(wilcox.nx1[j]==max(wilcox.nx1)){
				assigntable = rbind(assigntable, c(paste(inside$csv_Number[j], "->", inout$csv_Number[i]), 2))
         		      assignments <- rbind(assignments, data.frame('Incoming'=inside$csv_Number[j], 'Outgoing'=inout$csv_Number[i], 
					'TypeOfCollision'="Nx1", 'Duration'=Duration, 'Forced'=0))
			 	needed$csv_Number[needed$csv_Number==inside$csv_Number[j]] <- inout$csv_Number[i]
      			collisions$n_frames[inside$csv_Number[j]] <- nrow(subset(needed,needed$csv_Number==inside$csv_Number[j]))
        			collisions$n_frames[inout$csv_Number[i]] <- nrow(subset(needed,needed$csv_Number==inout$csv_Number[i]))
				collisions$Origin[inout$csv_Number[i]] <- collisions$Origin[inside$csv_Number[j]]
				collisions$New1[inside$csv_Number[j]] <- collisions$New1[inout$csv_Number[i]]
        			collisions$End[inside$csv_Number[j]] <- subset(collisions, collisions$csv_Number==inout$csv_Number[i])$End	
        			collisions$Start[inout$csv_Number[i]] <- subset(collisions, collisions$csv_Number==inside$csv_Number[j])$Start
				collisions[inside$csv_Number[j],] <- c(inside$csv_Number[j],NA,NA,NA,NA,1,0,NA,NA)
				inside <- inside[-c(j),]
				new=1
								 }
                                }
		     }
		  }

		if(nrow(inside)>=3)
		  {
		    if(((max(wilcox.nx1[1:nrow(inside)])!=0)
                  && !is.na(min(sort(wilcox.nx1,decreasing=TRUE)[2]/max(wilcox.nx1[1:nrow(inside)]), max(wilcox.nx1[1:nrow(inside)])/sort(wilcox.nx1,decreasing=TRUE)[2]) < 0.001)
                  && if(!is.na(min(sort(wilcox.nx1,decreasing=TRUE)[2]/max(wilcox.nx1[1:nrow(inside)]), max(wilcox.nx1[1:nrow(inside)])/sort(wilcox.nx1,decreasing=TRUE)[2]) < 0.001))
                       {min(sort(wilcox.nx1,decreasing=TRUE)[2]/max(wilcox.nx1[1:nrow(inside)]), max(wilcox.nx1[1:nrow(inside)])/sort(wilcox.nx1,decreasing=TRUE)[2]) < 0.001})
                  && 
                    (
                       (max(wilcox.nx1[1:nrow(inside)])>retest.max)
                       && if(!is.na(min(min(max(wilcox.nx1[1:nrow(inside)]),retest.max)/max(max(wilcox.nx1[1:nrow(inside)]),retest.max), max(max(wilcox.nx1[1:nrow(inside)]),retest.max)/min(max(wilcox.nx1[1:nrow(inside)]),retest.max)) < 0.001))
                       {min(min(max(wilcox.nx1[1:nrow(inside)]),retest.max)/max(max(wilcox.nx1[1:nrow(inside)]),retest.max), max(max(wilcox.nx1[1:nrow(inside)]),retest.max)/min(max(wilcox.nx1[1:nrow(inside)]),retest.max)) < 0.001}    )
			  )
		     {
			 for(j in 1:nrow(inside)){
			 if(wilcox.nx1[j]==max(wilcox.nx1)){
				assigntable = rbind(assigntable, c(paste(inside$csv_Number[j], "->", inout$csv_Number[i]),nrow(inside)))
         		      assignments <- rbind(assignments, data.frame('Incoming'=inside$csv_Number[j], 'Outgoing'=inout$csv_Number[i], 
					'TypeOfCollision'="Nx1", 'Duration'=Duration, 'Forced'=0))
				needed$csv_Number[needed$csv_Number==inside$csv_Number[j]] <- inout$csv_Number[i]
        			collisions$n_frames[inside$csv_Number[j]] <- nrow(subset(needed,needed$csv_Number==inside$csv_Number[j]))
        			collisions$n_frames[inout$csv_Number[i]] <- nrow(subset(needed,needed$csv_Number==inout$csv_Number[i]))
				collisions$Origin[inout$csv_Number[i]] <- collisions$Origin[inside$csv_Number[j]]
				collisions$New1[inside$csv_Number[j]] <- collisions$New1[inout$csv_Number[i]]
        			collisions$End[inside$csv_Number[j]] <- subset(collisions, collisions$csv_Number==inout$csv_Number[i])$End	
        			collisions$Start[inout$csv_Number[i]] <- subset(collisions, collisions$csv_Number==inside$csv_Number[j])$Start
				collisions[inside$csv_Number[j],] <- c(inside$csv_Number[j],NA,NA,NA,NA,1,0,NA,NA)
				inside <- inside[-c(j),]
				new=1
								}
                               		 }
		     }
		  }

	    } #max.frames>=10
	  } #mehr als 1 in border
	} #larve verl�sst border
  } #inout zeile loop
} #if inout>0
} #keine 2er bundles gehen rein!
} #while
} #for loop origin

assigntable <- as.data.frame(assigntable)
n_nx1_assign <- n_nx1_assign + nrow(assigntable)
n_assignments <- n_assignments + nrow(assigntable)

#numbers only for approximation of quality, 2x1 includes dependend 1x1 cases:
n_2x1_assign <- n_2x1_assign + nrow(subset(assigntable,assigntable$incoming==2))+ n_dependend
n_3x1_assign <- n_3x1_assign + nrow(subset(assigntable,assigntable$incoming==3))
n_4x1_assign <- n_4x1_assign + nrow(subset(assigntable,assigntable$incoming==4))
n_5x1_assign <- n_5x1_assign + nrow(subset(assigntable,assigntable$incoming==5))


#####################################################
##### ALL ASSIGNMENTS ###############################
#####################################################
ordered <- order(assignments$Outgoing)
assignments <- assignments[ordered,]

setwd(log.path)
write.csv(assignments, row.names =FALSE, "assignments.csv")

#####################################################
##### MERGE TRACKS ##################################
#####################################################
setwd(box.path)
collisionframe <- rep(0,78)
for(i in 1:nrow(assignments)){
	IN <- read.csv(paste(assignments$Incoming[i], ".csv", sep=""),header=FALSE)
	OUT <- read.csv(paste(assignments$Outgoing[i], ".csv", sep=""),header=FALSE)
	Duration = assignments$Duration[i]
	NAs <- as.data.frame(matrix(data=NA,nrow=Duration,ncol=78))

	NEW <- rbind(IN, NAs)
	NEW <- rbind(NEW, OUT)

	unlink(paste(assignments$Incoming[i], ".csv", sep=""), recursive = FALSE, force = FALSE)
	unlink(paste(assignments$Outgoing[i], ".csv", sep=""), recursive = FALSE, force = FALSE)		

	write.table(NEW, row.names=FALSE, col.names=FALSE, sep=",", paste(assignments$Outgoing[i], ".csv", sep=""))
}


if(FALSE){
########################
########################
#Get the number of identified individual larvae (tracks >= half video)

halfvideo <- halfvideo + nrow(subset(collisions,collisions$n_frames>=1440))

setwd(log.path)
stdout <- scan("stdout.log", character(0), sep = "\n") # separate each line		#load stdout.log
line <- grep("Second Pass: Updated Larvae on Dish guess:", stdout, fixed=T) 		#get number of line with guess
guess <- rbind(guess,as.integer(unlist(strsplit(stdout[line],"guess: "))[2]))		#get guess
allguess <- sum(guess)

line <- grep("RESOLVED_BY_MODEL:", stdout, fixed=T)[1] 						#get number of line with guess
model2x2 <- rbind(model2x2,as.integer(unlist(strsplit(stdout[line],"RESOLVED_BY_MODEL:"))[2]))			#get guess
allmodel <- sum(model2x2)

justNA <- subset(collisions, is.na(collisions$Origin)==TRUE)
single <- rbind(single, sum(subset(collisions, collisions$Bundle==1)$Bundle))
allsingle <- sum(single)

#halfvideo: number of identified individuals
#allguess: number of guessed individuals (should be equal for every approach)
#allmodel: number of solved collisions by the model
#allsingle: number of single larvae tracks that do not originate from border or odor cup

######################################################
### NUMBERS FOR COLLISIONS (2-bundle-collisions) #####
######################################################
#number of 2-larvae-bundles
n.two.bundle <- nrow(subset(collisions, collisions$Bundle==2))
#number of 2-bundle -> 3-bundle
two.to.three <- subset(subset(collisions, collisions$Bundle==2), is.na(subset(collisions, collisions$Bundle==2)$New2))
n.two.to.three <- nrow(subset(subset(collisions, collisions$Bundle==2), is.na(subset(collisions, collisions$Bundle==2)$New2)))
#number of 2-bundles with noch previous single-larvae data without '2->3'-cases
n.two.at.start <- 0
for(i in collisions$csv_Number)
  {
     
     if(collisions$Bundle[i]==2 && nrow(subset(collisions, collisions$New1==i))==0 
       && nrow(subset(two.to.three, two.to.three$csv_Number==i))==0)
       {
          n.two.at.start <- n.two.at.start + 1
       }
  }

#number of all '1+1 -> 2 -> 1+1'-collisions
gesamt <- gesamt + n.two.bundle - n.two.to.three - n.two.at.start
}


# }	
#####END OF LOOP THROUGH EVERY BOX#####
#Sys.time() - start

#Results are given by:
#	'assigntable' for Nx1
#     'wilcox.p' for 2x2
#	'wilcox3' for 3x3
#	'assignments' gives all assignments with type and duration of collision
