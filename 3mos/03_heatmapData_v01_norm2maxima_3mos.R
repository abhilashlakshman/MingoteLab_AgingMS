rm(list=ls())

library(plotly)

# IMPORT ORIGINAL IMAGE SIZES FROM TWO SOURCES (CONFOCAL AND EPI)
org.im.det.czi <- read.delim("~/OneDrive - Indiana University/ShaferLab/Others/Jackie/INTRSECT Ant.-Post._Medial-Lateral/Medial-Lateral/3 mos_age/Analysed/imageDetails_3-VTA czi files.txt", header = T)
org.im.det.lif <- read.delim("~/OneDrive - Indiana University/ShaferLab/Others/Jackie/INTRSECT Ant.-Post._Medial-Lateral/Medial-Lateral/3 mos_age/Analysed/imageDetails_2-VTA lif files.txt", header = T)

# MERGE THE TWO ORIGINAL IMAGE SIZE FILES
org.im.det <- rbind(org.im.det.czi, org.im.det.lif)

# IMPORT IMAGE SIZES FROM FILES THAT WERE USED FOR COUNTING
counted.im.det <- read.delim("~/OneDrive - Indiana University/ShaferLab/Others/Jackie/INTRSECT Ant.-Post._Medial-Lateral/Medial-Lateral/3 mos_age/Analysed/imageDetails_1-VTA_counted.txt", header = T)
counted.im.det$org.px.size <- org.im.det$px.size

# IMPORT THE BREGMA DETAILS FILE
pre.bregma.det <- read.delim("~/OneDrive - Indiana University/ShaferLab/Others/Jackie/INTRSECT Ant.-Post._Medial-Lateral/Medial-Lateral/3 mos_age/3mos_bregma.txt", header = T)
bregma <- data.frame(
  "bregma" = c("-3.7",	"-3.6",	"-3.5",	"-3.4",	"-3.3",	"-3.2")
)
bregma.det <- cbind(bregma, pre.bregma.det)

# IMPORT THE COUNTS FILES' NAMES
setwd("~/OneDrive - Indiana University/ShaferLab/Others/Jackie/INTRSECT Ant.-Post._Medial-Lateral/Medial-Lateral/3 mos_age/1-VTA_counted tifs/output_counts/")
files <- list.files()
vta.file.ls <- files[grep('_vta.tif', files)]

# SEPARATE EACH ANIMAL AND BREGMA INTO A DIFFERENT A HIERARCHICAL LIST
animal.ls <- list()
animal.ids <- colnames(bregma.det)[-1]
file.names.by.animals <- list()
for (i in 1:length(animal.ids)) {
  file.names.by.animals[[i]] <- vta.file.ls[grep(animal.ids[i], vta.file.ls)]
}
names(file.names.by.animals) <- animal.ids
file.ind <- list()
for (i in 1:length(file.names.by.animals)) {
  per.animal.bregma <- list()
  for (j in 1:length(bregma.det$bregma)) {
    if (is.na(bregma.det[j,i+1])) {
      per.animal.bregma[[j]] <- NA
    } else {
      x <- file.names.by.animals[[i]]
      y <- x[grep(bregma.det[j,i+1], x, fixed = T)]
      if (length(y) > 1) {
        y <- y[1]
      }
      dat <- read.delim(y, header = T)
      nam.split.to.match.with.image.det <- strsplit(y, split = ".tif.txt")[[1]][1]
      sub.counted.im.det <- counted.im.det[grep(nam.split.to.match.with.image.det, counted.im.det$file.name, fixed = F),]
      
      if (dim(sub.counted.im.det)[1] == 0) {
        dat$X.corr <- NA
        dat$Y.corr <- NA
      } else {
        if (sub.counted.im.det$px.size == 1) {
          dat$X.corr <- dat$X * sub.counted.im.det$org.px.size
          dat$Y.corr <- dat$Y * sub.counted.im.det$org.px.size
        } else {
          dat$X.corr <- (dat$X * (sub.counted.im.det$px.size / sub.counted.im.det$org.px.size)) / 25400
          dat$Y.corr <- (dat$Y * (sub.counted.im.det$px.size / sub.counted.im.det$org.px.size)) / 25400
        }
      }
      
      per.animal.bregma[[j]] <- dat
    }
  }
  names(per.animal.bregma) <- bregma.det$bregma
  file.ind[[i]] <- per.animal.bregma
}
names(file.ind) <- animal.ids

# COMPUTE DISTANCE FROM THE CENTER TO IDENTIFY MEDIAL-LATERAL DISTRIBUTION
final.out.da <- list()
final.out.daglu <- list()
final.out.TH <- list() ########################################################
for (i in 1:length(file.ind)) {
  final.out.sub1.da <- list()
  final.out.sub1.daglu <- list()
  final.out.sub1.TH <- list() #################################################
  x <- file.ind[[i]]
  for (j in 1:length(x)) {
    if (length(x[[j]]) == 1) {
      final.out.sub2 <- NA
      final.out.sub1.da[[j]] <- final.out.sub2
      final.out.sub1.daglu[[j]] <- final.out.sub2
      final.out.sub1.TH[[j]] <- final.out.sub2 ################################
    } else {
      pre.pre.df <- x[[j]]
      
      pre.df <- rbind(
        subset(pre.pre.df, pre.pre.df$Counter == 2),
        subset(pre.pre.df, pre.pre.df$Counter == 3),
        subset(pre.pre.df, pre.pre.df$Counter == 1), ##########################
        subset(pre.pre.df, pre.pre.df$Counter == 7), ##########################
        subset(pre.pre.df, pre.pre.df$Counter == 50),
        subset(pre.pre.df, pre.pre.df$Counter == 51)
      )
      
      df <- rbind(
        subset(pre.df, pre.df$Counter == 2),
        subset(pre.df, pre.df$Counter == 3),
        subset(pre.df, pre.df$Counter == 1), ##################################
        subset(pre.df, pre.df$Counter == 7) ###################################
      )
      
      new.dat <- rbind(
        subset(pre.df, pre.df$Counter == 50),
        subset(pre.df, pre.df$Counter == 51)
      )
      
      if (abs(new.dat$X.corr[1] - new.dat$X.corr[2]) < abs(new.dat$Y.corr[1] - new.dat$Y.corr[2])) {
        m <- ((new.dat[2,"Y.corr"] - new.dat[1,"Y.corr"])/(new.dat[2,"X.corr"] - new.dat[1,"X.corr"]))
        c <- new.dat[1,"Y.corr"] - m * (new.dat[1,"X.corr"])
      } else {
        m <- ((new.dat[2,"X.corr"] - new.dat[1,"X.corr"])/(new.dat[2,"Y.corr"] - new.dat[1,"Y.corr"]))
        c <- new.dat[1,"X.corr"] - m * (new.dat[1,"Y.corr"])
      }
      
      out.diff.vta <- matrix(NA, nrow = length(df[,1]), ncol = 1)
      out.coord.vta <- matrix(NA, nrow = length(df[,1]), ncol = 1)
      colnames(out.coord.vta) <- c("X.corr")
      
      
      for (ii in 1:length(df[,1])) {
        out.diff.vta[ii,1] <- df[ii,"X.corr"] - ((df[ii,"Y.corr"] - c)/(m))
        out.coord.vta[ii,"X.corr"] <- ((df[ii,"Y.corr"] - c)/(m))
      }
      
      df$diff <- out.diff.vta
      
      vta.out.dist.cent.DA <- matrix(NA, nrow = 1000, ncol = 1)
      vta.out.dist.cent.DAGLU <- matrix(NA, nrow = 1000, ncol = 1)
      vta.out.dist.cent.allTH <- matrix(NA, nrow = 1000, ncol = 1) ############
      
      for (k in 1:length(out.diff.vta)) {
        if (df$Counter[k] == 2) {
          vta.out.dist.cent.DA[k,1] <- abs(out.diff.vta[k])
        } else if (df$Counter[k] == 3) {
          vta.out.dist.cent.DAGLU[k,1] <- abs(out.diff.vta[k])
        }
      }
      
      vta.out.dist.cent.allTH[1:length(out.diff.vta),1] <- abs(out.diff.vta) #######
      
      final.out.sub1.da[[j]] <- data.frame(
        "DA" = vta.out.dist.cent.DA
      )
      final.out.sub1.daglu[[j]] <- data.frame(
        "DAGLU" = vta.out.dist.cent.DAGLU
      )
      final.out.sub1.TH[[j]] <- data.frame( ###################################
        "allTH" = vta.out.dist.cent.allTH #####################################
      ) #######################################################################
    }
  }
  names(final.out.sub1.da) <- bregma.det$bregma
  names(final.out.sub1.daglu) <- bregma.det$bregma
  names(final.out.sub1.TH) <- bregma.det$bregma ###############################
  final.out.da[[i]] <- final.out.sub1.da
  final.out.daglu[[i]] <- final.out.sub1.daglu
  final.out.TH[[i]] <- final.out.sub1.TH ######################################
}
names(final.out.da) <- animal.ids
names(final.out.daglu) <- animal.ids
names(final.out.TH) <- animal.ids #############################################


# CREATE ANIMAL-WISE HEATMAP DATA
out.ls.da <- list()
out.ls.daglu <- list()
out.ls.TH <- list() ###########################################################

for (i in 1:length(final.out.da)) { 
  da <- final.out.da[[i]]
  daglu <- final.out.daglu[[i]]
  th <- final.out.TH[[i]] #####################################################
  
  da.mat <- matrix(NA, nrow = length(seq(0.05, 1, by = 0.1)), ncol = length(bregma.det$bregma))
  rownames(da.mat) <- seq(0.05, 1, by = 0.1)
  colnames(da.mat) <- bregma.det$bregma
  
  daglu.mat <- matrix(NA, nrow = length(seq(0.05, 1, by = 0.1)), ncol = length(bregma.det$bregma))
  rownames(daglu.mat) <- seq(0.05, 1, by = 0.1)
  colnames(daglu.mat) <- bregma.det$bregma
  
  TH.mat <- matrix(NA, nrow = length(seq(0.05, 1, by = 0.1)), ncol = length(bregma.det$bregma)) #######
  rownames(TH.mat) <- seq(0.05, 1, by = 0.1) ##################################
  colnames(TH.mat) <- bregma.det$bregma #######################################
  
  for (j in 1:length(da)) {
    x.da <- da[[j]]/max(th[[j]], na.rm = T) ###################################
    x.daglu <- daglu[[j]]/max(th[[j]], na.rm = T) #############################
    x.th <- th[[j]]/max(th[[j]], na.rm = T) ###################################
    
    hist.x.da <- hist(x.da[,1], breaks = seq(0, 1, by = 0.1), plot = F)
    hist.x.daglu <- hist(x.daglu[,1], breaks = seq(0, 1, by = 0.1), plot = F)
    hist.x.TH <- hist(x.th[,1], breaks = seq(0, 1, by = 0.1), plot = F) #######
    
    da.mat[,j] <- hist.x.da$counts
    daglu.mat[,j] <- hist.x.daglu$counts
    TH.mat[,j] <- hist.x.TH$counts ############################################
  }
  
  out.ls.da[[i]] <- da.mat
  out.ls.daglu[[i]] <- daglu.mat
  out.ls.TH[[i]] <- TH.mat ####################################################
}

names(out.ls.da) <- animal.ids
names(out.ls.daglu) <- animal.ids
names(out.ls.TH) <- animal.ids ################################################

# EXPORT RDS FILES FOR ANIMAL-WISE HEATMAPS
setwd("~/OneDrive - Indiana University/ShaferLab/Others/Jackie/INTRSECT Ant.-Post._Medial-Lateral/Medial-Lateral/HeatmapData/IndData/")
saveRDS(out.ls.da, "Normalized_DA_3mos.rds")
saveRDS(out.ls.daglu, "Normalized_DAGLU_3mos.rds")
saveRDS(out.ls.TH, "Normalized_allTH_3mos.rds") ###############################