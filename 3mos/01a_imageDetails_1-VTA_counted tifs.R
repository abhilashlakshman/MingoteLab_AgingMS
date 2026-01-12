rm(list=ls())

setwd("~/OneDrive - Indiana University/ShaferLab/Others/Jackie/INTRSECT Ant.-Post._Medial-Lateral/Medial-Lateral/3 mos_age/1-VTA_counted tifs/output_image/")
file.ls <- list.files()

file.names <- c()

for (i in 1:length(file.ls)) {
  x <- strsplit(file.ls[[i]], split = ".txt")
  file.names[i] <- x[[1]]
}

out.dat <- as.data.frame(matrix(NA, nrow = length(file.ls), ncol = 5))
colnames(out.dat) <- c("file.name", "width.um", "height.um", "width.px", "height.px")

for (i in 1:length(file.ls)) {
  out.dat$file.name[i] <- file.names[i]
  
  df <- read.delim(file.ls[[i]], header = F, skip = 3, sep = " ")
  
  sub.df.width <- subset(df, df$V1 == "Width:", select = c("V3", "V5"))
  sub.df.height <- subset(df, df$V1 == "Height:", select = c("V3", "V5"))
  
  width.info <- as.numeric(sub.df.width[1,1])
  height.info <- as.numeric(sub.df.height[1,1])
  
  width.info.px <- as.numeric(strsplit(sub.df.width$V5, split = "\\(|\\)")[[1]][2])
  height.info.px <- as.numeric(strsplit(sub.df.height$V5, split = "\\(|\\)")[[1]][2])
  
  out.dat$width.um[i] <- width.info
  out.dat$height.um[i] <- height.info
  out.dat$width.px[i] <- width.info.px
  out.dat$height.px[i] <- height.info.px
}

out.dat$px.size <- out.dat$width.um/out.dat$width.px

write.table(out.dat, "~/OneDrive - Indiana University/ShaferLab/Others/Jackie/INTRSECT Ant.-Post._Medial-Lateral/Medial-Lateral/3 mos_age/Analysed/imageDetails_1-VTA_counted.txt",
            sep = "\t", row.names = F, quote = F)
