# Mingling and differentiation based on spatstat. Updated on 13.02.2021.

rm(list = ls())										

dataPath <- "/Users/arng0001/Dropbox/Gwaith/Cyhoeddi/IndexCorrelations/Data/"
resultPath <- "/Users/arng0001/Dropbox/Gwaith/Cyhoeddi/IndexCorrelations/Results/"

# install.packages("spatstat", dep = T)
library(spatstat)
library(ggplot2)
# install.packages("metafor", dep = T)
library(metafor)

# Strange behaviour for Knysna72_116.txt
filenames <- c("Knysna72_116.txt", "Ja.txt", "Jb.txt", "Jc.txt", "Jd.txt", "JSb.txt", "XSa.txt", "XSb.txt", "XSc.txt", 
               "XSd.txt", "Da.txt", "Db.txt", "Leye.txt")
XiaohongNamesZH <- c()
for (z in 1 : 12) {
  if(z < 10)
    XiaohongNamesZH <- c(XiaohongNamesZH, paste("ZH0", z, "_2013.csv", sep = ""))
  else
    XiaohongNamesZH <- c(XiaohongNamesZH, paste("ZH", z, "_2013.csv", sep = ""))
}

filenames <- c(filenames, XiaohongNamesZH)
XiaohongNamesYLK <- c()
for (z in 1 : 12) {
  if(z < 10)
    XiaohongNamesYLK <- c(XiaohongNamesYLK, paste("YLK0", z, "_2013.csv", sep = ""))
  else
    XiaohongNamesYLK <- c(XiaohongNamesYLK, paste("YLK", z, "_2013.csv", sep = ""))
}

filenames <- c(filenames, XiaohongNamesYLK)

# HongxiangNamesGuangxi <- c()
# for (z in 1 : 25) 
#     HongxiangNamesGuangxi <- c(HongxiangNamesGuangxi, paste("plot", z, ".txt", sep = ""))
# 
# filenames <- c(filenames, HongxiangNamesGuangxi)

# ArneNames <- c("DurangoData.txt", "Bov14btrans.txt", "PenyrAlltGanol.dat", "Soe55b.dat", "Pufferzone1.dat", "15B1.dat", "Walsdorf1.dat", "BialowiezaH.txt", "P1RG.txt")
# ArneNames <- c("DurangoData.txt", "PenyrAlltGanol.dat", "Pufferzone1.dat")

# myData <- read.table(paste(dataPath, ArneNames[3], sep = ""), header = TRUE)
# names(myData)
# sort(table(myData$Species))
# sort(table(myData$species))
# filenames <- c(filenames, ArneNames)

xmax <- c(116, 100, 100, 100, 100, 100, 70, 70, 61, 60, 90, 100, 202, rep(100, length(XiaohongNamesZH)), rep(100, length(XiaohongNamesYLK)))#, rep(100, length(HongxiangNamesGuangxi)))#, 50, 103.68, 100)
ymax <- c(116, 100, 100, 100, 100, 50, 70, 70, 50, 60, 110, 80, 86, rep(100, length(XiaohongNamesZH)), rep(100, length(XiaohongNamesYLK)))#, rep(100, length(HongxiangNamesGuangxi)))#, 50, 102, 67)
rm(XiaohongNamesYLK, XiaohongNamesZH)#, ArneNames) # HongxiangNamesGuangxi

# filenames <- filenames[-18] # ZH06 has a missing dbh
# xmax <- xmax[-18] # ZH06 has a missing dbh
# ymax <- ymax[-18] # ZH06 has a missing dbh
  
getEuclideanDistance <- function(xmax, ymax, x1, y1, x2, y2) { # Illian et al. (2008), p. 184
  dx <- abs(x1 - x2)
  dy <- abs(y1 - y2)
  # dx <- min(dx, xmax - dx)
  # dy <- min(dy, ymax - dy)
  dz <- sqrt(dx * dx + dy * dy)
  return(dz)
}

calcMingling <- function(mark, neighbours, k) {
  mingling <- NA
  for (i in 1 : length(mark)) {
    msum <- 0
    for (j in 1 : k) {
      if(mark[i] !=  mark[neighbours[i, j]]) 
        msum <- msum + 1
    }
    mingling[i] <- msum / k
  }
  return(mingling)
}


calcMingling2 <- function(mark, neighbours, k, ns) {
  mingling <- NA
  for (i in 1 : length(mark)) {
    msum <- 0
    n <- c()
    n[1] <- mark[i]
    for (j in 1 : k) {
      if(mark[i] !=  mark[neighbours[i, j]]) 
        msum <- msum + 1
      n[j + 1] <- mark[neighbours[i, j]]
    }
    m <- min(k + 1, ns)
    mingling[i] <- msum * length(unique(n))/ (k * m)
  }
  return(mingling)
}

calcExpectedMinglingAllSpecies <- function(type) {
  ta <- table(type)
  s <- length(ta)
  ka <- length(type)
  swm <- 0
  for (i in 1 : s)
    swm <- swm + ta[[i]] * (ka - ta[[i]]) / (ka * (ka - 1))
  return(swm)   
}

calcExpectedMinglingOneSpecies <- function(species.vector, species) {
  ka <- length(species.vector)
  swm <- (ka - length(species.vector[species.vector == species]))/(ka - 1)
  return(swm)   
}

xdiff <- function(d1, d2) {
  stopifnot(all(d1 > 0))
  stopifnot(all(d2 > 0))
  d <- d1 / d2
  if(d > 1)
    d <- 1 / d
  return(d)	
}

calcDiff <- function(mark, neighbours, k) {
  td <- NA
  for (i in 1 : length(mark)) {
    tsum <- 0
    for (j in 1 : k) {
      tsum <- tsum + xdiff(mark[i], mark[neighbours[i, j]])
    }
    td[i] <- 1 - tsum / k
  }
  return(td)
}

ddiff <- function(m1, m2) {
  stopifnot(all(m1 > 0))
  stopifnot(all(m2 > 0))
  m <- (m1 > m2) * (1 - m2 / m1) - ((m1 < m2) * (1 - m1 / m2))
  return(m)
}

calcDD <- function(mark, neighbours, k) {
  td <- NA
  for (i in 1 : length(mark)) {
    tsum <- 0
    for (j in 1 : k) {
      tsum <- tsum + ddiff(mark[i], mark[neighbours[i, j]])
    }
    td[i] <- tsum / k
  }
  return(td)
}

calcExpectedDiffOneSpecies <- function(mark, svector, species) {
  r <- c()
  r[1] <- 0
  s <- c()
  s[length(mark)] <- 0
  mark <- sort(mark)
  for (i in 2 : length(mark)) 
    r[i]  <-  r[i - 1] + mark[i - 1]
  for (i in 1 : (length(mark) - 1)) 
    s[length(mark) - i] <- s[length(mark) + 1 - i] + 1 / mark[length(mark) + 1 - i]
  ET <- 0
  for (i in 1 : length(mark)) {
    if(svector[i] == species)
      ET <- ET + r[i] / mark[i] + s[i] * mark[i]
  }
  ET <- 1 - 1 / ((length(mark) - 1) * length(mark[svector == species])) * ET
  return(ET)   
}

calcDominance <- function(mark, neighbours, k) {
  dm <- NA
  for (i in 1 : length(mark)) {
    dsum <- 0
    for (j in 1 : k) {
      if(mark[i] > mark[neighbours[i, j]])
        dsum <- dsum + 1
    }
    dm[i] <- dsum / k
  }
  return(dm)
}

calcDominance2 <- function(mark, neighbours, k) {
  dm <- NA
  for (i in 1 : length(mark)) {
    dsum <- 0
    tsum <- 0 # pi * (mark[i] / 200)^2
    for (j in 1 : k) {
      if(mark[i] > mark[neighbours[i, j]])
        dsum <- dsum + mark[neighbours[i, j]] # pi * (mark[neighbours[i, j]] / 200)^2
      tsum <- tsum + mark[neighbours[i, j]] # pi * (mark[neighbours[i, j]] / 200)^2
    }
    dm[i] <- dsum / tsum
  }
  return(dm)
}

calcTrigDiff <- function(mark, neighbours, k, alpha) {
  td <- NA
  for (i in 1 : length(mark)) {
    tsum <- 0
    for (j in 1 : k) {
      tsum <- tsum + mark[i]^(2 * alpha) / (mark[i]^(2 * alpha) + mark[neighbours[i, j]]^(2 * alpha))
      # tsum <- tsum + 0.5 * (tanh(alpha * log(mark[neighbours[i, j]] / mark[i])) + 1)
    }
    td[i] <- tsum / k 
    # td[i] <- 1 - tsum / k # Considering the compliment facilitates intuitive interpretation
  }
  return(td)
}

calcDCI <- function(mark, neighbours, k) {
  dci <- NA
  # mark <- rank(mark, ties.method = "average") / length(mark)
  for (i in 1 : length(mark)) {
    dsum <- md <- 0
    for (j in 1 : k) {
      dsum <- dsum + (mark[i] * mark[neighbours[i, j]])
      md <- md + mark[neighbours[i, j]]
    }
    # dci[i] <- mark[i] * dsum / ((mark[i] + dsum) / (k + 1))^2
    md <- (md + mark[i]) / (k + 1)
    dci[i] <- dsum / (k * md^2)
    # dci[i] <- mark[i] * dsum / mean(mark)^2
  }
  return(dci)
}

calcVarIndex <- function(mark, neighbours, k) {
  td <- NA
  for (i in 1 : length(mark)) {
    tsum <- 0
    marks <- c()
    marks <- c(marks, mark[i])
    for (j in 1 : k) {
      tsum <- tsum + 0.5 * (mark[i] - mark[neighbours[i, j]])^2
      marks <- c(marks, mark[neighbours[i, j]])
    }
    varm <- var(marks)
    if(is.na(varm))
      varm <- 1
    td[i] <- tsum / (k * varm)
  }
  return(td)
}

# varm <- NA
# if(is.na(varm))
#   varm <- 1


which(filenames == "YLK12_2013.csv")
filenames[37] # 5

# z <- 1
for (z in 1 : length(filenames)) {
# for (z in 25 : length(filenames)) {
  # if((z < 14) | (z > 37))
  if(z < 14)
    myData <- read.table(paste(dataPath, filenames[z], sep = ""), header = TRUE)
  else myData <- read.csv(paste(dataPath, filenames[z], sep = ""), header = TRUE)
  # names(myData)
  # head(myData)
  range(myData$x)
  range(myData$y)
  if(z == 1)
    myData$species <- myData$speciesCode
  myData <- myData[c("dbh", "species", "x", "y")] 
  # names(myData)
  if(z == 8) # 7
    myData$x <- myData$x - 70
  if(z == 9) # 8
    myData$x <- myData$x - 2
  if(z == 13) {
    myData$x <- myData$x - 199.261
    myData$y <- myData$y - 199.4
  }
  # if(z > 37)
  #   myData <- myData[myData$dbh >= 5,]
  # sort(table(myData$species))
  numberSpecies <- length(table(myData$species))
  if(numberSpecies <= 1) stop("Number ", z, " has only one species!\n")
  species <- as.numeric(names(sort(table(myData$species), decreasing = T)[sort(table(myData$species), decreasing = T) > 20]))

  xwindow <- owin(c(0, xmax[z]), c(0, ymax[z]))
  myDataP <- ppp(myData$x, myData$y, window = xwindow, marks = myData$species)
  k <- seq(1, 50, 1) # 40

  corMT <- corMT2 <- corMDD <- corMDD2 <- corMU <- corMU2 <- corMUU <- corMUU2 <- corPU <- mdist <- corMU3 <- corMU32 <- 
    corMdci <- corM2dci <- corMvi <- corM2vi <- c()
  # i <- 1
  for (i in 1 : length(k)) {
    neighbours <- nnwhich(myDataP, k = 1 : k[i])
    dist <- nndist(myDataP, k = k[i])
    c <-  bdist.points(myDataP)
    Fi <- (xmax[z] - 2 * dist) * (ymax[z] - 2 * dist) / 10000
    RF <- (dist < c) * 1 / Fi
    mi <- calcMingling(myData$species, as.matrix(neighbours), k[i])
    meanMingling <- c()
    for (j in 1 : length(species))
      meanMingling[j] <- sum(mi[myData$species == species[j]] * RF[myData$species == species[j]]) / sum(RF[myData$species == species[j]])
    mi2 <- calcMingling2(myData$species, as.matrix(neighbours), k[i], numberSpecies)
    meanMingling2 <- c()
    for (j in 1 : length(species))
      meanMingling2[j] <- sum(mi2[myData$species == species[j]] * RF[myData$species == species[j]]) / sum(RF[myData$species == species[j]])
    ti <- calcDiff(myData$dbh, as.matrix(neighbours), k[i])
    meanDifferentiation <- c()
    for (j in 1 : length(species))
      meanDifferentiation[j] <- sum(ti[myData$species == species[j]] * RF[myData$species == species[j]]) / sum(RF[myData$species == species[j]])
    dd <- calcDD(myData$dbh, as.matrix(neighbours), k[i])
    meanDD <- c()
    for (j in 1 : length(species))
      meanDD[j] <- sum(dd[myData$species == species[j]] * RF[myData$species == species[j]]) / sum(RF[myData$species == species[j]])
    ui <- calcDominance(myData$dbh, as.matrix(neighbours), k[i])
    meanU <- c()
    for (j in 1 : length(species))
      meanU[j] <- sum(ui[myData$species == species[j]] * RF[myData$species == species[j]]) / sum(RF[myData$species == species[j]])
    ui2 <- calcDominance2(myData$dbh, as.matrix(neighbours), k[i])
    meanU2 <- c()
    for (j in 1 : length(species))
      meanU2[j] <- sum(ui2[myData$species == species[j]] * RF[myData$species == species[j]]) / sum(RF[myData$species == species[j]])
    # mean(ui2)
    # hist(ui2)
    ui3 <- calcTrigDiff(myData$dbh, as.matrix(neighbours),  k[i], 0.5) # alpha = 0.5 so that exponent = 1
    meanU3 <- c()
    for (j in 1 : length(species))
      meanU3[j] <- sum(ui3[myData$species == species[j]] * RF[myData$species == species[j]]) / sum(RF[myData$species == species[j]])
    
    # dci <- calcDCI(myData$dbh, as.matrix(neighbours), k[i])
    # dci <- (dci - min(dci)) / (max(dci) - min(dci))
    # meanDci <- c()
    # for (j in 1 : length(species))
    #   meanDci[j] <- sum(dci[myData$species == species[j]] * RF[myData$species == species[j]]) / sum(RF[myData$species == species[j]])
    # 
    # vi <- calcVarIndex(myData$dbh, as.matrix(neighbours), k[i]) 
    # vi <- (vi - min(vi)) / (max(vi) - min(vi))
    # meanVi <- c()
    # for (j in 1 : length(species))
    #   meanVi[j] <- sum(vi[myData$species == species[j]] * RF[myData$species == species[j]]) / sum(RF[myData$species == species[j]])
   

    corMT[i] <- cor(meanDifferentiation, meanMingling)
    corMT2[i] <- cor(meanDifferentiation, meanMingling2)
    corMDD[i] <- cor(meanDD, meanMingling)
    corMDD2[i] <- cor(meanDD, meanMingling2)
    corMU[i] <- cor(meanU, meanMingling)
    corMU2[i] <- cor(meanU, meanMingling2)
    corMUU[i] <- cor(meanU2, meanMingling)
    corMUU2[i] <- cor(meanU2, meanMingling2)
    corMU3[i] <- cor(meanU3, meanMingling)
    corMU32[i] <- cor(meanU3, meanMingling2)
    # corMdci[i] <- cor(meanDci, meanMingling)
    # corM2dci[i] <- cor(meanDci, meanMingling2)
    # corMvi[i] <- cor(meanVi, meanMingling)
    # corM2vi[i] <- cor(meanVi, meanMingling2)
    mdist[i] <- sum(dist[myData$species == species[j]] * RF[myData$species == species[j]], na.rm = T) / sum(RF[myData$species == species[j]], na.rm = T)

    expMingling <- c()
    for (j in 1 : length(species))
      expMingling[j] <- calcExpectedMinglingOneSpecies(myData$species, species[j])
    expDifferentiation <- c()
    for (j in 1 : length(species))
      expDifferentiation[j] <- calcExpectedDiffOneSpecies(myData$dbh, myData$species, species[j])
    Psi <- 1 - meanMingling / expMingling
    Ups <- 1 - meanDifferentiation / expDifferentiation
    corPU[i] <- cor(Ups, Psi)

    # if(i == 5) {
    #   data1 <- data.frame(meanDifferentiation, meanMingling)
    #   data2 <- data.frame(meanDifferentiation, meanMingling2)
    #   data3 <- data.frame(Ups, Psi)
    #   data4 <- data.frame(meanDD, meanMingling)
    #   data5 <- data.frame(meanDD, meanMingling2)
    #   data6 <- data.frame(meanU, meanMingling)
    #   data7 <- data.frame(meanU, meanMingling2)
    #   data8 <- data.frame(meanU2, meanMingling)
    #   data9 <- data.frame(meanU2, meanMingling2)
    #   data10 <- data.frame(meanU3, meanMingling)
    #   data11 <- data.frame(meanU3, meanMingling2)
    #   p <- ggplot(data1, aes(y = meanMingling, x = meanDifferentiation)) + geom_point(color = "blue") + labs(title = "", x = "", y = "") + theme(axis.text = element_text(color = "black", size = 16),
    #          plot.margin = unit(c(0, 0, 0, 0), "mm")) + xlim(-0.5, 1) + ylim(-0.5, 1) + geom_smooth(aes(y = meanMingling, x = meanDifferentiation), color = "blue", method = "lm", fill = alpha("blue", 0.2)) +
    #          geom_point(data2, mapping = aes(y = meanMingling2, x = meanDifferentiation), color = "blue") + geom_smooth(aes(y = meanMingling2, x = meanDifferentiation), color = 'blue', method = "lm", fill = alpha("blue", 0.1)) +
    #          geom_point(data3, mapping = aes(y = Psi, x = Ups), color = "orange") + geom_smooth(aes(y = Psi, x = Ups), color = 'orange', method = "lm", fill = alpha("orange", 0.2)) +
    #          geom_point(data4, mapping = aes(y = meanMingling, x = meanDD), color = "darkorchid1") + geom_smooth(aes(y = meanMingling, x = meanDD), color = 'darkorchid1', method = "lm", fill = alpha("darkorchid1", 0.2)) +
    #          geom_point(data5, mapping = aes(y = meanMingling2, x = meanDD), color = "darkorchid1") + geom_smooth(aes(y = meanMingling2, x = meanDD), color = 'darkorchid1', method = "lm", fill = alpha("darkorchid1", 0.1)) +
    #          geom_point(data6, mapping = aes(y = meanMingling, x = meanU), color = "darkgreen") + geom_smooth(aes(y = meanMingling, x = meanU), color = 'darkgreen', method = "lm", fill = alpha("darkgreen", 0.2)) +
    #          geom_point(data7, mapping = aes(y = meanMingling2, x = meanU), color = "darkgreen") + geom_smooth(aes(y = meanMingling2, x = meanU), color = 'darkgreen', method = "lm", fill = alpha("darkgreen", 0.1)) +
    #          geom_point(data8, mapping = aes(y = meanMingling, x = meanU2), color = "black") + geom_smooth(aes(y = meanMingling, x = meanU2), color = 'black', method = "lm", fill = alpha("bisque4", 0.2)) +
    #          geom_point(data9, mapping = aes(y = meanMingling2, x = meanU2), color = "black") + geom_smooth(aes(y = meanMingling2, x = meanU2), color = 'black', method = "lm", fill = alpha("bisque4", 0.1)) +
    #          geom_point(data10, mapping = aes(y = meanMingling, x = meanU3), color = "cyan") + geom_smooth(aes(y = meanMingling, x = meanU3), color = 'cyan', method = "lm", fill = alpha("cyan", 0.2)) +
    #          geom_point(data11, mapping = aes(y = meanMingling2, x = meanU3), color = "cyan") + geom_smooth(aes(y = meanMingling2, x = meanU3), color = 'cyan', method = "lm", fill = alpha("cyan", 0.1)) 
    #     suppressMessages(ggsave(filename = paste(resultPath, "MinglingDifferentiation", substring(filenames[z], 1, nchar(filenames[z]) - 4), ".pdf", sep = ""),
    #          plot = p, # or an explicit ggplot object name,
    #          units = "mm", # other options c("in", "cm", "mm"),
    #          dpi = 400))
    # 
    #   rm(data1, data2, data3, p)
    # }
    
  }

  cat("Stand: ", z, " of ", length(filenames), "Name ", substring(filenames[z], 1, nchar(filenames[z]) - 4), "\n")
  # data1 <- data.frame(mdist, corMT)
  # data2 <- data.frame(mdist, corMT2)
  # data3 <- data.frame(mdist, corPU)
  # data4 <- data.frame(mdist, corMDD)
  # data5 <- data.frame(mdist, corMDD2)
  # data6 <- data.frame(mdist, corMU)
  # data7 <- data.frame(mdist, corMU2)
  # p <- ggplot(data1, aes(y = corMT, x = mdist)) + geom_point(color = "blue") + labs(title = "", x = "", y = "") + theme(axis.text = element_text(color = "black", size = 16),
  #    plot.margin = unit(c(0, 0, 0, 0), "mm")) + xlim(0, 15) + ylim(-0.5, 1) + geom_smooth(aes(y = corMT, x = mdist), color = 'blue', method = "loess", fill = alpha("blue", 0.2)) +
  #    geom_point(data2, mapping = aes(y = corMT2, x = mdist), color = "blue4") + geom_smooth(aes(y = corMT2, x = mdist), color = 'blue4', method = "loess", fill = alpha("blue4", 0.2)) +
  #    geom_point(data3, mapping = aes(y = corPU, x = mdist), color = "orange") + geom_smooth(aes(y = corPU, x = mdist), color = 'orange', method = "loess", fill = alpha("orange", 0.2)) +
  #    geom_point(data4, mapping = aes(y = corMDD, x = mdist), color = "darkred") + geom_smooth(aes(y = corMDD, x = mdist), color = 'darkred', method = "loess", fill = alpha("darkred", 0.2)) +
  #    geom_point(data5, mapping = aes(y = corMDD2, x = mdist), color = "red1") + geom_smooth(aes(y = corMDD2, x = mdist), color = 'red1', method = "loess", fill = alpha("red1", 0.2)) +
  #    geom_point(data6, mapping = aes(y = corMU, x = mdist), color = "darkgreen") + geom_smooth(aes(y = corMU, x = mdist), color = 'darkgreen', method = "loess", fill = alpha("darkgreen", 0.2)) +
  #    geom_point(data7, mapping = aes(y = corMU2, x = mdist), color = "darkolivegreen") + geom_smooth(aes(y = corMU2, x = mdist), color = 'darkolivegreen', method = "loess", fill = alpha("darkolivegreen", 0.2)) +
  #    geom_vline(xintercept = mdist[5], linetype = "dashed", color = "black", size = 0.2) +   
  #    geom_vline(xintercept = mdist[10], linetype = "dashed", color = "black", size = 0.2)    
  # suppressMessages(ggsave(filename = paste(resultPath, "MinglingDifferentiationCorr", substring(filenames[z], 1, nchar(filenames[z]) - 4), ".pdf", sep = ""),
  #        plot = p, # or an explicit ggplot object name,
  #        units = "mm", # other options c("in", "cm", "mm"), 
  #        dpi = 400))
  
  allData <- data.frame(dist = mdist, corr = corMT)
  allData <- rbind(allData, data.frame(dist = mdist, corr = corMT2))
  allData <- rbind(allData, data.frame(dist = mdist, corr = corPU))
  allData <- rbind(allData, data.frame(dist = mdist, corr = corMDD))
  allData <- rbind(allData, data.frame(dist = mdist, corr = corMDD2))
  allData <- rbind(allData, data.frame(dist = mdist, corr = corMU))
  allData <- rbind(allData, data.frame(dist = mdist, corr = corMU2))
  allData <- rbind(allData, data.frame(dist = mdist, corr = corMUU))
  allData <- rbind(allData, data.frame(dist = mdist, corr = corMUU2))
  # allData <- rbind(allData, data.frame(dist = mdist, corr = corMdci))
  # allData <- rbind(allData, data.frame(dist = mdist, corr = corM2dci))
  # allData <- rbind(allData, data.frame(dist = mdist, corr = corMvi))
  # allData <- rbind(allData, data.frame(dist = mdist, corr = corM2vi))
  funcs <- matrix(nrow = length(mdist), ncol = 11)
  funcs[, 1] <- corMT
  funcs[, 2] <- corMT2
  funcs[, 3] <- corPU
  funcs[, 4] <- corMDD
  funcs[, 5] <- corMDD2
  funcs[, 6] <- corMU
  funcs[, 7] <- corMU2
  funcs[, 8] <- corMUU
  funcs[, 9] <- corMUU2
  funcs[, 10] <- corMU3
  funcs[, 11] <- corMU32
  # funcs[, 12] <- corMdci
  # funcs[, 13] <- corM2dci
  # funcs[, 14] <- corMvi
  # funcs[, 15] <- corM2vi
  lo <- apply(funcs, MARGIN = 1, FUN = quantile, probs = 0.025) #, na.rm = T)
  hi <- apply(funcs, MARGIN = 1, FUN = quantile, probs = 0.975) #, na.rm = T)
  # me <- apply(funcs, MARGIN = 1, FUN = quantile, probs = 0.5)
  lo <- ksmooth(x = mdist, y = lo, kernel = "normal", bandwidth = 3)
  hi <- ksmooth(x = mdist, y = hi, kernel = "normal", bandwidth = 3)
  # me1 <- ksmooth(x = mdist, y = me, kernel = "normal", bandwidth = 3)
  # me2 <- me1
  # me2$y <- apply(cbind(lo$y, hi$y), 1, median) # median(lo$x, hi$x)
  # ?ksmooth
  pdf(file = paste(resultPath, "Funnel", substring(filenames[z], 1, nchar(filenames[z]) - 4), ".pdf", sep = ""))
  par(mar = c(2, 2.8, 0.5, 0.5)) 
  plot(allData$corr, allData$dist, main = "", xlab = "", ylab = "", pch = 16, cex = 0.5, col = "white", lwd = 2, axes = FALSE, xlim = c(-1, 1), ylim = c(0, 15))#, yaxt = 'n')
  mt <- ksmooth(x = mdist, y = corMT, kernel = "normal", bandwidth = 2)
  mt2 <- ksmooth(x = mdist, y = corMT2, kernel = "normal", bandwidth = 2)
  pu <- ksmooth(x = mdist, y = corPU, kernel = "normal", bandwidth = 2)
  mdd <- ksmooth(x = mdist, y = corMDD, kernel = "normal", bandwidth = 2)
  mdd2 <- ksmooth(x = mdist, y = corMDD2, kernel = "normal", bandwidth = 2)
  mu <- ksmooth(x = mdist, y = corMU, kernel = "normal", bandwidth = 2)
  mu2 <- ksmooth(x = mdist, y = corMU2, kernel = "normal", bandwidth = 2)
  muu <- ksmooth(x = mdist, y = corMUU, kernel = "normal", bandwidth = 2)
  muu2 <- ksmooth(x = mdist, y = corMUU2, kernel = "normal", bandwidth = 2)
  mu3 <- ksmooth(x = mdist, y = corMU3, kernel = "normal", bandwidth = 2)
  mu32 <- ksmooth(x = mdist, y = corMU32, kernel = "normal", bandwidth = 2)
  # mdci <- ksmooth(x = mdist, y = corMdci, kernel = "normal", bandwidth = 2)
  # m2dci <- ksmooth(x = mdist, y = corM2dci, kernel = "normal", bandwidth = 2)
  # mvi <- ksmooth(x = mdist[-1], y = corMvi[-1], kernel = "normal", bandwidth = 2)
  # m2vi <- ksmooth(x = mdist[-1], y = corM2vi[-1], kernel = "normal", bandwidth = 2)
  polygon(c(lo$y, rev(hi$y)), c(lo$x, rev(hi$x)), col = "lightgray", border = "lightgray")
  points(corMT, mdist, pch = 16, col = "blue", cex = 0.5)
  points(corMT2, mdist, pch = 16, col = "blue", cex = 0.5)
  points(corPU, mdist, pch = 16, col = "orange", cex = 0.5)
  points(corMDD, mdist, pch = 16, col = "purple", cex = 0.5)
  points(corMDD2, mdist, pch = 16, col = "purple", cex = 0.5)
  points(corMU, mdist, pch = 16, col = "green", cex = 0.5)
  points(corMU2, mdist, pch = 16, col = "green", cex = 0.5)
  points(corMUU, mdist, pch = 16, col = "black", cex = 0.5)
  points(corMUU2, mdist, pch = 16, col = "black", cex = 0.5)
  points(corMU3, mdist, pch = 16, col = "cyan", cex = 0.5)
  points(corMU32, mdist, pch = 16, col = "cyan", cex = 0.5)
  # points(corMdci, mdist, pch = 16, col = "yellow", cex = 0.5)
  # points(corM2dci, mdist, pch = 16, col = "yellow", cex = 0.5)
  # points(corMvi, mdist, pch = 16, col = "red", cex = 0.5)
  # points(corM2vi, mdist, pch = 16, col = "red", cex = 0.5)
  # lines(lo$y, lo$x, lwd = 4, lty = 1,  col = "black")
  # lines(hi$y, hi$x, lwd = 4, lty = 1, col = "black")
  lines(mt$y, mt$x, lwd = 1, lty = 1, col = "blue")
  lines(mt2$y, mt2$x, lwd = 1, lty = 2, col = "blue")
  lines(pu$y, pu$x, lwd = 1, lty = 1, col = "orange")
  lines(mdd$y, mdd$x, lwd = 1, lty = 1, col = "purple")
  lines(mdd2$y, mdd2$x, lwd = 1, lty = 2, col = "purple")
  lines(mu$y, mu$x, lwd = 1, lty = 1, col = "green")
  lines(mu2$y, mu2$x, lwd = 1, lty = 2, col = "green")
  lines(muu$y, muu$x, lwd = 1, lty = 1, col = "black")
  lines(muu2$y, muu2$x, lwd = 1, lty = 2, col = "black")
  lines(mu3$y, mu3$x, lwd = 1, lty = 1, col = "cyan")
  lines(mu32$y, mu32$x, lwd = 1, lty = 2, col = "cyan")
  # lines(mdci$y, mdci$x, lwd = 1, lty = 1, col = "yellow")
  # lines(m2dci$y, m2dci$x, lwd = 1, lty = 2, col = "yellow")
  # lines(mvi$y, mvi$x, lwd = 1, lty = 1, col = "red")
  # lines(m2vi$y, m2vi$x, lwd = 1, lty = 2, col = "red")
  # lines(me2$y, me2$x, lwd = 4, lty = 1, col = "black")
  abline(h = mdist[5], lty = 2, lwd = 1)
  abline(h = mdist[10], lty = 2, lwd = 1)
  axis(1, lwd = 2, cex.axis = 1.7)
  axis(2, las = 1, lwd = 2, cex.axis = 1.7)
  box(lwd = 2)
  dev.off()
  # rm(data1, data2, data3, data4, data5, data6, data7, p, allData, funcs, lo, hi, me, me1, me2, corMT, corMT2, corPU, corMDD,
  #    corMDD2, corMU, corMU2, corMUU, corMUU2, corMU3, corMU32, mt, mt2, 
  #    pu, mdd, mdd2, mu, mu2, muu, muu2, mu3, mu23)
  rm(allData, funcs, lo, hi, corMT, corMT2, corPU, corMDD,
     corMDD2, corMU, corMU2, corMUU, corMUU2, corMU3, corMU32, mt, mt2,
     pu, mdd, mdd2, mu, mu2, muu, muu2, mu3, mu32) # me, me1, me2, corMdci, corM2dci, corMvi, corM2vi, , mdci, m2dci, mvi, m2vi
}  



