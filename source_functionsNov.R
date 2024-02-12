print(getwd())

library(BiocManager)
options(repos = BiocManager::repositories())
# if(!require("devtools"))
#   install.packages("devtools")
# library(rsconnect)
# if(!require(devtools)) {
#   BiocManager::install("devtools")
#   library(devtools)
# }
# 
# if(!require(rsconnect)) {
#   BiocManager::install("rsconnect")
#   library(rsconnect)
# }


if(!require(minfiData)) {
  BiocManager::install("minfiData")
  library(minfiData)
}
message("minfidata done")
if(!require(sva)){
  BiocManager::install("sva")
  library(sva)
}
message("sva done")
if(!require(bumphunter)){
  BiocManager::install("bumphunter")
  library(bumphunter)
}
message("bumphunter done")
# if(!require(lumi)){
#   BiocManager::install("lumi")
#   library(lumi)
# }
# message("lumi done")
# 
# if(!require(minfi)){
#   BiocManager::install("minfi")
#   library(minfi)
# }
# 
# message("minfi done")
 if(!require(IlluminaHumanMethylationEPICmanifest)){
   BiocManager::install("IlluminaHumanMethylationEPICmanifest")
   library(IlluminaHumanMethylationEPICmanifest)
 }
# message("illumina1 done")
 if(!require(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)){
   BiocManager::install("IlluminaHumanMethylationEPICanno.ilm10b4.hg19")
   library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
 }
#message("manifest loaded")
if(!require(shiny)){
  install.packages("shiny")
  library(shiny)
}
message("shiny done")

if(!require(DT)){
  install.packages("DT")
  library(DT)
}

message("DT done")

if(!require(shinyWidgets)){
  install.packages("shinyWidgets")
  library(shinyWidgets)
}
message("shinyWidgets done")


if(!require(shinydashboard)){
  install.packages("shinydashboard")
  library(shinydashboard)
}
message("shinydashboard done")


if(!require(ggplot2)){
  install.packages("ggplot2")
  library(ggplot2)
}
message("ggplot2 done")


if(!require(ggpubr)){
  install.packages("ggpubr")
  library(ggpubr)
}
message("ggpubr done")

if(!require(ggnewscale)){
  install.packages("ggnewscale")
  library(ggnewscale)
}
message("ggnewscale done")

if(!require(foreach)){
  install.packages("foreach")
  library(foreach)
}
message("foreach done")

if(!require(ggrepel)){
  install.packages("ggrepel")
  library(ggrepel)
}
message("ggrepel done")

if(!require(pheatmap)){
  install.packages("pheatmap")
  library(pheatmap)
}
message("pheatmap done")

if(!require(shinycssloaders)){
  install.packages("shinycssloaders")
  library(shinycssloaders)
}
message("shinycssloaders done")

if(!require(shinybusy)){
  install.packages("shinybusy")
  library(shinybusy)
}
message("shinybusy done")

if(!require(waiter)){
  install.packages("waiter")
  library(waiter)
}
message("waiter done")

if(!require(tinytex)){
  install.packages('tinytex')
  library(tinytex)
}


library(gridExtra)
library(grid)




message("packages loaded")
#load(file = "https://github.com/hackingjpr/Idat-Shiny/blob/main/ATRT.v3.abs.chun.Rdata")
 # load(file = "./ATRT.v3.abs.chun.Rdata")
 readRDS(file = "./atrt34.model.rds") -> ATRT
# atrt.meth.os.meta.n8.extract -> ATRT
#load(file = "https://github.com/hackingjpr/Idat-Shiny/blob/main/ECRT.v3.abs.chun.Rdata")
 # load(file = "./ECRT.v3.abs.chun.Rdata")
readRDS(file = "./ecrt32.model.rds") -> ECRT
# ecrt.meth.os.meta.n32.extract -> ECRT
#load(file = "https://github.com/hackingjpr/Idat-Shiny/blob/main/ALL.v3.abs.chun.Rdata")
 # load(file = "./ALL.v3.abs.chun.Rdata")
 readRDS(file = "./mrt13.model.rds") -> ALL
 # all.meth.os.meta.n13.extract -> ALL

beta2m <- function (beta) {
  m <- log2(beta/(1 - beta))
  return(m)
}



process_idats <- function(basenames){
  
  #check file sizes, epic will be over 1e7
  # idx.epic <-
  #   file.info(paste0(basenames, "_Red.idat"))$size > 1e7
  
  # Split into mixed or single test e.g. 450k & Epic or one or the other
  file.sizes <- file.info(paste0(basenames, "_Red.idat"))$size
  
  length(which(file.sizes> 1e7))>0 & length(which(file.sizes< 1e7))>0 -> mixed.test
  length(which(file.sizes> 1e7))>0 | length(which(file.sizes< 1e7))>0 -> single.test
  
  if(mixed.test) {
    #  as before where you process them seperately then combine
    idx.epic <-
      file.info(paste0(basenames, "_Red.idat"))$size > 1e7
    m450k <- read.metharray(basenames[!idx.epic],
                            verbose = TRUE,
                            force = TRUE) 
    mEpic <- read.metharray(basenames[idx.epic],
                            verbose = TRUE,
                            force = TRUE) 
    comb.rgSet <- combineArrays(m450k, 
                                mEpic,
                                outType = "IlluminaHumanMethylation450k")
    comb.mSet <- preprocessNoob(comb.rgSet, 
                                dyeMethod = "single")
    comb.detP <- detectionP(comb.rgSet)
    comb.gmrSet <- mapToGenome(ratioConvert(comb.mSet))
    comb.beta <- getBeta(comb.gmrSet)
    comb.MRT.pd <- pData(comb.gmrSet)
  }
  
  if((!mixed.test)&single.test){
    # only need the one
    single <- read.metharray(basenames,
                             verbose = TRUE,
                             force = TRUE)
    comb.mSet <- preprocessNoob(single, dyeMethod = "single")
    comb.detP <- detectionP(single)
    comb.gmrSet <- mapToGenome(ratioConvert(comb.mSet))
    comb.beta <- getBeta(comb.gmrSet)
    comb.MRT.pd <- pData(comb.gmrSet)
  }
  
  return(list(betas = comb.beta, pvals = comb.detP))
}

get_basenames <- function(dir){
  list.files(path = dir, pattern = "_Red.idat", full.names = T) -> temp.files
  return(gsub("_Red.idat","", temp.files))
}

extract.metagene <- function(index, weights, exp.matrix, scaling) {
  exp.matrix[index, ] -> temp.exp
  apply(temp.exp, 2, function(x) {
    mean(x * weights)
  }) -> raw.metagene
  round(raw.metagene, digits = 3)
  as.numeric(scale(raw.metagene, center = scaling[1], scale = scaling[2])) -> Risk_Value
  round(Risk_Value, digits = 3)
  names(Risk_Value) <- colnames(exp.matrix)
  df <- as.data.frame(Risk_Value)
  return(df)
}

digits = 0:9
createRandString<- function() {
  v = c(sample(LETTERS, 5, replace = TRUE),
        sample(digits, 4, replace = TRUE),
        sample(LETTERS, 1, replace = TRUE))
  return(paste0(v,collapse = ""))
}

generate_figure_highlight_mrt <- function(new.sample.meta.score, indexRow) {
  if(is.null(indexRow)){indexRow=1}
  #temp.df <- readRDS(file = "https://github.com/hackingjpr/Idat-Shiny/blob/main/temp.df.rds")
   temp.df <- readRDS(file = "./mrt13.dist.rds")
  #df.cat.atrt <- readRDS(file = "https://github.com/hackingjpr/Idat-Shiny/blob/main/df.cat.atrt.rds")
   # df.cat.atrt <- readRDS(file = "./df.cat.atrt.rds")
  #comb.SDb.atrt <- readRDS(file = "https://github.com/hackingjpr/Idat-Shiny/blob/main/comb.SDb.atrt.rds")
   # comb.SDb.atrt <- readRDS(file = "./comb.SDb.atrt.rds")
  ggplot(aes(x = 1:nrow(temp.df), y = mrt13), data = temp.df) +
    geom_line() +
    scale_shape_manual(values = c(1, 4,  3)) +
    scale_color_manual(values = c('#E69F00', '#999999', "white")) +
    ylab("MRT-13") +
    theme_minimal() +
    theme(
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      legend.position = "none"
    ) -> b
  
  df.lines.hor <-
    foreach(i = 1:length(new.sample.meta.score),
            .combine = rbind) %do% {
              data.frame(
                x = 0,
                xend = max(which(
                  temp.df$mrt13 < new.sample.meta.score[i]
                )),
                y = new.sample.meta.score[i],
                yend = new.sample.meta.score[i]
              )
            }
  df.lines.hor$labels <- names(new.sample.meta.score)
  
  df.lines.ver <-
    foreach(i = 1:length(new.sample.meta.score),
            .combine = rbind) %do% {
              data.frame(
                x = max(which(
                  temp.df$mrt13 < new.sample.meta.score[i]
                )),
                xend = max(which(
                  temp.df$mrt13 < new.sample.meta.score[i]
                )),
                y = new.sample.meta.score[i],
                yend = min(temp.df$mrt13)
              )
            }
  df.lines.ver$perc <-
    paste0(round(df.lines.ver$xend / length(temp.df$mrt13) * 100), "th")
  
  df.lines.ver$colour <- factor(ifelse(1:nrow(df.lines.ver)==indexRow,"highlight","no.highlight"), levels = c("highlight","no.highlight"))
  df.lines.hor$colour <- factor(ifelse(1:nrow(df.lines.hor)==indexRow,"highlight","no.highlight"), levels = c("highlight","no.highlight"))
  
  b <- b +
    geom_segment(
      aes(
        x = x,
        y = y,
        xend = xend,
        yend = yend,
        colour = as.character(colour)
      ),
      #colour = "red",
      linetype = "dashed",
      data = df.lines.hor
    ) +
    geom_segment(
      aes(
        x = x,
        y = y,
        xend = xend,
        yend = yend,
        colour = colour
      ),
      #colour = "red",
      linetype = "dashed",
      data = df.lines.ver
    ) +
    #geom_text_repel(aes(x = x+10, y = y+0.1, label = labels), direction = "y", data = df.lines.hor)
    geom_text(aes(
      x = x + 10,
      y = y + 0.1,
      label = labels,
      colour = colour
    ), data = df.lines.hor) 
  #scale_color_discrete("red", "lightgrey")
  
  df <- data.frame()
  c <- ggplot() + theme_void()
  
  
  #ggarrange(a,a2,ggarrange(
  # c,b, ncol = 2, nrow = 1, widths = c(0.015,1)),ncol=1,nrow=3)
  ggarrange(c,
            b,
            ncol = 2,
            nrow = 1,
            widths = c(0.015, 1))
  #return(data.frame(perc = df.lines.ver[,5], row.names=rownames(df.lines.ver)))
}

generate_figure_highlight_ecrt <- function(new.sample.meta.score, indexRow) {
  if(is.null(indexRow)){indexRow=1}
  #temp.df <- readRDS(file = "https://github.com/hackingjpr/Idat-Shiny/blob/main/temp.df.rds")
  temp.df <- readRDS(file = "./ecrt32.dist.rds")
  #df.cat.atrt <- readRDS(file = "https://github.com/hackingjpr/Idat-Shiny/blob/main/df.cat.atrt.rds")
  # df.cat.atrt <- readRDS(file = "./df.cat.atrt.rds")
  #comb.SDb.atrt <- readRDS(file = "https://github.com/hackingjpr/Idat-Shiny/blob/main/comb.SDb.atrt.rds")
  # comb.SDb.atrt <- readRDS(file = "./comb.SDb.atrt.rds")
  ggplot(aes(x = 1:nrow(temp.df), y = ecrt32), data = temp.df) +
    geom_line() +
    scale_shape_manual(values = c(1, 4,  3)) +
    scale_color_manual(values = c('#E69F00', '#999999', "white")) +
    ylab("ECRT-32") +
    theme_minimal() +
    theme(
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      legend.position = "none"
    ) -> b
  
  df.lines.hor <-
    foreach(i = 1:length(new.sample.meta.score),
            .combine = rbind) %do% {
              data.frame(
                x = 0,
                xend = max(which(
                  temp.df$ecrt32 < new.sample.meta.score[i]
                )),
                y = new.sample.meta.score[i],
                yend = new.sample.meta.score[i]
              )
            }
  df.lines.hor$labels <- names(new.sample.meta.score)
  
  df.lines.ver <-
    foreach(i = 1:length(new.sample.meta.score),
            .combine = rbind) %do% {
              data.frame(
                x = max(which(
                  temp.df$ecrt32 < new.sample.meta.score[i]
                )),
                xend = max(which(
                  temp.df$ecrt32 < new.sample.meta.score[i]
                )),
                y = new.sample.meta.score[i],
                yend = min(temp.df$ecrt32)
              )
            }
  df.lines.ver$perc <-
    paste0(round(df.lines.ver$xend / length(temp.df$ecrt32) * 100), "th")
  
  df.lines.ver$colour <- factor(ifelse(1:nrow(df.lines.ver)==indexRow,"highlight","no.highlight"), levels = c("highlight","no.highlight"))
  df.lines.hor$colour <- factor(ifelse(1:nrow(df.lines.hor)==indexRow,"highlight","no.highlight"), levels = c("highlight","no.highlight"))
  
  b <- b +
    geom_segment(
      aes(
        x = x,
        y = y,
        xend = xend,
        yend = yend,
        colour = as.character(colour)
      ),
      #colour = "red",
      linetype = "dashed",
      data = df.lines.hor
    ) +
    geom_segment(
      aes(
        x = x,
        y = y,
        xend = xend,
        yend = yend,
        colour = colour
      ),
      #colour = "red",
      linetype = "dashed",
      data = df.lines.ver
    ) +
    #geom_text_repel(aes(x = x+10, y = y+0.1, label = labels), direction = "y", data = df.lines.hor)
    geom_text(aes(
      x = x + 10,
      y = y + 0.1,
      label = labels,
      colour = colour
    ), data = df.lines.hor) 
  #scale_color_discrete("red", "lightgrey")
  
  df <- data.frame()
  c <- ggplot() + theme_void()
  
  
  #ggarrange(a,a2,ggarrange(
  # c,b, ncol = 2, nrow = 1, widths = c(0.015,1)),ncol=1,nrow=3)
  ggarrange(c,
            b,
            ncol = 2,
            nrow = 1,
            widths = c(0.015, 1))
  #return(data.frame(perc = df.lines.ver[,5], row.names=rownames(df.lines.ver)))
}

generate_figure_highlight_atrt <- function(new.sample.meta.score, indexRow) {
  if(is.null(indexRow)){indexRow=1}
  #temp.df <- readRDS(file = "https://github.com/hackingjpr/Idat-Shiny/blob/main/temp.df.rds")
  temp.df <- readRDS(file = "./atrt34.dist.rds")
  #df.cat.atrt <- readRDS(file = "https://github.com/hackingjpr/Idat-Shiny/blob/main/df.cat.atrt.rds")
  # df.cat.atrt <- readRDS(file = "./df.cat.atrt.rds")
  #comb.SDb.atrt <- readRDS(file = "https://github.com/hackingjpr/Idat-Shiny/blob/main/comb.SDb.atrt.rds")
  # comb.SDb.atrt <- readRDS(file = "./comb.SDb.atrt.rds")
  ggplot(aes(x = 1:nrow(temp.df), y = atrt34), data = temp.df) +
    geom_line() +
    scale_shape_manual(values = c(1, 4,  3)) +
    scale_color_manual(values = c('#E69F00', '#999999', "white")) +
    ylab("ATRT-34") +
    theme_minimal() +
    theme(
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      legend.position = "none"
    ) -> b
  
  df.lines.hor <-
    foreach(i = 1:length(new.sample.meta.score),
            .combine = rbind) %do% {
              data.frame(
                x = 0,
                xend = max(which(
                  temp.df$atrt34 < new.sample.meta.score[i]
                )),
                y = new.sample.meta.score[i],
                yend = new.sample.meta.score[i]
              )
            }
  df.lines.hor$labels <- names(new.sample.meta.score)
  
  df.lines.ver <-
    foreach(i = 1:length(new.sample.meta.score),
            .combine = rbind) %do% {
              data.frame(
                x = max(which(
                  temp.df$atrt34 < new.sample.meta.score[i]
                )),
                xend = max(which(
                  temp.df$atrt34 < new.sample.meta.score[i]
                )),
                y = new.sample.meta.score[i],
                yend = min(temp.df$atrt34)
              )
            }
  df.lines.ver$perc <-
    paste0(round(df.lines.ver$xend / length(temp.df$atrt34) * 100), "th")
  
  df.lines.ver$colour <- factor(ifelse(1:nrow(df.lines.ver)==indexRow,"highlight","no.highlight"), levels = c("highlight","no.highlight"))
  df.lines.hor$colour <- factor(ifelse(1:nrow(df.lines.hor)==indexRow,"highlight","no.highlight"), levels = c("highlight","no.highlight"))
  
  b <- b +
    geom_segment(
      aes(
        x = x,
        y = y,
        xend = xend,
        yend = yend,
        colour = as.character(colour)
      ),
      #colour = "red",
      linetype = "dashed",
      data = df.lines.hor
    ) +
    geom_segment(
      aes(
        x = x,
        y = y,
        xend = xend,
        yend = yend,
        colour = colour
      ),
      #colour = "red",
      linetype = "dashed",
      data = df.lines.ver
    ) +
    #geom_text_repel(aes(x = x+10, y = y+0.1, label = labels), direction = "y", data = df.lines.hor)
    geom_text(aes(
      x = x + 10,
      y = y + 0.1,
      label = labels,
      colour = colour
    ), data = df.lines.hor) 
  #scale_color_discrete("red", "lightgrey")
  
  df <- data.frame()
  c <- ggplot() + theme_void()
  
  
  #ggarrange(a,a2,ggarrange(
  # c,b, ncol = 2, nrow = 1, widths = c(0.015,1)),ncol=1,nrow=3)
  ggarrange(c,
            b,
            ncol = 2,
            nrow = 1,
            widths = c(0.015, 1))
  #return(data.frame(perc = df.lines.ver[,5], row.names=rownames(df.lines.ver)))
}

#generate_figure_highlight(c("sampleA" = 0.2, "sampleB" = 0.3, "sampleC" = 0.7), 1)

generate_figure_percentage_atrt <- function(new.sample.meta.score, indexRow){
  temp.df <- readRDS(file =  "./atrt34.dist.rds")

  
  df.lines.hor <-foreach(i = 1:length(new.sample.meta.score), .combine = rbind)%do%{
    data.frame(x=0,
               xend=max(which(temp.df$atrt34 <new.sample.meta.score[i])),
               y=new.sample.meta.score[i],
               yend=new.sample.meta.score[i])}
  df.lines.hor$labels <- names(new.sample.meta.score)
  
  df.lines.ver <-foreach(i = 1:length(new.sample.meta.score), .combine = rbind)%do%{
    data.frame(x=max(which(temp.df$atrt34 < new.sample.meta.score[i])),
               xend=max(which(temp.df$atrt34 < new.sample.meta.score[i])),
               y=new.sample.meta.score[i],
               yend=min(temp.df$atrt34))}
  df.lines.ver$perc <- paste0(round(df.lines.ver$xend/ length(temp.df$atrt34)*100),"th")
  
  return(df.lines.ver$perc[indexRow])
}

generate_figure_percentage_mrt <- function(new.sample.meta.score, indexRow){
  temp.df <- readRDS(file =  "./mrt13.dist.rds")
  
  
  df.lines.hor <-foreach(i = 1:length(new.sample.meta.score), .combine = rbind)%do%{
    data.frame(x=0,
               xend=max(which(temp.df$mrt13 <new.sample.meta.score[i])),
               y=new.sample.meta.score[i],
               yend=new.sample.meta.score[i])}
  df.lines.hor$labels <- names(new.sample.meta.score)
  
  df.lines.ver <-foreach(i = 1:length(new.sample.meta.score), .combine = rbind)%do%{
    data.frame(x=max(which(temp.df$mrt13 < new.sample.meta.score[i])),
               xend=max(which(temp.df$mrt13 < new.sample.meta.score[i])),
               y=new.sample.meta.score[i],
               yend=min(temp.df$mrt13))}
  df.lines.ver$perc <- paste0(round(df.lines.ver$xend/ length(temp.df$mrt13)*100),"th")
  
  return(df.lines.ver$perc[indexRow])
}

generate_figure_percentage_ecrt <- function(new.sample.meta.score, indexRow){
  temp.df <- readRDS(file =  "./ecrt32.dist.rds")
  
  
  df.lines.hor <-foreach(i = 1:length(new.sample.meta.score), .combine = rbind)%do%{
    data.frame(x=0,
               xend=max(which(temp.df$ecrt32 <new.sample.meta.score[i])),
               y=new.sample.meta.score[i],
               yend=new.sample.meta.score[i])}
  df.lines.hor$labels <- names(new.sample.meta.score)
  
  df.lines.ver <-foreach(i = 1:length(new.sample.meta.score), .combine = rbind)%do%{
    data.frame(x=max(which(temp.df$ecrt32 < new.sample.meta.score[i])),
               xend=max(which(temp.df$ecrt32 < new.sample.meta.score[i])),
               y=new.sample.meta.score[i],
               yend=min(temp.df$ecrt32))}
  df.lines.ver$perc <- paste0(round(df.lines.ver$xend/ length(temp.df$ecrt32)*100),"th")
  
  return(df.lines.ver$perc[indexRow])
}

# df.cat.atrt <- readRDS(file = "./df.cat.atrt.rds")
# comb.SDb.atrt <- readRDS(file = "./comb.SDb.atrt.rds")

# ggplot(aes(x=1:nrow(temp.df)), data = temp.df) +
#   geom_point(aes(y=time.os, color = factor(status.os), shape = factor(status.os)), size = 2, show.legend = FALSE) +
#   scale_shape_manual(values=c(1, 4,  3)) +
#   scale_color_manual(values=c('#E69F00','#999999',"white")) +
#   ylab("Time in days") +
#   # new_scale_color() +
#   # geom_point(aes(y=time.fruh*30, color = factor(status.fruh), shape = factor(status.fruh)), size = 4, show.legend = FALSE) +
#   # scale_color_manual(values=c('red','blue',"white")) +
#   theme_minimal() +
#   theme(axis.title.x=element_blank(),
#         axis.text.x=element_blank(),
#         axis.ticks.x=element_blank()) -> a
# 
# ggplot(aes(x=1:nrow(temp.df)), data = temp.df) +
#   # geom_point(aes(y=time.os, color = factor(status.os), shape = factor(status.os)), size = 4, show.legend = FALSE) +
#   # scale_shape_manual(values=c(1, 4,  3)) +
#   # scale_color_manual(values=c('#E69F00','#999999',"white")) +
#   # new_scale_color() +
#   geom_point(aes(y=time.fruh*30, color = factor(status.fruh), shape = factor(status.fruh)), size = 2, show.legend = FALSE) +
#   scale_color_manual(values=c('#E69F00','#999999',"white")) +
#   scale_shape_manual(values=c(1, 4,  3)) +
#   ylab("Time in days") +
#   theme_minimal() +
#   theme(axis.title.x=element_blank(),
#         axis.text.x=element_blank(),
#         axis.ticks.x=element_blank()) -> a2
# 
# ggplot(aes(x=1:nrow(temp.df),y=V1), data = temp.df) +
#   geom_line() +
#   scale_shape_manual(values=c(1, 4,  3)) +
#   scale_color_manual(values=c('#E69F00','#999999',"white")) +
#   ylab("ATRT-8") +
#   theme_minimal() +
#   theme(axis.title.x=element_blank(),
#         axis.text.x=element_blank(),
#         axis.ticks.x=element_blank()) -> b
