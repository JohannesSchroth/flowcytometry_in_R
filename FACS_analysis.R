####Script for standard analysis of FACS data, later to be included in more generic functions
####16-October-2018

##User Input----

#Selecting Directories
#location of data
data_dir <- '~/Desktop/WHRI Summer Project 2018/Raw Data + FlowJo + Prism/p-akt/p-AKT 11 10 18 N=2/'
#directory of gating template .csv file
gating_dir <- 'extdata/gating_template/G.csv'
#location for output files to be saved
save_dir <- '~/Desktop/'


#markers used in staining panel
markers_used <- c('FSC-A','SSC-A', 'CD27-FITC', 'CD8-PerCP', 'CD4-Pecf594', 'CD45RA-Bv605', 'p-akt-PE')

#loading required packages
library(cytofkit)
library(openCyto)
library(flowViz)
library(ggcyto)
library(gridExtra)
library(grid)
library(ggplot2)


#Improting and formatting compensation data----
setwd(data_dir)

#import compensation controls
compensation_data <- read.flowSet(pattern='Compensation')
#double check all samples loaded
compensation_data

#import marker data.frame
markers <- read.delim('/Library/Frameworks/R.framework/markers.txt', header = F)
#select specific markers used from generic marker data.frame
matches <- cbind(as.character(markers$V1[match(markers_used, markers$V2)]), 
                 as.character(markers$V2[match(markers_used, markers$V2)]))

matches

#set colnames to actual marker used
colnames(compensation_data) <- matches[,2][match(colnames(compensation_data), matches[,1])]
colnames(compensation_data)

#bind markers used by row, matching to channel
chanels <- rbind(as.vector(markers[,1]),as.vector(markers[,2]))
chanels
#name samples by their marker
sampleNames(compensation_data) <- gsub(pattern='(.*)(Unstained)(.*)', replacement = 'UN', x=sampleNames(compensation_data))
#need to make this line more generic, some sample names have (,1f,) instead
sampleNames(compensation_data) <- gsub(pattern='(.*\\s)(.*\\s)(.*)(,2f,)(\\d*)(.*)', replacement = '\\2\\3/\\5-A', x=sampleNames(compensation_data))
compensation_data
sampleNames(compensation_data) <- chanels[2,][match(sampleNames(compensation_data), chanels[1,])]

#naming metadata
pData(compensation_data)$name <- sampleNames(compensation_data)

#Calculating Spillover Matrix and plotting transformed data scatter plot----

#match sample names to colnames for matrix order
matc <- na.omit(match(pData(compensation_data)$name, colnames(compensation_data)))

#apply function to order compensation data in specific order
compensation_ctrls <- fsApply(compensation_data, function(ff) ff[,matc])

#extract unstained samples
unstained_ctrl <- pData(compensation_ctrls)[grep(pattern='UN', x=sampleNames(compensation_ctrls), fixed = T),]

#calculate spillover matrix
spillover_matrix <- spillover(x = compensation_ctrls,
                          unstained=unstained_ctrl,
                          stain_match='ordered',
                          method='mean')

#calculate compensation matrix from spillover, should be inverse? but exactly same here
compensation_matrix <- compensation(spillover_matrix)

#apply compensation to compensation samples, could do an overlay thing to compare uncompensated w/ compensated
compensated_ctrls <- compensate(x=compensation_ctrls, compensation_matrix)

class(compensation_ctrls)

tData <- transform(compensated_ctrls, transformList(colnames(compensated_ctrls), 'linearTransform'))

print(splom(compensated_ctrls[[1]]))

#Organise Sample Data----

setwd(data_dir)

#Import raw data only
raw_data <- read.flowSet(pattern = '.fcs')
raw_data <- raw_data[grep(pattern = '(^[[:upper:]]{2}\\d{2})', x=sampleNames(raw_data))]
sampleNames(raw_data)

#Set column names
colnames(raw_data) <- matches[,2][match(colnames(raw_data), matches[,1])]
colnames(raw_data)

#apply compensation to data
raw_data_compensated <- compensate(raw_data, compensation_matrix)

#select variables excluding FSC and SSC
vars <- colnames(raw_data_compensated)
vars <- vars[-grep('FSC', vars)]
vars <- vars[-grep('SSC', vars)]



#transfor data on vars channels
for (i in 1:length(raw_data_compensated)) {
  d <- transform(raw_data_compensated[[i]], estimateLogicle(raw_data_compensated[[i]], channels=vars))
  processed_data[[i]] <- flowSet(d) 
}

processed_data

#Gating----

#Import Gating template 
gating_template <- gatingTemplate(system.file(gating_dir, package = 'openCyto'))

#Create gating set
data_gs <- GatingSet(processed_data)

#apply gating hierarchy to data
gating(gating_template, dara_gs, mc.cores=2, parallel_type='multicore')






