####Script for standard analysis of FACS data, later to be included in more generic functions
####Stared on 16-October-2018

##User Input----


#Defining Directories

#location of data
data_dir <- '~/Desktop/WHRI Summer Project 2018/Raw Data + FlowJo + Prism/p-akt/p-AKT 11 10 18 N=2/'
#directory of gating template .csv file
gating_dir <- 'extdata/gating_template/IC p-AKT.csv'
#location for output files to be saved
save_dir <- '~/Desktop/'
#location of markers and channel name data frame containing all names/channels
marker_dir <- '/Library/Frameworks/R.framework/markers.txt'

#markers used in staining panel
markers_used <- c('FSC-A','SSC-A', 'CD27-FITC', 'CD8-PerCP', 'CD4-Pecf594', 'CD45RA-Bv605', 'p-akt-PE')

hide_pop <- c('CD4+', 'CD8+', 'CD4+CD8+', 'CD4-CD8-', 
              'CD4+CD8-/CD45RA+', 'CD4+CD8-/CD27+', 'CD4-CD8+/CD45RA+', 'CD4-CD8+/CD27+')


##Improting and formatting compensation data----

#loading required packages
packages <- c('cytofkit', 'openCyto', 'flowViz', 'ggcyto', 'gridExtra', 'grid', 'ggplot2')
lapply(packages, library, character.only=TRUE)

setwd(data_dir)

#import compensation controls
compensation_data <- read.flowSet(pattern='Compensation')

#import marker data.frame
markers <- read.delim(marker_dir, header = F)

#select specific markers used from generic marker data.frame
matches <- cbind(as.character(markers$V1[match(markers_used, markers$V2)]), 
                 as.character(markers$V2[match(markers_used, markers$V2)]))

#set colnames to marker used
colnames(compensation_data) <- matches[,2][match(colnames(compensation_data), matches[,1])]

#bind markers used by row, matching to channel
chanels <- rbind(as.vector(markers[,1]),as.vector(markers[,2]))

#name samples by their marker
sampleNames(compensation_data) <- gsub(pattern='(.*)(Unstained)(.*)', replacement = 'UN', x=sampleNames(compensation_data))
#need to make this line more generic, some sample names have (,1f,) instead
sampleNames(compensation_data) <- gsub(pattern='(.*\\s)(.*\\s)(.*)(,\\d{1}f,)(\\d*)(.*)', replacement = '\\2\\3/\\5-A', x=sampleNames(compensation_data))

sampleNames(compensation_data) <- chanels[2,][match(sampleNames(compensation_data), chanels[1,])]

#naming metadata
pData(compensation_data)$name <- sampleNames(compensation_data)


##Calculating Spillover Matrix and plotting transformed data scatter plot----

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
class(compensated_ctrls)

for (i in 1:length(compensated_ctrls)) {
  compensated_ctrls[[i]] <- transform(compensated_ctrls[[i]], estimateLogicle(compensated_ctrls[[i]], channels=colnames(compensated_ctrls)))
}

print(splom(compensated_ctrls[[1]]))

##Organise Sample Data----

#Import raw data only
raw_data <- read.flowSet(pattern = '.fcs')
raw_data <- raw_data[grep(pattern = '(^[[:upper:]]{2}\\d{2})', x=sampleNames(raw_data))]

#Set column names
matches <- cbind(as.character(markers$V1[match(markers_used, markers$V2)]), 
                 as.character(markers$V2[match(markers_used, markers$V2)]))

colnames(raw_data) <- matches[,2][match(colnames(raw_data), matches[,1])]

#apply compensation to data
raw_data_compensated <- compensate(raw_data, compensation_matrix)

#select variables excluding FSC and SSC
vars <- colnames(raw_data_compensated)
vars <- vars[-grep('FSC', vars)]
vars <- vars[-grep('SSC', vars)]

#transfor data on vars channels
for (i in 1:length(raw_data_compensated)) {
  raw_data_compensated[[i]] <- transform(raw_data_compensated[[i]], estimateLogicle(raw_data_compensated[[i]], channels=vars))
}

processed_data <- raw_data_compensated

colnames(processed_data) <- as.character(markers$V3[match(colnames(processed_data), markers$V2)])

##Gating----

#Import Gating template 
gating_template <- gatingTemplate(system.file(gating_dir, package = 'openCyto'))

#Create gating set
data_gs <- GatingSet(processed_data)

#apply gating hierarchy to data
gating(gating_template, data_gs, mc.cores=2, parallel_type='multicore')

#hide unwanted populations
lapply(hide_pop, function(thisNode)setNode(data_gs, thisNode, FALSE))

#print gating hierarchy pdf to desktop
pdf(file='/Users/johannesschroth/Desktop/Gating.pdf', height = 6, width = 8.48 )
plot(data_gs)
for (i in 1:length(data_gs)) {plotGate(data_gs[[i]], path=2)}
dev.off()



data_gs[[1]]

for (i in 1:length(data_gs)){
  CD4_data <- as.matrix(exprs(getData(data_gs[[i]], 'CD4+CD8-')))
}

CD4_data

str(data)
