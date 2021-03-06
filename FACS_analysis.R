####Script for standard analysis of FACS data, later to be included in more generic functions
####Stared on 16-October-2018

##User Input----

#Defining Directories

#location of data
data_dir <- '~/Desktop/'
#directory of gating template .csv file
gating_dir <- '~/Desktop/G.csv'
#location for output files to be saved
save_dir <- '~/Desktop/'
#location of markers and channel name data frame containing all names/channels
marker_dir <- '~/Desktop/markers.txt'
#markers used in staining panel
markers_used <- c('FSC-A', 'FSC-H', 'FSC-W', 'SSC-A', 'SSC-H', 'SSC-W', 'CD27-FITC', 'CD8-PerCP', 'CD4-Pecf594', 'CD45RA-Bv605', 'CD220-APC')
#populations to hide in gating hierarchy plot
hide_pop <- c('CD4+', 'CD8+', 'CD4+CD8+', 'CD4-CD8-', 
              'CD4+CD8-/CD45RA+', 'CD4+CD8-/CD27+', 
              'CD4-CD8+/CD45RA+', 'CD4-CD8+/CD27+')


##Improting and formatting compensation data----

#loading required packages
packages_1 <- c('cytofkit', 'openCyto', 'flowViz', 'ggcyto', 'gridExtra', 'grid', 'ggplot2')
lapply(packages_1, library, character.only=TRUE)

#set working directory
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
sampleNames(compensation_data) <- gsub(pattern='(.*\\s)(.*_)(.*)(,\\d{1}f,)(\\d*)(.*)', replacement = '\\3/\\5-A', x=sampleNames(compensation_data))
sampleNames(compensation_data) <- chanels[2,][match(sampleNames(compensation_data), chanels[1,])]

compensation_data

#naming metadata
pData(compensation_data)$name <- sampleNames(compensation_data)

##Calculating Spillover Matrix and plotting transformed data scatter plot----

#match sample names to colnames for matrix order
matc <- na.omit(match(pData(compensation_data)$name, colnames(compensation_data)))

#apply function to order compensation data in specific order
compensation_ctrls <- fsApply(compensation_data, function(ff) ff[,matc])

#extract unstained samples
unstained_ctrl <- grep(pattern='UN', x=sampleNames(compensation_ctrls), fixed = T)

#calculate spillover matrix
spillover_matrix <- spillover(x = compensation_ctrls,
                          unstained=unstained_ctrl,
                          stain_match='ordered',
                          method='mean')

#calculate compensation matrix from spillover, should be inverse? but exactly same here
compensation_matrix <- compensation(spillover_matrix)

#apply compensation to compensation samples, could do an overlay thing to compare uncompensated w/ compensated
compensated_ctrls <- compensate(x=compensation_ctrls, compensation_matrix)

for (i in 1:length(compensated_ctrls)) {
  compensated_ctrls[[i]] <- transform(compensated_ctrls[[i]], estimateLogicle(compensated_ctrls[[i]], channels=colnames(compensated_ctrls)))
}

comp_mat_plot <- splom(compensated_ctrls[[1]])
comp_mat_plot
##Organise Sample Data----

#Import raw data only
raw_data <- read.flowSet(pattern = '.fcs')
raw_data <- raw_data[grep(pattern = '(^[[:upper:]]{2}\\d{2})', x=sampleNames(raw_data))]

#Set column names
matches <- cbind(as.character(markers$V1[match(markers_used, markers$V2)]), 
                 as.character(markers$V2[match(markers_used, markers$V2)]))

colnames(raw_data) <- matches[,2][match(colnames(raw_data), matches[,1])]
raw_data
#apply compensation to data
raw_data_compensated <- compensate(raw_data, compensation_matrix)

#select variables excluding FSC and SSC
vars <- colnames(raw_data_compensated)
vars <- vars[-grep('FSC', vars)]
vars <- vars[-grep('SSC', vars)]

#transfor data on vars channels, using logicle transformation, try others?
for (i in 1:length(raw_data_compensated)) {
  raw_data_compensated[[i]] <- transform(raw_data_compensated[[i]], estimateLogicle(raw_data_compensated[[i]], channels=vars))
}
processed_data <- raw_data_compensated
colnames(processed_data) <- as.character(markers$V3[match(colnames(processed_data), markers$V2)])

##Quality Asessment----

raw_data_compensated[[1]]

library(Biobase)

processed_data[[1]]

for (i in 1:length(processed_data)){
  data <- AnnotatedDataFrame(as.data.frame(exprs(processed_data[[i]])))
}

#Need to find way to annotate data frame names


##Gating----

#Import Gating template
gating_template <- gatingTemplate(gating_dir)

#Create gating set
data_gs <- GatingSet(processed_data)

#apply gating hierarchy to data
gating(gating_template, data_gs, mc.cores=2, parallel_type='multicore')

#hide unwanted populations
lapply(hide_pop, function(thisNode)setNode(data_gs, thisNode, FALSE))

#print gating hierarchy pdf to desktop
pdf(file='~/Desktop/Gating.pdf', height = 6, width = 8.48 )
print(comp_mat_plot)
plot(data_gs)
for (i in 1:length(data_gs)) {plotGate(data_gs[[i]], path=2)}
dev.off()

data_gs[[17]]

#can now extract gated cell pop data from gating set, for further analysis
#may be easier to keep using gating set for most analyses, but Rtsne requires data frame
#example of CD4 data extraction below

singlets_1 <- as.data.frame(exprs(getData(data_gs[[17]], 'singlets')))


#save data as .csv file to not have to re-gate data every time
setwd(save_dir)
write.csv(singlets_1, 'singlets_1')

