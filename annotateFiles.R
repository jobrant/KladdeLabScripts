setwd("~/Desktop/KLADDE_LAB/FOR_ANQI/ANQI_ANNOS/")

library(tidyverse)
library(openxlsx)

# Processing ATAC-Seq files is listed first. Scripts for processing metilene
# data are farther down.

# Get file names for ATAC-Seq differential analysis Excel files.
atac_files <- list.files(pattern = "^M.*\\.xlsx$")

# Generates a TSV peak file for HOMER for each sheet in each Excel file.
processAtacFiles <- function(x){
        for (i in 1:length(atac_files)){
                sheets <- getSheetNames(x[i])
                temp <- read.xlsx(x[i], sheet = sheets[1])
                temp$Unique_ID <- paste(temp$Chrom, temp$Start, sep = ".")
                temp$Strand <- 0
                temp <- temp[,c(10, 1, 2, 3, 11, 4:9)]
                name <- substr(atac_files[i], 0, (nchar(atac_files[i]) - 5))
                write.table(temp, file = paste(paste("TSVs", name, sep = "/"), 
                                               sheets[1], "tsv", sep = "."), 
                            quote = F, sep = '\t', row.names = F, col.names = T)
                temp2 <- read.xlsx(x[i], sheet = sheets[2])
                temp2$Unique_ID <- paste(temp2$Chrom, temp2$Start, sep = ".")
                temp2$Strand <- 0
                temp2 <- temp2[,c(10, 1, 2, 3, 11, 4:9)]
                write.table(temp2, file = paste(paste("TSVs", name, sep = "/"), 
                                                sheets[2], "tsv", sep = "."), 
                            quote = F, sep = '\t', row.names = F, col.names = T)
        }
}
# Now you can transfer the files just generated to hpc to run homer to annotate
# Then download the annotated files.

# This will combine the annotation data with the original data, yielding a 
# complete dataset in Excel format
anno_files <- readLines("./TSVs/filenames.txt")
merge_atacAnnos <- function(x){
        for (i in 1:length(anno_files)){
                data <- read.table(paste(paste("TSVs", x[i], sep = '/'), "tsv", 
                                         sep = '.'), header = T, sep = '\t')
                annos <- read_tsv(paste(paste("TSVs", x[i], sep = '/'), 
                                        "annotated.tsv", sep = "."))
                colnames(annos)[1] <- "Unique_ID"
                annos <- annos[,c(1, 8:19)]
                annos <- as.data.frame(annos)
                combined <- merge(data, annos, by = "Unique_ID", all = T)
                name <- strsplit(x[i], "[.]")
                write.xlsx(combined, file = paste(x[i], "annotated.xlsx", 
                                                  sep = '.'), 
                           sheetName = name[[1]][5])
        }
}

# This will combine the atac-seq contrasts into one Excel workbook with
# multiple worksheets (as in the original data)
atac_filenames <- readLines("./TSVs/atac_filenames.txt")
mergeSheets <- function(x){
        for (i in 1:length(atac_filenames)){
                name_split <- strsplit(x[i], "[.]")
                temp1 <- read.xlsx(paste(paste(paste(x[i], "Genes", sep = "."), 
                                               name_split[[1]][1], sep = "-"), 
                                         "annotated.xlsx", sep = "."))
                temp2 <- read.xlsx(paste(paste(paste(x[i], "Genes", sep = "."), 
                                               name_split[[1]][3], sep = "-"), 
                                         "annotated.xlsx", sep = "."))
                wb <- createWorkbook(paste(x[i], "annotated", "xlsx", sep = '.'))
                addWorksheet(wb, paste("Genes", name_split[[1]][1], sep = '-'))
                addWorksheet(wb, paste("Genes", name_split[[1]][3], sep = '-'))
                writeData(wb, 1, temp1)
                writeData(wb, 2, temp2)
                saveWorkbook(wb, file = (paste(x[i], "annotated", "xlsx", 
                                               sep = '.')), overwrite = T)
        }
}



# ********************* Metilene processing below *********** #

# Get file names for Metilene differential methylation analysis Excel files.
meti_files <- list.files(pattern = "^G.*\\.xlsx$|^H.*\\.xlsx$")

# Generates a TSV peak file for HOMER for each sheet in each Excel file.
processMetiFiles <- function(x){
        for (i in 1:length(meti_files)){
                temp <- read.xlsx(x[i], sheet = 1)
                temp$Unique_ID <- paste(temp$chr, temp$start, sep = ".")
                temp$Strand <- 0
                temp <- temp[,c(15, 1, 2, 3, 16, 4:14)]
                name <- substr(meti_files[i], 0, (nchar(meti_files[i]) - 5))
                write.table(temp, file = paste(paste("TSVs", name, sep = "/"), 
                                               "tsv", sep = "."), quote = F, 
                            sep = '\t', row.names = F, col.names = T)
        }
}

# Now you can transfer the files just generated to hpc to run homer to annotate

# This will combine the annotation data with the original data, yielding a 
# complete dataset in Excel format
meti_annoFiles <- readLines("./TSVs/meti_files.txt")
merge_metiAnnos <- function(x){
        for (i in 1:length(meti_annoFiles)){
                data <- read.table(paste(paste("TSVs", x[i], sep = '/'), "tsv", 
                                         sep = '.'), header = T, sep = '\t', 
                                   comment.char = "")
                annos <- read_tsv(paste(paste("TSVs", x[i], sep = '/'), 
                                        "annotated.tsv", sep = "."))
                colnames(annos)[1] <- "Unique_ID"
                annos <- annos[,c(1, 8:19)]
                annos <- as.data.frame(annos)
                combined <- merge(data, annos, by = "Unique_ID", all = T)
                name <- strsplit(x[i], "[.]")
                write.xlsx(combined, file = paste(x[i], "annotated.xlsx", 
                                                  sep = '.'), 
                           sheetName = name[[1]][1])
        }
}

