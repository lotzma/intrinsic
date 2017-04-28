## R functions for preparing the data set for analysis

library(nnls)
library(mice)

# Extract exposure of mutation catalog M to signatures P using non-negative least squares

extract <- function(M, P) {
  # Should include bootstrapping in this procedure?
  f <- function(y) nnls(P,y)$x
  X <- apply(M, 2, f)
  return(X)
}
    
# Prepare a data set. Loads the columns with names listed in 'header' into a dataframe
df <- read.csv('data/ncomms12064-s3.csv', header=TRUE, skip=1)
  
# Include BRAF/NRAS/NF1 information
tcga_meta <- read.csv('data/brafnrasnf2.csv', header=FALSE)
    
indices <- match(tcga_meta$V1, df$BARCODE)
indices <- indices[!is.na(indices)]
df <- df[indices,]
indices <- match(df$BARCODE, tcga_meta$V1)
indices <- indices[!is.na(indices)]
df$Cohort <- tcga_meta[indices, c('V2')]

# Label 'other'

# Remove 'OTHER' columns as this consists of mostly NA
df <- df[, !(names(df) %in% c("OTHER"))]
# Drop rows with missing data (may have to eventually replace with clever imputation)
#df <- df[complete.cases(df),]
  
## Include the exposures to signatures into the data

# Extract rows containing mutations
mutations <- df[, grep(".in.", colnames(df))]
alphabet <- colnames(mutations)
# Read in the signatures and sort the coordinates in the right order
all_sigs <- read.csv('data/signatures_probabilities.txt', sep="\t")
all_sigs$alpha <- gsub(' ', '', sub('>','.',paste(all_sigs[,1],'.in.',all_sigs[,2])))
all_sigs <- all_sigs[match(alphabet, all_sigs$alpha),]
# Signature matrix
P <- data.matrix(all_sigs[,-c(1:3,34:41)])
# Extract exposures to signatures
M <- t(data.matrix(mutations))
E <- t(extract(M, P[,c(1,5,7)]))
# Determine the number of SNVs attributed to each signture
df$sig1 <- round(rowSums(outer(E[,1],P[,1])))
df$sig5 <- round(rowSums(outer(E[,2],P[,2])))
df$sig7 <- round(rowSums(outer(E[,3],P[,3])))
#df$sig1 <- E[,1]
#df$sig5 <- E[,2]
#df$sig7 <- E[,3]
        
labels <- c('age_at_diagnosis', 'breslow', 'gender', 'clark_level', 'ulceration',
           'tissue_source_site', 'body_area', 'tissue_type', 'Cohort',
           'V60L', 'D84E', 'V92M', 'R142H', 'R151C', 'I155T', 'R160W', 'R163Q', 'D294H', 
           'other.rscore', 'rgeno', 'sig1', 'sig5', 'sig7', 'totalNonUV', 'totalSNV')
data <- df[, labels]
data <- data[!is.na(data$age_at_diagnosis),]
#data <- data[complete.cases(data),]