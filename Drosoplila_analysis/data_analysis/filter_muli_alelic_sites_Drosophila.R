# Save bi-alelic cites only

library(data.table)

path <- getwd()

# Select the datafile
count.file <- "counts_Drosophila.cf"
full.file <- file.path(path, count.file)
# Read datatable from the data file 
dt <- fread(full.file)

# Add regex to match mono- and bi-alelic sites only (e.g. 5,0,0,0 or 0,0,1,1)
re<-"(?:.*(?:\\b(?:0)\\b)){2}"

# Filter out the enries that don't match regex
filtered <- dt[grepl(re, Erecta_Dark) & grepl(re, Erecta_Light) & grepl(re, Santomea) & grepl(re, Yakuba) & grepl(re, Melanogaster)& grepl(re, Simulans)]
re.0 <- "0,0,0,0"
filtered1 <- filtered[!(grepl(re.0, Erecta_Dark) | grepl(re.0, Erecta_Light) | grepl(re.0, Santomea) | grepl(re.0, Yakuba) | grepl(re.0, Melanogaster) | grepl(re.0, Simulans))]
filtered1[,CHROM:=as.character(CHROM)]
filtered1[,CHROM:="X"]

# Save result
filename <- file.path(path, "counts_Drosophila_filtered.cf")
writeLines(c(paste0("COUNTSFILE NPOP ", ncol(filtered1)-2, " NSITES ", nrow(filtered1)), "CHROM POS Erecta_Dark Erecta_Light Santomea Yakuba Melanogaster Simulans"), filename)
fwrite(filtered1, file=filename, append=TRUE, sep="\t")