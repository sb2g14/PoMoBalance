library(data.table)
library(Biostrings)

cwd <- getwd()
path <- file.path(cwd,"combine_for_counts")

file.name <- "Drosophila_aligned_raw.fasta"
fasta.file <- file.path(path, file.name)

fasta <- readDNAStringSet(fasta.file)
vec <- c()
for (p in 1:length(names(fasta))){
    name <- names(fasta[p])
    if (grepl("_R_", name)){
        vec <- c(vec, substr(name,4,nchar(name)))
    } else {
        vec <- c(vec, name)
    }
}
names(fasta) <- vec
fasta.output <- file.path(path,"Drosophila_aligned.fasta")
writeXStringSet(fasta, fasta.output, append=FALSE,
                compress=FALSE, compression_level=NA, format="fasta")

