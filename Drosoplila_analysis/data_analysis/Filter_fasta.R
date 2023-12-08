library(data.table)
library(Biostrings)

cwd <- getwd()
path <- file.path(cwd,"combine_for_counts")

file.name <- "Drosophila_aligned.fasta"
fasta.file <- file.path(path, file.name)
fasta <- readDNAStringSet(fasta.file)
freq <- alphabetFrequency(fasta)
missing <- freq[,ncol(freq)-2] + freq[,ncol(freq)-3]
n.letters <- width(fasta)[1]
filter <- missing/n.letters
ind <- which(filter < 0.5)
fasta.filtered <- fasta[ind]

fasta.output <- file.path(path, "Drosophila_filtered_test.fasta")
writeXStringSet(fasta.filtered, fasta.output, append=FALSE,
            compress=FALSE, compression_level=NA, format="fasta")