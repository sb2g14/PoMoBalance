library(data.table)
library(Biostrings)

cwd <- getwd()
path <- file.path(cwd,"combine_for_counts")

seq_type <- c("Melanogaster", "Santomea", "Yakuba", "Simulants", "Erecta_Dark", "Erecta_Light")
for (p in 1:length(seq_type)){
    type <- seq_type[p]
    file.name <- paste0(type, ".txt")
    fasta.file <- file.path(path, file.name)

    fasta <- readDNAStringSet(fasta.file)
    n_ind <- length(fasta)
    names(fasta) <- paste0(type, "-", seq(n_ind))

    fasta.output <- file.path(path, paste0("data_renamed/", type, "_renamed.txt"))
    writeXStringSet(fasta, fasta.output, append=FALSE,
                compress=FALSE, compression_level=NA, format="fasta")
}

