### Process data from BLAST

- Download aligned sequence from BLAST
    
- Name "Melanogaster" or "Si,mulans", etc, and copy to folder `combine_for_counts`
    
- Run `process_data_for_counts.R` to rename individuals in fasta file
    
- Combine files in `combine_for_counts/data_renamed` with `cat *.txt > ../Drosophila.fasta`
    
- Align with `mafft --reorder --adjustdirection --anysymbol --leavegappyregion --auto /Users/Documents/Tools/PoMoBalance/Drosoplila_analysis/data_analysis/combine_for_counts/Drosophila.fasta > /Users/Documents/Tools/PoMoBalance/Drosoplila_analysis/data_analysis/combine_for_counts//Drosophila_aligned_raw.fasta`
    
- To remove `_R_` from the fasta names run `Drosophila/rename_after_align.R`
    
- To filter fasta data run `Filter_fasta.R`
    
- To create counts file run using `cflib` package run `generate_counts_file.sh`
    
- To filter counts file keeping only non-zero bi-alelic and less sites run `filter_muli_alelic_sites_Drosophila.R` (there is a bug, it doesn't include the first line for some reason, add `COUNTSFILE NPOP POP_NUM NSITES SITES_NUM` manually)
    
- Convert to PoMo with `counts_to_pomo_states_converter.R`