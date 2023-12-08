#!/bin/bash

# Move generated VCFs to a folder
mkdir vcf
mv p1.vcf p2.vcf p3.vcf p4.vcf vcf/
# zip FASTA and VCFs with bi-allelic sites only
gzip anc.fa
cd vcf/
bgzip -c p1.vcf > p1_t.vcf.gz
tabix p1_t.vcf.gz
bgzip -c p2.vcf > p2_t.vcf.gz
tabix p2_t.vcf.gz
bgzip -c p3.vcf > p3_t.vcf.gz
tabix p3_t.vcf.gz
bgzip -c p4.vcf > p4_t.vcf.gz
tabix p4_t.vcf.gz
bcftools view -m2 -M2 -v snps p1_t.vcf.gz > p1_bi.vcf
bcftools view -m2 -M2 -v snps p2_t.vcf.gz > p2_bi.vcf
bcftools view -m2 -M2 -v snps p3_t.vcf.gz > p3_bi.vcf
bcftools view -m2 -M2 -v snps p4_t.vcf.gz > p4_bi.vcf
bgzip -c p1_bi.vcf > p1.vcf.gz
tabix p1.vcf.gz
bgzip -c p2_bi.vcf > p2.vcf.gz
tabix p2.vcf.gz
bgzip -c p3_bi.vcf > p3.vcf.gz
tabix p3.vcf.gz
bgzip -c p4_bi.vcf > p4.vcf.gz
tabix p4.vcf.gz

# Run cflib to generate count file from FASTA and VCF
cd /Users/sb442/Documents/Tools/PoMoBalance/SLiM/data_generation/vcf/
vcf_files=(p1.vcf.gz p2.vcf.gz p3.vcf.gz p4.vcf.gz)
python3 FastaVCFToCounts.py /Users/Documents/Tools/PoMoBalance/SLiM/data_generation/anc.fa.gz "${vcf_files[@]}" /Users/Documents/Tools/PoMoBalance/SLiM/data_generation/counts.cf --merge --ploidy 2
