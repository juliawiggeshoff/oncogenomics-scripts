#!/bin/bash
dir="${PWD}/misc/"

mkdir -p ${dir}

## Download input files: https://www.cancerhotspots.org/#/download
## spreadhseet = Hotspot Results V2
## MAF = V2 Mutational Data (MAF)

## save snvs and indels, separately, without headers

tsv_indels_hotspots="${dir}/hotspots_v2_hg19_indels.spreadsheet.tsv"
tsv_snvs_hotspots="${dir}/hotspots_v2_hg19_snvs.spreadsheet.tsv"
maf_hotspots="${dir}/cancerhotspots.v2.hg19.maf"

py="mk_cancerhotspots.py"

## 1. Select hotspots spreadsheet for just mutations in lung samples, column with info changes between indels and snvs 

cat ${tsv_snvs_hotspots} | awk -F"\t" '$39 ~ /lung/' > ${dir}/lung_cancer.hotspots_v2_hg19_snvs.spreadsheet.tsv
cat ${tsv_indels_hotspots} | awk -F"\t" '$37 ~ /lung/' > ${dir}/lung_cancer.hotspots_v2_hg19_indels.spreadsheet.tsv

## 2. Select MAF file with hotspots for just mutations in lung samples; simplify the output

tail -n +3 ${maf_hotspots} | cut -f1,5,6,7,35,136,145 | grep lung | awk '{ print "chr"$2, $3, $4, $1, $6, $5 }' OFS="\t" > ${dir}/lung_cancer.simplified.cancerhotspots.v2.hg19.maf

## 3. Run python script for indels and snvs to generate a simplified tsv with the spreadsheet

python3 ${py} --tsv ${dir}/lung_cancer.hotspots_v2_hg19_snvs.spreadsheet.tsv --variant snv --output ${dir}/maf_lookup_table.lung_cancer.hotspots_v2_hg19_snvs.spreadsheet.tsv
python3 ${py} --tsv ${dir}/lung_cancer.hotspots_v2_hg19_indels.spreadsheet.tsv --variant indel --output ${dir}/maf_lookup_table.lung_cancer.hotspots_v2_hg19_indels.spreadsheet.tsv

## 4. Concatenate, sort and output unique lines 

cat ${dir}/maf_lookup_table.lung_cancer.hotspots_v2_hg19_snvs.spreadsheet.tsv ${dir}/maf_lookup_table.lung_cancer.hotspots_v2_hg19_indels.spreadsheet.tsv | \
    sort -k1,1V -k2,2n | uniq > ${dir}/maf_lookup_table.lung_cancer.hotspots_v2_hg19_indels_snvs.spreadsheet.tsv

## 5. Run python script again to compare maf with modified tsv
python3 ${py} --tsv ${dir}/maf_lookup_table.lung_cancer.hotspots_v2_hg19_indels_snvs.spreadsheet.tsv \
    --maf ${dir}/lung_cancer.simplified.cancerhotspots.v2.hg19.maf --output ${dir}/unsort.lung_cancer.hotspots_v2_hg19_indels_snvs.bed

## 6. Sort bed file and add extra columns before liftover, otherwise browser tool truncates the file

sort -k1,1V -k2,2n -k5,5V ${dir}/unsort.lung_cancer.hotspots_v2_hg19_indels_snvs.bed | \
    uniq | awk -F'\t' -v OFS='\t' '{print $0, 1, 1, $6}' > ${dir}/sort.lung_cancer.hotspots_v2_hg19_indels_snvs.bed

## 7. run liftover here: https://genome.ucsc.edu/cgi-bin/hgLiftOver

## 8. Keep columns 1-5 and 9 from the lifted-over bed file

# cut -f1-5,9 ${dir}/LONG.lung_cancer.hotspots_v2_liftOver_hg38_indels_snvs.bed > ${dir}/lung_cancer.hotspots_v2.liftOver_hg38.bed

