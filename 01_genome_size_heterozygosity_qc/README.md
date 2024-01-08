
# Organzing genomic reads (stLRF/NGS) of the two enset landraces: landrace Mazia and Bedadeti

```sh
prog=~/apps/jellyfish/
wgs_dir=home/data/wgs

# Download SRA files of landrace Bedadeti (SRR1515268, SRR1515269)  and convert it to fastq format

 esearch -db sra -query PRJNA432894 | efetch -format runinfo | cut -d ',' -f 1 | grep SRR | grep -E 'SRR1515268|SRR1515269' | xargs -n 1 -P 12 fastq-dump --split-files --gzip --skip-technical

zcat $wgs_dir/trim_out/SRR1515269_1.fastq.gz $wgs_dir/trim_out/SRR1515268_1.fastq.gz > $wgs_dir/bedadeti_1.fastq.gz
zcat $wgs_dir/trim_out/SRR1515269_2.fastq.gz $wgs_dir/trim_out/SRR1515268_2.fastq.gz > $wgs_dir/bedadeti_2.fastq.gz

rm $wgs_dir/trim_out/SRR151526?_?.fastq.gz

# For landrace Mazia, remove stLRF adapter sequences using 
# Run [stLFR2Supernova (v2.1.1)](https://github.com/BGI-Qingdao/stlfr2supernova_pipeline) and merge them with Mazia NGS reads

## Create output directory
mkdir stLFR_assemble
cd stLFR_assemble

## set variables
data_dir=/home/data
out_path=output
out_path_step1=output_step1
max_memory=252G
threads=64
apps_dir=/home/apps/stLFRdenovo-v1.0.5

## Run stLFRdenovo
$apps_dir/stLFRdenovo_pipeline.sh \
-f $data_dir/stLFR_raw_data/CL200149738_L01_read_1.fq.gz  $data_dir/stLFR_raw_data/CL200149738_L01_read_1.fq.gz \
-o $out_path \
-m $max_memory \
-t $threads 

# Combine the adapter trimmed stLFR reads 100-bp raw paired-end reads

zcat stLFR_assemble/output/read.1.fq.gz.clean.gz $data_dir/WGS_data/DP8400009105TR_L01_520_1.fq.gz > $wgs_dir/mazia.1.fq.gz
zcat stLFR_assemble/output/read.2.fq.gz.clean.gz $data_dir/WGS_data/DP8400009105TR_L01_520_2.fq.gz > $wgs_dir/mazia.2.fq.gz 

ls -lht $wgs_dir/*.?.fq.gz | awk '{print $9}' | sed 's/...fq.gz//g' | uniq > accession_list

cat accession_list
mazia 
bedadeti
```

# Trimming adapter sequences and poor quality reads, and generating k-mer frequency tables and GenomeScope statistics. 


```bash
# Run jellyfish and genomescope2
for K in 17 21 27 31;
do
# This will generate kmer count and kmer histogram for esimated genome size of 580M and with estimated 43X coverage.
    for L in $(cat accession_list)
        do
		
		# Trim adapter and poor quality reads 
		trim_galore --cores 40 --quality 30 --fastqc -out jungle_seed --paired ${i}.1.fq.gz ${i}.2.fq.gz
		
        # generate Kmer tabel 
        ${prog}/jellyfish-linux count -C -m $k -s 580M --bf-size 25G -t 16 <(zcat ${K}.${L}_1_val_1.fq.gz) <(zcat ${K}${L}_2_val_2.fq.gz) -o ${L}.${K}_mer.reads.jf
        
        # generate Kmer histogram
        ${prog}/jellyfish-linux histo -t 16 ${L}.${K}_mer.reads.jf > ${L}.${K}_mer.histo
        rm -rf ${L}.${K}_mer.reads.jf

        # Run genomescope 
        genomescope2 \
        - ${L}${K}_mer.histo \
        -o genomescope_${i}_k17 \
        -p 2 \
        -k ${K} \
        --testing \
        --fitted_hist
        done
done
````

# Assess assembly contiguity and completness  
```bash
# Input files: 
# mazia.fna 
# bedadeti.fna

for i in *.fna
do 

BASENAME=$(basename $i | sed -e "s/.fna/g")

# contigs contiguity 
~/apps/gfastats/build/bin/gfastats \
-b s \
-f ${i}.fna  \
-s s \
-t \
--seq-report \
--stats \
--sort descending > ${i}.gfastats.summmary.txt


# BUSCO genes 

busco  \
-i ${i}.fna \
-m genome \
-l embryophyta_odb10 \
-c 16 \
-o busco.${i}.out \
--long \
-f

done 
````
