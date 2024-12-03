# Structural identification and annotation repeatitive elements
```bash
#Installations :
conda create -n repeatmasking -c bioconda RepeatMasker repeatmodeler trf 

# Set path to assembled genome
genome=path_to_assembled_genome

# RepeatProteinMask
RepeatProteinMask $genome -noLowSimple -pvalue 0.0001 
cat ${genome}.annot | awk -v OFS='\t' '{print $4,$5,$6,$7,$8,$9,$10,$11}' | grep -v SeqID > repeatproteinmask.bed
bedtools getfasta -fi ${genome}.fna -bed repeatproteinmask.bed -fo RM.out.fasta  

# RepeatMasker
RepeatMasker -species Viridiplantae $genome -pa 8 

RM2Bed.py ${genome}.fna.out
bedtools getfasta -fi ${genome}.fna -bed ${genome}.rm.bed -fo RM.out.fasta  

# RepeatModeler 
BuildDatabase -name your_species_name ${genome}.fna
RepeatModeler -database your_species_name -pa 32 -LTRStruct  

# Miniature Inverted-repeat Transposable Elements (MITE)-Hunter
~/apps/MITE-hunter/MITE_Hunter_manager.pl -i ${genome}.fna -g ${genome}_MH -c 32 -n 5 -S 12345678 –P 1

## MITE_Hunter_manager TEs:  ${genome}_MH_Step8_singlet.fa 

# Tandem Repeat Finder (TRF)
trf ${genome}.fna 2 5 7 80 10 50 2000 -m -h 

# Convert TRF to fasta:
# source:TRFdat_to_bed.py #source: https://raw.githubusercontent.com/hdashnow/TandemRepeatFinder_scripts/master/TRFdat_to_bed.py

./TRFdat_to_bed.py --dat ${genome}.fna.2.7.7.80.10.50.500.dat --bed ${genome}.trf.bed 

bedtools getfasta -fi ${genome}.fna -bed ${genome}.trf.bed  > ${genome}.trf.fasta

# PILER as described in https://www.drive5.com/piler/. 
pals -self ${genome}.fna hit.gff
piler -trs hit.gff -out trs.gff
mkdir fams
piler -trs2fasta trs.gff -seq ${genome}.fasta -path fams
mkdir aligned_fams
cd fams
for fam in *
do
   muscle -in $fam -out ../aligned_fams/$fam -maxiters 1 -diags1
done

cd ..
mkdir cons
cd aligned_fams

for fam in *
do
   piler -cons $fam -out ../cons/$fam -label $fam
done
cd ../cons
cat * > ../${genome}.piler_library.fasta

# sequence redundancy removal - CD-HIT
# contatenate into one file  
cat ${genome}.RPM.out.fasta ${genome}.RM.fasta ${genome}.RMod.out.fasta ${genome}.EDTA.out.fasta ${genome}.trf.fasta ${genome}_MH_Step8_singlet.fa ${genome}.piler_library.fasta > redundant_repeat_lib.fasta 

cd-hit -i redundant_repeat_lib.fasta -d 0 -o redundant_repeat_lib.non_redun.fasta -c 0.80 -n 5 -G 1 -g 0 -M 0 -T 24

# Excluding highly conserved protein coding genes or multi-member gene families
# Here BLASTP search was used (e-value threshold of 1e-10) against monocots protein database downloaded from EnsemblPlants (https://plants.ensembl.org/index.html) and Banana Genome Hub.

# Extract the fasta files the predicted repeats from the various tools

rep_lib=redundant_repeat_lib.non_redun

makeblastdb -in blastdb/musa.arab_T_ensemblemonocot.enseteg.fa -input_type prot 
### BLASTP
blastx 
-db blastdb/musa.arab_T_ensemblemonocot.enseteg.fa \
-query $rep_lib \
-outfmt '6 qseqid qlen slen qstart qend sstart send stitle pident length evalue bitscore' \
-evalue 1e-10 \
-max_hsps 1 \
-num_threads 24 \
-out $rep_lib.blastx.out

# convert blast hits to bed file and exclude the bed file interval from fasta 
cut -f 1,4,5 $rep_lib.blastx.out | awk '{if ($2>$3) print $1,$3,$2,".",".","-"; else print $1,$2,$3,".",".","+";}' OFS='\t' | awk '{a=$2-1;print $1,a,$3,$4,$5,$6;}' OFS='\t'| bedtools sort > $rep_lib.protein_cleared.draft.bed
bedtools merge -i $rep_lib.protein_cleared.draft.bed -s -c 6 -o distinct > $rep_lib.protein_cleared.final.bed
# generate the index file 
samtools faidx $rep_lib
# reformate the index file to genome file 
awk -v OFS='\t' {'print $1,$2'} $rep_lib.fai > $rep_lib.txt
# final bed file to exclude  <br />
bedtools complement -i $rep_lib.protein_cleared.final.bed -g $rep_lib.txt > $rep_lib.final.bed
# final protein coding gene free final de novo repeats 
bedtools getfasta -fi ${rep_lib}.fasta -bed $rep_lib.final.bed > $rep_file.final

# Extensive de-novo Transposable Element Annotator (EDTA) 
# installation: https://github.com/oushujun/EDTA 

EDTA.pl \
-genome ${genome}.fna \
--species other \
--step all \
--overwrite 0 \
--sensitive 1 \
--anno 1 \
--evaluate 1 \
--threads 48 \
--debug 1

# Final combined repeat identification and repeatitive sequence masking 
cat $rep_file.final ${genome}.fna.mod.EDTA.TElib.fa > ${genome}_repeat_lib.fasta

```

# Annotation of structurally identified repeatitive elements

```bash
RepeatMasker -nolow -no_is -norna -engine ncbi -lib ${genome}_repeat_lib.fasta ${genome}.fna -pa 12 -gff -poly -u -xsmall

RM2Bed.py ${genome}.fna.out
bedtools getfasta -fi ${genome}.fna -bed ${genome}.rm.bed -fo ${genome}.RM.final.out.fasta  

### use the hard masked coordinate bed file for soft masking in bedtools 
# bedtools maskfasta \
# -fi ${genome}.fna \
# -bed ${genome}.rm.bed -fo ${genome}.softmasked.1.fna -soft

# Annotate RepeatMasker output repeats using TEsorter
TEsorter \
${genome}.RM.final.out.fasta \
--seq-type nucl \
--processors 48 \
--pass2-rule 80-80-80 \
--prefix ${genome}_RM_tesorter.out
```
