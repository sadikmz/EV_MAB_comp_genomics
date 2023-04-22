#!/bin/bash

for i in mazia bedadeti
do 
~/apps/gfastats/build/bin/gfastats \
-b s \
-f EV_${i}.assembly.fna  \
-s s \
-t \
--seq-report \
--stats \
--sort descending > ${i}.gfastats.summmary.txt

done 

~/apps/gfastats/build/bin/gfastats \
-f purged.fa  \
-t \
--stats \
--cmd > summary_assembly_stat/assembly_${i}.hap1.stat 


~/apps/gfastats/build/bin/gfastats \
-b s \
-f ~/analysis/annotation/mazia/mazia_v1.500.softmasked.1.fna  \
-s s \
-t \
--seq-report \
--stats 

#!/bin/bash
for i in mazia bedadeti
do 
~/apps/gfastats/build/bin/gfastats \
-b s \
-f /home/u1866313/Novogene/annotation/mazia/qc/genomescope2/${i}.unmasked.1.fna  \
-s s \
-t \
--seq-report \
--stats 
done 


~/apps/gfastats/build/bin/gfastats \
-f /home/u1866313/mnt/Shared258/Enset/pacbio/raw.data/epo.fasta/EV_wild_epo.pacbio_clr.fasta \
-t > epo.pacbio_clr.reads.stat.txt 

~/apps/gfastats/build/bin/gfastats \
-f ../annotation/decontamination/epo.flye.v2.scaffolds.fna \
-t > epo.flye_v2.stat.tx

~/apps/gfastats/build/bin/gfastats \
-f EV_wild_c.adapter_seq_discarded.fasta \
-t > wild_c.reads.stat.txt 

~/apps/gfastats/build/bin/gfastats \
-f wild_b.adapt.discarded.fasta \
-t > wil_d.reads.stat.txt 

~/apps/gfastats/build/bin/gfastats \
-f mazia_hifi.adapt.discarded.fasta \
-t > mazia.reads.stat.txt 

/apps/gfastats/build/bin/gfastats -f assembled_fasta/mazia_s33_k61_adapt_discarded.fasta --nstar-report assembly_nx.stat 

        ~/apps/gfastats/build/bin/gfastats \
        -f assembled_fasta/${i}.hap2.fasta  \
        -t \
        --nstar--report \
        --cmd > summary_assembly_stat/assembly_${i}.hap2.nxstat.stat


