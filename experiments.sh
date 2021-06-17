#!/bin/bash
wget https://pstools.s3.us-east-2.amazonaws.com/hifiasm
chmod 777 hifiasm
for i in HG02109 HG02080 HG01928 HG01891 HG01243 HG01175 HG01358 HG01109;
do mkdir $i; mkdir analyses; cd $i; mkdir PacBio_HiFi; cd PacBio_HiFi;
aws s3 cp --recursive --no-sign-request s3://human-pangenomics/working/HPRC/${i}/raw_data/PacBio_HiFi/ ./;
for j in *.bam; do samtools bam2fq -@64 $j > $j.fastq; rm $j; done
rm -r *.bam;
~/final_exps/hifiasm -o ../../analyses/${i}.asm -t64 -z40 *.fastq > ${i}.asm.log 2>&1;
rm -r *.fastq;
awk '/^S/{print ">"$2;print $3}' ../../analyses/${i}.asm.r_utg.gfa > ../../analyses/${i}.asm.r_utg.fa;
cd ..;
 mkdir hic; cd hic; 
aws s3 cp --recursive --no-sign-request s3://human-pangenomics/working/HPRC/${i}/raw_data/hic/ ./;
wget https://pstools.s3.us-east-2.amazonaws.com/pstools_1;
mv pstools_1 pstools; 
chmod 777 pstools;
./pstools hic_mapping_unitig -t64 ../../analyses/${i}.asm.r_utg.fa <(zcat *R1*.fastq.gz) <(zcat *R2*.fastq.gz);
./pstools resolve_haplotypes -t64 hic_name_connection.output ../../analyses/${i}.asm.r_utg.gfa ./; 
./pstools hic_mapping_haplo -t64 pred_haplotypes.fa <(zcat *R1*.fastq.gz) <(zcat *R2*.fastq.gz) -o scaff_connections.txt;
./pstools haplotype_scaffold -t64 scaff_connections.txt pred_haplotypes.fa ./ 
 cp pred_hap1.fa ../../analyses/${i}.hap1.fa;
 cat pred_hap2.fa pred_broken_nodes.fa > ../../analyses/${i}.hap2.fa;
 cd ..  
 rm -rf hic; 
 rm -rf PacBio_HiFi; 
 cd .. 
 done 
