#!/bin/sh                                                                                                                                                                                                 

###################################################################                                                                                                                                       
# pipeline for running snp-calling pipeline on re-sequencing data #                                                                                                                                       
#                                                                 #                                                                                                                                       
# usage:  snp_calling_April_13_2018.sh  $basedir $num_of_processors     #                                                                                                                                       
#                                                                 #                                                                                                                                       
###################################################################                                                                                                                                       
cd $1
CPU=$2

#hisat2-build ~/eric/TAIR10_chr_all.fas ~/eric/TAIR10_chr_all
#bwa/bowtie2

bwa index ~/path-to-cowpea-genome/genome.fasta #this should be edited to the path of the cowpea genome file

for file in `dir -d *_1_cln.fq.gz ` ; do #change the regular expression (regex) to something that matches your trimmed data                                                                                                                                                               

    file2=`echo "$file" |sed 's/_1_cln.fq.gz/_2_cln.fq.gz/'` #might need to change the regex to matche your trimmed read data
    samfile=`echo "$file" | sed 's/_1_cln.fq.gz/.sam/'`                                                                                                                                                       

    bwa mem -t $CPU genome $file $file2 #change "genome" to the prefix of the reference cow pea genome {no .fasta}
done                                                                                                                                                                                                     

ls *.sam | parallel -j $CPU samtools view -bS -o {.}.bam {}                                                                                                                                              
ls *.bam | parallel -j $CPU samtools sort -o {.}.sort.bam {}

#all of these steps are derived from the Broad GATK best practices: https://gatk.broadinstitute.org/hc/en-us/articles/360035894711-About-the-GATK-Best-Practices
#install Picard: https://broadinstitute.github.io/picard/
ls *.sort.bam | parallel -j $CPU java -Djava.io.tmpdir=`pwd`~/tmp -jar ~/snp_calling/picard-tools-1.119/MarkDuplicates.jar INPUT={} OUTPUT={.}.md.bam METRICS_FILE={.}.metrics REMOVE_DUPLICATES=false ASSUME_SORTED=true VALIDATION_STRINGENCY=SILENT                                                                                                                                                               
ls *.md.bam | parallel -j $CPU java -Djava.io.tmpdir=`pwd`~/tmp -jar ~/snp_calling/picard-tools-1.119/AddOrReplaceReadGroups.jar INPUT={} OUTPUT={}.rg.bam SORT_ORDER=coordinate RGID={} RGLB=1 RGPL=illumina RGPU=run RGSM={} RGCN=tom360 RGDS={}                                                                                                                                                                   
ls *.rg.bam | parallel -j $CPU samtools index {}

#Indel realignment - install GATK https://gatk.broadinstitute.org/hc/en-us
ls *.rg.bam | parallel -j $CPU java -Djava.io.tmpdir=`pwd`~/tmp -jar ~/snp_calling/GenomeAnalysisTK.jar -T RealignerTargetCreator -R ~/eric/TAIR10_chr_all.fas -U ALLOW_N_CIGAR_READS -nt 1 -I {}  -o {.}.intervals                                                                                                                             
ls *rg.bam | parallel -j $CPU java -Djava.io.tmpdir=`pwd`~/tmp -jar ~/snp_calling/GenomeAnalysisTK.jar -T IndelRealigner -R ~/eric/TAIR10_chr_all.fas -U ALLOW_N_CIGAR_READS -I {}  -o {.}.ralgn.bam -targetIntervals {.}.intervals

#calling SNPs
ls *.ralgn.bam |parallel -j $CPU java -jar ~/snp_calling/GenomeAnalysisTK.jar -T HaplotypeCaller -R ~/eric/TAIR10_chr_all.fas -U ALLOW_N_CIGAR_READS -I {} -o {.}.g.vcf --emitRefConfidence GVCF #substitute ~/eric/TAIR10_chr_all.fas with the cow pea ref genome

#CombineGVCFs step
ls *.vcf > gvcf.list
java -Djava.io.tmpdir=`pwd`~/tmp -jar ~/snp_calling/GenomeAnalysisTK.jar -T CombineGVCFs -R ~/eric/TAIR10_chr_all.fas -U ALLOW_N_CIGAR_READS --variant gvcf.list -o gvcf.vcf


