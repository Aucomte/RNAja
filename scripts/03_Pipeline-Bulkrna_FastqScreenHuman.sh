#!/bin/bash

GENERAL=/Users/COHEN/Documents/BulkRNA/AnalysesBioinformatiques/
GENOME=/Users/COHEN/Documents/BulkRNA/AnalysesBioinformatiques/02_Genome/PlasmoDB-46_Pfalciparum3D7_AnnotatedTranscripts.fasta # rajouter a l'avenir les genes var référencés et enlever les genes vars de 3D7. 
READS=/Users/COHEN/Documents/BulkRNA/AnalysesBioinformatiques/01_RawData/RNA-seq
FASTQC=/usr/local/bin/fastqc
TRIMM=/Users/COHEN/Documents/Logiciels/trimmomatic-0.39/trimmomatic-0.39.jar
STAR=/Users/COHEN/Documents/Logiciels/star-2.7.3a/bin/macosx_x86_64/
PATH_DATA_TRIM=/Users/COHEN/Documents/BulkRNA/AnalysesBioinformatiques/04_Trimm



#Trimming avec trimmomatic

echo "########################### Mapping Human ##########################"



cd $GENERAL

if [ ! -d "06_Humancov" ];then #Creation du fichier
	mkdir "06_Humancov"
cd $READS
	for fm in *_r1.fastq.gz;
      do
      f1=$fm
      f2=${fm/_r1.fastq.gz/}_r2.fastq.gz
      f3=${fm%%_*}

      fastq_screen --conf /Users/COHEN/Documents/Logiciels/FastQ-Screen-0.14.1/fastq_screen.conf --aligner BOWTIE2 $f1 $f2 --outdir /Users/COHEN/Documents/BulkRNA/AnalysesBioinformatiques/06_Humancov/
  done 
fi 
