#!/bin/bash

GENERAL=/Users/COHEN/Documents/BulkRNA/AnalysesBioinformatiques/
GENOME=/Users/COHEN/Documents/BulkRNA/AnalysesBioinformatiques/02_Genome/PlasmoDB-46_Pfalciparum3D7_AnnotatedTranscripts.fasta # rajouter a l'avenir les genes var référencés et enlever les genes vars de 3D7. 
READS=//Users/COHEN/Documents/BulkRNA/AnalysesBioinformatiques/01_RawData
FASTQC=/usr/local/bin/fastqc
TRIMM=/Users/COHEN/Documents/Logiciels/trimmomatic-0.39/trimmomatic-0.39.jar
STAR=/Users/COHEN/Documents/Logiciels/star-2.7.3a/bin/macosx_x86_64/
PATH_DATA_TRIM=/Users/COHEN/Documents/BulkRNA/AnalysesBioinformatiques/04_Trimm


#---------------------------------------------------------------------------------------------------------------------------------------------
#FastQC
cd $GENERAL
echo "########################## FastQC ##########################"
if [ ! -d "03_FastQC" ];then #Creation du fichier
	mkdir -p 03_FastQC #-p créer tous les dossiers dans l'aborescence
	for fq in $READS/*_R1.fastq.gz
      do
      f2=${fq%%_R1.fastq.gz}"_R2.fastq.gz"
      f3=${fq%%_*}
      echo "We do FastQC for:" $f3"."
      fastqc -o FastQC -t 4 $f1 $f2 #-o repertoire output pour l'échantillon -t est le nombre de taches que lui dit de faire en mm temps 1 coeur= 1 tache
     done     
fi


#--------------------------------------------------------------------------------------------------------------------------------------------

#Trimming avec trimmomatic

echo "########################### Trimming des séquences ##########################"
if [ ! -d "04_Trimm" ];then #Creation du fichier
	mkdir "04_Trimm"
cd $READS
	for ft in *_R1.fastq.gz;
      do
      f1=$ft
      f2=${ft/_R1.fastq.gz/}_R2.fastq.gz
      f3=${ft%%_*}

		java -jar /Users/COHEN/Documents/Logiciels/trimmomatic-0.39/trimmomatic-0.39.jar PE -trimlog Resultslog.txt $f1 $f2 $f3"_1_paired.fq.gz" $f3"_1_unpaired.fq.gz" $f3"_2_paired.fq.gz" $f3"_2_unpaired.fq.gz"  ILLUMINACLIP:/Users/COHEN/Documents/Logiciels/trimmomatic-0.39/adapters/TruSeq3-PE-2.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:50
    mv $f3"_1_paired.fq.gz" $f3"_1_unpaired.fq.gz" $f3"_2_paired.fq.gz" $f3"_2_unpaired.fq.gz" /Users/COHEN/Documents/BulkRNA/AnalysesBioinformatiques/04_Trimm
	done 
fi


#--------------------------------------------------------------------------------------------------------------------------------------------

## /!\ DO NOT RUN /!\
#Indexing for STAR
cd $GENERAL
if [ ! -d "05_Index_STAR" ];then
    echo "########################### Indexing for STAR ##########################"
    cd /Users/COHEN/Documents/BulkRNA/AnalysesBioinformatiques/02_Genome/P0060623_1612/
    mkdir "05_Index_STAR"
    cd /Users/COHEN/Documents/BulkRNA/AnalysesBioinformatiques/
    /Users/COHEN/Documents/Logiciels/STAR-2-7-5c/bin/MacOSX_x86_64/STAR \
      --runMode genomeGenerate \
     --runThreadN 4 \
      --genomeDir /Users/COHEN/Documents/BulkRNA/AnalysesBioinformatiques/02_Genome/05_Index_STAR \
      --genomeFastaFiles /Users/COHEN/Documents/BulkRNA/AnalysesBioinformatiques/02_Genome/plasmodb-46_pfalciparum3d7_genome.fasta \
      --sjdbGTFfile /Users/COHEN/Documents/BulkRNA/AnalysesBioinformatiques/02_Genome/plasmodb-46_pfalciparum3d7_gffread.gtf\
      
fi

cd $GENERAL
#Mapping with STAR
if [ ! -d "05-Mapping_STAR" ];then
  echo "########################### Mapping with STAR ##########################"
  mkdir "05-Mapping_STAR"
  cd $PATH_DATA_TRIM

  for f1 in *_1_paired.fq.gz
    do    f2=${f1%%_1_paired.fq.gz}"_2_paired.fq.gz"
      echo $f2
      f3=${f1%_1_paired.fq.gz}
      echo $f3
      echo "We do the mapping for:" $f3"."

      /Users/COHEN/Documents/Logiciels/STAR-2-7-5c/bin/MacOSX_x86_64/STAR --runThreadN 4 \
         --genomeDir  /Users/COHEN/Documents/BulkRNA/AnalysesBioinformatiques/02_Genome/05_Index_STAR \
         --readFilesIn $f1 $f2 \
         --readFilesCommand gunzip -c \
         --outFileNamePrefix /Users/COHEN/Documents/BulkRNA/AnalysesBioinformatiques/05-Mapping_BWA/"/"$f3 \
         --outSAMtype BAM SortedByCoordinate \
         --sjdbGTFfile /Users/COHEN/Documents/BulkRNA/AnalysesBioinformatiques/02_Genome/plasmodb-46_pfalciparum3d7_gffread.gtf \
         --sjdbGTFtagExonParentTranscript Parent \
         --runMode alignReads \
         --outFilterMultimapNmax 1\
         --outFilterMatchNmin 35 \
         --twopassMode Basic
  done
fi

cd 