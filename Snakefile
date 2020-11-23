# Importation des librairies glob, sys, et re
import glob,sys,re
from pathlib import Path

## appel du fichier de config
configfile: "configPT.yaml"

## Definition des variables "constantes" à partir du fichier config (attention aux arborescences dans le yaml)
FASTA_REF= config['FASTA']
GFF_REF= config['GFF']
GTF_REF= config['GTF']
READS_DIR = config['READS']
OUTDIR = Path(config["OUTDIR"]).resolve().as_posix()

FASTQC_DIR = f"{OUTDIR}/01_FastQC/"
FASTQSCREEN_DIR = f"{OUTDIR}/02_FastQScreen/"
FASTP_DIR = f"{OUTDIR}/03_FastP/"
FASTQC_AF_DIR = f"{OUTDIR}/04_FastQC_afterFT/"
MULTIQC_DIR = f"{OUTDIR}/05_MULTIQC/"
INDEX_DIR = f"{OUTDIR}/06_Index_HISAT/"
INDEX_STAR_DIR = f"{OUTDIR}/07_index_HISAT/"
MAPPINGH_DIR = f"{OUTDIR}/08_mapping_HISAT/"
MAPPINGS_DIR = f"{OUTDIR}/09_mapping_STAR/"
BAM_DIR = f"{OUTDIR}/12_BAMfiles/"
STATS_DIR = f"{OUTDIR}/14_Stats_BAM/"
FILTER_DIR = f"{OUTDIR}/13_Filtered_BAM/"
DUPL_DIR = f"{OUTDIR}/15_MarkedDupli_BAM/"
HTSEQ_DIR = f"{OUTDIR}/16_HTseq-count/"
STRINGTIE_DIR = f"{OUTDIR}/17_Stringtie/"
name_hisat_index = config['name_hisat_index']

## Software path

fastQScreen_conf = config['FastQScreen_conf']

#Wild cards
READS, = glob_wildcards(f"{READS_DIR}{{reads}}.fastq.gz")
MAPPING = ["STAR", "HISAT2"]
COUNT = ["HTSEQ", "STINGTIE"]

# Regle finale pour vérifier la présence des outputs et si ils sont présents ils ne lancent pas le job
rule final:
     input:
          out_fastqc = expand(f"{FASTQC_DIR}{{reads}}_fastqc.html", reads = READS),
          out_fastqs = expand(f"{FASTQSCREEN_DIR}{{reads}}_screen.html", reads = READS),
          out_fastp = expand(f"{FASTP_DIR}{{reads}}_report_fastp.html", reads = READS),
          #out_fastqc_af = expand(f"{FASTQC_AF_DIR}{{reads}}_r1_paired_fastqc.html", reads = READS),
          out_multiqc = expand(f"{MULTIQC_DIR}multiqc_report.html"),
          out_hisat_map = expand(f"{MAPPINGH_DIR}{{reads}}_HISAT.bam" , reads = READS),
          out_star_map = expand(f"{MAPPINGS_DIR}{{reads}}Aligned.sortedByCoord.out.bam" , reads = READS),
          out_bam_sort = expand(f"{BAM_DIR}{{reads}}_STAR_sort.bam", reads = READS),
          out_bam_filt = expand(f"{FILTER_DIR}{{reads}}_STAR_sort_mapped.bam", reads = READS),
          out_bam_index = expand(f"{BAM_DIR}{{reads}}_STAR_sort.bam.bai", reads = READS),
          out_bam_stats = expand(f"{STATS_DIR}{{reads}}_STAR_sort_flagstat.txt", reads = READS),
          out_bam_md = expand(f"{DUPL_DIR}{{reads}}_HISAT_sort_mapped_markdup_metric.txt", reads = READS),
          out_htseq = expand(f"{HTSEQ_DIR}HISAT/{{reads}}.txt", reads = READS),
          #out_gtf_list = expand(f'{STRINGTIE_DIR}STRINGTIE_MERGE/STAR_gtf_list.txt'),
          #out_stringtie = expand(f'{STRINGTIE_DIR}STRINGTIE_MERGE/STAR_stringtie_merged.gtf'),
          out_stringtiecc = expand(f'{STRINGTIE_DIR}STAR/{{reads}}/{{reads}}.gtf', reads=READS),
          #out_gtf_list = expand(f'{STRINGTIE_DIR}HISAT_gtf_list.txt'),
          out_strg_list = expand(f'{STRINGTIE_DIR}HISAT_Stringtie_list.txt'),
          out_csv = expand(f'{STRINGTIE_DIR}HISAT_transcript_count_matrix.csv')


# ----------------------------------- QUALITY -------------------------------------------------------------------
# FastQC before filtering
rule fastqc:
     """ 
          FastQC - Control Quality of reads
     """
     threads:4
     params:
        out_fastqc=directory(f"{FASTQC_DIR}"),
     input:
        r1 = f"{READS_DIR}{{reads}}.fastq.gz"
     output:
        out_fastqc_r1 = f"{FASTQC_DIR}{{reads}}_fastqc.html",
        out_fastqc_r1_zip = f"{FASTQC_DIR}{{reads}}_fastqc.zip",
     message:
          """
               input=
                    r1 = {input.r1}
               threads={threads}
          """
     conda:
        "envs/quality.yaml"
     shell:
          """
               fastqc -o {params.out_fastqc} -t 4 {input.r1}
          """

# FastQScreen agaisnt 3D7 and human
rule fastqscreen:
     """ 
          FastQScreen - Control Quality of reads
     """
     threads:4
     params:
        out_fastqs=directory(f"{FASTQSCREEN_DIR}"),
        fastqs_conf = f"{fastQScreen_conf}"
     input:
        fastqc= rules.fastqc.output.out_fastqc_r1_zip,
        r1 = f"{READS_DIR}{{reads}}.fastq.gz",
     output:
        html_r1= f"{FASTQSCREEN_DIR}{{reads}}_screen.html",
        txt_r1=f"{FASTQSCREEN_DIR}{{reads}}_screen.txt",
     message:
          """
               Execute {rule}
               input=
                    r1 = {input.r1}
               threads={threads}
          """
     conda:
        "envs/quality.yaml"
     shell:
          """
               fastq_screen --conf {params.fastqs_conf} --aligner BOWTIE2 {input.r1} --outdir {params.out_fastqs}
          """

# Filtering and trimming with FastP
rule fastp:
     """ 
          FastP - Filtering and trimming before mapping
     """
     params:
        out_fastp=directory(f"{FASTP_DIR}"),
     input:
        fastqc = rules.fastqc.output.out_fastqc_r1_zip,
        r1 = f"{READS_DIR}{{reads}}.fastq.gz",
     output:
        out_r1 = f"{FASTP_DIR}{{reads}}_fastp.fastq.gz",
        out_report= f"{FASTP_DIR}{{reads}}_report_fastp.html",
        out_json= f"{FASTP_DIR}{{reads}}_fastp.json"

     message:
        """
            Execute {rule}
               input=
                    r1 = {input.r1}
        """
     conda:
        "envs/quality.yaml"
     shell:
          """
               fastp -i {input.r1} -o {output.out_r1} -h {output.out_report} -j {output.out_json}  
          """

# FastQC after filtering
rule fastqc_af:
     """ 
          FastQC - Control Quality of reads after FastP
     """
     params:
        out_fastqc_af = directory(f"{FASTQC_AF_DIR}"),
     input:
        r1 = rules.fastp.output.out_r1,
     output:
        out_fastqc_af_r1 = f"{FASTQC_AF_DIR}{{reads}}_paired_fastqc.html",
        out_fastqc_af_r1_zip = f"{FASTQC_AF_DIR}{{reads}}_paired_fastqc.zip",
     message:
          """
            Execute {rule}
               input=
                    r1 = {input.r1}
          """
     conda:
        "envs/quality.yaml"
     shell:
          """
               fastqc -o {params.out_fastqc_af} -t 4 {input.r1}
          """

# FastQScreen agaisnt 3D7 and human
rule multi_qc:
    """
    MultiQC
    """
    params:
        out_multi = directory(f"{MULTIQC_DIR}"),
    input:
        expand(f"{FASTQSCREEN_DIR}{{reads}}_screen.html", reads = READS),
        expand(f"{FASTQSCREEN_DIR}{{reads}}_screen.txt", reads = READS),
        expand(f"{FASTQC_DIR}{{reads}}_fastqc.html", reads=READS),
        expand(f"{FASTQC_DIR}{{reads}}_fastqc.zip", reads=READS),
        expand(f"{FASTQC_AF_DIR}{{reads}}_paired_fastqc.html", reads=READS),
        expand(f"{FASTQC_AF_DIR}{{reads}}_paired_fastqc.zip", reads=READS),
        expand(f"{FASTP_DIR}{{reads}}_report_fastp.html", reads=READS),
        expand(f"{FASTP_DIR}{{reads}}_fastp.json",reads=READS)
    output:
        f"{MULTIQC_DIR}multiqc_report.html",
    conda:
        "envs/quality.yaml"
    shell:
        """
             multiQC {input} -o {params.out_multi}
        """

# ----------------------------------- MAPPING -------------------------------------------------------------------
# Indexing of 3D7 genome (environ 1 min)
rule hisat2_index:
    """                                                                                                                                                                                       
    Make index with HISAT2 for 3D7                                                                                                                                                 
    """
    input:
        fasta = f"{FASTA_REF}"
    output:
        index = f"{INDEX_DIR}{name_hisat_index}",
    message:
        """                                                                                                                                                                                   
        Execute {rule}                                                                                                                                                                    
        input:                                                                                                                                                                            
            fasta : {input.fasta}                                                                                                                                                       
        output:                                                                                                                                                                           
            index: {output.index}                                                                                                                     
        """
    log:
        output = f"{INDEX_DIR}LOG/{name_hisat_index}_HISAT-INDEX.o",
        error = f"{INDEX_DIR}LOG/{name_hisat_index}_HISAT-INDEX.e",
    conda:
        "envs/hisat2.yaml"
    shell:
        """
        ln -s {input.fasta} {output.index} 1>{log.output} 2>{log.error}
        hisat2-build {input.fasta} {output.index} --quiet
        """

# Mapping with HISAT2 on Pf 3D7 (max 2 min par échantillon, total 10 min)
rule hisat2_map:
    """
    Map reads on the genome with HISAT2 on paired end
    """
    threads: 4
    input:
        reference = rules.hisat2_index.output.index,
        r1 = rules.fastp.output.out_r1,
    output:
        summary = f"{MAPPINGH_DIR}{{reads}}_HISAT_summary.txt",
        bam = f"{MAPPINGH_DIR}{{reads}}_HISAT.bam",
    message:
        """                                                                                                                                                                            
        Execute {rule}                                                                                                                                                                   
        input:                                                                                                                                                                            
            reference : {input.reference}
            R1 : {input.r1}                                                                                                                                                   
        output:                                                                                                                                                                           
            bam: {output.bam}                                                                                                 
        """
    conda:
        "envs/hisat2.yaml"
    shell:
        """
        hisat2 -p {threads} -x {input.reference} -u {input.r1} --summary-file {output.summary} | samtools view -S -b >{output.bam} 
        """


#Indexing for RNA-STAR (environ 2 min)
rule star_index:
    """
    Make index with STAR
    """
    threads: 4
    params:
        dir = directory(f"{INDEX_STAR_DIR}")
    input:
        fasta = f"{FASTA_REF}",
        gtf = f"{GTF_REF}",
    output:
        index = f"{INDEX_STAR_DIR}Log.out",
    message:
        """
        Execute {rule}
        input:
            fasta : {input.fasta}
            gtf : {input.gtf}
        output:
            index: {output.index}
        """
    conda:
        "envs/Star.yaml"
    shell:
        """
        ln -s {input.fasta} {output.index}
        STAR --runMode genomeGenerate --runThreadN {threads} --genomeDir {params.dir} --genomeFastaFiles {input.fasta} --sjdbGTFfile {input.gtf}
        """

#Mapping with RNA-STAR on Pf 3D7 (environ 10 min par échantillon, total 1h)
rule star_mapping:
    """
    Map reads on the genome with STAR on paired end
    """
    threads: 4
    params:
        idx = directory(f"{INDEX_STAR_DIR}"),
        bam = f"{MAPPINGS_DIR}{{reads}}"
    input:
        idx = rules.star_index.output.index,
        gtf = f"{GTF_REF}",
        r1 = rules.fastp.output.out_r1
    output:
        bam = f"{MAPPINGS_DIR}{{reads}}Aligned.sortedByCoord.out.bam"
    message:
        """
        Execute {rule}
        input:
            idx : {params.idx}
            gtf : {input.gtf}
            reads : 
                r1: {input.r1}
        output:
            bam: {output.bam}
        """
    conda:
        "envs/Star.yaml"
    shell:
        """
        STAR --runThreadN 4 --genomeDir {params.idx} --readFilesIn {input.r1} --readFilesCommand gunzip -c --outFileNamePrefix {params.bam} --outSAMtype BAM SortedByCoordinate --sjdbGTFfile {input.gtf} 
        """

#Sort the BAM files with samtools
rule bam_sort:
    """
    Sort the BAM files     
    """
    threads: 4
    input:
        star =f'{MAPPINGS_DIR}{{reads}}Aligned.sortedByCoord.out.bam',
        hisat = f"{MAPPINGH_DIR}{{reads}}_HISAT.bam"
    output:
        out_star = f"{BAM_DIR}{{reads}}_STAR_sort.bam",
        out_hisat = f"{BAM_DIR}{{reads}}_HISAT_sort.bam"
    message:
        """
        Execute {rule}
        input:
            star : {input.star}
            hisat : {input.hisat}
        output:
            bam: {output.out_star}, {output.out_hisat}
        """
    conda:
        "envs/samtools.yaml"
    shell:
        """
        samtools sort  {input.star} > {output.out_star};
        samtools sort  {input.hisat} > {output.out_hisat}
        """

#Filtering of the BAM files with samtools
rule bam_filter:
    """
    Filtering BAM files      
    """
    threads: 4
    input:
        star_bam = f"{BAM_DIR}{{reads}}_STAR_sort.bam",
        hisat_bam = f"{BAM_DIR}{{reads}}_HISAT_sort.bam"
    output:
        mapped_star = f"{FILTER_DIR}{{reads}}_STAR_sort_mapped.bam",
        unmapped_star = f"{FILTER_DIR}{{reads}}_STAR_sort_unmapped.bam",
        mapped_hisat = f"{FILTER_DIR}{{reads}}_HISAT_sort_mapped.bam",
        unmapped_hisat = f"{FILTER_DIR}{{reads}}_HISAT_sort_unmapped.bam"
    message:
        """
        Execute {rule}
            input: {input}
        
        """
    conda:
        "envs/samtools.yaml"
    shell:
        """
        samtools view -b -F 4 {input.star_bam} > {output.mapped_star}
        samtools view -b -f 4 {input.star_bam} > {output.unmapped_star}
        samtools view -b -F 4 {input.hisat_bam} > {output.mapped_hisat}
        samtools view -b -f 4 {input.hisat_bam} > {output.unmapped_hisat}
        """

# Identify the duplicates in BAM with Picard tools
rule bam_markdup:
    """
    Mark duplicates in BAM files      
    """
    threads: 4
    input:
         bam_star=f"{FILTER_DIR}{{reads}}_STAR_sort_mapped.bam",
         bam_hisat = f"{FILTER_DIR}{{reads}}_HISAT_sort_mapped.bam",
    output:
        star = f"{DUPL_DIR}{{reads}}_STAR_sort_mapped_markdup.bam",
        star_txt = f"{DUPL_DIR}{{reads}}_STAR_sort_mapped_markdup_metric.txt",
        hisat = f"{DUPL_DIR}{{reads}}_HISAT_sort_mapped_markdup.bam",
        hisat_txt = f"{DUPL_DIR}{{reads}}_HISAT_sort_mapped_markdup_metric.txt"
    message:
        """
        Execute {rule}
        input: {input}
        """
    conda:
        "envs/picard.yaml"
    shell:
        """
        java -jar picard.jar MarkDuplicates I={input.bam_star} O={output.star} M={output.star_txt}
        java -jar picard.jar MarkDuplicates I={input.bam_hisat} O={output.hisat} M={output.hisat_txt}
        """

#Index the BAM files with samtools
rule bam_index:
    """
    Index BAM files      
    """
    threads: 4
    input:
        bam_star = f"{BAM_DIR}{{reads}}_STAR_sort.bam",
        bam_hisat = f"{BAM_DIR}{{reads}}_HISAT_sort.bam",
        bam_star_filt = f"{FILTER_DIR}{{reads}}_STAR_sort_mapped.bam",
        bam_hisat_filt =f"{FILTER_DIR}{{reads}}_HISAT_sort_mapped.bam",
        bam_star_filt_md = f"{DUPL_DIR}{{reads}}_STAR_sort_mapped_markdup.bam",
        bam_hisat_filt_md = f"{DUPL_DIR}{{reads}}_HISAT_sort_mapped_markdup.bam",
    output:
        bai = f"{BAM_DIR}{{reads}}_STAR_sort.bam.bai"
    message:
        """
        Execute {rule}
            input: {input}
        """
    conda:
        "envs/samtools.yaml"
    shell:
        """
        samtools index {input.bam_star}
        samtools index {input.bam_hisat}
        samtools index {input.bam_star_filt}
        samtools index {input.bam_hisat_filt}
        samtools index {input.bam_star_filt_md}
        samtools index {input.bam_hisat_filt_md}
        """

#Stats on the BAM files with samtools
rule bam_stats:
    """
    Statistics data on the BAM files      
    """
    threads: 4
    input:
        bam_star = f"{BAM_DIR}{{reads}}_STAR_sort.bam",
        bam_hisat = f"{BAM_DIR}{{reads}}_HISAT_sort.bam",
        bam_star_filt = f"{FILTER_DIR}{{reads}}_STAR_sort_mapped.bam",
        bam_hisat_filt =f"{FILTER_DIR}{{reads}}_HISAT_sort_mapped.bam",
        bam_star_filt_md = f"{DUPL_DIR}{{reads}}_STAR_sort_mapped_markdup.bam",
        bam_hisat_filt_md = f"{DUPL_DIR}{{reads}}_HISAT_sort_mapped_markdup.bam",
    output:
        flagstat_star = f"{STATS_DIR}{{reads}}_STAR_sort_flagstat.txt",
        flagstat_hisat = f"{STATS_DIR}{{reads}}_HISAT_sort_flagstat.txt",
        flagstat_star_filt= f"{STATS_DIR}{{reads}}_STAR_sort_mapped_flagstat.txt",
        flagstat_hisat_filt = f"{STATS_DIR}{{reads}}_HISAT_sort_mapped_flagstat.txt",
        flagstat_star_filt_md= f"{STATS_DIR}{{reads}}_STAR_sort_mapped_markdup_flagstat.txt",
        flagstat_hisat_filt_md= f"{STATS_DIR}{{reads}}_HISAT_sort_mapped_markdup_flagstat.txt",
    message:
        """
        Execute {rule}
            input: {input}
        """
    conda:
        "envs/samtools.yaml"
    shell:
        """
        samtools flagstat {input.bam_star} > {output.flagstat_star}
        samtools flagstat {input.bam_hisat} > {output.flagstat_hisat}
        samtools flagstat {input.bam_star_filt} > {output.flagstat_star_filt}
        samtools flagstat {input.bam_star_filt_md} > {output.flagstat_star_filt_md}
        samtools flagstat {input.bam_hisat_filt} > {output.flagstat_hisat_filt}
        samtools flagstat {input.bam_hisat_filt_md} > {output.flagstat_hisat_filt_md}
        """

# ----------------------------------- COUNT -------------------------------------------------------------------
#Create the count table from the BAM files with HTseq count
rule htseq_counts :
    """
    Create the count table with HTseq-count      
    """
    threads: 4
    input:
         bam_index = rules.bam_index.output.bai,
         bam_star = f"{DUPL_DIR}{{reads}}_STAR_sort_mapped_markdup.bam",
         bam_hisat = f"{DUPL_DIR}{{reads}}_HISAT_sort_mapped_markdup.bam",
         gtf = f"{GTF_REF}"
    output:
         count_table_star = f"{HTSEQ_DIR}STAR/{{reads}}.txt",
         count_table_hisat = f"{HTSEQ_DIR}HISAT/{{reads}}.txt"
    message:
        """
        Execute {rule}
        input: 
            bam: {input.bam_star}
            gtf: {input.gtf}
        output: {output.count_table_star}
        """
    conda:
        "envs/htseq.yaml"
    shell:
        """
        htseq-count -t exon -i gene_id -f bam {input.bam_star} {input.gtf} > {output.count_table_star}
        htseq-count -t exon -i gene_id -f bam {input.bam_hisat} {input.gtf} > {output.count_table_hisat}
        """

#Create the count table from the BAM files with stringtie
rule stringtie :
    """
    Create the count table with HTseq-count
    """
    threads: 4
    input:
         bam_index = rules.bam_index.output.bai,
         bam_star = f"{DUPL_DIR}{{reads}}_STAR_sort_mapped_markdup.bam",
         bam_hisat = f"{DUPL_DIR}{{reads}}_HISAT_sort_mapped_markdup.bam",
         gtf = f"{GTF_REF}"
    output:
         star_gtf = f'{STRINGTIE_DIR}STAR/{{reads}}/{{reads}}.gtf',
         hisat_gtf= f'{STRINGTIE_DIR}HISAT/{{reads}}/{{reads}}.gtf',
         star_tsv = f'{STRINGTIE_DIR}STAR/{{reads}}/{{reads}}.tsv',
         hisat_tsv = f'{STRINGTIE_DIR}HISAT/{{reads}}/{{reads}}.tsv'
    message:
        """
        Execute {rule}
        input:
            bam: {input.bam_star}, {input.bam_hisat}
            gtf: {input.gtf}
        output: {output.star_tsv}
        """
    conda:
        "envs/stringtie.yaml"
    shell:
        """
        stringtie -p 8 -G {input.gtf} -e -B -o {output.star_gtf} -A {output.star_tsv} {input.bam_star}
        stringtie -p 8 -G {input.gtf} -e -B -o {output.hisat_gtf} -A {output.hisat_tsv} {input.bam_hisat}
        """

# Create a list of GTF created by Stringtie
rule stringtie_gtf_list:
    threads : 1
    input:
        list_gtf_star = expand(f'{STRINGTIE_DIR}STAR/{{reads}}/{{reads}}.gtf', reads = READS),
        list_gtf_hisat = expand(f'{STRINGTIE_DIR}HISAT/{{reads}}/{{reads}}.gtf', reads = READS),
    output :
        gtf_star = f'{STRINGTIE_DIR}STAR_gtf_list.txt',
        gtf_hisat = f'{STRINGTIE_DIR}HISAT_gtf_list.txt',
    shell:
        """
        find {input.list_gtf_star} >> {output.gtf_star}
        find {input.list_gtf_hisat} >> {output.gtf_hisat}
        """

# Create a list with the ID sample and the path of the stringtie GTF
rule list_for_prepDE:
    input: rules.stringtie_gtf_list.output.gtf_star, rules.stringtie_gtf_list.output.gtf_star
    output:
        star_list= f'{STRINGTIE_DIR}STAR_Stringtie_list.txt',
        hisat_list= f'{STRINGTIE_DIR}HISAT_Stringtie_list.txt'
    run:
        star_list_name_reads= open(output.star_list, "w")
        hisat_list_name_reads= open(output.hisat_list, "w")
        star_gtf = open(rules.stringtie_gtf_list.output.gtf_star, "r")
        hisat_gtf = open(rules.stringtie_gtf_list.output.gtf_hisat, "r")
        for i in range(len(READS)):
            name= READS[i]
            x=re.split("-", name)
            ID = x[0]+"_"+x[2]
            star_list_name_reads.write(ID + " " + star_gtf.readline())
            hisat_list_name_reads.write(ID + " " + hisat_gtf.readline())

# Convert the stringtie GTF to a count table
rule prepDE_stringtie_table:
    threads : 1
    input :
        mergelist_star = rules.list_for_prepDE.output.star_list,
        mergelist_hisat = rules.list_for_prepDE.output.hisat_list,
    output :
        gcsv_star = f'{STRINGTIE_DIR}STAR_gene_count_matrix.csv',
        tcsv_star = f'{STRINGTIE_DIR}STAR_transcript_count_matrix.csv',
        gcsv_hisat = f'{STRINGTIE_DIR}HISAT_gene_count_matrix.csv',
        tcsv_hisat = f'{STRINGTIE_DIR}HISAT_transcript_count_matrix.csv',
    conda:
        "envs/stringtie.yaml"
    shell:
        """
        python2 prepDE.py -i {input.mergelist_star} -t {output.tcsv_star} -g {output.gcsv_star}
        python2 prepDE.py -i {input.mergelist_hisat} -t {output.tcsv_hisat} -g {output.gcsv_hisat}
       """

# pseudo count:

rule Kalisto:




