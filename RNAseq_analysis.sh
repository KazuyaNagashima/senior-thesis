#####set data#########
filename=("ERR3153937")
ref_gene=("susScr11")
######################

###PAIRE END###

dir="/home/ohganelab/Desktop/nagashima"
cd ${dir}
mkdir ${dir}/${filename}

#move file to a new directory

mv ${filename}.sra ${dir}/${filename}
cd ${dir}/${filename}

#use SRA Toolkit

fastq-dump --split-files ${filename}.sra --gzip

#quality check by fastp

fastp -i ${filename}_1.fastq.gz -I ${filename}_2.fastq.gz -o fastp_${filename}_1.fastq.gz -O fastp_${filename}_2.fastq.gz -h report_fastp.html -j report_fastp.json -q 25 -n 10 -w 14 --length_required 20 --length_limit 500

#extract files

gunzip fastp_${filename}_1.fastq.gz
gunzip fastp_${filename}_2.fastq.gz

#mapping by HISAT2

hisat2 -p 14 -x /home/ohganelab/analyzetools/reference/Sus_scrofa/Ensembl/Sscrofa11.1/Sequence/Hisat2Index/susScr11 -1 fastp_${filename}_1.fastq -2 fastp_${filename}_2.fastq -S ${filename}.sam

#quality check by samtools (MAPQ>4)

samtools view -q 4 ${filename}.sam > ${filename}.qc.sam

#from sam to bam

cat ${filename}.sam | perl -e 'while(<>){print $_ if(/\@SQ||\@PG/);}' >${filename}.header.sam
cat ${filename}.header.sam ${filename}.qc.sam > ${filename}.qc.hd.sam

samtools view -b -S ${filename}.qc.hd.sam > ${filename}.bam
samtools sort ${filename}.bam > ${filename}_sorted.bam
samtools index ${filename}_sorted.bam

#assemble by stringtie
stringtie ${filename}_sorted.bam -p 14 -G /home/ohganelab/analyzetools/reference/Sus_scrofa/Ensembl/Sscrofa11.1/Annotation/Genes/genes.gff3 -o ${filename}.gtf -l ${filename} -A ${filename}.abund.tab

#organize directory
mkdir stringtie
mv ${filename}.gtf ${filename}.abund.tab stringtie
mkdir fastp_report
mkdir trush
mv report_fastp.html report_fastp.json fastp_report
mv ${filename}.sam ${filename}.qc.sam ${filename}.qc.hd.sam ${filename}.header.sam ${filename}.bam trush
