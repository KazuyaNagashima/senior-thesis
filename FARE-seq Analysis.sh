#####set data#########
filename=("fibro-3min-170728")
ref_gene=("susScr11")
ngs_ref=("sscr11")
######################

dir="/home/ohganelab/Desktop/nagashima"
cd ${dir}
mkdir ${dir}/${filename[k]}

#move file to a new directory

mv ${filename[k]}_S1_L001_R1_001.fastq.gz ${dir}/${filename[k]}/${filename[k]}_S1_L001_R1_001.fastq.gz
mv ${filename[k]}_S1_L001_R2_001.fastq.gz ${dir}/${filename[k]}/${filename[k]}_S1_L001_R2_001.fastq.gz

cd ${dir}/${filename[k]}

#quality check with fastp

fastp -i ${filename[k]}_S1_L001_R1_001.fastq.gz -I ${filename[k]}_S1_L001_R2_001.fastq.gz -o fastp_${filename[k]}_S1_L001_R1_001.fastq.gz -O fastp_${filename[k]}_S1_L001_R2_001.fastq.gz -h report_fastp.html -j report_fastp.json -q 25 -n 10 -w 14 --length_required 20 --length_limit 500

#extract files

gunzip fastp_${filename[k]}_S1_L001_R1_001.fastq.gz
gunzip fastp_${filename[k]}_S1_L001_R2_001.fastq.gz

#mapping with bowtie2

bowtie2 -p 18 -x /home/ohganelab/analyzetools/reference/Sus_scrofa/Ensembl/Sscrofa11.1/Sequence/Bowtie2Index/genome -1 fastp_${filename[k]}_S1_L001_R1_001.fastq -2 fastp_${filename[k]}_S1_L001_R2_001.fastq -S ${filename[k]}.sam

#remove multi-mapped reads

cat ${filename[k]}.sam | perl -e 'while(<>){print $_ if(/\@SQ||\@PG/);}' >${filename[k]}.header.sam
grep -v "XS" ${filename[k]}.sam >${filename[k]}.uq.sam

#quality check with samtools (MAPQ>4)

samtools view -q 4 ${filename[k]}.uq.sam > ${filename[k]}.uq.qc.sam

#from sam to bam,bed,bedgraph

cat ${filename[k]}.header.sam ${filename[k]}.uq.qc.sam >${filename[k]}.uq.qc.hd.sam
samtools view -b -S ${filename[k]}.uq.qc.hd.sam > ${filename[k]}.bam
samtools sort ${filename[k]}.bam > ${filename[k]}_sorted.bam
samtools index ${filename[k]}_sorted.bam
bamToBed -i ${filename[k]}_sorted.bam > ${filename[k]}_sorted.bed
genomeCoverageBed -ibam ${filename[k]}_sorted.bam -bg -trackline -trackopts 'name="Mytrack" visibility=2 color=255.0.0' > ${filename[k]}_sorted.bedGraph

#validating library construction with picard

java -jar /home/ohganelab/analyzetools/picard/build/libs/picard.jar CollectInsertSizeMetrics I=${filename[k]}_sorted.bam O=${filename[k]}_insert_size_metrics.txt H=${filename[k]}_insert_size_metrics.pdf

#organize directory

mkdir fastp_report
mkdir trush
mkdir picard_report
mv report_fastp.html report_fastp.json fastp_report
mv ${filename[k]}.sam ${filename[k]}.uq.sam ${filename[k]}.uq.qc.sam ${filename[k]}.header.sam ${filename[k]}.uq.qc.hd.sam ${filename[k]}.bam  trush
mv ${filename[k]}_insert_size_metrics.txt ${filename[k]}_insert_size_metrics.pdf picard_report

rm -r trush
mkdir trush

#peak calling(macs2)

macs2 callpeak -t ${filename[k]}_sorted.bam --name ${filename[k]} \--nomodel --nolambda --keep-dup all --call-summits -B
mkdir macs2data

#see tss region with ngsplot

ngs.plot.r -G ${ngs_ref} -R tss -C ${filename[k]}_sorted.bam -T ${filename[k]} -O FAIRE.tss
ngs.plot.r -G ${ngs_ref} -R genebody -C ${filename[k]}_sorted.bam -T ${filename[k]} -O FAIRE.genebody

mkdir ngsplotdata
mkdir ${dir}/${filename[k]}/ngsplotdata/tss
mkdir ${dir}/${filename[k]}/ngsplotdata/genebody
mv FAIRE.tss.avgprof.pdf FAIRE.tss.heatmap.pdf ${filename[k]}_sorted.bam.cnt FAIRE.tss.zip  ${dir}/${filename[k]}/ngsplotdata/tss
mv FAIRE.genebody.avgprof.pdf FAIRE.genebody.heatmap.pdf FAIRE.genebody.zip  ${dir}/${filename[k]}/ngsplotdata/genebody

#peaks annotation with homer

cat ${filename[k]}_peaks.narrowPeak |perl -e 'while(<>){chomp; @array=split/\t/; print "chr$array[0]\t$array[1]\t$array[2]\t$array[4]\t$array[3]\n"}'>${filename[k]}.mspeaks.bed 
bed2pos.pl ${filename[k]}.mspeaks.bed >${filename[k]}.mspeaks.hb
annotatePeaks.pl ${filename[k]}.mspeaks.hb ${ref_gene} >${filename[k]}_peaks.annotated.txt

mv ${filename[k]}_control_lambda.bdg ${filename[k]}_summits.bed ${filename[k]}_peaks.narrowPeak ${filename[k]}_treat_pileup.bdg ${filename[k]}_peaks.xls macs2data

mkdir peak_annotation
mv ${filename[k]}_peaks.annotated.txt peak_annotation

mv fastp_${filename[k]}_S1_L001_R1_001.fastq fastp_${filename[k]}_S1_L001_R2_001.fastq ${filename[k]}.mspeaks.bed ${filename[k]}.mspeaks.hb trush
