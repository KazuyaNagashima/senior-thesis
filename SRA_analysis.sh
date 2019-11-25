#####set data#########
filename=("ERR3154128")
ref_gene=("susScr11")
######################

###SINGLE END###

dir="/home/ohganelab/Desktop/nagashima"
cd ${dir}
mkdir ${dir}/${filename}

#move file to a new directory

mv ${filename}.sra ${dir}/${filename}
cd ${dir}/${filename}

#use SRA Toolkit
fastq-dump ${filename}.sra --gzip

#quality check by fastp

fastp -i ${filename}.fastq.gz -o fastp_${filename}.fastq.gz -h report_fastp.html -j report_fastp.json -q 25 -n 10 -w 14 --length_required 20 --length_limit 500

#extract files

gunzip fastp_${filename}.fastq.gz

#mapping by bowtie2

bowtie2 -p 18 -x /home/ohganelab/analyzetools/reference/Sus_scrofa/Ensembl/Sscrofa11.1/Sequence/Bowtie2Index/genome -U fastp_${filename}.fastq -S ${filename}.sam

#remove multi-mapped reads

cat ${filename}.sam | perl -e 'while(<>){print $_ if(/\@SQ||\@PG/);}' >${filename}.header.sam
grep -v "XS" ${filename}.sam >${filename}.uq.sam

#quality check by samtools (MAPQ>4)

samtools view -q 4 ${filename}.uq.sam > ${filename}.uq.qc.sam

#from sam to bam,bed,bedgraph

cat ${filename}.header.sam ${filename}.uq.qc.sam >${filename}.uq.qc.hd.sam
samtools view -b -S ${filename}.uq.qc.hd.sam > ${filename}.bam
samtools sort ${filename}.bam > ${filename}_sorted.bam
samtools index ${filename}_sorted.bam
bamToBed -i ${filename}_sorted.bam > ${filename}_sorted.bed
genomeCoverageBed -ibam ${filename}_sorted.bam -bg -trackline -trackopts 'name="Mytrack" visibility=2 color=255.0.0' > ${filename}_sorted.bedGraph

#validating library construction by picard

java -jar /home/ohganelab/analyzetools/picard/build/libs/picard.jar CollectInsertSizeMetrics I=${filename}_sorted.bam O=${filename}_insert_size_metrics.txt H=${filename}_insert_size_metrics.pdf

#organize directory

mkdir fastp_report
mkdir trush
mkdir picard_report
mv report_fastp.html report_fastp.json fastp_report
mv ${filename}.sam ${filename}.uq.sam ${filename}.uq.qc.sam ${filename}.header.sam ${filename}.uq.qc.hd.sam ${filename}.bam  trush
mv ${filename}_insert_size_metrics.txt ${filename}_insert_size_metrics.pdf picard_report

#peak calling(macs2)

macs2 callpeak -t ${filename}_sorted.bam --name ${filename} \--nomodel --nolambda --keep-dup all --call-summits -B
mkdir macs2data

#see tss region by ngsplot

ngs.plot.r -G ${ngs_ref} -R tss -C ${filename}_sorted.bam -T ${filename} -O FAIRE.tss
ngs.plot.r -G ${ngs_ref} -R genebody -C ${filename}_sorted.bam -T ${filename} -O FAIRE.genebody

mkdir ngsplotdata
mkdir ${dir}/${filename}/ngsplotdata/tss
mkdir ${dir}/${filename}/ngsplotdata/genebody
mv FAIRE.tss.avgprof.pdf FAIRE.tss.heatmap.pdf ${filename}_sorted.bam.cnt FAIRE.tss.zip  ${dir}/${filename}/ngsplotdata/tss
mv FAIRE.genebody.avgprof.pdf FAIRE.genebody.heatmap.pdf FAIRE.genebody.zip  ${dir}/${filename}/ngsplotdata/genebody

#peaks annotation by homer

cat ${filename}_peaks.narrowPeak |perl -e 'while(<>){chomp; @array=split/\t/; print "chr$array[0]\t$array[1]\t$array[2]\t$array[4]\t$array[3]\n"}'>${filename}.mspeaks.bed
bed2pos.pl ${filename}.mspeaks.bed >${filename}.mspeaks.hb
annotatePeaks.pl ${filename}.mspeaks.hb ${ref_gene} >${filename}_peaks.annotated.txt

mv ${filename}_control_lambda.bdg ${filename}_summits.bed ${filename}_peaks.narrowPeak ${filename}_treat_pileup.bdg ${filename}_peaks.xls macs2data

mkdir peak_annotation
mv ${filename}_peaks.annotated.txt peak_annotation

mv fastp_${filename}.fastq ${filename}.mspeaks.bed ${filename}.mspeaks.hb trush
