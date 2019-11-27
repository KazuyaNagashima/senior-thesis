####set data########
filename1=("fibro3min-190515")
filename2=("fibro-3min-170624")
filename3=("fibro-3min-170728")
ngs_ref=("susScr11")
#2_homer.sh
date="19.11.20"
######################

# import data in directory

dir="/home/ohganelab/Desktop/nagashima"
cd ${dir}/comparison
#mkdir ${date}
cd ${date}
mkdir ${dir}/comparison/${date}/ngsplot/
cd ngsplot

cp ${dir}/${filename1}/${filename1}_sorted.bam ${dir}/comparison/${date}/ngsplot
cp ${dir}/${filename2}/${filename2}_sorted.bam ${dir}/comparison/${date}/ngsplot
cp ${dir}/${filename3}/${filename3}_sorted.bam ${dir}/comparison/${date}/ngsplot

echo ${filename1}_sorted.bam    -1  "${filename1}">list.txt
echo ${filename2}_sorted.bam    -1  "${filename2}">>list.txt
echo ${filename3}_sorted.bam  -1    "${filename3}">>list.txt

ngs.plot.r -G ${ngs_ref} -D ensemble -R genebody -C list.txt -O comparison_genebody -L 4000 -FL 200
ngs.plot.r -G ${ngs_ref} -D ensemble -R tss -C list.txt -O comparison_tss -L 4000 -FL 200
ngs.plot.r -G ${ngs_ref} -D ensemble -R genebody -C list.txt -O k-means -L 4000 -FL 200 -GO km -KNC 10 -NRS 200

mkdir genebody
mv comparison_genebody.avgprof.pdf comparison_genebody.heatmap.pdf comparison_genebody.zip genebody
mkdir tss
mv comparison_tss.avgprof.pdf comparison_tss.heatmap.pdf comparison_tss.zip tss
mkdir k-means
mv k-means.avgprof.pdf k-means.heatmap.pdf k-means.zip k-means
mkdir trush
mv ${filename1}_sorted.bam ${filename1}_sorted.bam.bai ${filename1}_sorted.bam.cnt ${filename2}_sorted.bam ${filename2}_sorted.bam.bai ${filename2}_sorted.bam.cnt ${filename3}_sorted.bam ${filename3}_sorted.bam.bai ${filename3}_sorted.bam.cnt list.txt trush
