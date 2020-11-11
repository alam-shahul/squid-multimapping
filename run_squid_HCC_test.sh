#!/bin/bash

argv=("$@")

OutputDir="" # downloaded reference and RNA-seq folder, the same as output folder
DataHCC=""
threads=2
seed=666
ExeDir=$(pwd)

for ((i=0; i<${#argv}; i+=2)); do
	if [[ ${argv[$i]} == "--fusioncatcher" ]]; then
		FusionCatcherExe=${argv[(($i+1))]}
	elif [[ ${argv[$i]} == "--fusioncatcherdata" ]]; then
		FusionCatcherDataDir=${argv[(($i+1))]}
	elif [[ ${argv[$i]} == "--chimeriscan" ]]; then
		ChimerascanDir=${argv[(($i+1))]}
	elif [[ ${argv[$i]} == "--jaffa" ]]; then
		JaffaDir=${argv[(($i+1))]}
	elif [[ ${argv[$i]} == "--defuse" ]]; then
		DefuseDir=${argv[(($i+1))]}
	elif [[ ${argv[$i]} == "--picard" ]]; then
		PicardDir=${argv[(($i+1))]}
	elif [[ ${argv[$i]} == "--seed" ]]; then
		seed=${argv[(($i+1))]}
	elif [[ ${argv[$i]} == "--datadownload" ]]; then
		OutputDir=${argv[(($i+1))]}
	elif [[ ${argv[$i]} == "--dataHCC" ]]; then
		DataHCC=${argv[(($i+1))]}
	elif [[ ${argv[$i]} == "--threads" ]]; then
		threads=${argv[(($i+1))]}
	fi
done

GenomeFasta=${OutputDir}/Ensemble75/Homo_sapiens.GRCh37.75.dna_sm.primary_assembly.fa
AnnotationGTF=${OutputDir}/Ensemble75/Homo_sapiens.GRCh37.75.gtf

# preparing STAR index
echo "Generating STAR index"
mkdir -p ${OutputDir}/Ensemble75/STARIndex
STAR --runThreadN ${threads} --runMode genomeGenerate --genomeDir ${OutputDir}/Ensemble75/STARIndex --genomeFastaFiles ${GenomeFasta} --sjdbGTFfile ${AnnotationGTF} --genomeSAsparseD 3 --genomeSAindexNbases 12

echo $OutputDir
mkdir -p ${OutputDir}
# 
cd $OutputDir
# align reads using STAR
echo "STAR alignment"
mkdir -p ${OutputDir}/${seed}/StarAlign/HCC1954
STAR  --runThreadN ${threads} --genomeDir ${OutputDir}/Ensemble75/STARIndex --readFilesIn ${OutputDir}/RNAseq/RNAHCC1954_1.fastq ${OutputDir}/RNAseq/RNAHCC1954_2.fastq --outFileNamePrefix ${OutputDir}/${seed}/StarAlign/HCC1954/ --outSAMtype BAM Unsorted --outSAMstrandField intronMotif --chimSegmentMin 15 --chimOutType Junctions SeparateSAMold --outReadsUnmapped Fastx --outMultimapperOrder Random --runRNGseed ${seed}
samtools view -Shb ${OutputDir}/${seed}/StarAlign/HCC1954/Chimeric.out.sam -o ${OutputDir}/${seed}/StarAlign/HCC1954/Chimeric.out.bam
samtools sort -o ${OutputDir}/${seed}/StarAlign/HCC1954/Aligned.sortedByCoord.out.bam ${OutputDir}/${seed}/StarAlign/HCC1954/Aligned.out.bam

# SQUID
echo "predicting TSV with SQUID"
mkdir -p ${OutputDir}/${seed}/TSVprediction/squid_1954
squid -b ${OutputDir}/${seed}/StarAlign/HCC1954/Aligned.sortedByCoord.out.bam -c ${OutputDir}/${seed}/StarAlign/HCC1954/Chimeric.out.bam -o ${OutputDir}/${seed}/TSVprediction/squid_1954/squidtsv
awk 'BEGIN{FS="\t";} {if(substr($0,0,1)=="#" || !($1~"M" || $1~"G" || $4~"M" || $4~"G")) print $0;}' ${OutputDir}/${seed}/TSVprediction/squid_1954/squidtsv_sv.txt > ${OutputDir}/${seed}/TSVprediction/squid_1954/squidtsv_sv_final.txt
 awk 'BEGIN{FS="\t";} {if(substr($0,0,1)=="#" || !($1~"M" || $1~"G" || $4~"M" || $4~"G")) print $0;}' ${OutputDir}/TSVprediction/squid_1395/squidtsv_sv.txt > ${OutputDir}/TSVprediction/squid_1395/squidtsv_sv_final.txt
