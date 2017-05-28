#!/bin/bash

#### This will take tophat2 mapped reads and call HTseq to quantify within 
# The gff3 of exons only used here is from gerstein et al, 2014
if [[ "$#" -lt 1 ]]
then
    echo "USEAGE: SE or PE read mapping:" 1>&2
    echo "This program accepts trimGalore'd fastq files and aligns the reads to ce10Mito with bowtie2 via tophat2" 1>&2
    echo "$(basename $0) [BAMDir]"  1>&2
    echo "   [BAMDir]: directory of directories, each containing an accepted_hits.bam" 1>&2
    exit 1
fi

gtf="/Users/acd13/Desktop/ucscRefseqGenes26Oct2014/ucscRefSeqGenes_26Oct2014.gtf"

BAM_DIR=$(echo $1 | sed 's:/$::g')

for DIR in $(find "${BAM_DIR}" -type d)
do
	if [[ -e "${DIR}/accepted_hits.bam" ]]
	then
		cd "${DIR}"
		justDirName=$( pwd | perl -lane '$name=(split/\//)[-1]; print $name;')
		echo $justDirName
		bam="${DIR}/accepted_hits.bam"
		outoutFile="${justDirName}_htseqCountsOverUcscRefSeqGenesOnly.txt"
		htseq-count -f bam -r pos -i gene_id -s no -a 30 -t exon -m union "${bam}" "${gtf}" > "${outoutFile}"
	fi
done
