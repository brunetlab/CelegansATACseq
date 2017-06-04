# split each state in each stage into a separate bed file

hiHMModel=$1
outDir=$2
mkdir -p $outDir

cut -f 4 $hiHMModel | sort -u > compare_hiHMMToMyChromHMM.tmp
# the result of that is 
#10_Rep1
#11_Rep2
#12_Het1
#13_Het2
#14_Low1
#15_Low2
#16_Low3
#17_Unmap
#1_Pro
#2_Enh1
#3_Enh2
#4_Egn1
#5_Egn2
#6_Egn3
#7_Egn4
#8_Egn5
#9_Egn6

while read hiState
do
    echo "Starting on ${hiState}"
    outFile="${outDir}/${hiState}.bed"
    rm -f ${outFile}
    grep ${hiState} ${hiHMModel} | perl -lane 'print join("\t", ("chr".$F[0], $F[1], $F[2], $F[3]));' >> ${outFile}
done < compare_hiHMMToMyChromHMM.tmp

rm compare_hiHMMToMyChromHMM.tmp

# then combine specific states together
mkdir -p $outDir/merged
cat "${outDir}/10_Rep1.bed" "${outDir}/11_Rep2.bed" | sortBed -i - | mergeBed -i - > "${outDir}/merged/mergedRep.bed"
cat "${outDir}/12_Het1.bed" "${outDir}/13_Het2.bed" | sortBed -i - | mergeBed -i - > "${outDir}/merged/mergedHet.bed"
cat "${outDir}/16_Low3.bed" "${outDir}/15_Low2.bed" "${outDir}/14_Low1.bed" | \
sortBed -i - | \
mergeBed -i - > "${outDir}/merged/mergedLow.bed"

# since I didn't include unmappable, those regions should be in the low bin, so I'll also merge that in
cat "${outDir}/16_Low3.bed" "${outDir}/15_Low2.bed" "${outDir}/14_Low1.bed" "${outDir}/17_Unmap.bed" | \
sortBed -i - | \
mergeBed -i - > "${outDir}/merged/mergedLowUnmap.bed"
cat "${outDir}/2_Enh1.bed" "${outDir}/3_Enh2.bed" | sortBed -i - | mergeBed -i - > "${outDir}/merged/mergedEnh.bed"
cat "${outDir}/4_Egn1.bed" "${outDir}/5_Egn2.bed" | sortBed -i - | mergeBed -i - > "${outDir}/merged/merged5pTx.bed"
cat "${outDir}/9_Egn6.bed" "${outDir}/8_Egn5.bed" "${outDir}/7_Egn4.bed" | \
sortBed -i - | \
mergeBed -i - > "${outDir}/merged/merged3pTx.bed"

# I don't differentiate b/t gene positions, so I'll also merge all of those
cat "${outDir}/merged/merged3pTx.bed" "${outDir}/merged/merged5pTx.bed" "${outDir}/6_Egn3.bed" | \
cut -f 1-3 | \
sortBed -i - | \mergeBed -i - > "${outDir}/merged/mergedAllTx.bed"
