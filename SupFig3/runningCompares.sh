cd /Users/acd13/Desktop/ChromHMM
mkdir -p comparingToHiHMM

# first do it for EE
# I selected and drug these into the terminal, and I'll just write them to a file for ease

rm -f eeChromHMMStatesToCompareTo.txt
for f in /Users/acd13/Desktop/ChromHMM/ce10_ws220_EE_L3_YA/8marks_EE_L3_YA_1e-3_19states/reordered/beds/EE/EENotInActiveHetRep.bed /Users/acd13/Desktop/ChromHMM/ce10_ws220_EE_L3_YA/8marks_EE_L3_YA_1e-3_19states/reordered/beds/EE/EEActiveEnhancers_notRepEnh13.bed /Users/acd13/Desktop/ChromHMM/ce10_ws220_EE_L3_YA/8marks_EE_L3_YA_1e-3_19states/reordered/beds/EE/allEERepMerged_notRepEnh13.bed /Users/acd13/Desktop/ChromHMM/ce10_ws220_EE_L3_YA/8marks_EE_L3_YA_1e-3_19states/reordered/beds/EE/allEEHetChromMerged.bed /Users/acd13/Desktop/ChromHMM/ce10_ws220_EE_L3_YA/8marks_EE_L3_YA_1e-3_19states/reordered/beds/EE/allEETSSMerged.bed /Users/acd13/Desktop/ChromHMM/ce10_ws220_EE_L3_YA/8marks_EE_L3_YA_1e-3_19states/reordered/beds/EE/allEETxMerged.bed /Users/acd13/Desktop/ChromHMM/ce10_ws220_EE_L3_YA/8marks_EE_L3_YA_1e-3_19states/reordered/beds/EE/EE_E19.bed /Users/acd13/Desktop/ChromHMM/ce10_ws220_EE_L3_YA/8marks_EE_L3_YA_1e-3_19states/reordered/beds/EE/EE_E18.bed /Users/acd13/Desktop/ChromHMM/ce10_ws220_EE_L3_YA/8marks_EE_L3_YA_1e-3_19states/reordered/beds/EE/EE_E17.bed /Users/acd13/Desktop/ChromHMM/ce10_ws220_EE_L3_YA/8marks_EE_L3_YA_1e-3_19states/reordered/beds/EE/EE_E16.bed /Users/acd13/Desktop/ChromHMM/ce10_ws220_EE_L3_YA/8marks_EE_L3_YA_1e-3_19states/reordered/beds/EE/EE_E15.bed /Users/acd13/Desktop/ChromHMM/ce10_ws220_EE_L3_YA/8marks_EE_L3_YA_1e-3_19states/reordered/beds/EE/EE_E14.bed /Users/acd13/Desktop/ChromHMM/ce10_ws220_EE_L3_YA/8marks_EE_L3_YA_1e-3_19states/reordered/beds/EE/EE_RepEnh.bed /Users/acd13/Desktop/ChromHMM/ce10_ws220_EE_L3_YA/8marks_EE_L3_YA_1e-3_19states/reordered/beds/EE/EE_E13.bed /Users/acd13/Desktop/ChromHMM/ce10_ws220_EE_L3_YA/8marks_EE_L3_YA_1e-3_19states/reordered/beds/EE/EE_E12.bed /Users/acd13/Desktop/ChromHMM/ce10_ws220_EE_L3_YA/8marks_EE_L3_YA_1e-3_19states/reordered/beds/EE/EE_E11.bed /Users/acd13/Desktop/ChromHMM/ce10_ws220_EE_L3_YA/8marks_EE_L3_YA_1e-3_19states/reordered/beds/EE/EE_E10.bed /Users/acd13/Desktop/ChromHMM/ce10_ws220_EE_L3_YA/8marks_EE_L3_YA_1e-3_19states/reordered/beds/EE/EE_E9.bed /Users/acd13/Desktop/ChromHMM/ce10_ws220_EE_L3_YA/8marks_EE_L3_YA_1e-3_19states/reordered/beds/EE/EE_E8.bed /Users/acd13/Desktop/ChromHMM/ce10_ws220_EE_L3_YA/8marks_EE_L3_YA_1e-3_19states/reordered/beds/EE/EE_E7.bed /Users/acd13/Desktop/ChromHMM/ce10_ws220_EE_L3_YA/8marks_EE_L3_YA_1e-3_19states/reordered/beds/EE/EE_E6.bed /Users/acd13/Desktop/ChromHMM/ce10_ws220_EE_L3_YA/8marks_EE_L3_YA_1e-3_19states/reordered/beds/EE/EE_E5.bed /Users/acd13/Desktop/ChromHMM/ce10_ws220_EE_L3_YA/8marks_EE_L3_YA_1e-3_19states/reordered/beds/EE/EE_E4.bed /Users/acd13/Desktop/ChromHMM/ce10_ws220_EE_L3_YA/8marks_EE_L3_YA_1e-3_19states/reordered/beds/EE/EE_E3.bed /Users/acd13/Desktop/ChromHMM/ce10_ws220_EE_L3_YA/8marks_EE_L3_YA_1e-3_19states/reordered/beds/EE/EE_E2.bed /Users/acd13/Desktop/ChromHMM/ce10_ws220_EE_L3_YA/8marks_EE_L3_YA_1e-3_19states/reordered/beds/EE/EE_E1.bed
do
    echo $f >> eeChromHMMStatesToCompareTo.txt
done
bash compare_hiHMMToMyChromHMM.sh hiHMM_EE_fromModEncode.bed eeChromHMMStatesToCompareTo.txt comparingToHiHMM/EE

# and now for L3
rm -f l3ChromHMMStatesToCompareTo.txt
for f in /Users/acd13/Desktop/ChromHMM/ce10_ws220_EE_L3_YA/8marks_EE_L3_YA_1e-3_19states/reordered/beds/L3/L3NotInActiveHetRep.bed /Users/acd13/Desktop/ChromHMM/ce10_ws220_EE_L3_YA/8marks_EE_L3_YA_1e-3_19states/reordered/beds/L3/L3_allActiveEnhancers_notRepEnh13.bed /Users/acd13/Desktop/ChromHMM/ce10_ws220_EE_L3_YA/8marks_EE_L3_YA_1e-3_19states/reordered/beds/L3/allL3RepMerged_notRepEnh13.bed /Users/acd13/Desktop/ChromHMM/ce10_ws220_EE_L3_YA/8marks_EE_L3_YA_1e-3_19states/reordered/beds/L3/allL3HetChromMerged.bed /Users/acd13/Desktop/ChromHMM/ce10_ws220_EE_L3_YA/8marks_EE_L3_YA_1e-3_19states/reordered/beds/L3/allL3TSSMerged.bed /Users/acd13/Desktop/ChromHMM/ce10_ws220_EE_L3_YA/8marks_EE_L3_YA_1e-3_19states/reordered/beds/L3/allL3TxMerged.bed /Users/acd13/Desktop/ChromHMM/ce10_ws220_EE_L3_YA/8marks_EE_L3_YA_1e-3_19states/reordered/beds/L3/L3_E19.bed /Users/acd13/Desktop/ChromHMM/ce10_ws220_EE_L3_YA/8marks_EE_L3_YA_1e-3_19states/reordered/beds/L3/L3_E18.bed /Users/acd13/Desktop/ChromHMM/ce10_ws220_EE_L3_YA/8marks_EE_L3_YA_1e-3_19states/reordered/beds/L3/L3_E17.bed /Users/acd13/Desktop/ChromHMM/ce10_ws220_EE_L3_YA/8marks_EE_L3_YA_1e-3_19states/reordered/beds/L3/L3_E16.bed /Users/acd13/Desktop/ChromHMM/ce10_ws220_EE_L3_YA/8marks_EE_L3_YA_1e-3_19states/reordered/beds/L3/L3_E15.bed /Users/acd13/Desktop/ChromHMM/ce10_ws220_EE_L3_YA/8marks_EE_L3_YA_1e-3_19states/reordered/beds/L3/L3_E14.bed /Users/acd13/Desktop/ChromHMM/ce10_ws220_EE_L3_YA/8marks_EE_L3_YA_1e-3_19states/reordered/beds/L3/L3_repressedEnh.bed /Users/acd13/Desktop/ChromHMM/ce10_ws220_EE_L3_YA/8marks_EE_L3_YA_1e-3_19states/reordered/beds/L3/L3_E13.bed /Users/acd13/Desktop/ChromHMM/ce10_ws220_EE_L3_YA/8marks_EE_L3_YA_1e-3_19states/reordered/beds/L3/L3_E12.bed /Users/acd13/Desktop/ChromHMM/ce10_ws220_EE_L3_YA/8marks_EE_L3_YA_1e-3_19states/reordered/beds/L3/L3_E11.bed /Users/acd13/Desktop/ChromHMM/ce10_ws220_EE_L3_YA/8marks_EE_L3_YA_1e-3_19states/reordered/beds/L3/L3_E10.bed /Users/acd13/Desktop/ChromHMM/ce10_ws220_EE_L3_YA/8marks_EE_L3_YA_1e-3_19states/reordered/beds/L3/L3_E9.bed /Users/acd13/Desktop/ChromHMM/ce10_ws220_EE_L3_YA/8marks_EE_L3_YA_1e-3_19states/reordered/beds/L3/L3_E8.bed /Users/acd13/Desktop/ChromHMM/ce10_ws220_EE_L3_YA/8marks_EE_L3_YA_1e-3_19states/reordered/beds/L3/L3_E7.bed /Users/acd13/Desktop/ChromHMM/ce10_ws220_EE_L3_YA/8marks_EE_L3_YA_1e-3_19states/reordered/beds/L3/L3_E6.bed /Users/acd13/Desktop/ChromHMM/ce10_ws220_EE_L3_YA/8marks_EE_L3_YA_1e-3_19states/reordered/beds/L3/L3_E5.bed /Users/acd13/Desktop/ChromHMM/ce10_ws220_EE_L3_YA/8marks_EE_L3_YA_1e-3_19states/reordered/beds/L3/L3_E4.bed /Users/acd13/Desktop/ChromHMM/ce10_ws220_EE_L3_YA/8marks_EE_L3_YA_1e-3_19states/reordered/beds/L3/L3_E3.bed /Users/acd13/Desktop/ChromHMM/ce10_ws220_EE_L3_YA/8marks_EE_L3_YA_1e-3_19states/reordered/beds/L3/L3_E2.bed /Users/acd13/Desktop/ChromHMM/ce10_ws220_EE_L3_YA/8marks_EE_L3_YA_1e-3_19states/reordered/beds/L3/L3_E1.bed 
do
    echo $f >> l3ChromHMMStatesToCompareTo.txt
done
bash compare_hiHMMToMyChromHMM.sh hiHMM_L3_fromModEncode.bed l3ChromHMMStatesToCompareTo.txt comparingToHiHMM/L3