source('/Users/acd13/Desktop/ATAC/visualizing/gviz/gvizSupportFunctions_noLog.R')
setwd('/Users/acd13/Desktop/ATAC/visualizing/gviz/differentialExamples/')


# EE > L3
xol1Locus <- 'X:8,038,000-8,045,000'
plotEE_v_L3_withDiffPeaks(xol1Locus, 10, 'up', withChromStates=TRUE, plotName="xol1_EEvL3_10bpMeanEnrich_withChromHMM_noLog")
plotEE_v_L3_withDiffPeaks(xol1Locus, 10, 'up', withChromStates=FALSE, plotName="xol1_EEvL3_10bpMeanEnrich_noLog")
plotAllStagesAndMetaPeaks(xol1Locus, 10, plotName="xol1_allStages_10bpMeanEnrich_noLog")
plotEE_v_L3_withDiffPeaks(xol1Locus, 25, 'up', withChromStates=FALSE, plotName="xol1_EEvL3_25bpMeanEnrich_noLog")
plotAllStagesAndMetaPeaks(xol1Locus, 25, plotName="xol1_allStages_25bpMeanEnrich_noLog")

sdz21Locus <- "III:5,654,638-5,660,251"
plotAllStagesAndMetaPeaks(sdz21Locus, 25, revStrand=TRUE, plotName="sdz21_allStages_25bpMeanEnrich_noLog")
plotEE_v_L3_withDiffPeaks(sdz21Locus, 25, revStrand=TRUE, 'up', withChromStates=TRUE, plotName="sdz21_EEvL3_25bpMeanEnrich_withChromHMM_noLog")
plotEE_v_L3_withDiffPeaks(sdz21Locus, 25, revStrand=TRUE, 'up', withChromStates=FALSE, plotName="sdz21_EEvL3_25bpMeanEnrich_noLog")
plotAllStagesAndMetaPeaks(sdz21Locus, 25, revStrand=TRUE, plotName="sdz21_allStages_25bpMeanEnrich_noLog")

cav1 <- 'IV:9,770,419-9,774,512'
plotAllStagesAndMetaPeaks(cav1, 10, plotName="cav1_allStages_10bpMeanEnrich_noLog")
plotEE_v_L3_withDiffPeaks(cav1, 10, 'up', withChromStates=TRUE, plotName="cav1_EEvL3_10bpMeanEnrich_withChromHMM_noLog")
plotEE_v_L3_withDiffPeaks(cav1, 10, 'up', withChromStates=FALSE, plotName="cav1_EEvL3_10bpMeanEnrich_noLog")
plotAllStagesAndMetaPeaks(cav1, 10, plotName="cav1_allStages_10bpMeanEnrich_noLog")


#ergo1 <- 'V:1,002,893-1,018,480'
#plotEE_v_L3_withDiffPeaks(ergo1, 25, 'up')

#srz64 <- 'II:1,288,635-1,297,075'
#plotEE_v_L3_withDiffPeaks(srz64, 25, 'up')
# not enough known about the gene
#plotAllStagesAndMetaPeaks(srz64, 25)

#ceh13 <- 'III:7,543,624-7,562,314'
#plotEE_v_L3_withDiffPeaks(ceh13, 25, 'up')

#nob1Php3 <- 'III:12,071,738-12,101,388'
#plotEE_v_L3_withDiffPeaks(nob1Php3, 25, 'up')

#mex3 <-'I:120,771-137,676'
#plotEE_v_L3_withDiffPeaks(mex3, 25, 'up')

#pop1Plus <- 'I:2,815,611-2,850,137'
#plotEE_v_L3_withDiffPeaks(pop1Plus, 25, 'up')

#mom2 <- 'V:8,353,017-8,359,248'
#plotEE_v_L3_withDiffPeaks(mom2, 25, 'up')

#mom5 <- 'I:9,961,992-9,969,845'
#plotEE_v_L3_withDiffPeaks(mom5, 25, 'up')

pes2.1 <- 'I:11,347,805-11,363,784'
plotEE_v_L3_withDiffPeaks(pes2.1, 25, 'up')
plotAllStagesAndMetaPeaks(pes2.1, 25)

#cfz2 <- 'V:3,433,290-3,459,911'
#plotEE_v_L3_withDiffPeaks(cfz2, 25, 'up', revStrand=TRUE)

#elt1 <- 'IV:9,601,862-9,625,436'
#plotEE_v_L3_withDiffPeaks(elt1, 25, 'up')

#ceh36 <- 'X:14,181,834-14,197,629'
#plotEE_v_L3_withDiffPeaks(ceh36, 25, 'up')
#plotAllStagesAndMetaPeaks(ceh36, 25)


#mls2 <- 'X:4,755,135-4,766,999'
#plotEE_v_L3_withDiffPeaks(mls2, 25, 'up')




# these are examples that are higher in L3 than EE
daf12 <- 'X:10,629,613-10,667,976'
plotEE_v_L3_withDiffPeaks(daf12, 25, 'down', withChromStates=TRUE, plotName="daf12_EEvL3_25bpMeanEnrich_withChromHMM_noLog")
plotEE_v_L3_withDiffPeaks(daf12, 25, 'down', withChromStates=FALSE, plotName="daf12_EEvL3_25bpMeanEnrich_noLog")
plotAllStagesAndMetaPeaks(daf12, 25, plotName="daf12_allStages_25bpMeanEnrich_noLog")




hlh11Locus <- "III:9,617,000-9,629,000"
plotEE_v_L3_withDiffPeaks(hlh11Locus, 25,'down', withChromStates=TRUE, plotName="hlh11_EEvL3_25bpMeanEnrich_withChromHMM_noLog")


lin41 <- 'I:9,334,200-9,346,087'
plotEE_v_L3_withDiffPeaks(lin41, 25, 'down', withChromStates=TRUE, revStrand=TRUE, plotName="lin41_EEvL3_25bpMeanEnrich_withChromHMM_noLog")
plotEE_v_L3_withDiffPeaks(lin41, 50, 'down', withChromStates=TRUE, revStrand=TRUE, plotName="lin41_EEvL3_50bpMeanEnrich_withChromHMM_noLog")


mlk1 <- 'V:5,056,492-5,071,676'
plotEE_v_L3_withDiffPeaks(mlk1, 25, 'down', withChromStates=TRUE, revStrand=TRUE, plotName="mlk1_EEvL3_25bpMeanEnrich_withChromHMM_noLog")



tax6 <- 'IV:10,491,927-10,509,794'
plotEE_v_L3_withDiffPeaks(tax6, 25, 'down', withChromStates=TRUE, revStrand=TRUE, plotName="tax6_EEvL3_25bpMeanEnrich_withChromHMM_noLog")


kin1 <- 'I:14,916,548-14,950,960'
plotEE_v_L3_withDiffPeaks(kin1, 25, 'down', withChromStates=TRUE, plotName="kin1_EEvL3_25bpMeanEnrich_withChromHMM_noLog")


frm1 <- 'I:14,890,765-14,908,000'
plotEE_v_L3_withDiffPeaks(frm1, 25, 'down', withChromStates=TRUE, revStrand=TRUE, plotName="frm1_EEvL3_25bpMeanEnrich_withChromHMM_noLog")
#


din1 <- 'II:11,591,659-11,656,458'
plotEE_v_L3_withDiffPeaks(din1, 25, 'down', withChromStates=TRUE)

gpa16 <- 'I:894,449-908,179'
plotEE_v_L3_withDiffPeaks(gpa16, 25, 'down', withChromStates=TRUE)

#npr23 <- 'I:1,055,794-1,067,773'
#plotEE_v_L3_withDiffPeaks(npr23, 25, 'down', withChromStates=TRUE)

#kin29 <- 'X:2,841,921-2,854,471'
#plotEE_v_L3_withDiffPeaks(kin29, 25, 'down', withChromStates=TRUE)

#unc23 <- 'V:8,932,462-8,943,553'
#plotEE_v_L3_withDiffPeaks(unc23, 25, 'down', withChromStates=TRUE)

pde6 <- 'I:871,484-884,716'
plotEE_v_L3_withDiffPeaks(pde6, 25, 'down', withChromStates=TRUE)

#gsa1 <- 'I:1,063,814-1,079,793'
#plotEE_v_L3_withDiffPeaks(gsa1, 25, 'down', withChromStates=TRUE)

#pck1 <- 'III:61,035-69,053'
#plotEE_v_L3_withDiffPeaks(pck1, 25, 'down', withChromStates=TRUE)

#crh1 <- 'III:11,646,429-11,704,885'
#plotEE_v_L3_withDiffPeaks(crh1, 25, 'down', withChromStates=TRUE)

#lim9 <- 'I:9,164,110-9,228,032'
#plotEE_v_L3_withDiffPeaks(lim9, 25, 'down', withChromStates=TRUE)

gsy1 <- 'II:12,865,602-12,880,437'
plotEE_v_L3_withDiffPeaks(gsy1, 25, 'down', withChromStates=TRUE)

#hrg1 <- 'X:3,999,201-4,018,037'
#plotEE_v_L3_withDiffPeaks(hrg1, 25, 'down', withChromStates=TRUE)

sma2 <- 'III:8,747,691-8,758,988'
plotEE_v_L3_withDiffPeaks(sma2, 25, 'down', withChromStates=TRUE)

unc43 <- 'IV:10,321,095-10,351,796'
plotEE_v_L3_withDiffPeaks(unc43, 25, 'down', withChromStates=TRUE)

#col10 <- 'V:9,162,274-9,167,284'
#plotEE_v_L3_withDiffPeaks(col10, 10, 'down', withChromStates=TRUE)

#unc14 <- 'I:7,121,227-7,132,224'
#plotEE_v_L3_withDiffPeaks(unc14, 25, 'down', withChromStates=TRUE)

wdr23 <- 'I:7,594,503-7,604,195'
plotEE_v_L3_withDiffPeaks(wdr23, 25, 'down', withChromStates=TRUE)

#sma10 <- 'IV:273,378-283,615'
#plotEE_v_L3_withDiffPeaks(sma10, 25, 'down', withChromStates=TRUE)

#nhr143 <- 'V:8,971,165-8,976,904'
#plotEE_v_L3_withDiffPeaks(nhr143, 25, 'down', withChromStates=TRUE)

unc39 <- 'V:14,362,627-14,376,010'
plotEE_v_L3_withDiffPeaks(unc39, 25, 'down', withChromStates=TRUE)

# Others I tried
gap3Locus <- "I:1,987,659-2,034,549"
plotEE_v_L3_withDiffPeaks(gap3Locus, 50,'down', withChromStates=TRUE) #, plotName="gap3_EEvL3_25bpMeanEnrich")

gap3Locus2 <- "I:2,015,659-2,034,549"
plotEE_v_L3_withDiffPeaks(gap3Locus2, 50,'down', withChromStates=TRUE) #, plotName="gap3_EEvL3_25bpMeanEnrich")

#dpy5 <- 'I:5,431,021-5,434,034'
#plotEE_v_L3_withDiffPeaks(dpy5, 10, 'down', withChromStates=TRUE)

snap1 <- 'V:8,126,282-8,137,373'
plotEE_v_L3_withDiffPeaks(snap1, 25, 'down', withChromStates=TRUE)

grh1 <- 'I:1,256,612-1,281,507'
plotEE_v_L3_withDiffPeaks(grh1, 25, 'down', withChromStates=TRUE)

#tba2 <- 'I:12,962,668-12,968,865'
#plotEE_v_L3_withDiffPeaks(tba2, 25, 'down', withChromStates=TRUE)

soc2 <- 'IV:5,117,852-5,154,947'
plotEE_v_L3_withDiffPeaks(soc2, 25, 'down', withChromStates=TRUE)

#pept1 <- 'X:6,454,582-6,463,974'
#plotEE_v_L3_withDiffPeaks(pept1, 25, 'down', withChromStates=TRUE)

#pyp1 <- 'IV:9,993,443-10,000,904'
#plotEE_v_L3_withDiffPeaks(pyp1, 25, 'down', withChromStates=TRUE)

#soc1 <- 'V:4,637,728-4,648,410'
#plotEE_v_L3_withDiffPeaks(soc1, 25, 'down', withChromStates=TRUE)

#cct6 <- 'III:5,852,521-5,857,004'
#plotEE_v_L3_withDiffPeaks(cct6, 25, 'down', withChromStates=TRUE)

#rps12 <- 'III:5,678,774-5,681,016'
#plotEE_v_L3_withDiffPeaks(rps12, 10, 'down', withChromStates=TRUE)

#mtr4 <- 'IV:9,823,678-9,832,596'
#plotEE_v_L3_withDiffPeaks(mtr4, 25, 'down', withChromStates=TRUE)

#act5 <- 'III:13,597,796-13,612,409'
#plotEE_v_L3_withDiffPeaks(act5, 25, 'down', withChromStates=TRUE)



