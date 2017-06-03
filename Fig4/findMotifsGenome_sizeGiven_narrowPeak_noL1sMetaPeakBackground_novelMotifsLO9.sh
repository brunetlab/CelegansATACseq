#!/bin/bash

if [[ "$#" -lt 1 ]]
then
    echo "$(basename $0) [npDir]"  1>&2
    echo "   [npDir]: " 1>&2
    exit 1
fi

npDir=$(echo $1 | sed 's:/$::g')

cd ${npDir}
mkdir -p newCEMotifsLO10
for f in $(find ${npDir} -name '*.narrowPeak')
do
        justName=$(basename "${f}" | sed 's/\.narrowPeak//g')
        findMotifsGenome.pl $f ce10 newCEMotifsLO10/metaPeakBackground/$justName \
            -size given -len 6,8,10,12 -bits \
            -bg /srv/gsfs0/projects/brunet/Aaron/homerMotifs/ATAC/ATAC_consensusPeaks/allRepsAllButL1MetaPeaks_smartMerged300bpSummitDist.bed \
            -mcheck /srv/gsfs0/projects/brunet/Aaron/homerMotifs/motifs/hughesCePWMSForHomerLO9.motif \
            -mknown /srv/gsfs0/projects/brunet/Aaron/homerMotifs/motifs/hughesCePWMSForHomerLO9.motif
done
