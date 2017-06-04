# this is specifically to run on scg3 so I just hardcoded everything in
### Setup everything
module load bedtools/2.21.0 # to stay consistent with what I run on my own machine

regionSuffix='withoutProms'
FILE=$1 # the file is passed in

OUTPUT_DIR_BASE="/srv/gsfs0/projects/brunet/Aaron/distalNoncoding/shufflingDistalRegions/${regionSuffix}/"
REGION_TO_EXCLUDE='/srv/gsfs0/projects/brunet/Aaron/distalNoncoding/rawData/ExtendedPromAndExonsMerged.bed'

MAPPABLE_GENOME="/srv/gsfs0/projects/brunet/Aaron/distalNoncoding/rawData/ce10.K100.mappable.subtract_blackList_N2GdnaMacs2.1q5e-2peaks0.5Overlap.bed"
GENOME="/srv/gsfs0/projects/brunet/Aaron/distalNoncoding/rawData/ce10.chrom.sizes"

PREFIX=$(basename "${FILE}" | rev | cut -d'.' -f2- | rev)
OUTPUT_DIR="${OUTPUT_DIR_BASE}${PREFIX}"
[[ ! -d "${OUTPUT_DIR}" ]] && mkdir "${OUTPUT_DIR}"
TRIMMED_BED="${OUTPUT_DIR}/${PREFIX}_trimmed.bed"


# Trim the bed file to the minimal info (i.e. just position)
cut -f 1-3 "${FILE}" > "${TRIMMED_BED}"
	
for i in {1..10000}
do
	OUT_FILE="${OUTPUT_DIR}/${PREFIX}.${regionSuffix}.shuffled${i}.bed.gz"
	intersectBed -f 0.5 -v -a "${TRIMMED_BED}" -b "${REGION_TO_EXCLUDE}" | \
	shuffleBed -incl "${MAPPABLE_GENOME}" -excl "${REGION_TO_EXCLUDE}" -f 0.5 -i - -g "${GENOME}" -noOverlapping -seed "${i}" -maxTries 10000 | \
	gzip -c > \
	"${OUT_FILE}"
done

rm ${TRIMMED_BED}
