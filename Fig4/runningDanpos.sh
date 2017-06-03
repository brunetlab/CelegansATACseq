# get setup
cd /Users/acd13/Desktop/ATAC/Analysis/nucleosomeCalls

mkdir -p H3_EE_v_L3_fromDanpos
cd H3_EE_v_L3_fromDanpos
for dir in EE_H3 EE_H3_INPUT L3_H3 L3_H3_INPUT
do
	mkdir -p $dir
done

# link to the necessary files
for f in /Users/acd13/Desktop/ChIPs/EE/AB1791_H3/bams/AB1791_H3_Eemb_*_filter.filt.nodup.srt.ba*
do
	ln -s $f EE_H3/$(basename $f)
done

mv EE_H3/*INPUT* EE_H3_INPUT/
rm EE_H3/*INPUT*

for f in /Users/acd13/Desktop/ChIPs/L3/AB1791_H3/bams/AB1791_H3_L3_*_filter.filt.nodup.srt.ba*
do
	ln -s $f L3_H3/$(basename $f)
done

mv L3_H3/*INPUT* L3_H3_INPUT/
L3_H3/*INPUT*

# actually run the program
python /Users/acd13/Softwares/danpos-2.2.2/danpos.py dpos EE_H3:L3_H3 -b EE_H3:EE_H3_INPUT,L3_H3:L3_H3_INPUT > runningDanpos_H3.log 2>&1

for wig in $(find result/ -name '*.wig')
do
	python /Users/acd13/Dropbox/Scripts/pythonScripts/wigToBedGraph.py -w $wig > "${wig}.bg"
done
