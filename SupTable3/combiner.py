def get_key(l):
    return "_".join(l[:3])
def rename_gene_loc(s):
    if s == 'inExtProm' or s == 'inCoreProm':
        return 'promoter'
    elif s == 'inIntrons':
        return 'intron'    
    elif s == 'inExons':
        return 'exon'
    elif s == 'within0.5kbTES':
        return 'nearTTS'
    return s
    
def add_diff(d, fn, comp):
    with open(fn, 'r') as f:
        for l in f:
            frs = l.rstrip("\n").split("\t")
            k = get_key(frs[5:])
            d[k][comp + "_FC"] = "{0:.3f}".format(float(frs[3]))
            d[k][comp + "_FDR"] = "{0:.2e}".format(float(frs[4]))
    return d

def add_chromHMM(d, fn, stage):
    with open(fn, 'r') as f:
        for l in f:
            frs = l.rstrip("\n").split("\t")
            k = get_key(frs)
            d[k][stage + '_chromHMMState'] = rename_chromHMM(frs[4], stage)
    return d

def rename_chromHMM(s, stage):
    if s == stage+'NotInActiveHetRep':
        return stage+'_lowSignal'
    elif s in ['L3_allActiveEnhancers_notRepEnh13', 'EEActiveEnhancers_notRepEnh13', 'allYAEnh_notRepEnh13']:
        return stage+'_activeEnhancer'
    elif s in ['L3_repressedEnh', 'EE_RepEnh', 'YA_repEnh']:
        return stage+'_repressedEnhancer'
    elif s in ['allL3HetChromMerged', 'allEEHetChromMerged', 'allYAHetChrom']:
        return stage+'_heterochromatin'
    elif s in ['allL3RepMerged_notRepEnh13', 'allEERepMerged_notRepEnh13', 'allYARep_notRepEnh13']:
        return stage+'_H3K27me3Repressed'
    elif s in ['allL3TSSMerged', 'allEETSSMerged', 'allYATSS']:
        return stage+'_promoterTSS'
    elif s in ['allL3TxMerged', 'allEETxMerged', 'allYATxMerged']:
        return stage+'_transcribedGeneBody'
    

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='combine files')
    parser.add_argument("nearestTSS")
    parser.add_argument("EEChrom")
    parser.add_argument("L3Chrom")
    parser.add_argument("YAChrom")
    parser.add_argument("GeneLoc")
    parser.add_argument("EEvL3")
    parser.add_argument("L3vYA")

    args=parser.parse_args()
    
    out_header = ['chr', 'start', 'stop', 'peak', 'nearestTSS', 'distToTSS']
    all_peaks = {}
    with open(args.nearestTSS, 'r') as f:
        for l in f:
            frs = l.rstrip("\n").split("\t")
            k = get_key(frs)
            if k in all_peaks:
                frs[out_header.index('nearestTSS')] = frs[out_header.index('nearestTSS')] + ',' + all_peaks[k]['nearestTSS']
            all_peaks[k] = {}
            for i in xrange(len(frs)):
            	all_peaks[k][out_header[i]] = frs[i]
            	
    out_header.append('GeneBasedLocation')
    with open(args.GeneLoc, 'r') as f:
        for l in f:
            frs = l.rstrip("\n").split("\t")
            k = get_key(frs)
            all_peaks[k]['GeneBasedLocation'] = rename_gene_loc(frs[4])

    out_header += ['EEvL3_FC', 'EEvL3_FDR', 'L3vYA_FC', 'L3vYA_FDR']
    all_peaks = add_diff(all_peaks, args.EEvL3, 'EEvL3')
    all_peaks = add_diff(all_peaks, args.L3vYA, 'L3vYA')

    out_header += [stage + '_chromHMMState' for stage in ['EE', 'L3', 'YA']]
    all_peaks = add_chromHMM(all_peaks, args.EEChrom, 'EE')
    all_peaks = add_chromHMM(all_peaks, args.L3Chrom, 'L3')
    all_peaks = add_chromHMM(all_peaks, args.YAChrom, 'YA')

    with open ('fullyAnnotated_consensusAtacPeaks.bed', 'w') as f:
        f.write("\t".join(out_header) + "\n")
        for k in all_peaks:
            out = []
            for h in out_header:
                try:
                    out.append(all_peaks[k][h])
                except KeyError:
                    out.append('NA')
            f.write("\t".join(out) + "\n")
