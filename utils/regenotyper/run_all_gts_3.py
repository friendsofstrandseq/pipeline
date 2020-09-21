import numpy as np
import os, os.path 
import glob, pdb
import re

# Collect foldernames with results
dirs_pre = '/home/hoeps/PhD/projects/huminvs/mosaicatcher/analysis/results/sept18_U32/sv_probabilities/*'
c_dir = '/home/hoeps/PhD/projects/huminvs/mosaicatcher/tracks2/Mapping_Normalization/CN_track_plots/result/nonred_inversions_n32_CN.txt'
dirs = glob.glob(dirs_pre)
# Make a regex object that will search for HG, GM or NA - the typical starts of samplenames.
regex = re.compile(r'HG|GM|NA')



# Run every sample in regenotyper
for dir in dirs:
    print(dir)
    if True:
        # Extract 'real' samplenames by looking for HG, GM and NA and extracting that +7 characters
        samplename = dir[regex.search(dir).span()[0]:regex.search(dir).span()[0]+7]
        # Find the accodring bedfile
        print(samplename)
        #bedpath = glob.glob('/home/hoeps/PhD/projects/huminvs/mosaicatcher/bed_factory/audano_3/wgot/done_withnames/newnames/{}*'.format(samplename))[0]
        try:
            c_dir = glob.glob('/home/hoeps/PhD/projects/huminvs/mosaicatcher/tracks/CN_map_U24/{}*'.format(samplename[2:]))[0]
        except:
            c_dir = glob.glob('/home/hoeps/PhD/projects/huminvs/mosaicatcher/tracks/CN_map_U24/{}*'.format('00733'))[0]
        if not os.path.exists(c_dir):
            print('missing dir: ', c_dir)
        #else:
        #    print('exists', c_dir)
            # Run 
        os.system('Rscript regenotype.R -f {}/100000_fixed_norm.selected_j0.1_s0.1/probabilities.Rdata -o {}/res/ -c {}'.format(dir, dir, c_dir))
    if False:
        print('error')
