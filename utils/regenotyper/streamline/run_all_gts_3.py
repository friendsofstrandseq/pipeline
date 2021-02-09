import numpy as np
import os, os.path
import glob, pdb
import re
from pathlib import Path

cp_from_seneca = True
run_regenotyper = True
concat_results = True

### First: copy over the according files ###
path_to_results = '/home/hoeps/s/g/korbel2/StrandSeq/Test_WH/pipeline_7may/pipeline/sv_probabilities'
res_remote = f'{path_to_results}/*/100000_fixed_norm.selected_j0.1_s0.1/probabilities.Rdata'
debugfile = f'{path_to_results}/../counts/HG00733/manual_segments_counts.txt.debug'

res_local_raw = '/home/hoeps/PhD/projects/huminvs/mosaicatcher/analysis/results/review_peter4xx/'
res_local = f'{res_local_raw}sv_probabilities'

# There is a lot that could have gone wrong by now. Ask the user to verify that this is correct.
print(f'Copy from seneca: {cp_from_seneca}')
print(f'Path to seneca results: {path_to_results}')
print(f'Path to local results: {res_local}')
print(f'Paths correct? [y/n]')
yes = {'yes','y', 'ye'}
no = {'no','n'}

choice = input().lower()
if choice in yes:
   pass
elif choice in no:
   sys.exit("Please adjust paths and rerun the script.")
else:
   sys.stdout.write("Please respond with 'yes' or 'no'")

if cp_from_seneca:
    # cp the results
    print('Creating output folder')
    Path(f'{res_local}').mkdir(parents=True, exist_ok=True)

    print('Please make sure seneca is mounted')

    print(f'Copying {path_to_results} to {res_local}')
    os.system(f'''
    cp -r --parents {res_remote} {res_local}
    ''')

    os.system(f'''
    mv {res_local}{path_to_results}/* {res_local}/
    ''')
    # Also awk the debug file
    print('Also awking the debug file.')
    os.system(f'''
    awk '!seen[$1,$2,$3]++' {debugfile} > {res_local_raw}msc.debug
    ''')
    # clean up
    os.system(f'''
    sudo rm -r {res_local}/home
    ''')

### Second: run regenotyper ###
if run_regenotyper:
    # Collect foldernames with results
    dirs = glob.glob(f'{res_local}/*')
    print(dirs)
    # Make a regex object that will search for HG, GM or NA - the typical starts of samplenames.
    regex = re.compile(r"HG|GM|NA")

    # Run every sample in regenotyper
    for dir in dirs:
        if True:
            if dir.find("HG002") > -1:
                samplename = "HG002"
            else:
                # Extract 'real' samplenames by looking for HG, GM and NA and extracting that  7 characters
                print(dir)
                samplename = dir[regex.search(dir).span()[0]:regex.search(dir).span()[0]+7]
            # Find the accodring bedfile

            # bedpath = glob.glob('/home/hoeps/PhD/projects/huminvs/mosaicatcher/bed_factory/audano_3/wgot/done_withnames/newnames/{}*'.format(samplename))[0]
            try:
                bedpath = glob.glob(
                    "NOTHING/home/hoeps/PhD/projects/huminvs/mosaicatcher/david_list_15sep/david_list_26oct/individual_beds_gt/{}*".format(
                        samplename
                    )
                )[0]
            except:
                print(f"No bedfile found for {samplename}")
                bedpath = "none"
            try:
                c_dir = f'{res_local_raw}msc.debug'
            except:
                print("Missing CN annotation: {}".format(samplename))
                c_dir = glob.glob(
                    "/home/hoeps/PhD/projects/huminvs/mosaicatcher/tracks/tracks_hufsah_21sept_second/result/{}*".format(
                        "00733"
                    )
                )[0]
            if not os.path.exists(c_dir):
                print("missing dir: ", c_dir)
            # else:
            #    print('exists', c_dir)
            # Run
            # if ((dir.find('HG002')>-1) | (dir.find('HG03125')>-1) | (dir.find('HG03486')>-1) | (dir.find('NA12878')>-1)):
            print(dir)
            print(bedpath)
            if bedpath != "none":
                os.system(
                    f"Rscript ../regenotype.R \
                -f {dir}/100000_fixed_norm.selected_j0.1_s0.1/probabilities.Rdata \
                -b {bedpath} \
                -o {dir}/res/ \
                -c {c_dir}"
                )
            else:
                os.system(
                    f"/usr/bin/Rscript ../regenotype.R \
                    -f {dir}/100000_fixed_norm.selected_j0.1_s0.1/probabilities.Rdata \
                    -o {dir}/res/ \
                    -c {c_dir}"
                )
        if False:
            print("error")

    print('Regenotyper run finished.')

### Third: make genotypes ###

if concat_results:
    # Combine sv_calls_bulk.txt
    os.system(f'''
    awk 'FNR==1 && NR!=1 {{ while (/^chrom/) getline; }} 1 {{print}}' {res_local}/*/res/all/sv_calls_bulk.txt > all.txt
    ''')
