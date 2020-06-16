import argparse
import sys
import pysam
from pathlib import Path
import glob
import pdb

from collections import defaultdict
def counts(sample, input_bam, input_bed ,norm_count_output, mapping_counts, norm_plot_output):
    dictionary= defaultdict(lambda: defaultdict(tuple))
    mapping_counts_file = open(mapping_counts, 'r')
    for lines in mapping_counts_file:
        line= lines.strip().split("\t")
        chrom= line[0]
        interval_start= int(line[1])
        interval_end= int(line[2])
        reads_originated= int(line[3])
        reads_mapped = int(line[4])
        dictionary[chrom][(interval_start, interval_end)]=(reads_originated, reads_mapped)
    
    watson_count=0
    crick_count=0
    norm_counts_file= open(norm_count_output, 'w')
    # Write header to output file
    norm_counts_file.write('chrom'+ "\t"+ 'start'+ "\t" + 'end'+ "\t" + 'sample' + "\t" + 'cell' + "\t" + 'C'+ "\t" + 'W' +"\n")
    norm_plots_file= open(norm_plot_output, 'w')
    # Write header to output file
    norm_plots_file.write('chrom'+ "\t"+ 'start'+ "\t" + 'end'+ "\t" + 'sample' + "\t" + 'cell' + "\t" + 'C'+ "\t" + 'W' +"\n")
    
    # Get all bam file paths
    path = Path(input_bam)
    glob_path = path.glob('*.bam')
    
    # Iterate over each cell / bam file
    for file in glob_path:

        # Load the according file
        file_name = str(file).strip().split("/")[-1]
        cell=file_name.strip().split(".bam")[0]
        bam_file = pysam.AlignmentFile(file, "rb")

        # Also get the bed_file
        norm_dictionary= defaultdict(lambda: defaultdict(tuple))
        with open(input_bed, 'r') as bed_file:
            next(bed_file) #skipping the header in bed file
            # Iterate over each manual segment aka each bed file line
            for line in bed_file:
                
                line_r = line.strip().split("\t")
                chromosome= line_r[0]
                sub_dictionary= dictionary[chromosome]
                seg_start= int(line_r[1])
                seg_end= int(line_r[2])
                interval= int(seg_end-seg_start)
                bin_size=100
                sc_TRIP_bin=100000
                origin_count=0
                mapped_count=0
                # This super nice .fetch function seems to collect all bam entries that overlap 
                # with that segment.
                for read in bam_file.fetch(chromosome, seg_start, seg_end):
                    
                    # We want to exclude reads that match any of these 5 failing criteria
                    # We use the fact that 'any' is lazy and stops as soon as it finds a true
                    # value. So e.g. c4 only has to be tested if c1, c2 and c3 all returned false.
                    c1 = 'read.is_read2'
                    c2 = 'read.is_qcfail'
                    c3 = 'read.is_secondary'
                    c4 = 'read.is_duplicate'
                    c5 = 'read.mapq < 10'
                    c6 = 'read.pos < seg_start'
                    c7 = 'read.pos >= seg_end'
                    if any([eval(c1),eval(c2),eval(c3),eval(c4),eval(c5),eval(c6),eval(c7)]):
                        pass 
                    else:
                        if read.is_reverse:
                            watson_count+=1
                        elif not read.is_reverse:
                            crick_count+=1
                        else:
                            pass
                            
            
                if (seg_start,seg_end) in norm_dictionary[chromosome]:
                    normalizing_factor_counts= norm_dictionary[chromosome][(seg_start,seg_end)][0]
                    normalizing_factor_plots= norm_dictionary[chromosome][(seg_start,seg_end)][1]
                    mapped_count= norm_dictionary[chromosome][(seg_start,seg_end)][2]
                else:
                    #bins=0
                    for m in sub_dictionary:
                        if seg_start >= m[1]:
                            continue
                        else:
                            origin_count+=sub_dictionary[m][0]
                            mapped_count+=sub_dictionary[m][1]
                            #bins+=1
                        if seg_end<=m[1]:
                            break
                        else:
                            pass
                    if (mapped_count!=0):
                        normalizing_factor_counts= interval/mapped_count
                        normalizing_factor_plots= sc_TRIP_bin/mapped_count
                        
                    else:
                        normalizing_factor_counts=0
                        normalizing_factor_plots=0
                    norm_dictionary[chromosome][(seg_start,seg_end)]= (normalizing_factor_counts, normalizing_factor_plots, mapped_count)
                #normalising both watson and crick counts to make the heights comparable
                norm_crick_counts= crick_count *normalizing_factor_counts
                norm_crick_plots= crick_count * normalizing_factor_plots
                norm_watson_counts= watson_count*normalizing_factor_counts
                norm_watson_plots= watson_count*normalizing_factor_plots
                norm_counts_file.write(str(chromosome)+ "\t"+ str(seg_start)+ "\t" + str(seg_end) + "\t" + str(sample) + "\t" + str(cell) +"\t" +str(norm_crick_counts)+ "\t" + str(norm_watson_counts)+ "\n")
                norm_plots_file.write(str(chromosome)+ "\t"+ str(seg_start)+ "\t" + str(seg_end) + "\t" + str(sample) + "\t" + str(cell) +"\t" +str(norm_crick_plots)+ "\t" + str(norm_watson_plots)+ "\n")
                
                # Reset count for next segment.
                watson_count = 0
                crick_count = 0


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-s", "--sample", help="The sample name", required=True)
    parser.add_argument("-i", "--input_bam", help="The input bam file", required=True)
    parser.add_argument("-b", "--input_bed", type=str, help="The bed file with segments",required=True)
    parser.add_argument("-n", "--norm_count_output",type=str, help="The output file with normalised watson and crick counts for downstream procesing", required=True)
    parser.add_argument("-p", "--norm_plot_output",type=str, help="The output file with normalised watson and crick counts for plots", required=True)
    parser.add_argument("-m", "--mapping_counts",type=str, help="The mapping info from simulated data", required=True)
    args = parser.parse_args()
    counts(args.sample, args.input_bam, args.input_bed, args.norm_count_output, args.mapping_counts, args.norm_plot_output)
