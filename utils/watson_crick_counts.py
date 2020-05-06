import argparse
import sys
import pysam
from pathlib import Path
import glob
import pdb
def counts(sample, input_bam, input_bed ,counts_output):
    counts_file = open (counts_output, 'w')
    watson_count=0
    crick_count=0
    # Write header to output file
    counts_file.write('chrom'+ "\t"+ 'start'+ "\t" + 'end'+ "\t" + 'sample' + "\t" + 'cell' + "\t" + 'C'+ "\t" + 'W' +"\n")
    
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
        with open(input_bed, 'r') as bed_file:
            next(bed_file) #skipping the header in bed file
            # Iterate over each manual segment aka each bed file line
            for line in bed_file:
                
                line_r = line.strip().split("\t")
                chromosome= line_r[0]
                seg_start= int(line_r[1])
                seg_end= int(line_r[2])
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
                counts_file.write(str(chromosome)+ "\t"+ str(seg_start)+ "\t" + str(seg_end) + "\t" + str(sample) + "\t" + str(cell) + "\t" + str(crick_count)+ "\t" + str(watson_count)+"\n")
                
                # Reset count for next segment.
                watson_count = 0
                crick_count = 0


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-s", "--sample", help="The sample name", required=True)
    parser.add_argument("-i", "--input_bam", help="The input bam file", required=True)
    parser.add_argument("-b", "--input_bed", type=str, help="The bed file with segments",required=True)
    parser.add_argument("-c", "--counts_output", type=str, help="The output file with number of watson and reads in each segment in the bed file", required=True)
    args = parser.parse_args()
    counts(args.sample, args.input_bam, args.input_bed, args.counts_output)
