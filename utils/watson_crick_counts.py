import argparse
import sys
import pysam
from pathlib import Path
import glob
def counts(input_bam, input_bed ,counts_output):
    counts_file = open (counts_output, 'w')
    counts_file.write('chrom'+ "\t"+ 'start'+ "\t" + 'end'+ "\t" + 'sample' + "\t" + 'cell' + "\t" + 'C'+ "\t" + 'W' +"\n")
    path = Path(input_bam)
    glob_path = path.glob('*.bam')
    for file in glob_path:
        file_name = str(file).strip().split("/")[-1]
        cell=file_name.strip().split(".bam")[0]
        sample= cell.strip().split(".")[0]
        bam_file = pysam.AlignmentFile(file, "rb")
        watson_count =0
        crick_count= 0
        with open(input_bed, 'r') as bed_file:
            for line in bed_file:
                line_r = line.strip().split("\t")
                chromosome= line_r[0]
                seg_start= int(line_r[1])
                seg_end= int(line_r[2])
                for read in bam_file.fetch(chromosome, seg_start, seg_end):
                    if read.is_reverse:
                        crick_count+=1
                    elif not read.is_reverse:
                        watson_count+=1
                    else:
                        pass

                counts_file.write(str(chromosome)+ "\t"+ str(seg_start)+ "\t" + str(seg_end) + "\t" + str(sample) + "\t" + str(cell) + "\t" + str(crick_count)+ "\t" + str(watson_count)+"\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-i", "--input_bam", help="The input bam file", required=True)
    parser.add_argument("-b", "--input_bed", type=str, help="The bed file with segments",required=True)
    parser.add_argument("-c", "--counts_output", type=str, help="The output file with number of watson and reads in each segment in the bed file", required=True)
    args = parser.parse_args()
    counts(args.input_bam, args.input_bed, args.counts_output)


