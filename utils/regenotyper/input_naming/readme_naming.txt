Whoeps, 17th Feb 2021

There are three kinds of names that each Sample can have. 

- a) "StrandSeq long": e.g. HVWVWAFXY_GM12329Ax02_19s004318-1-1
- b) "StrandSeq short: e.g. GM12329A
- c) name in rna files or in snp freeze (sometimes called "RNA name"): e.g. NA12329
	The difference between b and c is usually just replacing "GM" with "NA" in some samples


# FILES # (Oh wow i only realize now how terrible and confusing this is. 
# TODO improve xD the first three files could easily be one file.)

namimg_sseq_to_rna: we translate between naming a) and c)

samplenames_sseq_to_rna: translate b) and c)

samplenames_sseq_to_sseq_short.txt: translate between a) and b)


samples_pav: samples for which we have snp data in the PAV Freeze4 .vcf file. 

samples_pg: samples (in form "RNA name, (c)") which are missing in the PAV freeze4 vcf and therefore we 
			take from the pangenie vcf which has everything. But currently
			there are no such samples anymore at least for me. This was a bigger
			problem in PAV freeze3.


# THINGS THAT HAPPEN IN THE SNAKEMAKE #

In general, we usually stick with the long StrandSeq names (a). We use (b) names to extract SNP info from the freeze4 vcf 

[line 28] SAMPLES = glob_wildcards("../../sv_probabilities/{sample}/100000_fixed_norm.selected_j0.1_s0.1/probabilities.Rdata")]:
	we grab SAMPLES, whith samplenames in long a) version.

[line 61]     SAMPLEDICT[row['sseq']] = row['pav']
	we make a dictionary that can translate from a) to c)

[line 67] print('Samples considered:') print(SAMPLES)
	this should print the StrandSeq long names of the samples.

[line 177] if PHASE: \n rule rephase_all_txt:
	This is the R script that does the rephasing for us. For reasons of which i am not entirely sure, it likes the sseq_short names. I have a suspicion that this should be rna names, but i need to investigate that more. For now let's assume it does indeed use sseq_short.   


[line287] rule get_pav_vcf_gz:
	here we pull out the samples of interest from the freeze4 vcf. this is where we need the b) names, because these are the names that the vcf uses. So here, the dict from [line 61] is used. 
