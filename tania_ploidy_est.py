#!/usr/bin/env python

import sys
import gzip
from argparse import ArgumentParser
import scipy.stats
import math
from scipy.stats import binom

'''
*Ploidy-estimation* 
Mixture model for ploidy estimation by Tobias Marschall MPI-Saarbruecken, implemented by Tania Christiansen @Korbel-EMBL. 
*Input-file*: Each line corresponds to a segment for which a ploidy-state is to be estimated. Each line contains the name of the segment and then the calculated W/W+C fractions per cell. 
*Script*: 
	1)Reads in each line of file, assigns the first element as segment and the following entries of that list (fractions) are bundled as a list. 
	2)For each ploidy in the range 1:max_ploidy a mixture model is applied for each segment, that is each line of document or fractions-list.
	3)The estimated ploidy is the hightest output_value. 
	4)Output is written to file as a table. Cols: chr_seg	lh-ploidy1 .... lh-ploidN	 est-ploidy
'''

class Mixture:
	def __init__(self, means, weights):
		assert len(means) == len(weights)
		self.means = means
		self.weights = weights
		self.stddevs = [0.5]*len(means)

	def fit_stddevs_meanprop(self,frac_list):
		'''Fit standard deviations so that they are proportional to the means'''
		n = 0
		v = 1.0
		while True:
			n += 1
			new_v = 0.0
			for frac_element in frac_list:
				p = self.posterior_assignment(frac_element)
				for i in range(len(p)):
					new_v += p[i] * abs(frac_element - self.means[i]) / self.weights[i]
			new_v = new_v/len(frac_list)
			diff = abs(v - new_v)
			v = new_v
			self.stddevs = [v*m for m in self.weights]
			if diff <= 1e-10: 
				break
	def posterior_assignment(self, frac_element):
		'''Get posterior distribution over classes for frac_element'''
		p = [w*scipy.stats.norm.pdf(frac_element, loc=m, scale=stddev) for m,stddev,w in zip(self.means,self.stddevs,self.weights)]
		s = sum(p)
		return [frac_element/s for frac_element in p]

	def log_likelihood(self, frac_list):
		log_likelihood = 0.0
		for frac_element in frac_list:
			p = sum(w*scipy.stats.norm.pdf(frac_element, loc=m, scale=stddev) for m,stddev,w in zip(self.means,self.stddevs,self.weights))
			log_likelihood += math.log(p)
		return log_likelihood, self.stddevs

def main():
	parser = ArgumentParser(prog='ploidy-estimator.py', description=__doc__)
	parser.add_argument('filename', metavar='COUNT', help='Gzipped, tab-separated table with counts')
	parser.add_argument('--max-ploidy', default=4, type=int,
		help='Maximum ploidy that is considered (default=4)')
	args = parser.parse_args()

	# output to terminal
	print('Processing file', args.filename, file=sys.stderr)
	# preparing output file
	col_fields = ['chr_seg'] + ['lh-ploidy{}'.format(ploidy) for ploidy in range(1,args.max_ploidy+1)]+ ['est_ploidy']
	output_colnames = ('\t'.join(col_fields))

	# 1)Reads in each line of file, assigns the first element as segment and the following entries of that list (fractions) are bundled as a list. 
	with gzip.open(args.filename, 'r') as FH_in, open('ploidy_estimation.txt','w') as FH_out:
		FH_out.write(output_colnames+'\n')
		for i, line in enumerate(FH_in):
			if i%10 == 0:
				print('Currently calculating line ', i)
			seg_line = list(x.decode() for x in line.split())
			output_list = []
			cur_seg = seg_line[0]
			fractions = [float(j) for j in list(seg_line[]1:)]
			output_list.append(cur_seg)
			# 2)For each ploidy in the range 1:max_ploidy a mixture model is applied for each segment, that is each line of document or fractions-list.
			for ploidy in range(1,args.max_ploidy+1):
				binom_dist = binom(ploidy,0.5)				
				means = [j/ploidy for j in range(ploidy+1)]
				weights = [binom_dist.pmf(i) for i in range(ploidy+1)]
				mixture = Mixture(means=means, weights=weights)
				mixture.fit_stddevs_meanprop(fractions)
				likelihood, out_stdev = mixture.log_likelihood(fractions)
				output_list.append(likelihood)
			# 3)The estimated ploidy is the hightest output_value. 
			output_list.append((np.argmax(output_list[1:])+1))
			output_line=('\t'.join([str(i) for i in output_list]))
			# 4)Output is written to file as a table. Cols: chr_seg	lh-ploidy1 .... lh-ploidN	 est-ploidy
			FH_out.write(output_line+'\n')

	FH_out.close()
	FH_in.close()
if __name__ == '__main__':
	main()




