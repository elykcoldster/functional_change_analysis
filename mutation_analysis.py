import numpy as np
import sys
import matplotlib.pyplot as plt

if __name__=="__main__":
	orig_file = open(sys.argv[1], 'r')
	mut_file = open(sys.argv[2], 'r')
	motif_file = open(sys.argv[3], 'r')
	seq_length = int(sys.argv[4])
	
	threshold = 10

	orig = {}

	total_seqs = {}

	for line in orig_file:
		chr_pos = line.split()[0]
		score = float(line.split()[1])
		orig[chr_pos] = {'score': score, 'data': {}, 'motif': None}

	for line in motif_file:
		chr_pos = line.split()[0]
		motif_start = int(line.split()[1])
		motif_end = int(line.split()[2])
		orig[chr_pos]['motif'] = [motif_start, motif_end]

	num_muts = 0
	for line in mut_file:
		num_muts += 1
		chrinfo = line.split()[0].split('|')
		chr_pos = chrinfo[0]
		motif_center = int((orig[chr_pos]['motif'][0] + orig[chr_pos]['motif'][1]) / 2)
		
		seq_num = int(chrinfo[1]) - 1
		seq_num = seq_num - motif_center
		
		if seq_num in total_seqs:
			total_seqs[seq_num] += 1
		else:
			total_seqs[seq_num] = 1
		score = float(line.split()[1])
		"""if score > orig[chr_pos] + threshold or score < orig[chr_pos] - threshold:
			seqs[seq_num] += 1"""
		if np.sign(score) == -np.sign(orig[chr_pos]['score']):
			if seq_num in orig[chr_pos]['data']:
				orig[chr_pos]['data'][seq_num] += 1
			else:
				orig[chr_pos]['data'][seq_num] = 1

	seqs = {}
	rnd = 20
	for k, v in orig.items():
		data = v['data']
		for kd, vd in data.items():
			rkd = np.round(kd/rnd)*rnd
			if rkd in seqs:
				seqs[rkd] += vd / total_seqs[kd]
			else:
				seqs[rkd] = vd / total_seqs[kd]
	for k, v in seqs.items():
		plt.bar(k, v, width=rnd)
	plt.title('Functional Change Rate at Positions Relative to Motif Center ($N$={0})'.format(num_muts))
	plt.xlabel('Relative Nucleotide Position')
	plt.ylabel('Functional Change Rate')
	plt.show()