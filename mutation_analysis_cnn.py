import numpy as np
import sys
import matplotlib.pyplot as plt

if __name__=="__main__":
	orig_file = open(sys.argv[1], 'r')
	mut_file = open(sys.argv[2], 'r')
	motif_file = open(sys.argv[3], 'r')
	seq_length = int(sys.argv[4])
	
	threshold = 10
	rnd = 1

	orig = {}

	total_seqs = {}

	for line in orig_file:
		chr_pos = line.split()[0][1:]
		score = float(line.split()[1])
		orig[chr_pos] = {'score': score, 'data': {}, 'motif': None}

	for line in motif_file:
		chr_pos = line.split()[0]
		motif_start = int(line.split()[1])
		motif_end = int(line.split()[2])
		if chr_pos not in orig:
			continue
		orig[chr_pos]['motif'] = [motif_start, motif_end]

	num_muts = 0
	for line in mut_file:
		num_muts += 1
		chrinfo = line.split()[0].split('|')
		chr_pos = chrinfo[0][1:]
		if chr_pos not in orig:
			continue
		motif_center = int((orig[chr_pos]['motif'][0] + orig[chr_pos]['motif'][1]) / 2)
		
		seq_num = int(chrinfo[1]) - 1
		seq_num = seq_num - motif_center
		
		tsbin = np.round(seq_num / rnd) * rnd

		if tsbin in total_seqs:
			total_seqs[tsbin] += 1
		else:
			total_seqs[tsbin] = 1
		score = float(line.split()[2])
		"""if score > orig[chr_pos] + threshold or score < orig[chr_pos] - threshold:
			seqs[seq_num] += 1"""
		if (score > 0.5 and orig[chr_pos]['score'] < 0.5) or (score < 0.5 and orig[chr_pos]['score'] > 0.5):
			if seq_num in orig[chr_pos]['data']:
				orig[chr_pos]['data'][seq_num] += 1
			else:
				orig[chr_pos]['data'][seq_num] = 1

	seqs = {}
	for k, v in orig.items():
		data = v['data']
		for kd, vd in data.items():
			rkd = np.round(kd/rnd)*rnd
			if rkd in seqs:
				seqs[rkd] += vd
			else:
				seqs[rkd] = vd
	for k, v in seqs.items():
		plt.bar(k, v/total_seqs[k], width=rnd)
	plt.title('Functional Change Rate at Positions Relative to Motif Center ($N$={0})'.format(num_muts))
	plt.xlabel('Relative Nucleotide Position')
	plt.ylabel('Functional Change Rate')
	plt.show()