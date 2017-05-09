import numpy as np

def get_gc_content(line):
	gc = 0
	for c in line:
		if c.lower() == 'c' or c.lower() == 'g':
			gc += 1
	return gc / len(line)

fa_file = open('F:/data/chip_seq_peaks.fa', 'r')

peak_length = 280

chr_lengths = {'chr1':249250621,
			'chr2':243199373,
			'chr3':198022430,
			'chr4':191154276,
			'chr5':180915260,
			'chr6':171115067,
			'chr7':159138663,
			'chrX':155270560,
			'chr8':146364022,
			'chr9':141213431,
			'chr10':135534747,
			'chr11':135006516,
			'chr12':133851895,
			'chr13':115169878,
			'chr14':107349540,
			'chr15':102531392,
			'chr16':90354753,
			'chr17':81195210,
			'chr18':78077248,
			'chr20':63025520,
			'chrY':59373566,
			'chr19':59128983,
			'chr22':51304566,
			'chr21':48129895}

out_file = open('F:/data/chip_seq_peaks_negative.fa', 'w')

ln1 = ""
randindex = 0
for i, line in enumerate(fa_file):
	if i % 2 == 0:
		chrm = line.split(':')[0][1:]
		length = chr_lengths[chrm]
		randindex = np.random.randint(1,length - peak_length)

		ln1 = '>' + chrm + ':' + str(randindex) + '-'
	else:
		line_len = len(line) - 1
		gc = get_gc_content(line)

		ln1 += str(randindex + line_len) + '\n'
		ln2 = ""
		for i in range(0, line_len):
			if np.random.rand() < gc:
				if np.random.rand() < 0.5:
					ln2 += 'G'
				else:
					ln2 += 'C'
			else:
				if np.random.rand() < 0.5:
					ln2 += 'A'
				else:
					ln2 += 'T'
		ln2 += '\n'
		out_file.write(ln1)
		out_file.write(ln2)