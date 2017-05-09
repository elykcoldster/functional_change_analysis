import numpy as np

def seq2num(seq, slength=280, flat=True, transpose=False):
	num_array = []
	for i,bp in enumerate(seq):
		if bp.lower() == 'a':
			num_array.append(np.array([1,0,0,0]))
		elif bp.lower() == 'g':
			num_array.append(np.array([0,1,0,0]))
		elif bp.lower() == 'c':
			num_array.append(np.array([0,0,1,0]))
		elif bp.lower() == 't':
			num_array.append(np.array([0,0,0,1]))
	num_array = np.asarray(num_array)
	if flat == True:
		num_array = num_array.flatten()
	if transpose == True:
		num_array = np.transpose(num_array)
	return num_array
