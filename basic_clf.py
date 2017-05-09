import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.neural_network import MLPClassifier
from sklearn.svm import SVC
from sklearn.ensemble import AdaBoostClassifier, RandomForestClassifier

def seq2num(seq):
	num_array = []
	for bp in seq:
		if bp.lower() == 'a':
			num_array.append(1)
		elif bp.lower() == 'g':
			num_array.append(2)
		elif bp.lower() == 'c':
			num_array.append(3)
		elif bp.lower() == 't':
			num_array.append(4)
	return np.asarray(num_array)

if __name__ == '__main__':
	seq_length = 280
	pos_file = open('chip_seq_training_peaks.fa', 'r')
	neg_file = open('chip_seq_training_peaks_neg.fa', 'r')

	pos_test_file = open('chip_seq_test_peaks.fa', 'r')
	neg_test_file = open('chip_seq_test_peaks_neg.fa', 'r')

	pos_seqs = []
	neg_seqs = []
	labels = []
	for i, line in enumerate(pos_file):
		if i % 2 == 1:
			numseq = seq2num(line)
			if len(numseq) == seq_length:
				pos_seqs.append(numseq)
				labels.append(1)

	for i, line in enumerate(neg_file):
		if i % 2 == 1:
			numseq = seq2num(line)
			if len(numseq) == seq_length:
				neg_seqs.append(numseq)		
				labels.append(0)

	pos_test_seqs = []
	pos_test_labels = []
	neg_test_seqs = []
	neg_test_labels = []
	for i, line in enumerate(pos_test_file):
		if i % 2 == 1:
			numseq = seq2num(line)
			if len(numseq) == seq_length:
				pos_test_seqs.append(numseq)
				pos_test_labels.append(1)

	for i, line in enumerate(neg_test_file):
		if i % 2 == 1:
			numseq = seq2num(line)
			if len(numseq) == seq_length:
				neg_test_seqs.append(numseq)
				neg_test_labels.append(0)

	train_seqs = np.asarray(pos_seqs + neg_seqs)
	test_seqs = np.asarray(pos_test_seqs + neg_test_seqs)
	test_labels = pos_test_labels + neg_test_labels

	pos_test_seqs = np.asarray(pos_test_seqs)
	neg_test_seqs = np.asarray(neg_test_seqs)

	print('Finished loading training data. Size:', train_seqs.shape)
	
	#clf = MLPClassifier(hidden_layer_sizes=(280,), max_iter=50, verbose=True)
	clf = SVC(probability=True)
	clf.fit(train_seqs, labels)
	print('Finished training. Testing...')

	train_acc = clf.score(train_seqs, labels)
	print('Training Accuracy:', train_acc)
	test_acc = clf.score(test_seqs, test_labels)
	print('Test Accuracy:', test_acc)

	pos_predictions = clf.predict_proba(pos_test_seqs)
	neg_predictions = clf.predict_proba(neg_test_seqs)
	
	pos_out_file = open('custom_clf_predict_pos.csv', 'w')
	neg_out_file = open('custom_clf_predict_neg.csv', 'w')

	for p in pos_predictions:
		pos_out_file.write(str(p[1]) + '\n')

	for p in neg_predictions:
		neg_out_file.write(str(p[1]) + '\n')