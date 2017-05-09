import numpy as np
from seqtools import seq2num
from keras.models import Sequential
from keras.layers import Dense, Activation, Dropout, Flatten
from keras.layers.pooling import MaxPooling1D
from keras.layers.convolutional import Conv1D

if __name__ == '__main__':
	seq_length = 280
	pos_file = open('cnn_pos_training_1.fa', 'r')
	neg_file = open('cnn_neg_training.fa', 'r')

	pos_test_file = open('cnn_pos_test_1.fa', 'r')
	neg_test_file = open('cnn_neg_test.fa', 'r')

	pos_seqs = []
	neg_seqs = []
	labels = []

	print('Loading training data...')
	for i, line in enumerate(pos_file):
		if i % 2 == 1:
			numseq = seq2num(line, seq_length, flat=False, transpose=False)
			if numseq.shape[0] == seq_length:
				pos_seqs.append(numseq)
				labels.append(1)

	for i, line in enumerate(neg_file):
		if i % 2 == 1:
			numseq = seq2num(line, seq_length, flat=False, transpose=False)
			if numseq.shape[0] == seq_length:
				neg_seqs.append(numseq)		
				labels.append(0)

	train_seqs = np.asarray(pos_seqs + neg_seqs)
	print('Finished loading training data. Loading test data...')

	pos_test_seqs = []
	pos_test_labels = []
	neg_test_seqs = []
	neg_test_labels = []
	for i, line in enumerate(pos_test_file):
		if i % 2 == 1:
			numseq = seq2num(line, seq_length, flat=False)
			if numseq.shape[0] == seq_length:
				pos_test_seqs.append(numseq)
				pos_test_labels.append(1)

	for i, line in enumerate(neg_test_file):
		if i % 2 == 1:
			numseq = seq2num(line, seq_length, flat=False)
			if numseq.shape[0] == seq_length:
				neg_test_seqs.append(numseq)
				neg_test_labels.append(0)

	test_seqs = np.asarray(pos_test_seqs + neg_test_seqs)
	test_labels = pos_test_labels + neg_test_labels

	pos_test_seqs = np.asarray(pos_test_seqs)
	neg_test_seqs = np.asarray(neg_test_seqs)

	print('Finished loading test data. Building model...')

	model = Sequential()
	model.add(Conv1D(120, 32, input_shape=(seq_length, 4)))
	model.add(Activation('relu'))
	model.add(MaxPooling1D(pool_size=13, strides=13))
	model.add(Dropout(0.2))
	model.add(Flatten())
	model.add(Dense(1))
	model.add(Activation('sigmoid'))

	model.compile(loss='binary_crossentropy', optimizer='rmsprop', metrics=['accuracy'])

	print('Model compiled. Training...')
	model.fit(train_seqs, labels, epochs=25, batch_size=32)

	print('Training complete. Predicting test data...')
	predictions = model.predict(test_seqs)

	acc = 0
	pos_out_file = open('keras_cnn_pos.csv', 'w')
	neg_out_file = open('keras_cnn_neg.csv', 'w')
	for i, p in enumerate(predictions):
		if i < len(pos_test_seqs):
			pos_out_file.write(str(p[0]) + '\n')
			if p[0] > 0.5:
				acc += 1
		elif i >= len(pos_test_seqs):
			neg_out_file.write(str(p[0]) + '\n')
			if p[0] < 0.5:
				acc += 1
	acc = acc / (len(pos_test_seqs) + len(neg_test_seqs))
	print('Test accuracy:', acc)

	print('Predicting motif data...')
	motif_file = open('../jund_motif_seqs.fa', 'r')
	motif_seqs = []
	motif_chr = []
	for i, line in enumerate(motif_file):
		if i % 2 == 0:
			motif_chr.append(line[:-1])
		if i % 2 == 1:
			numseq = seq2num(line, seq_length, flat=False)
			if numseq.shape[0] == seq_length:
				motif_seqs.append(numseq)
	motif_seqs = np.asarray(motif_seqs)

	motif_predictions = model.predict(motif_seqs)
	motif_out = open('jund_motif_seqs_cnn.csv', 'w')

	motif_acc = 0
	for i,p in enumerate(motif_predictions):
		motif_out.write(motif_chr[i] + '\t' + str(p[0]) + '\n')
		if p[0] > 0.5:
			motif_acc += 1
	motif_acc /= len(motif_predictions)
	print('Motif prediction accuracy: ', motif_acc)
	
	print('Predicting mutation data...')
	mut_file = open('../jund_motif_mut_seqs.fa', 'r')
	mut_seqs = []
	mut_chr = []
	for i, line in enumerate(mut_file):
		if i % 2 == 0:
			mut_chr.append(line[:-1])
		if i % 2 == 1:
			numseq = seq2num(line, seq_length, flat=False)
			if numseq.shape[0] == seq_length:
				mut_seqs.append(numseq)
	mut_seqs = np.asarray(mut_seqs)

	mut_predictions = model.predict(mut_seqs)
	mut_out = open('jund_motif_mut_seqs_cnn.csv', 'w')

	for i,p in enumerate(mut_predictions):
		mut_out.write(mut_chr[i] + '\t' + str(p[0]) + '\n')