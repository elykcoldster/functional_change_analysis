import sys

file = open(sys.argv[1], 'r')
out_file = open(sys.argv[2], 'w')

lines = []
seq = ""
for line in file:
	if ':' not in line:
		seq += line[:-1]
		if len(line[:-1]) < 50:
			lines.append(seq)
			seq=""
	else:
		if len(seq) > 0:
			lines.append(seq)
			print(seq)
			seq=""
		lines.append(line[:-1])

for l in lines:
	out_file.write(l + '\n')