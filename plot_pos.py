# import numpy as np
import matplotlib.pyplot as plt

def plot_head_pos(filename = 'body.dat'):
	head_x = []
	head_y = []
	
	with open(filename, 'r') as fin:
		for line in fin:
			xy_temp = line.split()[1:3]
			head_x.append(float(xy_temp[0]))
			head_y.append(float(xy_temp[1]))

	print(len(head_x), len(head_y))

	plt.plot(head_x, head_y)
	plt.show()


if __name__ == '__main__':
	plot_head_pos()


