import numpy as np
import numpy.lib.recfunctions as rfn
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.animation as animation

def arr_bounds(arr):
	return (np.amin(arr), np.amax(arr))

def read_body_data(filename):
		data_raw = np.genfromtxt(filename, delimiter = ' ', dtype = None)

		data_raw = data_raw[:,1:]

		n_tstep = data_raw.shape[0]
		n_seg = int(data_raw.shape[1] / 3)

		data = np.full(
			shape = (n_tstep, n_seg),
			fill_value = np.nan,
			dtype = np.dtype([ ('x','f8'), ('y','f8'), ('phi','f8') ]),
		)

		for s in range(n_seg):
			data[:, s]['x'] = data_raw[:, s*3]
			data[:, s]['y'] = data_raw[:, s*3 + 1]
			data[:, s]['phi'] = data_raw[:, s*3 + 2]

		return data


class Plotters(object):
	@staticmethod
	def plot_head_pos(filename : str = '../body.dat'):
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

	@staticmethod
	def plot_worm_anim(filename = '../body.dat'):
		"""
		https://towardsdatascience.com/animations-with-matplotlib-d96375c5442c
		credit to the above for info on how to use FuncAnimation
		"""
		# idk what this does tbh
		matplotlib.use("Agg")
		
		# read the data
		data = read_body_data(filename)
		# data = data[:500]
		
		# set up the figure object
		fig = plt.figure()
		ax = plt.axes(
			xlim = arr_bounds(data['x']),
			ylim = arr_bounds(data['y']),
		)
		# ax = plt.axes()

		print('> positional bounds:\t', arr_bounds(data['x']), arr_bounds(data['y']))

		
		# Set up formatting for the movie files
		Writer = animation.writers['ffmpeg']
		writer = Writer(fps=30, metadata=dict(artist='Me'), bitrate=1800)

		# this function gets called on each frame
		def anim_update(i, line):
			print(f'\t{i}\t/\t{data.shape[0]}', end = '\r')
			line.set_data(data[i]['x'], data[i]['y'])

			return line,
		
		# set up the base worm
		line, = plt.plot([], [], 'r-')

		print('> finished setup!')

		# make the animation
		line_ani = animation.FuncAnimation(
			fig, 
			anim_update,
			fargs = (line,),
			frames = data.shape[0],
			interval = 50, 
			blit=True,
		)

		print('> animation created, saving...')

		# save it
		line_ani.save('data/worm.mp4', writer = writer)

		print('\n\n> done saving!')

	@staticmethod
	def plot_act(filename = '../act.dat'):
		data_raw = np.genfromtxt(filename, delimiter = ' ', dtype = np.float).T

		print(data_raw.shape, data_raw.dtype)

		T = data_raw[0]
		V = data_raw[1:]

		print(V.shape, V.dtype)

		plt.ylim(-10.0, 50.0)

		for vv in V:
			plt.plot(T, vv)

		plt.show()


	@staticmethod
	def test_anim():
		matplotlib.use("Agg")

		def update_line(num, data, line):
			line.set_data(data[..., :num])
			return line,

		# Set up formatting for the movie files
		Writer = animation.writers['ffmpeg']
		writer = Writer(fps=15, metadata=dict(artist='Me'), bitrate=1800)


		fig1 = plt.figure()

		data = np.random.rand(2, 25)
		l, = plt.plot([], [], 'r-')
		plt.xlim(0, 1)
		plt.ylim(0, 1)
		plt.xlabel('x')
		plt.title('test')
		line_ani = animation.FuncAnimation(fig1, update_line, 25, fargs=(data, l),
										interval=50, blit=True)
		line_ani.save('lines.mp4', writer=writer)





if __name__ == '__main__':
	import fire

	fire.Fire(Plotters)


