import numpy as np
import numpy.lib.recfunctions as rfn
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.animation as animation

WORM_RADIUS = 80e-6

def arr_bounds(
		arr, 
		pad_frac : float = 0.0,
	):
	arr_min = np.amin(arr)
	arr_max = np.amax(arr)
	
	arr_range = arr_max - arr_min

	arr_min = arr_min - arr_range * pad_frac
	arr_max = arr_max + arr_range * pad_frac

	return (arr_min, arr_max)

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

def _plot_collision_boxes(ax, blocks, vecs):
	from matplotlib.patches import Rectangle
	from matplotlib.collections import PatchCollection

	print(blocks)
	print(vecs)

	plot_boxes = []

	for bl in blocks:
		plot_boxes.append(Rectangle(
			xy = bl[0], 
			width = bl[1][0] - bl[0][0], 
			height = bl[1][1] - bl[0][1],
			fill = True,
		))

	pc = PatchCollection(
		plot_boxes, 
		facecolor = 'red', 
		alpha = 0.5,
		edgecolor = 'red',
	)

	ax.add_collection(pc)

def read_coll_objs_file(objs_file : str):
	blocks = []
	vecs = []
	
	with open(objs_file, 'r') as fin:
		for row in fin:
			row_lst = row.strip().split()
			row_lst = [ float(x) for x in row_lst ]

			blocks.append([ row_lst[0:2], row_lst[2:4] ])
			vecs.append(row_lst[4:])

	return (np.array(blocks), np.array(vecs))


def body_data_split_DV(data):
	n_tstep = data.shape[0]
	n_seg = data.shape[1]

	worm_thickness = (
		WORM_RADIUS / 2.0 * abs(
			np.sin(np.arccos(
				((np.linspace(0,n_seg,n_seg)) - n_seg / 2.0) 
				/ ( n_seg / 2.0 + 0.2)
			))
		)
	)

	data_Dorsal = np.full(
		shape = (n_tstep, n_seg),
		fill_value = np.nan,
		dtype = np.dtype([ ('x','f8'), ('y','f8')]),
	)

	data_Ventral = np.full(
		shape = (n_tstep, n_seg),
		fill_value = np.nan,
		dtype = np.dtype([ ('x','f8'), ('y','f8')]),
	)

	for t in range(n_tstep):
		dX = worm_thickness * np.cos(data[t]['phi'])
		dY = worm_thickness * np.sin(data[t]['phi'])
		data_Dorsal[t]['x'] = data[t]['x'] + dX
		data_Dorsal[t]['y'] = data[t]['y'] + dY   
		data_Ventral[t]['x'] = data[t]['x'] - dX   
		data_Ventral[t]['y'] = data[t]['y'] - dY 

	return (data_Dorsal, data_Ventral)




class Plotters(object):
	@staticmethod
	def plot_head_pos(filename : str = 'data/run/body.dat'):
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
	def plot_worm_anim_old(filename = 'data/run/body.dat'):
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
	def plot_worm_anim(
			filename = 'data/run/body.dat',
			collision_objs_file = 'data/collision_objs.tsv',
			out_file = 'data/worm.mp4',
		):
		"""
		https://towardsdatascience.com/animations-with-matplotlib-d96375c5442c
		credit to the above for info on how to use FuncAnimation
		"""
		# idk what this does tbh
		matplotlib.use("Agg")
		
		# read the data
		data = read_body_data(filename)
		# data = data[:250]

		# process it
		data_D, data_V = body_data_split_DV(data)
		
		# set up the figure object
		arrbd_x = arr_bounds(data['x'])
		arrbd_y = arr_bounds(data['y'])
		
		figsize = np.array([
			arrbd_x[1] - arrbd_x[0],
			arrbd_y[1] - arrbd_y[0],
		])

		figsize = figsize * 6 / min(figsize)
		print(f'> figsize:\t{figsize}')
		fig, ax = plt.subplots(1, 1, figsize = figsize)
		
		plt.xlim(*arrbd_x)
		plt.ylim(*arrbd_y)

		# fix the scaling
		ax.axis('equal')

		# draw the blocks
		_plot_collision_boxes(ax, *read_coll_objs_file(collision_objs_file))

		print('> positional bounds:\t', arr_bounds(data['x']), arr_bounds(data['y']))

		
		# Set up formatting for the movie files
		Writer = animation.writers['ffmpeg']
		writer = Writer(fps=30, metadata=dict(artist='Me'), bitrate=1800)

		# this function gets called on each frame
		def anim_update(i, line_D, line_V):
			print(f'\t{i}\t/\t{data.shape[0]}', end = '\r')
			plt.title(f'frame\t{i}')
			line_D.set_data(data_D[i]['x'], data_D[i]['y'])
			line_V.set_data(data_V[i]['x'], data_V[i]['y'])

			return line_D,line_V
		
		# set up the base worm
		line_D, = ax.plot([], [], 'r-')
		line_V, = ax.plot([], [], 'b-')

		print('> finished setup!')

		# make the animation
		line_ani = animation.FuncAnimation(
			fig, 
			anim_update,
			fargs = (line_D, line_V),
			frames = data.shape[0],
			interval = 50, 
			blit=True,
		)

		print(f'> animation created, saving to file `{out_file}`')

		# save it
		line_ani.save(out_file, writer = writer)

		print('\n\n> done saving!')

	@staticmethod
	def plot_act(filename = 'data/run/act.dat'):
		data_raw = np.genfromtxt(filename, delimiter = ' ', dtype = np.float).T

		print(data_raw.shape, data_raw.dtype)

		T = data_raw[0]
		V = data_raw[1:]

		print(V.shape, V.dtype)

		plt.ylim(-50.0, 50.0)

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


