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

