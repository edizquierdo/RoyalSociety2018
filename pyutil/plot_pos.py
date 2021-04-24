import os
import sys
from typing import *

from math import degrees

import numpy as np
import numpy.lib.recfunctions as rfn
from nptyping import NDArray,StructuredType

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.patches import Patch,Rectangle,Wedge
from matplotlib.collections import PatchCollection

import pandas as pd
import json

from collision_object import (
	Path,
	CollisionType,CollisionObject,
	read_collobjs_tsv,
	BoundingBox,AxBounds,BOUNDS_TEMPLATE,
	get_bounds,get_bbox_ranges,pad_BoundingBox,
	_bounds_tuples_to_bbox,_combine_bounds,
)

# types
# ==================================================
Axes = TypeVar('Axes') # TODO: make this actually reference matplotlib.Axes
CoordsArr = np.dtype([ ('x','f8'), ('y','f8')])
CoordsRotArr = np.dtype([ ('x','f8'), ('y','f8'), ('phi','f8') ])
OptInt = Optional[int]

WORM_RADIUS = 80e-6


"""
########   #######  ##     ## ##    ## ########   ######
##     ## ##     ## ##     ## ###   ## ##     ## ##    ##
##     ## ##     ## ##     ## ####  ## ##     ## ##
########  ##     ## ##     ## ## ## ## ##     ##  ######
##     ## ##     ## ##     ## ##  #### ##     ##       ##
##     ## ##     ## ##     ## ##   ### ##     ## ##    ##
########   #######   #######  ##    ## ########   ######
"""

def arr_bounds(
		arr : NDArray, 
		pad_frac : float = 0.0,
	) -> AxBounds:
	"""return the bounds of `arr` padded by some fraction of the range
	
	[extended_summary]
	
	### Parameters:
	 - `arr : NDArray`   
	   input array
	 - `pad_frac : float`   
	   multiplied by range to determine padding
	   (defaults to `0.0`)
	
	### Returns:
	 - `AxBounds` 
	   padded bounds
	"""
	arr_min : float = np.amin(arr)
	arr_max : float = np.amax(arr)
	
	arr_range : float = arr_max - arr_min

	arr_min = arr_min - arr_range * pad_frac
	arr_max = arr_max + arr_range * pad_frac

	return (arr_min, arr_max)



def _get_fig_bounds(
		collobjs : List[CollisionObject],
		arrbd_x : AxBounds = None, 
		arrbd_y : AxBounds = None,
		figsize_scalar : float = 6.0,
	) -> NDArray[2, float]:

	collobjs_bounds : Dict[str,float] = get_bounds(collobjs)

	# set up the figure object
	if arrbd_x is None:
		# arrbd_x = arr_bounds(data['x'])
		arrbd_x = (collobjs_bounds['bound_min_x'], collobjs_bounds['bound_max_x'])
	else:
		arrbd_x = tuple(arrbd_x)

	if arrbd_y is None:
		# arrbd_y = arr_bounds(data['y'])
		arrbd_y = (collobjs_bounds['bound_min_y'], collobjs_bounds['bound_max_y'])
	else:
		arrbd_y = tuple(arrbd_y)
	
	print('> positional bounds:\t', arrbd_x, arrbd_y)

	figsize : NDArray[2, float] = np.array([
		arrbd_x[1] - arrbd_x[0],
		arrbd_y[1] - arrbd_y[0],
	])
	
	# print(f'> figsize:\t{figsize}')

	return figsize * figsize_scalar / max(figsize)


def _get_fig_bounds_box(
		bounds : BoundingBox,
		figsize_scalar : float = 6.0,
	) -> NDArray[2, float]:

	figsize : NDArray[2, float] = np.array(list(get_bbox_ranges(bounds)))

	return figsize * figsize_scalar / max(figsize)


"""
########  ########    ###    ########
##     ## ##         ## ##   ##     ##
##     ## ##        ##   ##  ##     ##
########  ######   ##     ## ##     ##
##   ##   ##       ######### ##     ##
##    ##  ##       ##     ## ##     ##
##     ## ######## ##     ## ########
"""


def read_body_data(filename : Path) -> NDArray[(Any,Any), CoordsRotArr]:
	"""reads given tsv file into a numpy array
	
	array is a 2-D structured array of type `CoordsRotArr`
	with `'x', 'y', 'phi'` fields for each segment
	so essentially 3-D, where first index is timestep, second is segment, and third/field is x/y/phi
	
	### Parameters:
	- `filename : Path`   
	filename to read
	
	### Returns:
	- `NDArray[Any, CoordsRotArr]` 
	"""
	# read in
	data_raw : NDArray = np.genfromtxt(filename, delimiter = ' ', dtype = None)

	# trim first variable (time)
	data_raw = data_raw[:,1:]

	# compute dims
	n_tstep = data_raw.shape[0]
	n_seg = int(data_raw.shape[1] / 3)

	# allocate new array
	# data : NDArray[(n_tstep, n_seg), CoordsRotArr] = np.full(
	data : NDArray[(n_tstep, n_seg)] = np.full(
		shape = (n_tstep, n_seg),
		fill_value = np.nan,
		dtype = CoordsRotArr,
	)

	# organize by x pos, y pos, and rotation (phi)
	for s in range(n_seg):
		data[:, s]['x'] = data_raw[:, s*3]
		data[:, s]['y'] = data_raw[:, s*3 + 1]
		data[:, s]['phi'] = data_raw[:, s*3 + 2]

	return data


def read_coll_objs_file(objs_file : str) -> Tuple[NDArray,NDArray]:
	"""reads an old blocks/vecs style collider file
	
	### Parameters:
	 - `objs_file : str`   
	
	### Returns:
	 - `Tuple[NDArray,NDArray]` 
	"""
	blocks : list = []
	vecs : list = []
	
	with open(objs_file, 'r') as fin:
		for row in fin:
			row_lst = row.strip().split()
			row_lst = [ float(x) for x in row_lst ]

			blocks.append([ row_lst[0:2], row_lst[2:4] ])
			vecs.append(row_lst[4:])

	return (np.array(blocks), np.array(vecs))


# 	data : NDArray[Any, CoordsRotArr],
# ) -> Tuple[NDArray[Any, CoordsArr], NDArray[Any, CoordsArr]]:
def body_data_split_DV(
		data : NDArray,
	) -> Tuple[NDArray, NDArray]:
	"""splits a body data file into arrays of dorsal and ventral points
	
	takes in a `CoordsRotArr` (produced by `read_body_data()`)
	and splits it into arrays of `CoordsArr` by the same method as `WormView.m`
	by Cohen et al [1]

	```matlab
	R = D/2.0*abs(sin(acos(((0:Nbar-1)-NSEG./2.0)./(NSEG/2.0 + 0.2))));`
	```
	
	### Parameters:
	 - `data : NDArray[Any, CoordsRotArr]`   
	   [description]
	
	### Returns:
	 - `Tuple[NDArray[Any, CoordsArr], NDArray[Any, CoordsArr]]` 
	   [description]

	### References
	 - [1] Boyle, J. H., Berri, S. & Cohen, N. Gait Modulation in C. elegans: An Integrated Neuromechanical Model. Front. Comput. Neurosci. 6, (2012).
	 	https://www.frontiersin.org/articles/10.3389/fncom.2012.00010/full

	"""
	n_tstep : int = data.shape[0]
	n_seg : int = data.shape[1]

	worm_thickness : float = (
		WORM_RADIUS / 2.0 * abs(
			np.sin(np.arccos(
				((np.linspace(0,n_seg,n_seg)) - n_seg / 2.0) 
				/ ( n_seg / 2.0 + 0.2)
			))
		)
	)

	data_Dorsal : NDArray[(n_tstep, n_seg), CoordsArr] = np.full(
		shape = (n_tstep, n_seg),
		fill_value = np.nan,
		dtype = CoordsArr,
	)

	data_Ventral : NDArray[(n_tstep, n_seg), CoordsArr] = np.full(
		shape = (n_tstep, n_seg),
		fill_value = np.nan,
		dtype = CoordsArr,
	)

	# OPTIMIZE: this bit can be vectorized
	for t in range(n_tstep):
		dX : float = worm_thickness * np.cos(data[t]['phi'])
		dY : float = worm_thickness * np.sin(data[t]['phi'])
		data_Dorsal[t]['x'] = data[t]['x'] + dX
		data_Dorsal[t]['y'] = data[t]['y'] + dY   
		data_Ventral[t]['x'] = data[t]['x'] - dX   
		data_Ventral[t]['y'] = data[t]['y'] - dY 

	return (data_Dorsal, data_Ventral)


"""
        ########  ##        #######  ########
        ##     ## ##       ##     ##    ##
        ##     ## ##       ##     ##    ##
        ########  ##       ##     ##    ##
        ##        ##       ##     ##    ##
        ##        ##       ##     ##    ##
####### ##        ########  #######     ##
"""

def _plot_collision_boxes(ax : Axes, blocks : list, vecs : list):
	"""plots old-stype collision boxes

	### Parameters:
	 - `ax : Axes`   
	 - `blocks : list`   
	 - `vecs : list`   
	"""

	print(blocks)
	print(vecs)

	plot_boxes : List[Path] = []

	for bl in blocks:
		plot_boxes.append(Rectangle(
			xy = bl[0], 
			width = bl[1][0] - bl[0][0], 
			height = bl[1][1] - bl[0][1],
			fill = True,
		))

	pc : PatchCollection = PatchCollection(
		plot_boxes, 
		facecolor = 'red', 
		alpha = 0.5,
		edgecolor = 'red',
	)

	ax.add_collection(pc)


def _plot_collobjs(ax : Axes, collobjs : Path):
	"""reads collision objects from a tsv file and plots them on `ax`
	
	### Parameters:
	 - `ax : Axes`   
	   matplotlib axes object
	 - `collobjs : Path`   
	   tsv file of collision objects (the kind where first entry is collider type)
	"""
	plot_objs : List[Patch] = []

	for obj in collobjs:
		if obj.coll_type == CollisionType.Box_Ax:
			plot_objs.append(Rectangle(
				xy = [obj['bound_min_x'], obj['bound_min_y']], 
				width = obj['bound_max_x'] - obj['bound_min_x'], 
				height = obj['bound_max_y'] - obj['bound_min_y'],
				fill = True,
			))
		elif obj.coll_type == CollisionType.Disc:
			plot_objs.append(Wedge(
				center = [ obj['centerpos_x'], obj['centerpos_y'] ],
				r = obj['radius_outer'],
				theta1 = degrees(obj['angle_min']),
				theta2 = degrees(obj['angle_max']),
				width = obj['radius_outer'] - obj['radius_inner'],

				fill = True,
			))

	pc : PatchCollection = PatchCollection(
		plot_objs, 
		facecolor = 'red', 
		alpha = 0.5,
		edgecolor = 'red',
	)

	ax.add_collection(pc)





def _plot_foodPos(ax : Axes, params : Path, fmt = 'bo'):
	with open(params, 'r') as fin:
		params_data : dict = json.load(fin)
		if "ChemoReceptor" in params_data:
			foodpos_x : float = float(params_data["ChemoReceptors"]["foodPos"]["x"])
			foodpos_y : float = float(params_data["ChemoReceptors"]["foodPos"]["y"])
	
			ax.plot(foodpos_x, foodpos_y, fmt)



"""
 ######  ######## ######## ##     ## ########
##    ## ##          ##    ##     ## ##     ##
##       ##          ##    ##     ## ##     ##
 ######  ######      ##    ##     ## ########
      ## ##          ##    ##     ## ##
##    ## ##          ##    ##     ## ##
 ######  ########    ##     #######  ##
"""


def _draw_setup(
		rootdir : Path = 'data/run/',
		bodydat : Path = 'body.dat',
		collobjs : Path = 'coll_objs.tsv',
		params : Optional[Path] = 'params.json',
		time_window : Tuple[OptInt,OptInt] = (None,None),
		figsize_scalar : Optional[float] = None,
		bounds : Optional[BoundingBox] = None,
		pad_frac : Optional[float] = None,
	) -> Tuple[
		plt.figure, Axes, 
		NDArray[(Any,Any), CoordsRotArr],
		BoundingBox,
	]:
	"""sets up figure for drawing the worm
	
	- calculates bounding boxes (with padding)
	- sets up a scaled figure object
	- plots objects and food
	
	### Parameters:
	 - `rootdir : Path`   
		prepended to all other paths
	   (defaults to `'data/run/'`)
	 - `bodydat : Path`   
	   tsv of body data
	   (defaults to `'body.dat'`)
	 - `collobjs : Path`   
	   special tsv of collision objects
	   (defaults to `'coll_objs.tsv'`)
	 - `params : Optional[Path]`   
	   json of parameters, namely food position
	   (defaults to `'params.json'`)
	 - `time_window : Tuple[OptInt,OptInt]`   
	   will only preserve data between these two timesteps (also affects bounds)
	   (defaults to `(None,None)`)
	 - `figsize_scalar : float`   
	   sclar for figure size
	   (defaults to `6.0`)
	 - `bounds : Optional[BoundingBox]`   
	   manually specify bounds. if `None`, auto-generated
	   (defaults to `None`)
	 - `pad_frac : float`   
	   if auto generating bounds, pad the box by this fraction of the range on all sides
	   (defaults to `0.0`)
	
	### Returns:
	 - `Tuple[plt.figure, Axes, NDArray[(int,int), CoordsRotArr],BoundingBox]` 
	   returns figure (fig and ax objects), data, and bounding box of worm and objects
	"""

	# setting defaults
	# ==============================
	if figsize_scalar is None:
		figsize_scalar = 6.0
	if pad_frac is None:
		pad_frac = 0.0

	# getting the data
	# ==============================

	# prepend directory to paths
	bodydat = rootdir + bodydat
	collobjs = rootdir + collobjs
	params = rootdir + params if params is not None else None

	# read worm body
	data : NDArray[(int,int), CoordsRotArr] = read_body_data(bodydat)
	print(f'> raw data stats: shape = {data.shape}, \t dtype = {data.dtype}')
	# trim
	data = data[ time_window[0] : time_window[1] ]

	# read collision objects
	lst_collision_objects : List[CollisionObject] = list()
	if os.path.isfile(collobjs):
		lst_collision_objects = read_collobjs_tsv(collobjs)
	else:
		print(f'  >> WARNING: could not find file, skipping: {collobjs}')
		

	# get bounding boxes for contents
	# ==============================

	# get bounds
	bounds_objs : BoundingBox = get_bounds(lst_collision_objects)
	bounds_worm : BoundingBox = _combine_bounds([
		_bounds_tuples_to_bbox(arr_bounds(data['x'][:,0]), arr_bounds(data['y'][:,0])),
		_bounds_tuples_to_bbox(arr_bounds(data['x'][:,-1]), arr_bounds(data['x'][:,-1])),
	])

	# if no bounds given, replace with auto generated ones
	if bounds is None:
		bounds = _combine_bounds([bounds_objs, bounds_worm])

		# pad bounds
		bounds = pad_BoundingBox(bounds, pad_frac)
	
	print('test',bounds)

	# set up figure things
	# ==============================

	# get the figure size from the bounding box
	figsize : NDArray[2, float] = _get_fig_bounds_box(bounds, figsize_scalar)

	# create the figure itself
	print(f'> figsize:\t{figsize}')
	fig, ax = plt.subplots(1, 1, figsize = figsize)

	# fix the scaling
	ax.axis('equal')

	# plot preliminaries
	# ==============================

	# plot collision objects
	if collobjs is not None:
		_plot_collobjs(ax, lst_collision_objects)
	
	# plot food position
	if params is not None:
		#  and os.path.isfile(params):
		_plot_foodPos(ax, params)

	return (
		fig,ax,
		data,
		bounds,
	)





class Plotters(object):
	"""
	##     ## ########    ###    ########
	##     ## ##         ## ##   ##     ##
	##     ## ##        ##   ##  ##     ##
	######### ######   ##     ## ##     ##
	##     ## ##       ######### ##     ##
	##     ## ##       ##     ## ##     ##
	##     ## ######## ##     ## ########
	"""
	@staticmethod
	def pos(
			# args passed down to `_draw_setup()`
			rootdir : Path = 'data/run/',
			bodydat : Path = 'body.dat',
			collobjs : Path = 'coll_objs.tsv',
			params : Optional[Path] = 'params.json',
			time_window : Tuple[OptInt,OptInt] = (None,None),
			figsize_scalar : Optional[float] = None,
			pad_frac : Optional[float] = None,
			# args specific to this plotter
			idx : int = 0,
			show : bool = True,
		):

		fig,ax,data,bounds = _draw_setup(
			rootdir = rootdir,
			bodydat = bodydat,
			collobjs = collobjs,
			params = params,
			time_window = time_window,
			figsize_scalar = figsize_scalar,
			pad_frac = figsize_scalar,
		)

		head_data : NDArray[data.shape[0], CoordsRotArr] = data[:,idx]

		print(head_data.shape, head_data.dtype)
		ax.plot(head_data['x'], head_data['y'])

		if show:
			plt.show()

	"""
	   ###    ##    ## #### ##     ##
	  ## ##   ###   ##  ##  ###   ###
	 ##   ##  ####  ##  ##  #### ####
	##     ## ## ## ##  ##  ## ### ##
	######### ##  ####  ##  ##     ##
	##     ## ##   ###  ##  ##     ##
	##     ## ##    ## #### ##     ##
	"""

	@staticmethod
	def anim(
			rootdir : Path = 'data/run/',
			bodydat : Path = 'body.dat',
			collobjs : Path = 'coll_objs.tsv',
			params : Optional[Path] = 'params.json',
			output : Path = 'worm.mp4',
			time_window : Tuple[OptInt,OptInt] = (None,None),
			figsize_scalar : float = 6.0,
			fps : int = 30, bitrate : int = 1800,
		):
		"""
		https://towardsdatascience.com/animations-with-matplotlib-d96375c5442c
		credit to the above for info on how to use FuncAnimation
		"""
		output = rootdir + output
		# idk what this does tbh
		matplotlib.use("Agg")
		
		fig,ax,data,bounds = _draw_setup(
			rootdir = rootdir,
			bodydat = bodydat,
			collobjs = collobjs,
			params = params,
			time_window = time_window,
			figsize_scalar = figsize_scalar,
			pad_frac = figsize_scalar,
		)

		data_D, data_V = body_data_split_DV(data)
		
		# Set up formatting for the movie files
		writer : animation.FFMpegWriter = animation.writers['ffmpeg'](
			fps = fps, 
			metadata = dict(artist='Me'),
			bitrate = bitrate,
		)

		# this function gets called on each frame
		def anim_update(i, line_D, line_V):
			print(f'\t{i}\t/\t{data.shape[0]}', end = '\r')
			plt.title(f'frame   {i}')
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
			blit = True,
		)

		print(f'> animation created, saving to file `{output}`')

		# save it
		line_ani.save(output, writer = writer)

		print('\n\n> done saving!')

	"""
	######## ########     ###    ##     ## ########
	##       ##     ##   ## ##   ###   ### ##
	##       ##     ##  ##   ##  #### #### ##
	######   ########  ##     ## ## ### ## ######
	##       ##   ##   ######### ##     ## ##
	##       ##    ##  ##     ## ##     ## ##
	##       ##     ## ##     ## ##     ## ########
	"""

	@staticmethod
	def single_frame(
			rootdir : Path = 'data/run/',
			bodydat : Path = 'body.dat',
			collobjs : Path = 'coll_objs.tsv',
			params : Optional[Path] = 'params.json',
			arrbd_x = None, arrbd_y = None,
			i_frame : int = 0,
			figsize_scalar : float = 10.0,
			show : bool = True,
		):
		
		fig,ax,data,bounds = _draw_setup(
			rootdir = rootdir,
			bodydat = bodydat,
			collobjs = collobjs,
			params = params,
			time_window = (i_frame, i_frame+1),
			figsize_scalar = figsize_scalar,
			pad_frac = figsize_scalar,
		)

		data_D, data_V = body_data_split_DV(data)
		
		# set up the base worm
		line_D, = ax.plot(data_D[0]['x'], data_D[0]['y'], 'r-')
		line_V, = ax.plot(data_V[0]['x'], data_V[0]['y'], 'b-')

		if show:
			plt.show()


if __name__ == '__main__':
	import fire

	fire.Fire(Plotters)


