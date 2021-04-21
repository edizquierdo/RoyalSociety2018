from typing import *

import numpy as np
from nptyping import NDArray

import matplotlib.pyplot as plt

import pandas as pd

def act(filename = 'data/run/act.dat', names : Union[str,List[str],None] = None):
	
	if isinstance(names,str):
		names = names.split(',')

	print(f'> raw names: {names}')
	
	# data_raw = np.genfromtxt(filename, delimiter = ' ', dtype = np.float).T
	data_raw = pd.read_csv(filename, sep = ' ').to_records(index=False)
	fields_raw : List[str] = list(data_raw.dtype.fields.keys())

	if names is None:
		names = fields_raw
	else:
		# pattern matching
		names_new = []
		for x in names:
			if x in fields_raw:
				names_new.append(x)
			elif x.endswith('*'):
				x_searchstr : str = x[:-1]
				for y in fields_raw:
					if y.startswith(x_searchstr):
						names_new.append(y)
		names = names_new

	# if 't' in names:
	# 	del names['t']

	print(data_raw.shape, fields_raw)

	T : NDArray = data_raw['t']
	V : Dict[str, NDArray] = {
		x : data_raw[x]
		for x in names
	}

	# plt.ylim(-50.0, 50.0)

	for v_name,v_arr in V.items():
		# print(v_name, v_arr.shape, v_arr.dtype)
		plt.plot(T, v_arr, label = v_name)

	plt.title(filename)
	plt.legend()
	plt.show()
