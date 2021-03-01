from typing import *

import yaml


class process_yaml(object):

	@staticmethod
	def load_raw_and_comma_split(file_in : str, file_out : str):
		with open(file_in, 'r') as fin:
			data_raw = yaml.safe_load(fin.read())

		data_split = dict()

		for key,val in data_raw.items():
			multiKey = key.split(',')

			if len(multiKey) > 1:
				for k in multiKey:
					data_split[k] = val
			else:
				data_split[key] = val

		with open(file_out, 'w') as fout:
			yaml.dump(data_split, fout)


	@staticmethod
	def group_yaml(
			file_in : str, file_out : str,
			grouping : Dict[str,str] = {
				'g' : 'g',
				
			}
		):
		with open(file_in, 'r') as fin:
			data_raw = yaml.safe_load(fin.read())

		data_split = dict()

		for key,val in data_raw.items():
			multiKey = key.split(',')

			if len(multiKey) > 1:
				for k in multiKey:
					data_split[k] = val
			else:
				data_split[key] = val

		with open(file_out, 'w') as fout:
			yaml.dump(data_split, fout)















if __name__ == '__main__':
	import fire
	fire.Fire(process_yaml)
