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
	def yaml_to_bestgen_vec(file_in : str, file_out : str = 'best.gen.dat'):
		with open(file_in, 'r') as fin:
			data = yaml.safe_load(fin.read())
		
		output = [
			float('nan'),
			# // --------------------------------
			# // Parameters for the Ventral Nerve Cord Unit
			# // --------------------------------
			# // Bias
			# phen(1) = MapSearchParameter(gen(1), -BiasRange, BiasRange);        // DB, VBa, VBp
			data['theta']['B'],
			# phen(2) = MapSearchParameter(gen(2), -BiasRange, BiasRange);        // DD, VDa, VDp
			data['theta']['D'],

			# // Time Constant
			# phen(3) = MapSearchParameter(gen(3), TauMin, TauMax);               // DB, VBa, VBp
			data['tau']['B'],
			# phen(4) = MapSearchParameter(gen(4), TauMin, TauMax);               // DD, VDa, VDp
			data['tau']['D'],

			# // Self connections
			# phen(5) = MapSearchParameter(gen(5), -SCRange, SCRange);            // DB, VBa, VBp
			data['w']['self_B'],
			# phen(6) = MapSearchParameter(gen(6), -SCRange, SCRange);            // DD, VDa, VDp
			data['w']['self_D'],

			# // Chemical synapses
			# phen(7) = MapSearchParameter(gen(7), -CSRange, CSRange);            // DB -> DD, VBa -> VDa, VBp -> VDp
			data['w']['B_to_D'],

			# phen(8) = MapSearchParameter(gen(8), -CSRange, CSRange);          // DB -> VDa, DB -> VDp, VBa -> DD /2, VBp -> DD /2
			data['w']['DB_VD'],

			# phen(9) = MapSearchParameter(gen(9), -CSRange, CSRange);          // DD -> VDa
			data['w']['DD_VD'],

			# // Gap junctions across class within unit
			# phen(10) = MapSearchParameter(gen(10), 0.0, ESRange);      // DD - VDa, DD - VDp
			data['g']['D_across'],

			# // Gap junctions per class
			# phen(11) = MapSearchParameter(gen(11), 0.0, ESRange);      // VD - VD, DD - DD
			data['g']['D_fwd'],
			# phen(12) = MapSearchParameter(gen(12), 0.0, ESRange);      // VB - VB, DB - DB
			data['g']['B_fwd'],

			# // Gap junctions across class, across neural unit
			# phen(13) = MapSearchParameter(gen(13), 0.0, ESRange);      // VB -> DB+1
			data['g']['B_across'],

			# // Stretch receptor
			# phen(14) = MapSearchParameter(gen(14), -SRmax, 0.0);        // B- class SR weight
			# REVIEW
			data['r']['B'],

			# // NMJ Weight
			# phen(15) = MapSearchParameter(gen(15), 0.0, NMJmax);       // DB, VBa, VBp
			data['w']['NMJ']['B'],
			# phen(16) = MapSearchParameter(gen(16), -NMJmax, 0.0);      // DD, VDa, VDp
			data['w']['NMJ']['D'],

			# // --------------------------------
			# // Parameters for the Head circuit
			# // --------------------------------
			# // Bias
			# phen(17) = MapSearchParameter(gen(17), -BiasRange, BiasRange);    // SMDD, SMDV
			data['theta']['SMD'],
			# phen(18) = MapSearchParameter(gen(18), -BiasRange, BiasRange);    // RMDD, RMDV
			data['theta']['RMD'],

			# // Time Constant
			# phen(19) = MapSearchParameter(gen(19), TauMin, TauMax);           // SMDD, SMDV
			data['tau']['SMD'],
			# phen(20) = MapSearchParameter(gen(20), TauMin, TauMax);           // RMDD, RMDV
			data['tau']['RMD'],

			# // Self connections
			# phen(21) = MapSearchParameter(gen(21), -SCRange, SCRange);      // SMDD, SMDV
			data['w']['self_SMD'],
			# phen(22) = MapSearchParameter(gen(22), 4.0, SCRange);           // RMDD, RMDV
			data['w']['self_RMD'],

			# // Chemical synapses
			# phen(23) = MapSearchParameter(gen(23), -HCSRange, HCSRange);      // SMDD -> SMDV, SMDV -> SMDD
			data['w']['cross_SMD'],
			# phen(24) = MapSearchParameter(gen(24), -HCSRange, HCSRange);      // SMDD -> RMDV, SMDV -> RMDD
			data['w']['SMD_to_RMD'],
			# phen(25) = MapSearchParameter(gen(25), -HCSRange, HCSRange);      // RMDD -> RMDV, RMDV -> RMDD
			data['w']['cross_RMD'],

			# // Gap junctions across class within unit
			# phen(26) = MapSearchParameter(gen(26), 0.0, ESRange);      // SMDD - RMDD, SMDV - RMDV
			data['g']['SMD_to_RMD'],
			# phen(27) = MapSearchParameter(gen(27), 0.0, ESRange);      // RMDV - RMDD
			data['g']['RMD_across'],

			# // SMD Stretch Receptor
			# phen(28) = MapSearchParameter(gen(28), -SRmax, 0.0);        // SMD- class SR weight
			data['r']['SMD'],

			# // NMJ Weight
			# phen(29) = MapSearchParameter(gen(29), 0.0, NMJmax);    // SMDD, SMDV
			data['w']['NMJ']['SMD'],
			# phen(30) = MapSearchParameter(gen(30), 0.0, NMJmax);    // RMDD, RMDV
			data['w']['NMJ']['RMD'],
		]

		with open(file_out, 'w') as fout:
			fout.write(' '.join([
				'%.7f' % x
				for x in output
			]))















if __name__ == '__main__':
	import fire
	fire.Fire(process_yaml)
