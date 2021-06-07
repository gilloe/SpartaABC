# -*- coding: utf-8 -*-
"""
Created on Mon Feb 24 10:43:59 2020

@author: gillo
"""

import os
import subprocess
import logging
from  configuration import get_sparta_config

logger = logging.getLogger(__name__)


def run_sparta_abc(exe_path,conf_path):
	args = exe_path + ' ' + conf_path
	tmp = subprocess.call(args,shell=True)
	return tmp



def create_sims(general_conf, simulations_conf):
	'''
	simulate in spartaABC according to configuration.
	'''

	logger.debug('Inside create_sims function.')
	res_dir = general_conf.results_path
	cwd = general_conf.pipeline_path

	logger.info('Writing conf files.')
	for model in general_conf.available_models:
		sparta_config = get_sparta_config()
		# edit relevant config parameters
		sparta_config["_numSimulations"] = str(simulations_conf.simulations_num)
		sparta_config["_numBurnIn"] = str(simulations_conf.burnin_num)

		sparta_config['_outputGoodParamsFile'] = os.path.join(res_dir,f"SpartaABC_data_name_id{model}.posterior_params")
		sparta_config['_outputAlignmnetsFile'] = os.path.join(res_dir,f"alignments_{model}.fasta")
		sparta_config['_alignments_output'] = str(simulations_conf.alignments_num)
		sparta_config["_inputRealMSAFile"] = os.path.join(res_dir,general_conf.msa_file_name)
		sparta_config["_inputTreeFileName"] = os.path.join(res_dir,general_conf.tree_file_name)
		sparta_config["_minIRVal"] = str(round(simulations_conf.min_ir,2))
		sparta_config["_maxIRVal"] = str(round(simulations_conf.max_ir,2))
		sparta_config["_modelType"] = model

		# write sparta config file.
		with open(os.path.join(res_dir,f'_{model}.conf'), "wt") as fout:
			for key in sparta_config:
				to_write = f'{key} {sparta_config[key]}\n'
				fout.write(to_write)

		logger.info(f"Wrote {res_dir}_{model}.conf")
			
	stat = "didn't run"
	# execute spartaABC c++ program for all models(SIM,RIM).
	for model in general_conf.available_models:
		stat = run_sparta_abc(exe_path=os.path.join(cwd,'SpartaABC'), conf_path=os.path.join(res_dir,f'_{model}.conf'))
		logger.info(f"ran {os.path.join(cwd,'SpartaABC')},\n conf_path={os.path.join(res_dir,f'_{model}.conf')},\n stat={stat}")

	logger.info(f'simulations are done - finish stats {stat}.')
	if general_conf.verbose:
		print(f'simulations are done - finish stats {stat}.')
	return stat



def create_sims_from_data(general_conf, simulations_conf):
	'''
	wrapper function to initialize simulations
	'''
	# res_dir = general_conf.results_path
	# if not res_dir.endswith('/'):
	# 	res_dir += '/'
	# if not data_dir.endswith('/'):
	# 	data_dir += '/'

	logger.info('simulations started.')
	if general_conf.verbose:
		print('simulations started.')
	create_sims(general_conf, simulations_conf)
