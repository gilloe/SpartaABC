# -*- coding: utf-8 -*-
"""
Created on Wed Sep 30 13:02:55 2020

@author: gillo
"""

import click
import os
import logging
import sys
import warnings

from validation import tree_validator, msa_validator
import summarize_results as sumres
import simulations as runs
import msa_bias_corrector as corrector
import inference as infer_abc
from configuration import general_config, simulations_config, correction_config, inference_config
# set to environment variable 'DEBUG' to 1 for full list of runtime parameters.

DEBUG = True if os.environ.get("DEBUG","0") == "0" else False


pipeline_path = os.path.dirname(os.path.abspath(__file__))
pipeline_path = pipeline_path.replace('\\','/')
if pipeline_path[-1]!='/':
	pipeline_path = pipeline_path + '/'

os.chdir(pipeline_path)



def pipeline(general_conf,simulations_conf,correction_conf,inference_conf):
	# setup the logging system
	res_dir = os.path.join(general_conf.results_path, '')
	if os.path.isdir(res_dir + "logs/"):
		print("using existing log directory")
	else:
		os.mkdir(res_dir + "logs/")
		print("creating log directory")
	log_dir = res_dir + "logs/"
	log_id = general_conf.pipeline_path.split(r"/")[-2]
	if general_conf.verbose!=2:
		if not sys.warnoptions:
			warnings.simplefilter("ignore")
	logging.basicConfig(filename=log_dir+log_id+'.log',level=logging.INFO,
					format='%(asctime)s %(name)s %(levelname)-8s %(message)s',
					datefmt='%Y-%m-%d %H:%M:%S')  
	logger = logging.getLogger(__name__)

	# user file validation:
	validator_t = tree_validator()
	validator_m = msa_validator()
	try:
		validator_t.validate_tree(res_dir, general_conf.tree_file_name)
		validator_m.validate_msa(res_dir, general_conf.msa_file_name, correction_conf.model_parameters["mode"])
	except Exception as message:
		print(message)
		return

	
	if general_conf.skip_config["sparta"]:
		#Run the SpartaABC c++ program through a wrapper python file.
		runs.create_sims_from_data(general_conf, simulations_conf) 
		if general_conf.clean_run:
			os.remove(os.path.join(res_dir, '_eq.conf'))
			os.remove(os.path.join(res_dir, '_dif.conf'))
	else:
		logger.info("Skipping Sparta.")
		#check if sparta params file exists
		for model in general_conf.available_models:
			if os.path.isfile(res_dir + f'SpartaABC_data_name_id{model}.posterior_params'):
				print("retrieved existing params files.")
				logger.info("retrieved existing params files.")
			else:
				print("Could not find SpartaABC params file.")
				logger.error("Could not find SpartaABC params file.\nPlease provide the params files or run the full pipeline")
				return
	

	if general_conf.skip_config["correct_bias"]:
		#Correct the alignment bias by realigning few simulations using MAFFT and using ML to learn the corection.
		corrector.apply_correction(general_conf, correction_conf)
	else:
		logger.info("Skipping msa bias correction.")
		#check if sparta params file exists
		for model in general_conf.available_models:
			if os.path.isfile(res_dir + f'SpartaABC_msa_corrected_id{model}.posterior_params'):
				print("retrieved corrected param files.")
				logger.info("retrieved corrected param files.")
			else:
				print("Could not find corrected params file.")
				logger.error("Could not find corrected params file.\nPlease provide the param files or run without the --skip-bc option")
				return	
	

	if general_conf.skip_config["inference"]:
		# infer abc model parameters from the corrected simulations.
		infer_abc.calc_stats(general_conf, inference_conf) # Inferring parameters from the C++ simulations
	else:
		logger.info("Skipping inference step.")
		#check if sparta params file exists
		if os.path.isfile(res_dir + 'msa_corrected_res.csv'):
			print("retrieved inference results.")
			logger.info("retrieved inference results.")
		else:
			print("Could not find inference results.")
			logger.error("Could not find inference results.\nPlease provide the param files or run without the --skip-i option")
			return	
	
	# return inference results summary.
	sumres.get_stats(general_conf,result_file_name='msa_corrected_res.csv')

	return



@click.command()
@click.option('--path', help='Path of the phylogenetic tree and MSA (output will be saved to this path too).', required=True)
@click.option('--msaf', help='MSA filename', required=True)
@click.option('--trf', help='Phylogenetic tree filename (Newick without bootstrap)', required=True)
@click.option('--ver', default=0 ,help='Verbosity (0/1/2) Default: 1')
@click.option('--minr', default=0.0 ,help='Minimal INDEL rate. Default: 0')
@click.option('--maxr', default=0.05 ,help='Maximal INDEL rate. Default: 0.05')
@click.option('--bn', default=100 ,help='epsilon num to use for ABC. Default: 100')
@click.option('--numalign', default=200 ,help='Number of alignments for MSA bias correction')
@click.option('--skip-s', default=False, is_flag=True ,help='Skip SpartaABC simulation', hidden=DEBUG)
@click.option('--skip-m', default=False, is_flag=True ,help='Skip Mafft alignment', hidden=DEBUG)
@click.option('--skip-i', default=False, is_flag=True ,help='Skip inference step', hidden=DEBUG)
@click.option('--skip-bc', default=False, is_flag=True ,help='Skip MSA bias correction', hidden=DEBUG)
@click.option('--filterp', default=(0.9,15) ,help='MSA bias correction filtering parameter. Default: 0.9 15', hidden=DEBUG)
@click.option('--nonclean', default=False, is_flag=True ,help='Do not clean files at runtime', hidden=DEBUG)
@click.option('--mode', type=click.Choice(['amino', 'nuc']), required=True, help='Specify type of alignment, proteins or DNA.')
@click.option('--submodel', type=click.Choice(['JC', 'GTR']), default='JC', help='Specify substitution model.')
@click.option('--freq', default=(0.25,0.25,0.25,0.25), help='Specify rate parameters.')
@click.option('--rates', default=(1.0,1.0,1.0,1.0,1.0), help='Specify rate parameters.')
@click.option('--inv-prop', default=0.25, help='Specify invariable sites proportion.')
@click.option('--gamma-shape', default=0.50, help='Specify shape parameter for the gamma distribution.')
@click.option('--gamma-cats', default=10, help='Specify number of categories to use in the discrete gamma approximation.')
@click.option('--nsim', default=100000, help='Specify number of simulations performed by SpartaABC.', hidden=DEBUG)
@click.option('--nburnin', default=10000, help='Specify number of burn in simulations performed by SpartaABC.', hidden=DEBUG)
def pipeline_click(path,msaf,trf,ver,minr,maxr, bn,
				   numalign,
				   skip_s, skip_m, skip_i, skip_bc,
				   filterp, nonclean, 
				   mode, submodel, freq, rates, inv_prop, gamma_shape, gamma_cats,
				   nsim, nburnin):
	# ('not' present because of previous configurations logic)
	skip_config = {
		"sparta": not skip_s and  not skip_m and not skip_bc and not skip_i,
		"mafft": not skip_m and not skip_i,
		"inference": not skip_i ,
		"correct_bias": not skip_bc and not skip_i
	}



	submodel_params_ = {
		"mode": mode,
		"submodel": submodel,
		"freq": freq,
		"rates": rates,
		"inv_prop": inv_prop,
		"gamma_shape": gamma_shape,
		"gamma_cats": gamma_cats
	}

	clean_run = not nonclean # clean run is true when nonclean is false ('not' present because of previous configurations logic)
	op_sys= 'linux'
	model_list = ['eq', 'dif'] # eq->SIM , dif->RIM
	lib = "msa_corrected"
	size_threshold = 1E6


	general_conf = general_config(pipeline_path, path , trf, msaf, model_list ,skip_config, clean_run, ver ,op_sys)
	simulations_conf = simulations_config(nsim, nburnin, minr, maxr, numalign)
	correction_conf = correction_config(submodel_params_, filterp)
	inference_conf = inference_config(lib, bn, size_threshold)

	pipeline(general_conf,simulations_conf,correction_conf,inference_conf)


if __name__ == '__main__':
	pipeline_click()

