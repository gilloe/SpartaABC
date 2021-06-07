# -*- coding: utf-8 -*-
"""
Created on Thu Jan 28 15:21:57 2021

@author: gillo
"""

import os
import pandas as pd

def remove_files(res_path,to_remove):
	for i in to_remove:
		os.remove(os.path.join(res_path, i))


def get_stats(general_conf, result_file_name):
	'''
	Organize the results of the inference step in a summary file
	'''
	res_dir = general_conf.results_path

	chosen_model_field_name='bayes_class'
	rl_template='m_@@@_RL'
	ir_template='m_@@@_IR'
	ai_template='m_@@@_AIR'
	dr_template='m_@@@_DR'
	ad_template='m_@@@_ADR'


	df = pd.read_csv(os.path.join(res_dir, result_file_name))

	if general_conf.clean_run:
		remove_files(res_dir, to_remove=[
			result_file_name,
			"SpartaABC_msa_corrected_iddif.posterior_params",
			"SpartaABC_msa_corrected_ideq.posterior_params"
		])

	
	chosen_model = df[chosen_model_field_name].values[0]  # either 'eq' for SIM or 'dif' for RIM
	stats = {'chosen model': chosen_model.replace('eq','SIM').replace('dif','RIM')}  # add model name
	
	# prepare SIM summary

	rl_field_name = rl_template.replace('@@@', 'eq')
	ir_field_name = ir_template.replace('@@@', 'eq')
	ai_field_name = ai_template.replace('@@@', 'eq')
	

	stats.update({'SIM RL': df[rl_field_name].values[0],
				  'SIM R_ID': df[ir_field_name].values[0]*2,
				  'SIM A_ID': df[ai_field_name].values[0],
				  })
	stats_eq = [('RL', df[rl_field_name].values[0]),
				('R_ID', df[ir_field_name].values[0]*2),
				('A_ID', df[ai_field_name].values[0])]


	# prepare RIM summary
	rl_field_name = rl_template.replace('@@@', 'dif')
	ir_field_name = ir_template.replace('@@@', 'dif')
	dr_field_name = dr_template.replace('@@@', 'dif')
	ai_field_name = ai_template.replace('@@@', 'dif')
	ad_field_name = ad_template.replace('@@@', 'dif')
	stats.update({'RIM RL': df[rl_field_name].values[0],
			  'RIM R_I': df[ir_field_name].values[0],
			  'RIM R_D': df[dr_field_name].values[0],
			  'RIM A_I': df[ai_field_name].values[0],
			  'RIM A_D': df[ad_field_name].values[0],
			  })
	stats_dif = [('RL', df[rl_field_name].values[0]),
				('R_I', df[ir_field_name].values[0]),
				('R_D', df[dr_field_name].values[0]),
				('A_I', df[ai_field_name].values[0]),
				('A_D', df[ad_field_name].values[0])]

	# if verbose!=0:
	# os.system('cls||clear')
	print("\033[1m" + "Results:" + "\033[0m")
	print(f"Chosen model is: \033[1m{stats['chosen model']}\033[0m")
	print('+------------------------+	+------------------------+')
	print('| RIM                    |	| SIM                    |')
	print('+========================+	+========================+')
	print(f'| {stats_dif[0][0]}  |     {stats_dif[0][1]:.2f}        |	| {stats_eq[0][0]}   |     {stats_eq[0][1]:.2f}       |')
	print('+-----+------------------+	+------+-----------------+')
	print(f'| {stats_dif[1][0]} |     {stats_dif[1][1]:.5f}      |	| {stats_eq[1][0]} |     {stats_eq[1][1]:.5f}     |')
	print('+-----+------------------+	+------+-----------------+')
	print(f'| {stats_dif[2][0]} |     {stats_dif[2][1]:.5f}      |	| {stats_eq[2][0]} |     {stats_eq[2][1]:.5f}     |')
	print('+-----+------------------+	+------+-----------------+')
	print(f'| {stats_dif[3][0]} |     {stats_dif[3][1]:.5f}      |')
	print('+-----+------------------+')
	print(f'| {stats_dif[4][0]} |     {stats_dif[4][1]:.5f}      |')
	print('+-----+------------------+')


	pd.DataFrame(stats.items()).sort_values(by=0).to_csv(os.path.join(res_dir, 'sum_res.csv'), header=False, index=False)
	return stats


