# -*- coding: utf-8 -*-
"""
Created on Thu Mar 19 12:38:33 2020

@author: gillo
"""
import logging
import pandas as pd
import os
import click

#pwd = 'D:/university/projects/abc/code/abc_nn/abc_nn/'
pwd = os.path.dirname(os.path.abspath(__file__))
pwd = pipeline_path.replace('\\','/')
if pwd[-1]!='/':
    pwd = pwd + '/'

os.chdir(pwd)



logging.basicConfig(filename=pwd+'log/agg_all_res_to_single_file.log',level=logging.INFO,
                    format='%(asctime)s %(levelname)-8s %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S') 

res_path = 'results_small/infer_abc_params/eggnog/'
suffix = '.csv'
output_file_path = 'small_scripts/infer_abc_params_single_folder/df_summary.csv'     

def arg_all(res_path,output_file_path,suffix='.csv'):
    files_list = [filename for filename in os.listdir(pwd+res_path) if filename.endswith(suffix)]
    df_summary = pd.DataFrame()
    for filename in files_list:
        logging.debug(f'started {filename}')
        df_tmp = pd.read_csv(pwd+res_path+filename)
        df_tmp['filename']=filename
        df_summary = df_summary.append(df_tmp, ignore_index = True)
    df_summary.to_csv(pwd+output_file_path)
    logging.info(f'Done on res_path={res_path},suffix={suffix}, output_file_path={output_file_path}')

#arg_all(res_path=res_path,output_file_path=output_file_path,suffix=suffix)

@click.command()
@click.option('--rp', default='results_small/infer_abc_params/', help='Results to aggregate path.')
@click.option('--ofp', default='small_scripts/infer_abc_params_single_folder/df_summary.csv', help='Path for the output file.')
@click.option('--sf', default='.csv' ,help='Suffix of the files to aggregate.') 
def arg_all_click(rp,ofp,sf):
    arg_all(res_path=rp,output_file_path=ofp,suffix=sf)

if __name__ == '__main__':
    arg_all_click()