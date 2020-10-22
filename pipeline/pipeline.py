# -*- coding: utf-8 -*-
"""
Created on Wed Sep 30 13:02:55 2020

@author: gillo
"""


#%%
import click
import os
import logging
# pipeline_path = "D:/university/projects/abc/code/abc_nn/abc_nn/pipeline"
pipeline_path = os.path.dirname(os.path.abspath(__file__))
pipeline_path = pipeline_path.replace('\\','/')
if pipeline_path[-1]!='/':
    pipeline_path = pipeline_path + '/'

os.chdir(pipeline_path)
# os.chdir(pipeline_path+"/small_scripts/run_sparta_abc_single_folder")
import run_sparta_abc_single_folder_pipeline as runs
# os.chdir(pipeline_path+"/small_scripts/infer_abc_params_single_folder")
import infer_abc_params_single_folder_pipeline as sinf
# os.chdir(pipeline_path)

#%%

# ow_flag = False # over write results flag, if True, will delete other results
# verbose = 1
# cwd = pipeline_path #os.getcwd()
# data_dir = res_dir = cwd+'results/example4/' #results path
# op_sys= 'windows'#'linux'
# data_name = '' #don't change. Old stuff.
# msa_filename = 'msa.f'
# tree_filename = 'tree.t'
# minIR=0.01357
# maxIR=0.1675
# verbose = 1
# stat_eq,stat_dif = runs.create_sims_from_data(data_name=data_name,ow_flag=False,
#                       verbose=verbose,res_dir=res_dir,
#                       data_dir=data_dir,
#                       msa_filename=msa_filename,
#                       tree_filename=tree_filename,
#                       minIR=minIR,maxIR=maxIR,
#                       cwd=cwd, op_sys=op_sys) # Running spartaABC C++ code to get simulations (through another python script)
# #%%

# lib = 'data_name' #'Artropoda_ENOG41034WH_dummy'#'PF00747'   
# models_list=['ideq','iddif']
# bayes_combo_flag = False
# num_epochs = 200
# batch_size=4096
# b_num_top=100

# sinf.calc_stats(csv_out_path=res_dir,lib='data_name',path='' ,
#            verbose = verbose , models_list=['ideq','iddif'],
#            num_epochs = num_epochs,batch_size=batch_size,
#            b_num_top=b_num_top) # Inferring parameters from the C++ simulations

#%%


def pipeline(res_dir,msa_filename,tree_filename,pipeline_path='/bioseq/spartaabc/pipeline/',
             minIR=0,maxIR=0.05,op_sys='linux',verbose=0,
             num_epochs=200,batch_size=4096,b_num_top=100):
    
    # logging.basicConfig(filename=res_dir+'log_pipeline.log',
    #                 level=logging.DEBUG,
    #                 format='%(asctime)s %(levelname)-8s %(message)s',
    #                 datefmt='%Y-%m-%d %H:%M:%S')
    
    # logging.info('Started pipeline. Trying running create_sims_from_data.')
    
    stat_eq,stat_dif = runs.create_sims_from_data(data_name='',ow_flag=False,
                          verbose=verbose,res_dir=res_dir,
                          data_dir=res_dir,
                          msa_filename=msa_filename,
                          tree_filename=tree_filename,
                          minIR=minIR,maxIR=maxIR,
                          cwd=pipeline_path, op_sys=op_sys) # Running spartaABC C++ code to get simulations (through another python script)
    
    # logging.info(f'Done running create_sims_from_data with status eq:{stat_eq}, dif:{stat_dif}.')
    # logging.info('Starting running calc_stats')

    sinf.calc_stats(csv_out_path=res_dir,lib='data_name',path='' ,
               verbose = verbose , models_list=['ideq','iddif'],
               num_epochs = num_epochs,batch_size=batch_size,
               b_num_top=b_num_top) # Inferring parameters from the C++ simulations
    
    # logging.info('Done running calc_stats')

    return

#%% Debug example
#pipeline('/bioseq/data/results/spartaabc/1911/','msa.fas','tree.newick',verbose=1)

"""
#%% Example using the function without command line 
    

# ow_flag = False # over write results flag, if True, will delete other results
verbose = 1
res_dir = pipeline_path+'results/example4/' #results path
op_sys= 'linux'
data_name = '' #don't change. Old stuff.
msa_filename = 'msa.f'
tree_filename = 'tree.t'
minIR=0.01357
maxIR=0.1675
verbose = 1
num_epochs = 200
batch_size=4096
b_num_top=100

pipeline(pipeline_path=pipeline_path,res_dir=res_dir,msa_filename=msa_filename,tree_filename=tree_filename,
              minIR=minIR,maxIR=maxIR,op_sys=op_sys,verbose=verbose,
              num_epochs=num_epochs,batch_size=batch_size,b_num_top=b_num_top)
"""
# #%% click command

# ## python pipeline.py --path D:/university/projects/abc/code/abc_nn/abc_nn/pipeline/results/example4/ --msaf msa.f --trf tree.t

# @click.command()
# @click.option('--path', default=pipeline_path+'results/example4/', help='Path of the input data (output will be saved to this path too).')
# @click.option('--ops', default='linux', help='Opearating system (Linux and Windows are supported)')
# @click.option('--msaf', default='' ,help='MSA filename')
# @click.option('--trf', default='r', help='Phylogenetic tree filename')
# @click.option('--ver', default=1 ,help='Verbosity (0/1)')
# @click.option('--minr', default=0.0 ,help='Minimal INDEL rate')
# @click.option('--maxr', default=0.05 ,help='Maximal INDEL rate')
# @click.option('--ne', default=200, help='NN num epochs')
# @click.option('--bs', default=4096 ,help='NN batch size')
# @click.option('--bn', default=100 ,help='epsilon num to use for ABC')

# def pipeline_click(path,ops,msaf,trf,verb,minr,maxr,ne,bs,bn):
#     res_dir = path #results path
#     op_sys= ops#'linux'
#     data_name = '' #don't change. Old stuff.
#     msa_filename = msaf
#     tree_filename = trf
#     minIR = minr
#     maxIR = maxr
#     verbose = verb
#     num_epochs = ne
#     batch_size = bs
#     b_num_top = bn
#     pipeline(res_dir=res_dir,msa_filename=msa_filename,tree_filename=tree_filename,
#                  minIR=minIR,maxIR=maxIR,op_sys=op_sys,verbose=verbose,
#                  num_epochs=num_epochs,batch_size=batch_size,b_num_top=b_num_top)
#     return

    
# #%% main
        
# if __name__ == '__main__':
#     pipeline_click()

