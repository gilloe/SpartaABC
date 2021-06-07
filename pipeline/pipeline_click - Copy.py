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


#%%


def pipeline(res_dir,msa_filename,tree_filename,pipeline_path='/bioseq/spartaabc/pipeline/',
             minIR=0,maxIR=0.05,op_sys='linux',verbose=0,
             num_epochs=200,batch_size=4096,b_num_top=100):
    
    res_dir = os.path.join(res_dir, '')
    
    import sys
    if verbose!=2:
        if not sys.warnoptions:
            import warnings
            warnings.simplefilter("ignore")
        
    
    import run_sparta_abc_single_folder_pipeline as runs
    # os.chdir(pipeline_path+"/small_scripts/infer_abc_params_single_folder")
    import infer_abc_params_single_folder_pipeline as sinf
    # os.chdir(pipeline_path)
    import summarize_results as sumres
    #print('inside pipeline')
    #click.echo('inside pipeline')
    
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
    
    sumres.get_stats_v2(results_file_path=res_dir+'data_name_res.csv', minIR=minIR, maxIR=maxIR, minAI=0, maxAI=2, msa_path=res_dir+msa_filename,verbose=verbose)
    
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

## python pipeline.py --path D:/university/projects/abc/code/pipe_line_to_github/SpartaABC/pipeline/results/example_before_run --msaf msa.f --trf tree.t --ops win
# verbose = 1
# res_dir = 'D:/university/projects/abc/code/pipe_line_to_github/SpartaABC/pipeline/results/example_before_run' #results path
# op_sys= 'win'
# data_name = '' #don't change. Old stuff.
# msa_filename = 'msa.f'
# tree_filename = 'tree.t'
# minIR = 0
# maxIR = 0.05
# verbose = 1
# num_epochs = 200
# batch_size = 4096
# b_num_top = 100

# pipeline(pipeline_path=pipeline_path,res_dir=res_dir,msa_filename=msa_filename,tree_filename=tree_filename,
#               minIR=minIR,maxIR=maxIR,op_sys=op_sys,verbose=verbose,
#               num_epochs=num_epochs,batch_size=batch_size,b_num_top=b_num_top)




# #%% click command

# python pipeline_click.py --path D:/university/projects/abc/code/pipe_line_to_github/SpartaABC/pipeline/results/example_after_run/ --msaf msa.f --trf tree.t --ops win --ver 1



@click.command()
@click.option('--path', default=pipeline_path+'results/example4/', help='Path of the phylogenetic tree and MSA (output will be saved to this path too).')
@click.option('--ops', default='linux', help='Opearating system (linux / win are supported)')
@click.option('--msaf', default='' ,help='MSA filename')
@click.option('--trf', default='r', help='Phylogenetic tree filename (Newick without bootstrap)')
@click.option('--ver', default=1 ,help='Verbosity (0/1/2) Default: 1')
@click.option('--minr', default=0.0 ,help='Minimal INDEL rate. Default: 0')
@click.option('--maxr', default=0.05 ,help='Maximal INDEL rate. Default: 0.05')
@click.option('--ne', default=200, help='NN num epochs. Default 200')
@click.option('--bs', default=4096 ,help='NN batch size. Default 4096')
@click.option('--bn', default=100 ,help='epsilon num to use for ABC. Default: 100')

def pipeline_click(path,ops,msaf,trf,ver,minr,maxr,ne,bs,bn):
    #click.echo('inside pipeline_click')
    verbose = ver
    res_dir = path #results path
    op_sys= ops
    data_name = '' #don't change. Old stuff.
    msa_filename = msaf
    tree_filename = trf
    minIR=minr
    maxIR=maxr
    verbose = ver
    num_epochs = ne
    batch_size=bs
    b_num_top=bn
    
    pipeline(pipeline_path=pipeline_path,res_dir=res_dir,msa_filename=msa_filename,tree_filename=tree_filename,
                  minIR=minIR,maxIR=maxIR,op_sys=op_sys,verbose=verbose,
                  num_epochs=num_epochs,batch_size=batch_size,b_num_top=b_num_top)

    
#%% main
        
if __name__ == '__main__':
    pipeline_click()

