# -*- coding: utf-8 -*-
"""
Created on Mon Feb 24 10:43:59 2020

@author: gillo
"""

import os
import subprocess
import shutil
import logging

import click

from Bio import AlignIO
# from multiprocessing.dummy import Pool as ThreadPool 



#%%
#'sparta_abc/sparta_abc_iddif.exe' sparta_abc/sparta_conf_file_windows_template.conf')
#exe_path  = 'sparta_abc_iddif.exe'
#conf_path = 'sparta_abc/sparta_conf_file_windows_template.conf'


def phylip2fasta(data_path,dry_run_flag=False, verbose=1):
    if dry_run_flag:
        if verbose:
            print('dry run: entered phylip2fast')
            logging.info(data_path+' dry run: entered phylip2fast')
        return
#    for subdir, dirs, files in os.walk(data_path):
#        if subdir == data_path:
#            continue
#        AlignIO.convert(os.path.join(subdir,'ref_msa.aa.phy'), 'phylip-relaxed', os.path.join(subdir,'ref_msa.aa.fasta'),'fasta')
    AlignIO.convert(os.path.join(data_path,'ref_msa.aa.phy'), 'phylip-relaxed', os.path.join(data_path,'ref_msa.aa.fasta'),'fasta')
    logging.info(data_path + " Done convert phylip to fasta")
    if verbose:
        print("Done convert phylip to fasta")
    return

def run_sparta_abc(exe_path,conf_path):
    args = exe_path + ' ' + conf_path
    tmp = subprocess.call(args,shell=True)
    return tmp


# '_minIRVal 0'
# '_maxIRVal 0.05'
# minIR = 0
# maxIR = 0.05
# f'_minIRVal {round(minIR,2)}'
# f'_maxIRVal {round(maxIR,2)}'


def create_sims(data_name,verbose=1,res_dir='results/',
                data_dir='data',msa_filename='ref_msa.aa.fasta',
                tree_filename='RAxML_tree.tree',
                minIR = 0, maxIR = 0.05,
                cwd='/groups/pupko/gilloe/spartaABC/code/abc_nn/',
                op_sys='linux'):
    logging.info('Inside create_sims function.')
    # cwd = '/groups/pupko/gilloe/spartaABC/code/abc_nn/' #os.getcwd()
    # out_dir = '/results/' + data_name
    out_dir = res_dir + data_name
#    data_path = '\\data\\' + data_name
    # os.mkdir(cwd + out_dir) # In the pipeline version, assuming it exists
    #copy config file to this path
#    shutil.copyfile(src=cwd+'/sparta_abc/sparta_conf_file_windows_template.conf',dst=cwd+out_dir)
#    args_cp = 'copy ' + cwd +'\sparta_abc\sparta_conf_file_windows_template.conf ' + cwd+'\\results\\' + data_name
#    tmp = subprocess.call(args_cp,shell=True)
    # if not os.path.isfile(cwd+data_dir+data_name+'/ref_msa.aa.fasta'):
    #     phylip2fasta(data_path=cwd+data_dir+data_name,dry_run_flag=False, verbose=verbose)    
    logging.info('Writing conf files.')
    with open(cwd+'sparta_abc/sparta_conf_file_windows_template.conf', "rt") as fin:
        with open(out_dir+'_ideq.conf', "wt") as fout:
            for line in fin:
                if line.startswith('_outputGoodParamsFile'):
                    line_out = f'_outputGoodParamsFile {out_dir}SpartaABC_data_name_ideq.posterior_params\n'
                else:
                    line_out = line.replace('results/', res_dir).replace('data/', data_dir).replace('data_name/', data_name).replace('model_name','ideq').replace('ref_msa.aa.fasta',msa_filename).replace('RAxML_tree.tree',tree_filename).replace('_minIRVal 0',f'_minIRVal {round(minIR,2)}').replace('_maxIRVal 0.05',f'_maxIRVal {round(maxIR,2)}')
                fout.write(line_out)
                # print(line)
                # print(line.replace('results/', res_dir).replace('data/', data_dir).replace('data_name/', data_name).replace('model_name','ideq').replace('ref_msa.aa.fasta',msa_filename).replace('RAxML_tree.tree',tree_filename).replace('_minIRVal 0',f'_minIRVal {round(minIR,2)}').replace('_maxIRVal 0.05',f'_maxIRVal {round(maxIR,2)}'))
    logging.info(f"Wrote {out_dir+'_ideq.conf'}. Based on {cwd+'sparta_abc/sparta_conf_file_windows_template.conf'}")
    with open(cwd+'sparta_abc/sparta_conf_file_windows_template.conf', "rt") as fin:
        with open(out_dir+'_iddif.conf', "wt") as fout:
            for line in fin:
                if line.startswith('_outputGoodParamsFile'):
                    line_out = f'_outputGoodParamsFile {out_dir}SpartaABC_data_name_iddif.posterior_params\n'
                else:
                    line_out = line.replace('results/', res_dir).replace('data/', data_dir).replace('data_name/', data_name).replace('model_name','iddif').replace('ref_msa.aa.fasta',msa_filename).replace('RAxML_tree.tree',tree_filename).replace('_minIRVal 0',f'_minIRVal {round(minIR,2)}').replace('_maxIRVal 0.05',f'_maxIRVal {round(maxIR,2)}')
                fout.write(line_out)
    logging.info(f"Wrote {out_dir+'_iddif.conf'}. Based on {cwd+'sparta_abc/sparta_conf_file_windows_template.conf'}")
            
    stat_eq = stat_dif = "didn't run"
    if op_sys=='linux':
        logging.info(f"linux op system")
        stat_eq = run_sparta_abc(exe_path=cwd+'sparta_abc_ideq', conf_path=res_dir+data_name+data_name+'_ideq.conf')
        logging.info(f"ran {cwd+'sparta_abc_ideq'}, conf_path={res_dir+data_name+data_name+'_ideq.conf'}")
        stat_dif = run_sparta_abc(exe_path=cwd+'sparta_abc_iddif', conf_path=res_dir+data_name+data_name+'_iddif.conf')
        logging.info(f"ran {cwd+'sparta_abc_iddif'}, conf_path={res_dir+data_name+data_name+'_iddif.conf'}")
    else: #windows
        logging.info(f"windows op system")
        stat_eq = run_sparta_abc(exe_path=cwd+'sparta_abc_ideq.exe', conf_path=res_dir+data_name+data_name+'_ideq.conf')
        logging.info(f"ran {cwd+'sparta_abc_ideq.exe'}, conf_path={res_dir+data_name+data_name+'_ideq.conf'}")
        stat_dif = run_sparta_abc(exe_path=cwd+'sparta_abc_iddif.exe', conf_path=res_dir+data_name+data_name+'_iddif.conf')
        logging.info(f"ran {cwd+'sparta_abc_iddif.exe'}, conf_path={res_dir+data_name+data_name+'_iddif.conf'}")
    logging.info(f'{data_name} simulations are done - finish stats {stat_eq},{stat_dif}.')
    if verbose:
        print(f'{data_name} simulations are done - finish stats {stat_eq},{stat_dif}.')
    return (stat_eq,stat_dif)



#%% Editing confing file

# ow_flag = False # over write results flag, if True, will delete other results
# verbose = 1
# cwd = '/groups/pupko/gilloe/spartaABC/code/abc_nn' #os.getcwd()
# data_name = 'brassicales_ENOG410BSBN' #'EMGT00050000000002'#'EMGT00050000000053'
#data_path = '\\data\\' + data_name
#create in results file lib with data name
#out_dir = '\\results\\' + data_name

#%%

def create_sims_from_data(data_name,ow_flag=False,verbose=1,
                          res_dir='results',data_dir='data',
                          msa_filename='ref_msa.aa.fasta',
                          tree_filename='RAxML_tree.tree',
                          minIR = 0, maxIR = 0.05,
                          cwd='/groups/pupko/gilloe/spartaABC/code/abc_nn/',
                          op_sys='linux'):
    ##%% Configuring log
    if not res_dir.endswith('/'):
        res_dir += '/'
    if not data_dir.endswith('/'):
        data_dir += '/'
    logging.basicConfig(filename=res_dir+'log_run_spartaabc_on_single_folder.log',
                        level=logging.DEBUG,
                        format='%(asctime)s %(levelname)-8s %(message)s',
                        datefmt='%Y-%m-%d %H:%M:%S')
    logging.info(data_name +' started.')
    stat_eq = stat_dif = "didn't run"
    if verbose:
        print(data_name +' started.')
    try:
        stat_eq,stat_dif = create_sims(data_name=data_name,verbose=verbose,
                    res_dir=res_dir,data_dir=data_dir,
                    msa_filename=msa_filename,
                    tree_filename=tree_filename,
                    minIR=minIR,maxIR=maxIR,
                    cwd=cwd,op_sys=op_sys)
    except:
        if ow_flag:
            # shutil.rmtree(os.getcwd() + '/results/' + data_name)
            shutil.rmtree(os.getcwd() + '/'+res_dir + data_name)
            stat_eq,stat_dif = create_sims(data_name=data_name,verbose=verbose,
                        res_dir=res_dir,data_dir=data_dir,
                        msa_filename=msa_filename,
                        tree_filename=tree_filename,
                        minIR=minIR,maxIR=maxIR,
                        cwd=cwd,op_sys=op_sys)
        else:
            logging.info(data_name +' skipped.')
            if verbose:
                print(data_name +' skipped.')
    return (stat_eq,stat_dif)

#%%

# data_name = 'PF00747'   
# create_sims_from_data(data_name=data_name,ow_flag=False,verbose=1)

# ow_flag = False # over write results flag, if True, will delete other results
# verbose = 1
# cwd = 'D:/university/projects/abc/code/abc_nn/abc_nn/pipeline/' #os.getcwd()
# data_dir = res_dir = cwd+'results/example4/'
# op_sys= 'windows'#'linux'
# data_name = ''
# msa_filename = 'msa.f'
# tree_filename = 'tree.t'
# minIR=0.01357
# maxIR=0.1675
# verbose = 1
# stat_eq,stat_dif = create_sims_from_data(data_name=data_name,ow_flag=False,
#                       verbose=verbose,res_dir=res_dir,
#                       data_dir=data_dir,
#                       msa_filename=msa_filename,
#                       tree_filename=tree_filename,
#                       minIR=minIR,maxIR=maxIR,
#                       cwd=cwd, op_sys=op_sys)

#%%

# python small_scripts\run_sparta_abc_single_folder\run_sparta_abc_single_folder_pipeline.py --dird D:/university/projects/abc/code/abc_nn/abc_nn/pipeline/results/example2/ --msf msa.f --trf tree.t --ops windows --verb 1 --minr 0.02 --maxr 0.235


# dird = 'D:/university/projects/abc/code/abc_nn/abc_nn/pipeline/results/example2/' 
# msf = 'msa.f' 
# trf = 'tree.t' 
# ops = 'windows' 
# verb = 1

@click.command()
# @click.option('--owf', default=False, help='Overwrite flag')
# @click.option('--dirr', default='results', help='Results directory path.')
@click.option('--dird', default='data' ,help='Data (full) path')
@click.option('--msf', default='ref_msa.aa.fasta' ,help='MSA file name')
@click.option('--minr', default=0.0 ,help='Minimum indel rate')
@click.option('--maxr', default=0.05 ,help='Maximum indel rate')
@click.option('--trf', default='RAxML_tree.tree' ,help='Rooted phylogenetic tree file name')
@click.option('--ops', default='linux' ,help='Operation system: linux or windows')
@click.option('--verb', default=1, help='Verbosity (0/1)')

# @click.option('--count', default=1, help='number of greetings')
 
def create_sims_from_data_click(dird,msf,minr,maxr,trf,ops,verb):
    ow_flag = False # over write results flag, if True, will delete other results
    try:
        verbose = int(verb)
    except:
        verbose = 1
    data_dir = res_dir = dird
    cwd = 'D:/university/projects/abc/code/abc_nn/abc_nn/pipeline/' #pipeline location path
    op_sys= ops#'linux'
    data_name = ''
    msa_filename = msf
    tree_filename = trf
    stat_eq,stat_dif = create_sims_from_data(data_name=data_name,ow_flag=False,
                      verbose=verbose,res_dir=res_dir,
                      data_dir=data_dir,
                      msa_filename=msa_filename,
                      tree_filename=tree_filename,
                      minIR=minr,maxIR=maxr,
                      cwd=cwd, op_sys=op_sys)

# @click.command()
# @click.option('--dirp', default='data', help='Data directory name.')
# @click.option('--dirr', default='results', help='Results directory path.')
# @click.option('--dird', default='data' ,help='Data directory parent path')

# # @click.option('--count', default=1, help='number of greetings')
 
# def create_sims_from_data_click(dirp,dirr,dird):
#     create_sims_from_data(data_name=dirp,ow_flag=False,verbose=1,
#                           res_dir=dirr,data_dir=dird)

if __name__ == '__main__':
    create_sims_from_data_click()
