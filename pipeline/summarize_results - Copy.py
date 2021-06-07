# -*- coding: utf-8 -*-
"""
Created on Thu Jan 28 15:21:57 2021

@author: gillo
"""

import sys
import os
import shutil
import Bio.SeqUtils
import pandas as pd


def get_RL_bounds(msa_path, min_normalized_factor=0.8, max_normalized_factor=1.1):
    minIR = float('inf')
    maxIR = -float('inf')

    with open(msa_path) as f:
        f.readline()  # skip first header
        sequence = ''
        for line in f:
            line = line.rstrip()
            if not line.startswith('>'):
                sequence += line.rstrip()  # accumulate sequence
            else:
                # a new header was found
                seq_length = len(sequence.replace('-', ''))
                sequence = ''
                if seq_length < minIR:
                    minIR = seq_length
                if seq_length > maxIR:
                    maxIR = seq_length

    # don't forget last record!!
    if sequence != '':
        seq_length = len(sequence.replace('-', ''))
        if seq_length < minIR:
            minIR = seq_length
        if seq_length > maxIR:
            maxIR = seq_length

    return minIR * min_normalized_factor, maxIR * max_normalized_factor


def verify_results_validity(stats, df, ai_field_name, ir_field_name, rl_field_name,
                            dr_field_name, ad_field_name, minIR, maxIR,
                            minAI, maxAI, msa_path):

    is_lasso = True

    # logger.info('Calculating legal RL range...')
    minRL, maxRL = get_RL_bounds(msa_path)
    # logger.info(f'minRL={minRL}, maxRL={maxRL}')

    if (not minIR <= stats['IR'] <= maxIR) \
        or (not minAI <= stats['AI'] <= maxAI) \
        or (not minRL <= stats['RL'] <= maxRL)\
        or (stats['model'] == 'dif'
            and (not minIR <= stats['DR'] <= maxIR
                 or not minAI <= stats['AD'] <= maxAI)):
        is_lasso = False
        ir_field_name = 'm' + ir_field_name[1:]  # l_eq_IR -> m_eq_IR / l_dif_IR -> m_dif_IR
        stats['IR'] = df[ir_field_name].values[0]

        ai_field_name = 'm' + ai_field_name[1:]  # l_eq_AIR -> m_eq_AIR / l_dif_AIR -> m_dif_AIR
        stats['AI'] = df[ai_field_name].values[0]

        rl_field_name = 'm' + rl_field_name[1:]  # l_eq_RL -> m_eq_RL
        stats['RL'] = df[rl_field_name].values[0]

        if stats['model'] == 'dif':
            dr_field_name = 'm' + dr_field_name[1:]  # l_dif_DR -> m_dif_DR
            stats['DR'] = df[dr_field_name].values[0]

            ad_field_name = 'm' + ad_field_name[1:]  # l_dif_AD -> m_dif_AD
            stats['AD'] = df[ad_field_name].values[0]

    return is_lasso


# def get_stats(results_file_path, minIR, maxIR, minAI, maxAI, msa_path, chosen_model_field_name='nn_c_class', rl_template='l_@@@_RL',
#               ir_template='l_@@@_IR', ai_template='l_@@@_AIR', dr='l_@@@_DR', ad='l_@@@_ADR'):

#     df = pd.read_csv(results_file_path)

#     chosen_model = df[chosen_model_field_name].values[0]  # either 'eq' for SIM or 'dif' for RIM
#     stats = {'model': chosen_model}  # add model name

#     rl_field_name = rl_template.replace('@@@', df[chosen_model_field_name].values[0])
#     ir_field_name = ir_template.replace('@@@', df[chosen_model_field_name].values[0])
#     ai_field_name = ai_template.replace('@@@', df[chosen_model_field_name].values[0])
#     stats.update({'RL': df[rl_field_name].values[0],
#                   'IR': df[ir_field_name].values[0],
#                   'AI': df[ai_field_name].values[0]})

#     dr_field_name = ad_field_name = None
#     if chosen_model == 'dif':  # RIM model
#         dr_field_name = dr.replace('@@@', df[chosen_model_field_name].values[0])
#         ad_field_name = ad.replace('@@@', df[chosen_model_field_name].values[0])
#         stats.update({'DR': df[dr_field_name].values[0],
#                       'AD': df[ad_field_name].values[0]})

#     logger.info(f'Results BEFORE verification:\n{stats}')
#     is_lasso = verify_results_validity(stats, df, ai_field_name, ir_field_name, rl_field_name,
#                                     dr_field_name, ad_field_name, minIR, maxIR, minAI, maxAI, msa_path)
#     logger.info(f'Results AFTER verification:\n{stats}')

#     return stats, is_lasso

def get_stats_v2(results_file_path, minIR, maxIR, minAI, maxAI, msa_path, verbose=1, chosen_model_field_name='nn_c_class', rl_template='l_@@@_RL',
              ir_template='l_@@@_IR', ai_template='l_@@@_AIR', dr_template='l_@@@_DR', ad_template='l_@@@_ADR'):

    df = pd.read_csv(results_file_path)

    chosen_model = df[chosen_model_field_name].values[0]  # either 'eq' for SIM or 'dif' for RIM
    stats = {'chosen model': chosen_model.replace('eq','SIM').replace('dif','RIM')}  # add model name
    
    # SIM('eq') lasso / mean

    rl_field_name = rl_template.replace('@@@', 'eq')
    ir_field_name = ir_template.replace('@@@', 'eq')
    dr_field_name = ir_field_name
    ai_field_name = ai_template.replace('@@@', 'eq')
    ad_field_name = ai_field_name
    
    stats.update({'SIM RL': df[rl_field_name].values[0],
                  'SIM IR': df[ir_field_name].values[0],
                  'SIM AI': df[ai_field_name].values[0],
                  })
    stats_tmp = {'RL': df[rl_field_name].values[0],
                  'IR': df[ir_field_name].values[0],
                  'DR': df[ir_field_name].values[0],
                  'AI': df[ai_field_name].values[0],
                  'AD': df[ai_field_name].values[0],
                  'model':'eq',
                  }
    stats['SIM Lasso Flag'] = verify_results_validity(stats_tmp, df, ai_field_name, ir_field_name, rl_field_name,
                                    dr_field_name, ad_field_name, minIR, maxIR, minAI, maxAI, msa_path)
    if not stats['SIM Lasso Flag']:
        stats.update({'SIM RL': df['m'+rl_field_name[1:]].values[0],
                      'SIM IR': df['m'+ir_field_name[1:]].values[0],
                      'SIM AI': df['m'+ai_field_name[1:]].values[0],
                      })

    # RIM('dif') lasso / mean
    rl_field_name = rl_template.replace('@@@', 'dif')
    ir_field_name = ir_template.replace('@@@', 'dif')
    dr_field_name = dr_template.replace('@@@', 'dif')
    ai_field_name = ai_template.replace('@@@', 'dif')
    ad_field_name = ad_template.replace('@@@', 'dif')
    stats.update({'RIM RL': df[rl_field_name].values[0],
              'RIM IR': df[ir_field_name].values[0],
              'RIM DR': df[dr_field_name].values[0],
              'RIM AI': df[ai_field_name].values[0],
              'RIM AD': df[ad_field_name].values[0],
              })
    stats_tmp = {'RL': df[rl_field_name].values[0],
              'IR': df[ir_field_name].values[0],
              'DR': df[ir_field_name].values[0],
              'AI': df[ai_field_name].values[0],
              'AD': df[ai_field_name].values[0],
              'model':'dif'}
    stats['RIM Lasso Flag'] = verify_results_validity(stats_tmp, df, ai_field_name, ir_field_name, rl_field_name,
                                    dr_field_name, ad_field_name, minIR, maxIR, minAI, maxAI, msa_path)

    if not stats['RIM Lasso Flag']:
        stats.update({'SIM RL': df['m'+rl_field_name[1:]].values[0],
                      'SIM IR': df['m'+ir_field_name[1:]].values[0],
                      'SIM DR': df['m'+dr_field_name[1:]].values[0],
                      'SIM AI': df['m'+ai_field_name[1:]].values[0],
                      'SIM AD': df['m'+ad_field_name[1:]].values[0],
                      })
    if verbose!=0:
        print(stats)
    pd.DataFrame(stats.items()).to_csv(results_file_path[:-17]+'sum_res.csv', header=False, index=False)
    return stats



# res_dir =  'D:/university/projects/abc/code/pipe_line_to_github/SpartaABC/pipeline/results/example_after_run/'
# msa_filename = 'msa.f'
# minIR=0
# maxIR=0.05
# verbose=1
             
# stats = get_stats_v2(results_file_path=res_dir+'data_name_res.csv', minIR=minIR, maxIR=maxIR, minAI=0, maxAI=2, msa_path=res_dir+msa_filename,verbose=1)


