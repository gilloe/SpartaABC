import sys
import os
import logging
import shutil
import Bio.SeqUtils
import pandas as pd
from time import sleep
import CONSTANTS as CONSTS  # from /bioseq/sincopa/
from auxiliaries import *

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger('main')


def verify_fasta_format(fasta_path):
    logger.info('Validating FASTA format')
    Bio.SeqUtils.IUPACData.ambiguous_dna_letters += 'U-'
    legal_chars = set(Bio.SeqUtils.IUPACData.ambiguous_dna_letters.lower() +
                      Bio.SeqUtils.IUPACData.ambiguous_dna_letters +
                      Bio.SeqUtils.IUPACData.protein_letters.lower() +
                      Bio.SeqUtils.IUPACData.protein_letters)
    with open(fasta_path) as f:
        line_number = 0
        try:
            line = f.readline()
            line_number += 1
            if not line.startswith('>'):
                return f'Illegal <a href="https://www.ncbi.nlm.nih.gov/blast/fasta.shtml" target="_blank">FASTA format</a>. First line in MSA starts with "{line[0]}" instead of ">".'
            previous_line_was_header = True
            putative_end_of_file = False
            curated_content = f'>{line[1:]}'.replace("|", "_")
            for line in f:
                line_number += 1
                line = line.strip()
                if not line:
                    if not putative_end_of_file: # ignore trailing empty lines
                        putative_end_of_file = line_number
                    continue
                if putative_end_of_file:  # non empty line after empty line
                    return f'Illegal <a href="https://www.ncbi.nlm.nih.gov/blast/fasta.shtml" target="_blank">FASTA format</a>. Line {putative_end_of_file} in MSA is empty.'
                if line.startswith('>'):
                    if previous_line_was_header:
                        return f'Illegal <a href="https://www.ncbi.nlm.nih.gov/blast/fasta.shtml" target="_blank">FASTA format</a>. MSA contains an empty record. Both lines {line_number-1} and {line_number} start with ">".'
                    else:
                        previous_line_was_header = True
                        curated_content += f'>{line[1:]}\n'.replace("|", "_")
                        continue
                else:  # not a header
                    previous_line_was_header = False
                    for c in line:
                        if c not in legal_chars:
                            return f'Illegal <a href="https://www.ncbi.nlm.nih.gov/blast/fasta.shtml" target="_blank">FASTA format</a>. Line {line_number} in MSA contains illegal DNA character "{c}".'
                    curated_content += f'{line}\n'
        except UnicodeDecodeError as e:
            logger.info(e.args)
            line_number += 1  # the line that was failed to be read
            return f'Illegal <a href="https://www.ncbi.nlm.nih.gov/blast/fasta.shtml" target="_blank">FASTA format</a>. Line {line_number} in MSA contains one (or more) non <a href="https://en.wikipedia.org/wiki/ASCII" target="_blank">ascii</a> character(s).'
    # override the old file with the curated content
    with open(fasta_path, 'w') as f:
        f.write(curated_content)


def verify_newick_format(tree_path):
    logger.info('Validating NEWICK format')
    # TODO
    pass


def verify_msa_is_consistent_with_tree(msa_path, tree_path):
    logger.info('Validating MSA is consistent with the species tree')
    logger.info('TODO!!')
    return
    logger.info('Validating MSA is consistent with the species tree')

    tree_strains = get_tree_labels(tree_path)
    logger.info(f'Phylogenetic tree contains the following strains:\n{tree_strains}')

    with open(msa_path) as f:
        logger.info(f'Checking MSA...')
        for line in f:
            if line.startswith('>'):
                # make sure each header appears in tree
                strain = line.lstrip('>').rstrip('\n')
                logger.debug(f'Checking {strain}...')
                if strain not in tree_strains:
                    msg = f'{strain} species appears in the input MSA but not in the phylogenetic tree. ' \
                        f'Please make sure the phylogenetic tree you provide contains (at least) all ' \
                        f'the species in the provided MSA.'
                    logger.error(msg)
                    return msg
                else:
                    logger.info(f'{strain} appears in tree!')


def validate_input(msa_path, tree_path, error_path):
    logger.info('Validating input...')

    error_msg = verify_fasta_format(msa_path)
    if error_msg:
        fail(error_msg, error_path)

    error_msg = verify_newick_format(tree_path)
    if error_msg:
        fail(error_msg, error_path)

    error_msg = verify_msa_is_consistent_with_tree(msa_path, tree_path)
    if error_msg:
        fail(error_msg, error_path)


def fix_input(msa_path, tree_path, output_dir, tmp_dir):
    tree_name = os.path.split(tree_path)[-1]
    adjusted_tree = f'{output_dir}/{os.path.splitext(tree_name)[0]}_fixed{os.path.splitext(tree_name)[-1]}'
    fix_tree(msa_path, tree_path, tmp_dir, adjusted_tree)

    msa_name = os.path.split(msa_path)[-1]
    fixed_msa_path = f'{output_dir}/{os.path.splitext(msa_name)[0]}_fixed{os.path.splitext(msa_name)[-1]}'
    fix_msa(msa_path, fixed_msa_path)

    return fixed_msa_path, adjusted_tree


def main(msa_path, tree_path, minIR, maxIR, minAI, maxAI, output_dir_path, html_path):

    error_path = f'{output_dir_path}/error.txt'
    if html_path:
        run_number = initialize_html(CONSTS, output_dir_path, html_path)
        final_zip_path = f'{os.path.split(output_dir_path)[0]}/{CONSTS.WEBSERVER_NAME}_{run_number}'


    os.makedirs(output_dir_path, exist_ok=True)

    tmp_dir = f'{os.path.split(msa_path)[0]}/tmp'  # same folder as the input msa
    os.makedirs(tmp_dir, exist_ok=True)

    validate_input(msa_path, tree_path, error_path)

    results_file_path = os.path.join(output_dir_path, CONSTS.RESULT_CSV_NAME)
    if not os.path.exists(results_file_path):
        sys.path.insert(0, '/bioseq/spartaabc/pipeline')
        from pipeline import pipeline
        msa_name = os.path.split(args.input_msa_path)[-1]
        tree_name = os.path.split(args.input_tree_path)[-1]
        logger.error(f'''Calling:\npipeline('{output_dir_path}/', '{msa_name}', '{tree_name}', {minIR}, {maxIR})\n\n''')
        pipeline(output_dir_path+'/', msa_name, tree_name, minIR=minIR, maxIR=maxIR)
        sleep(CONSTS.RELOAD_INTERVAL)
    else:
        logger.error('#'*50+f'\nSkipping SpartaABC analysis!\nResults exist at: {results_file_path}\n'+'#'*50)
    assert os.path.exists(results_file_path), \
        f'SpartaABC run was finished but due to unknown reason the results file can not be found :-('

    if html_path:
        # shutil.make_archive(final_zip_path, 'zip', output_dir_path)
        finalize_html(html_path, error_path, run_number, CONSTS, results_file_path, minIR, maxIR, minAI, maxAI, msa_path)

    # except Exception as e:
    #     logger.info(f'SUCCEEDED = False')
    #     if html_path:
    #         error_msg = e.args[-1]
    #         if os.path.exists(error_path):  # in case any internal module crashed
    #             with open(error_path) as f:
    #                 error_msg = f.read()
    #         edit_failure_html(CONSTS, f'{error_msg}\n{sys.version}', html_path, run_number)
    #         add_closing_html_tags(html_path, CONSTS, run_number)


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


def finalize_html(html_path, error_path, run_number, CONSTS, results_file_path, minIR, maxIR, minAI, maxAI, msa_path):
    succeeded = not os.path.exists(error_path)
    logger.info(f'SUCCEEDED = {succeeded}')
    if succeeded:
        edit_success_html(CONSTS, html_path, results_file_path, minIR, maxIR, minAI, maxAI, msa_path)
    else:
        edit_failure_html(CONSTS, error_path, html_path, run_number)
    add_closing_html_tags(html_path, CONSTS, run_number)


def verify_results_validity(stats, df, ai_field_name, ir_field_name, rl_field_name,
                            dr_field_name, ad_field_name, minIR, maxIR,
                            minAI, maxAI, msa_path):

    non_lasso_flags = set()
    if not minIR <= stats['IR'] <= maxIR:
        ir_field_name = 'm' + ir_field_name[1:]  # l_eq_IR -> m_eq_IR / l_dif_IR -> m_dif_IR
        stats['IR'] = df[ir_field_name].values[0]
        non_lasso_flags.add('IR')

    if not minAI <= stats['AI'] <= maxAI:
        ai_field_name = 'm' + ai_field_name[1:]  # l_eq_AIR -> m_eq_AIR / l_dif_AIR -> m_dif_AIR
        stats['AI'] = df[ai_field_name].values[0]
        non_lasso_flags.add('AI')
        logger.info(type(stats['AI']))  # TODO
        print(type(df[ai_field_name].values[0]))
        print(type(str(df[ai_field_name].values[0])))
        # stats['AI'] = str(df[ai_field_name].values[0]) + " (the mean is used instead of Lasso regression's inference)"

    logger.info('Calculating legal RL range...')
    minRL, maxRL = get_RL_bounds(msa_path)
    logger.info(f'minRL={minRL}, maxRL={maxRL}')

    if not minRL <= stats['RL'] <= maxRL:
        rl_field_name = 'm' + rl_field_name[1:]  # l_eq_RL -> m_eq_RL
        stats['RL'] = df[rl_field_name].values[0]
        non_lasso_flags.add('RL')

    if stats['model'] == 'dif':  # RIM model, i.e., dr_field_name and ad_field_name are not None
        if not minIR <= stats['DR'] <= maxIR:
            dr_field_name = 'm' + dr_field_name[1:]  # l_dif_DR -> m_dif_DR
            stats['DR'] = df[dr_field_name].values[0]
            non_lasso_flags.add('DR')

        if not minAI <= stats['AD'] <= maxAI:
            ad_field_name = 'm' + ad_field_name[1:]  # l_dif_AD -> m_dif_AD
            stats['AD'] = df[ad_field_name].values[0]
            non_lasso_flags.add('AD')

    return non_lasso_flags


def get_stats(results_file_path, minIR, maxIR, minAI, maxAI, msa_path, chosen_model_field_name='nn_c_class', rl_template='l_@@@_RL',
              ir_template='l_@@@_IR', ai_template='l_@@@_AIR', dr='l_@@@_DR', ad='l_@@@_ADR'):

    df = pd.read_csv(results_file_path)

    chosen_model = df[chosen_model_field_name].values[0]  # either 'eq' for SIM or 'dif' for RIM
    stats = {'model': chosen_model}  # add model name

    rl_field_name = rl_template.replace('@@@', df[chosen_model_field_name].values[0])
    ir_field_name = ir_template.replace('@@@', df[chosen_model_field_name].values[0])
    ai_field_name = ai_template.replace('@@@', df[chosen_model_field_name].values[0])
    stats.update({'RL': df[rl_field_name].values[0],
                  'IR': df[ir_field_name].values[0],
                  'AI': df[ai_field_name].values[0]})

    dr_field_name = ad_field_name = None
    if chosen_model == 'dif':  # RIM model
        dr_field_name = dr.replace('@@@', df[chosen_model_field_name].values[0])
        ad_field_name = ad.replace('@@@', df[chosen_model_field_name].values[0])
        stats.update({'DR': df[dr_field_name].values[0],
                      'AD': df[ad_field_name].values[0]})

    logger.info(f'Results BEFORE verification:\n{stats}')
    non_lasso_flags_set = verify_results_validity(stats, df, ai_field_name, ir_field_name, rl_field_name,
                                                  dr_field_name, ad_field_name, minIR, maxIR, minAI, maxAI,
                                                  msa_path)
    logger.info(f'Results AFTER verification:\n{stats}')

    return stats, non_lasso_flags_set


def edit_success_html(CONSTS, html_path, results_file_path, minIR, maxIR, minAI, maxAI, msa_path):
    logger.info('Updating html from "running" to "finished"...')
    update_html(html_path, 'RUNNING', 'FINISHED')

    logger.info('Getting stats...')
    param2values, non_lasso_flags_set = get_stats(results_file_path, minIR, maxIR, minAI, maxAI, msa_path)
    non_lasso_txt = '*'

    additional_params = ''
    if 'DR' in param2values:
        additional_params += f'''<b>DR</b>: <font color='red'>{param2values['DR']:.4f}</font>{non_lasso_txt if 'DR' in non_lasso_flags_set else ''}<br>
                                <b>AD</b>: <font color='red'>{param2values['AD']:.2f}</font>{non_lasso_txt if 'AD' in non_lasso_flags_set else ''}<br>'''
    logger.info('Appending to html...')
    append_to_html(html_path, f'''
            <div class="container" style="{CONSTS.CONTAINER_STYLE}">
            <h3><b>
            Best parameters found are as follows:</a>
            </b></h3>
            <br>
            <p>
                <b>Model:</b> {'SIM' if param2values['model'] == 'eq' else 'RIM'}<br>
                <b>RL:</b> <font color='red'>{int(param2values['RL']+0.5)}</font>{non_lasso_txt if 'RL' in non_lasso_flags_set else ''}<br>
                <b>IR:</b> <font color='red'>{param2values['IR']:.4f}</font>{non_lasso_txt if 'IR' in non_lasso_flags_set else ''}<br>
                <b>AI:</b> <font color='red'>{param2values['AI']:.2f}</font>{non_lasso_txt if 'AI' in non_lasso_flags_set else ''}<br>
                {additional_params}
            </p>
            {non_lasso_txt+'Mean value is used instead of lasso regression (lasso inferred value is out of prior).' if non_lasso_flags_set else ''}
            </div>
''')


def edit_failure_html(CONSTS, error_msg, html_path, run_number):
    update_html(html_path, 'RUNNING', 'FAILED')
    append_to_html(html_path,
                   f'<div class="container" style="{CONSTS.CONTAINER_STYLE}" align="justify"><h3>\n'
                   f'<font color="red">{error_msg}</font></h3><br><br>'
                   f'Please make sure your input is OK and then try to re-run your job or '
                   f'<a href="mailto:{CONSTS.ADMIN_EMAIL}?subject={CONSTS.WEBSERVER_NAME}%20Run%20Number:%20{run_number}">'
                   f'contact us'
                   f'</a> '
                   f'for further information.<br>'
                   f'</div>\n')


def add_closing_html_tags(html_path, CONSTS, run_number):
    update_html(html_path, CONSTS.PROCESSING_MSG, '')  # remove "web server is now processing your request" message
    update_html(html_path, 'progress-bar-striped active', 'progress-bar-striped')  # stop_progress_bar

    append_to_html(html_path, f'''
            <hr>
                <h4 class=footer>
                    <p align='center'>Questions and comments are welcome! Please
                        <span class="admin_link"> 
                        <a href="mailto:{CONSTS.ADMIN_EMAIL}?subject={CONSTS.WEBSERVER_NAME.upper()}%20Run%20Number%20{run_number}">contact us</a> 
                        </span>
                    </p>
                </h4>
                <div id="bottom_links" align="center"><span class="bottom_link">
                <a href="{CONSTS.WEBSERVER_URL}" target="_blank">Home</a> 
                 | 
                <a href="{CONSTS.WEBSERVER_URL}/overview.php" target="_blank">Overview</a> 
            </span>
        </div>
        <br><br><br>
    </body>
</html>''')

    # have to be last thing that is done
    sleep(2*CONSTS.RELOAD_INTERVAL)
    update_html(html_path, CONSTS.RELOAD_TAGS, f'<!--{CONSTS.RELOAD_TAGS}-->')  # stop refresh



def initialize_html(CONSTS, output_dir_path, html_path):

    path_tokens = output_dir_path.split('/')
    # e.g., "/bioseq/data/results/spartaabc/12345678"
    run_number = path_tokens[path_tokens.index(CONSTS.WEBSERVER_NAME) + 1]

    update_html(html_path, 'QUEUED', 'RUNNING')

    return run_number


if __name__ == '__main__':
        from sys import argv
        print(f'Starting {argv[0]}. Executed command is:\n{" ".join(argv)}')

        import argparse
        parser = argparse.ArgumentParser()
        parser.add_argument('input_msa_path',
                            help='A path to a DNA MSA file to infer indel parameters.',
                            type=lambda path: path if os.path.exists(path) else parser.error(f'{path} does not exist!'))
        parser.add_argument('input_tree_path',
                            help='A path to a background species tree that contains all the species in the input MSA.',
                            type=lambda path: path if os.path.exists(path) else parser.error(f'{path} does not exist!'))
        parser.add_argument('output_dir_path',
                            help='A path to a folder in which the analysis results will be written.',
                            type=lambda path: path.rstrip('/'))
        parser.add_argument('minIR', type=lambda x: float(x),
                            help='Lower bound of the indel rate posterior.')
        parser.add_argument('maxIR', type=lambda x: float(x),
                            help='Upper bound of the indel rate posterior.')
        parser.add_argument('--minAI', default=1.001, type=lambda x: float(x),
                            help='Lower bound of the "a parameter" posterior.')
        parser.add_argument('--maxAI', default=2, type=lambda x: float(x),
                            help='Upper bound of the "a parameter" posterior.')
        parser.add_argument('--html_path', default=None,
                            help='A path to an html file that will be updated during the run.',
                            type=lambda path: path if os.path.exists(path) else parser.error(f'{path} does not exist!'))

        parser.add_argument('-v', '--verbose', help='Increase output verbosity', action='store_true')

        args = parser.parse_args()

        if args.verbose:
            logging.basicConfig(level=logging.DEBUG)
        else:
            logging.basicConfig(level=logging.INFO)

        main(args.input_msa_path, args.input_tree_path, args.minIR, args.maxIR,
             args.minAI, args.maxAI, args.output_dir_path, args.html_path)


