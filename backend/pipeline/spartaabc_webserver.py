import sys
import os
import logging
import shutil
import Bio.SeqUtils
import pandas as pd
from time import sleep
import CONSTANTS as CONSTS  # from /bioseq/sincopa/
from newick_validator import is_newick
from auxiliaries import *

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger('main')


def verify_fasta_format(fasta_path, sequence_type):
    logger.info('Validating FASTA format')
    legal_chars = set(Bio.SeqUtils.IUPACData.ambiguous_dna_letters+'U' if sequence_type == 'nuc'
                      else Bio.SeqUtils.IUPACData.protein_letters)
    legal_chars.add('-')
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
                        if c.upper() not in legal_chars:
                            return f'Illegal content detected in <a href="https://www.ncbi.nlm.nih.gov/blast/fasta.shtml" target="_blank">FASTA</a> file. Line {line_number} in MSA contains illegal {"DNA" if sequence_type == "nuc" else "amino acid"} character "{c}".'
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
    with open(tree_path) as f:
        if not is_newick(f.read()):
            return 'Invalid newick format.'


def verify_msa_is_consistent_with_tree(msa_path, tree_path):
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
                    msg = f'{strain} species appears in the input MSA but not in the phylogenetic tree. Please make ' \
                        f'sure the phylogenetic tree you provide contains all the species in the provided ' \
                        f'MSA. and re-submit your job'
                    logger.error(msg)
                    return msg
                else:
                    logger.info(f'{strain} appears in tree!')
                    tree_strains.remove(strain)

    if tree_strains:
        msg = f'There is one (or more) species that appears in the input tree but not in the MSA ' \
            f'(e.g., {tree_strains.pop()}). Please make sure the phylogenetic tree you provide does ' \
            f'not contain irrelevant species and re-submit your job.'
        logger.error(msg)
        return msg


def validate_input(msa_path, sequence_type, tree_path, error_path):
    logger.info('Validating input...')

    error_msg = verify_fasta_format(msa_path, sequence_type)
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


def main(msa_path, tree_path, sequence_type, submodel, freqs, rates, inv_prop,
         gamma_shape, gamma_cats, maxIR, output_dir_path, html_path):

    error_path = f'{output_dir_path}/error.txt'
    try:
        if html_path:
            run_number = initialize_html(CONSTS, output_dir_path, html_path)
            # final_zip_path = f'{os.path.split(output_dir_path)[0]}/{CONSTS.WEBSERVER_NAME}_{run_number}'

        os.makedirs(output_dir_path, exist_ok=True)

        tmp_dir = f'{os.path.split(msa_path)[0]}/tmp'  # same folder as the input msa
        os.makedirs(tmp_dir, exist_ok=True)

        validate_input(msa_path, sequence_type, tree_path, error_path)

        results_file_path = os.path.join(output_dir_path, CONSTS.RESULT_CSV_NAME)
        print(results_file_path)
        if not os.path.exists(results_file_path):
            msa_name = os.path.split(msa_path)[-1]
            tree_name = os.path.split(tree_path)[-1]
            cmd = f"python {CONSTS.MAIN_SCRIPT} --path {output_dir_path}/" \
                f" --msaf {msa_name} --trf {tree_name} --mode {sequence_type}" \
                f" --submodel {submodel} --maxr {maxIR}"
            if submodel == 'GTR':
                cmd += f' --freq {" ".join(freqs)}'
                cmd += f' --rates {" ".join(rates)}'
                cmd += f' --inv-prop {inv_prop}'
                cmd += f' --gamma-shape {gamma_shape}'
                cmd += f' --gamma-cats {gamma_cats}'
            logger.info(f'''Calling:\n{cmd}\n\n''')
            os.system(cmd)
            sleep(CONSTS.RELOAD_INTERVAL)
        else:
            logger.info('\n'+'#'*50+f'\nSkipping SpartaABC analysis!\nResults exist at:\n{results_file_path}\n'+'#'*50)
        assert os.path.exists(results_file_path), \
            f'SpartaABC run was finished but due to unknown reason. The result file can not be found :-('

        if html_path:
            # shutil.make_archive(final_zip_path, 'zip', output_dir_path)
            finalize_html(html_path, error_path, run_number, CONSTS, results_file_path)

    except Exception as e:
        logger.info(f'SUCCEEDED = False')
        if html_path:
            error_msg = e.args[-1]
            if os.path.exists(error_path):  # in case any internal module crashed
                with open(error_path) as f:
                    error_msg = f.read()
            edit_failure_html(CONSTS, f'{error_msg}', html_path, run_number)
            add_closing_html_tags(html_path, CONSTS, run_number)


def finalize_html(html_path, error_path, run_number, CONSTS, results_file_path):
    succeeded = not os.path.exists(error_path)
    logger.info(f'SUCCEEDED = {succeeded}')
    if succeeded:
        edit_success_html(CONSTS, html_path, results_file_path)
    else:
        edit_failure_html(CONSTS, error_path, html_path, run_number)
    add_closing_html_tags(html_path, CONSTS, run_number)


def get_stats(results_file_path):

    ''' example results file:
    RIM A_D,1.5554392000000004
    RIM A_I,1.5710081
    RIM RL,77.33
    RIM R_D,0.030431140000000006
    RIM R_I,0.031627869999999995
    SIM A_ID,1.5909845000000005
    SIM RL,82.99
    SIM R_ID,0.07235858
    chosen model,SIM
    '''
    with open(results_file_path) as f:
        model = None
        for line in f:
            tokens = line.rstrip().split(',')
            if tokens[-1].isalpha():
                model = tokens[-1]
        assert model, 'No model was found in results file!'

    stats = {'model': model}
    with open(results_file_path) as f:
        for line in f:
            tokens = line.rstrip().split(' ')
            if model in tokens[0]:
                key, value = tokens[-1].split(',')
                stats[key] = float(value)

    logger.info(stats)
    return stats




def edit_success_html(CONSTS, html_path, results_file_path):
    logger.info('Updating html from "running" to "finished"...')
    update_html(html_path, 'RUNNING', 'FINISHED')

    logger.info('Getting stats...')
    stats = get_stats(results_file_path)

    additional_params = ''
    if stats['model'] == 'RIM':
        additional_params += f'''<b>R_D</b>: <font color='red'>{stats['R_D']:.4f}</font><br>
                                <b>A_D</b>: <font color='red'>{stats['A_D']:.2f}</font><br>'''
    logger.info('Appending to html...')
    append_to_html(html_path, f'''
            <div class="container" style="{CONSTS.CONTAINER_STYLE}">
            <h3><b>
            Best configuration found is as follows:</a>
            </b></h3>
            <br>
            <p>
                <b>Model:</b> <font color='red'>{stats['model']}</font><br>
                <b>RL:</b> <font color='red'>{int(stats['RL']+0.5)}</font><br>
                <b>{'R_ID' if stats['model'] == 'SIM' else 'R_I'}:</b> <font color='red'>{2*stats['R_ID'] if stats['model'] == 'SIM' else stats['R_I']:.4f}</font><br>
                <b>{'A_ID' if stats['model'] == 'SIM' else 'A_I'}:</b> <font color='red'>{stats['A_ID'] if stats['model'] == 'SIM' else stats['A_I']:.4f}</font><br>
                {additional_params}
            </p>
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
                            help='A path to an MSA file to infer indel parameters.',
                            type=lambda path: path if os.path.exists(path) else parser.error(f'{path} does not exist!'))
        parser.add_argument('input_tree_path',
                            help='A path to a background species tree that contains all the species in the input MSA.',
                            type=lambda path: path if os.path.exists(path) else parser.error(f'{path} does not exist!'))
        parser.add_argument('mode', choices=['nuc', 'amino'],
                            help='Type of sequences that will be simulated in Sparta')
        parser.add_argument('submodel', choices=['JC', 'GTR'],
                            help='Type of substitution model used for simulations')
        parser.add_argument('maxIR', type=lambda x: float(x),
                            help='Upper bound of the indel rate posterior.')
        parser.add_argument('output_dir_path',
                            help='A path to a folder in which the analysis results will be written.',
                            type=lambda path: path.rstrip('/'))
        parser.add_argument('--freqs', default=None, type=lambda x: x.split(','),
                            help='Comma-separated floats representing the a,c,g,t frequencies, respectively')
        parser.add_argument('--rates', default=None, type=lambda x: x.split(','),
                            help='Comma-separated floats representing the rate parameters a,b,c,d,e, respectively')
        parser.add_argument('--inv-prop', default=None, type=float, help='Invariable sites proportion')
        parser.add_argument('--gamma-shape', default=None, type=float, help='Gamma distribution shape')
        parser.add_argument('--gamma-cats', default=None, type=int, help='Number of discrete gamma categories')
        parser.add_argument('--html_path', default=None,
                            help='A path to an html file that will be updated during the run.',
                            type=lambda path: path if os.path.exists(path) else parser.error(f'{path} does not exist!'))

        parser.add_argument('-v', '--verbose', help='Increase output verbosity', action='store_true')

        args = parser.parse_args()

        if args.verbose:
            logging.basicConfig(level=logging.DEBUG)
        else:
            logging.basicConfig(level=logging.INFO)

        main(args.input_msa_path, args.input_tree_path, args.mode, args.submodel, args.freqs,
             args.rates, args.inv_prop, args.gamma_shape, args.gamma_cats,
             args.maxIR, args.output_dir_path, args.html_path)


