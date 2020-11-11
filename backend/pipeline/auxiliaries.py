# from Bio import Phylo
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger('main')


def get_tree_labels(tree_to_prune_path):
    i = 0
    logger.info(f'Extracting leaf labels from {tree_to_prune_path}')
    while True:
        i += 1
        logger.info(f'Iteration #{i}')
        if i > 100:
            logger.info(f'Failed 100 times to extract leaf labels from {tree_to_prune_path}!')
            raise AssertionError(f'Failed 100 times to extract leaf labels from {tree_to_prune_path}!')
        try:
            tree = Phylo.read(tree_to_prune_path, 'newick')
            result = [leaf.name for leaf in tree.get_terminals()]
            break
        except:
            pass

    return result


def fail(error_msg, error_file_path):
    with open(error_file_path, 'w') as error_f:
        error_f.write(error_msg + '\n')
    raise Exception(error_msg)


def update_html(html_path, src, dst):
    # The initial file exists (generate by the cgi) so we can read and parse it.
    with open(html_path) as f:
        html_content = f.read()
    html_content = html_content.replace(src, dst)
    with open(html_path, 'w') as f:
        f.write(html_content)


def append_to_html(html_path, new_content):
    with open(html_path) as f:
        html_content = f.read()
    html_content += new_content
    with open(html_path, 'w') as f:
        f.write(html_content)


# def load_header2sequences_dict(fasta_path, get_length=False, upper_sequence=False):
#     header_to_sequence_dict = {}
#     seq_length = 0
#
#     with open(fasta_path) as f:
#         header = f.readline().lstrip('>').rstrip()
#         sequence = ''
#         for line in f:
#             line = line.rstrip()
#             if line.startswith('>'):
#                 seq_length = len(sequence)
#                 if upper_sequence:
#                     header_to_sequence_dict[header] = sequence.upper()
#                 else:
#                     # leave untouched
#                     header_to_sequence_dict[header] = sequence
#                 header = line.lstrip('>')
#                 sequence = ''
#             else:
#                 sequence += line
#
#         # don't forget last record!!
#         if sequence != '':
#
#             if upper_sequence:
#                 header_to_sequence_dict[header] = sequence.upper()
#             else:
#                 # leave untouched
#                 header_to_sequence_dict[header] = sequence
#
#     if get_length:
#         return header_to_sequence_dict, seq_length
#     else:
#         return header_to_sequence_dict



