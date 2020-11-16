#!/powerapps/share/centos7/python-anaconda3.6.5/bin/python

#############################################################################################################
# this file should be saved as part of the pipeline and the cgi should import it rather than copy it twice! #
#############################################################################################################

import os

# constants to use when sending e-mails using the server admin's email address.
ADMIN_EMAIL = 'TAU BioSequence <bioSequence@tauex.tau.ac.il>'
SMTP_SERVER = 'mxout.tau.ac.il'

OWNER_EMAIL = 'orenavram@gmail.com'

# general variables
SERVERS_RESULTS_DIR = '/bioseq/data/results'
SERVERS_LOGS_DIR = '/bioseq/data/logs'
RELOAD_INTERVAL = 5
RELOAD_TAGS = f'<META HTTP-EQUIV="REFRESH" CONTENT="{RELOAD_INTERVAL}"/>'

WEBSERVER_NAME = 'spartaabc'
WEBSERVER_URL = f'https://{WEBSERVER_NAME}.tau.ac.il'
WEBSERVER_TITLE = '<b>A web server to simulate sequences based on indel parameters inferred using an approximate Bayesian computation algorithm</b>'
WEBSERVER_TITLE = '<h2><i>SpartaABC</i><br>Inferring indel parameter values</h2>'

WEBSERVER_RESULTS_DIR = os.path.join(SERVERS_RESULTS_DIR, WEBSERVER_NAME)
WEBSERVER_LOGS_DIR = os.path.join(SERVERS_LOGS_DIR, WEBSERVER_NAME)
WEBSERVER_HTML_DIR = f'/data/www/html/{WEBSERVER_NAME}'

WEBSERVER_RESULTS_URL = os.path.join(WEBSERVER_URL, 'results')

Q_SUBMITTER_SCRIPT = '/bioseq/bioSequence_scripts_and_constants/q_submitter_power.py'
MAIN_SCRIPT = f'/bioseq/{WEBSERVER_NAME}/{WEBSERVER_NAME}_webserver.py'
REQUIRED_MODULES = ['python/python-anaconda3.6.5-orenavr2']

SUBMISSIONS_LOG = f'/bioseq/{WEBSERVER_NAME}/submissions_log.txt'
OUTPUT_DIR_NAME = 'outputs'
RESULT_WEBPAGE_NAME = 'result.html'
RESULT_CSV_NAME = 'data_name_res.csv'
EMAIL_FILE_NAME = 'email.txt'

# path to example data
EXAMPLE_MSA = os.path.join(WEBSERVER_HTML_DIR, 'example_msa.fas')
EXAMPLE_TREE = os.path.join(WEBSERVER_HTML_DIR, 'example_tree.newick')

CONTAINER_WIDTH = 'width: 850px'
CONTAINER_NO_MARGIN = 'margin: 0 auto'
CONTAINER_FONT = 'font-size: 20px'

CONTAINER_STYLE = f'{CONTAINER_WIDTH}; {CONTAINER_NO_MARGIN}; {CONTAINER_FONT}'

PROCESSING_MSG = f'<i>{WEBSERVER_NAME.upper()}</i> is now processing your request. This page will be automatically ' \
    f'updated every few seconds (until the job is done). You can also reload it manually. Once the job has finished, ' \
    f'several links to the output files will appear below. '

PROGRESS_BAR_ANCHOR = '''<!--progress_bar_anchor-->'''
PROGRESS_BAR_TAG = '''<div class="progress">
        <div class="progress-bar progress-bar-striped active" role="progressbar" style="width:100%">
        </div>
    </div>'''