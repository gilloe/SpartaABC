#!/powerapps/share/centos7/python-anaconda3.6.5/bin/python
import os
import shutil
import sys
import cgi
import cgitb
import subprocess
from time import time, ctime
from random import randint

if os.path.exists('/bioseq'):  # remote run
    sys.path.insert(0, '/bioseq/spartaabc')
    sys.path.insert(1, '/bioseq/bioSequence_scripts_and_constants')

import CONSTANTS as CONSTS
from auxiliaries import update_html, append_to_html  # from /bioseq/spartaabc/
from email_sender import send_email  # from /bioseq/bioSequence_scripts_and_constants/


def write_to_debug_file(cgi_debug_path_f, content):
    cgi_debug_path_f.write(f'{ctime()}: {content}\n')


def write_html_prefix(output_path, run_number):
    with open(output_path, 'w') as output_path_f:
        # html prefix. When creating a new one, copy it from the index page and replace the relevant values (i.e., web server name, etc...)
        output_path_f.write(f'''
<!DOCTYPE html>
<html lang="en">
    <head>
        <title>{CONSTS.WEBSERVER_NAME.upper()} Job #{run_number}</title>
        <meta http-equiv="cache-control" content="no-cache, must-revalidate, post-check=0, pre-check=0" />
        <meta http-equiv="cache-control" content="max-age=0" />
        <meta http-equiv="expires" content="0" />
        <meta http-equiv="expires" content="Tue, 01 Jan 1980 1:00:00 GMT" />
        <meta http-equiv="pragma" content="no-cache" />
        {CONSTS.RELOAD_TAGS}

        <link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.7/css/bootstrap.min.css">
        <link rel="stylesheet" href="https://gitcdn.github.io/bootstrap-toggle/2.2.2/css/bootstrap-toggle.min.css">
        <!-- <link rel="stylesheet" href="{CONSTS.WEBSERVER_URL}/css/general.css">
        <link rel="stylesheet" href="{CONSTS.WEBSERVER_URL}/css/nav.css"> -->
        <link rel="shortcut icon" type="image/x-icon" href="{CONSTS.WEBSERVER_URL}/SpartaABC_icon.gif" />

        <script src="https://ajax.googleapis.com/ajax/libs/jquery/3.3.1/jquery.min.js"></script>
        <script src="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.7/js/bootstrap.min.js"></script>
        <script src="https://code.jquery.com/jquery-1.10.2.js"></script>
    </head>
    <body>
        <nav role="navigation" class="navbar navbar-inverse navbar-fixed-top" id="nav">
            <div class="jumbotron" id="jumbo" style="margin-bottom: 10px; margin-top: 10px;">
                <div class="container">
                    <div class="row">
                        <div class="col-md-1">
                        </div>
                        <div class="col-md-2">
                            <img src="{CONSTS.WEBSERVER_URL}/SpartaABC_logo.gif" style="width:60%; float:right;">&nbsp;&nbsp;
                        </div>
                        <div class="col-md-8">
                            <div>
                            {CONSTS.WEBSERVER_TITLE}
                            </div>
                        </div>
                    </div>
                </div>
            </div>
        </nav>
        <div style="margin-top: 100px">
        </div>
        <br><br>
        <div class="container" style="{CONSTS.CONTAINER_STYLE}" align="justify">
            <br> 
            <br> 
            <br> 
            <H1 align=center>Job status: <font color='red'>QUEUED</font></h1>
            <br>
            {CONSTS.PROGRESS_BAR_ANCHOR}
            <br>
            {CONSTS.PROCESSING_MSG}A link to this page was sent to your email in case you wish to view these results at a later time without recalculating them. Please note that the results will be kept in the server for 3 months.
            <br><br>
        </div>\n\t\t''')


def exit_if_bot(form, wd):
    # email field should ALWAYS exist in the form (even if it's empty!!).
    # If it's not there, someone sent a request not via the website so they should be blocked.
    # confirm_email is a hidden field that only spamming bots might fill in...
    if 'email' not in form or ('confirm_email' in form and form['confirm_email'].value != ''):
        shutil.rmtree(wd)
        exit()


def peek_form(cgi_debug_path_f, form):
    # for debugging
    sorted_form_keys = sorted(form.keys())
    write_to_debug_file(cgi_debug_path_f,
                        f'These are the keys that the CGI received:\n{"; ".join(sorted_form_keys)}\n\n')
    write_to_debug_file(cgi_debug_path_f, 'Form values are:\n')
    for key in sorted_form_keys:
        if form[key].filename:
            write_to_debug_file(cgi_debug_path_f, f'100 first characters of {key} = {form[key].value[:100]}\n')
        else:
            write_to_debug_file(cgi_debug_path_f, f'{key} = {form[key]}\n')

    cgi_debug_path_f.flush()


def append_running_parameters_to_html(html_path, msa_name, tree_name, mode, substitution_model, maxIR,
                                      msa_is_provided_as_text, tree_is_provided_as_text, job_title, additional_parameters):

    content = f'<div class="container" style="{CONSTS.CONTAINER_STYLE}">\n'

    if job_title:
        content += '<div class="row" style="font-size: 20px;"><div class="col-md-12">\n'
        content += f'<h2><b>{job_title}</b></h2>\n'
        content += '</div></div>\n'

    content += '<div class="row" style="font-size: 20px;"><div class="col-md-12">\n'
    content += f'<b>Multiple sequence alignment:</b> {"raw text" if msa_is_provided_as_text else msa_name}\n'
    content += '</div></div>\n'

    content += '<div class="row" style="font-size: 20px;"><div class="col-md-12">\n'
    content += f'<b>Phylogenetic tree:</b> {"raw text" if tree_is_provided_as_text else tree_name}\n'
    content += '</div></div>\n'

    content += '<div class="row" style="font-size: 20px;"><div class="col-md-12">\n'
    content += f'<b>Sequence type:</b> {"DNA" if mode=="nuc" else "AA"}\n'
    content += '</div></div>\n'

    content += '<div class="row" style="font-size: 20px;"><div class="col-md-12">\n'
    content += f'<b>Substitution model:</b> {substitution_model}\n'
    content += '</div></div>\n'

    if substitution_model == 'GTR':
        for key in additional_parameters:
            content += '<div class="row" style="font-size: 20px;"><div class="col-md-12">\n'
            content += f'<b>{key}:</b> {additional_parameters[key]}\n'
            content += '</div></div>\n'

    content += '<div class="row" style="font-size: 20px;"><div class="col-md-12">\n'
    content += f'<b>Max indel rate:</b> {maxIR}\n'
    content += '</div></div>\n'

    content += '</div><br>\n'
    content += '\n\n<!--result-->\n\n\t\t'

    append_to_html(html_path, content)


def write_cmds_file(cmds_file, parameters, run_number):
    # the queue does not like very long commands so I use a dummy delimiter (!@#) to break the commands for q_submitter
    new_line_delimiter = '!@#'

    required_modules_as_str = ' '.join(CONSTS.REQUIRED_MODULES)
    with open(cmds_file, 'w') as f:
        f.write(f'module load {required_modules_as_str};')
        f.write(new_line_delimiter)
        f.write(f'python {CONSTS.WEBSERVER_WRAPPER} {parameters}\t{CONSTS.WEBSERVER_NAME}_{run_number}')
        f.write('\n')


def save_file_to_disk(cgi_debug_path_f, form, wd, form_field_name, file_name_on_disk=None):
    write_to_debug_file(cgi_debug_path_f, f'Saving FILE to disk')

    uploaded_file_name = form[form_field_name].filename
    write_to_debug_file(cgi_debug_path_f, f'Uploaded file name is:\n{uploaded_file_name}')

    data = form[form_field_name].value
    write_to_debug_file(cgi_debug_path_f, f'{uploaded_file_name} first 100 chars are: {data[:100]}\n')

    if file_name_on_disk is None:
        file_name_on_disk = uploaded_file_name
    data_path = os.path.join(f'{wd}/{file_name_on_disk}')

    with open(data_path, 'wb') as data_f:
        data_f.write(data)

    write_to_debug_file(cgi_debug_path_f, f'Uploaded data was saved to {data_path} successfully\n')

    return data_path, file_name_on_disk


def save_text_to_disk(cgi_debug_path_f, form, wd, form_field_name, file_name_on_disk):
    write_to_debug_file(cgi_debug_path_f, f'Saving TEXT to disk')

    data = form[form_field_name].value.rstrip()
    write_to_debug_file(cgi_debug_path_f, f'first 100 chars are: {data[:100]}\n')

    data_path = os.path.join(f'{wd}/{file_name_on_disk}')
    with open(data_path, 'w') as data_f:
        data_f.write(data)

    write_to_debug_file(cgi_debug_path_f, f'Uploaded data was saved to {data_path} successfully\n')

    return data_path


def run_cgi():
    # prints detailed error report on BROWSER when backend crashes
    # This line MUST appear (as is) BEFORE any error occurs to get a report about the exception!! otherwise you'll get a non-informatvie message like "internal server error"
    cgitb.enable()

    # print_hello_world() # for debugging
    form = cgi.FieldStorage()  # extract POSTed object

    # random_chars = "".join(choice(string.ascii_letters + string.digits) for x in range(20))
    # adding 20 random digits to prevent users guess a number and see data that are not their's
    run_number = str(round(time())) + str(randint(10 ** 9, 10 ** 10 - 1))
    # run_number = str(randint(10000, 10**5-1))
    # if form['example_page'].value == 'yes':
    #     run_number = 'example'

    results_url = os.path.join(CONSTS.WEBSERVER_RESULTS_URL, run_number)

    output_url = os.path.join(results_url, CONSTS.RESULT_WEBPAGE_NAME)

    wd = os.path.join(CONSTS.WEBSERVER_RESULTS_DIR, run_number)
    os.makedirs(wd)

    html_path = os.path.join(wd, CONSTS.RESULT_WEBPAGE_NAME)
    cgi_debug_path = os.path.join(wd, 'cgi_debug.txt')

    page_is_ready = os.path.exists(html_path)
    if not page_is_ready:
        write_html_prefix(html_path, run_number)  # html's prefix must be written BEFORE redirecting...

    print(f'Location: {output_url}')  # Redirects to the results url. MUST appear before any other print.
    print('Content-Type: text/html\n')  # For more details see https://www.w3.org/International/articles/http-charset/index#scripting
    sys.stdout.flush()  # must be flushed immediately!!!

    try:
        # first, make sure the submission is from a "real user"
        exit_if_bot(form, wd)

        # a file handler for cgi log
        cgi_debug_path_f = open(cgi_debug_path, 'w')

        write_to_debug_file(cgi_debug_path_f, f'{"#" * 100}\n{ctime()}: A new CGI request has been received!\n')

        # extract form's values:
        email = ''
        if 'email' in form and form['email'].value != '':
            email = form['email'].value.strip()

        # Send me a notification email every time there's a new request
        send_email(smtp_server=CONSTS.SMTP_SERVER,
                   sender=CONSTS.ADMIN_EMAIL,
                   receiver=f'{CONSTS.OWNER_EMAIL}',
                   subject=f'{CONSTS.WEBSERVER_NAME.upper()} - A new job has been submitted: {run_number}',
                   content=f'{os.path.join(CONSTS.WEBSERVER_URL, "results", run_number, "cgi_debug.txt")}\n'
                   f'{os.path.join(CONSTS.WEBSERVER_URL, "results", run_number, CONSTS.RESULT_WEBPAGE_NAME)}\n'
                   f'{email if email else "NO EMAIL"}')

        peek_form(cgi_debug_path_f, form)

        job_title = ''
        if form['job_title'].value != '':
            job_title = form['job_title'].value.strip()

        # minIR = 0
        # if 'IR_min_user_val' in form and form['IR_min_user_val'].value != '':
        #     minIR = float(form['IR_min_user_val'].value.strip())

        mode = 'nuc'
        if 'SeqType' in form and form['SeqType'].value != '':
            mode = form['SeqType'].value.strip()

        substitution_model = 'nuc'
        if 'SubType' in form and form['SubType'].value != '':
            substitution_model = form['SubType'].value.strip().upper()

        maxIR = 0.05
        if form['IR_max_user_val'].value != '':
            maxIR = float(form['IR_max_user_val'].value.strip())

        # until here it should roughly be the same in every CGI file!

        write_to_debug_file(cgi_debug_path_f, f'\n{"#" * 80}\nUploading data\n')

        msa_is_provided_as_text = tree_is_provided_as_text = False
        if form['msa_text'].value:
            msa_name = 'msa.fas'
            msa_path = save_text_to_disk(cgi_debug_path_f, form, wd,
                                         form_field_name='msa_text', file_name_on_disk=msa_name)
            msa_is_provided_as_text = True
        else:
            msa_path, msa_name = save_file_to_disk(cgi_debug_path_f, form, wd, form_field_name='msa_file')

        if form['tree_text'].value:
            tree_name = 'tree.newick'
            tree_path = save_text_to_disk(cgi_debug_path_f, form, wd,
                                          form_field_name='tree_text', file_name_on_disk=tree_name)
            tree_is_provided_as_text = True
        else:
            tree_path, tree_name = save_file_to_disk(cgi_debug_path_f, form, wd, form_field_name='tree_file')

        write_to_debug_file(cgi_debug_path_f, f'ls of {wd} yields:\n{os.listdir(wd)}\n')

        write_to_debug_file(cgi_debug_path_f, f'Extracting running parameters...\n')

        parameters = f'{msa_path} {tree_path} {mode} {substitution_model} {maxIR} {wd} --html_path {html_path}'

        additional_parameters = {}
        if substitution_model == 'GTR':
            freqs = []
            if form['freq_a'].value != '':
                freqs.append(form['freq_a'].value.strip())
            if form['freq_c'].value != '':
                freqs.append(form['freq_c'].value.strip())
            if form['freq_g'].value != '':
                freqs.append(form['freq_g'].value.strip())
            freqs.append(str(1-(sum(float(x) for x in freqs))))
            assert len(freqs) == 4, 'Not enough frequencies!'
            parameters += f' --freqs {",".join(freqs)}'
            additional_parameters['Nucleotide frequencies'] = ' ; '.join(f'{x}={float(y):.4f}' for x,y in zip('ACGT', freqs))

            rates = []
            if form['rate_a'].value != '':
                rates.append(form['rate_a'].value.strip())
            if form['rate_b'].value != '':
                rates.append(form['rate_b'].value.strip())
            if form['rate_c'].value != '':
                rates.append(form['rate_c'].value.strip())
            if form['rate_d'].value != '':
                rates.append(form['rate_d'].value.strip())
            if form['rate_e'].value != '':
                rates.append(form['rate_e'].value.strip())
            assert len(rates) == 5, 'Not enough rates!'
            parameters += f' --rates {",".join(rates)}'
            additional_parameters['Rate parameters'] = ' ; '.join(f'{x}={float(y):.4f}' for x,y in zip('abcde', rates))

            inv_prop = form['param_i'].value.strip()
            parameters += f' --inv-prop {inv_prop}'
            additional_parameters['Invariable sites proportion'] = inv_prop

            gamma_shape = form['param_g'].value.rstrip()
            parameters += f' --gamma-shape {gamma_shape}'
            additional_parameters['Gamma distribution shape'] = gamma_shape

            gamma_cats = form['param_cat'].value.strip()
            parameters += f' --gamma-cats {gamma_cats}'
            additional_parameters['# of discrete gamma categories'] = gamma_cats


        write_to_debug_file(cgi_debug_path_f, f'{ctime()}: Writing running parameters to html...\n')

        if not page_is_ready:
            append_running_parameters_to_html(html_path, msa_name, tree_name, mode, substitution_model, maxIR,
                                              msa_is_provided_as_text, tree_is_provided_as_text, job_title, additional_parameters)

        write_to_debug_file(cgi_debug_path_f, f'{ctime()}: Running parameters were written to html successfully.\n')
        cmds_file = os.path.join(wd, 'qsub.cmds')
        write_cmds_file(cmds_file, parameters, run_number)

        job_id_file = os.path.join(wd, 'job_id.txt')

        # a simple command when using shebang header (#!) in q_submitter_power.py
        submission_cmd = f'{CONSTS.Q_SUBMITTER_SCRIPT} {cmds_file} {wd} -q pupkowebr --verbose > {job_id_file}'

        if not page_is_ready:
            write_to_debug_file(cgi_debug_path_f, f'\nSUBMITTING JOB TO QUEUE:\n{submission_cmd}\n')
            subprocess.call(submission_cmd, shell=True)
        else:
            write_to_debug_file(cgi_debug_path_f, f'\nPage already exists! no need to run the analysis again\n')

        user_email_file = os.path.join(wd, CONSTS.EMAIL_FILE_NAME)
        if email != '':
            with open(user_email_file, 'w') as email_f:
                email_f.write(f'{email}\n')
            try:
                # Send the user a notification email for their submission
                notify_user(run_number, email, job_title,
                            'raw text' if msa_is_provided_as_text else msa_name,
                            'raw text' if tree_is_provided_as_text else tree_name,
                            minIR, maxIR)
            except:
                write_to_debug_file(cgi_debug_path_f, f'\nFailed sending notification to {email}\n')

        else:
            try:
                os.remove(user_email_file)  # for example mode
            except OSError:
                pass

        job_title_file = os.path.join(wd, 'job_title.txt')
        if job_title != '':
            with open(job_title_file, 'w') as job_f:
                job_f.write(f'{job_title}\n')
        else:
            try:
                os.remove(job_title_file)  # for example mode
            except OSError:
                pass

        write_to_debug_file(cgi_debug_path_f, f'\n\n{"#" * 50}\nCGI finished running!\n{"#" * 50}\n')

        cgi_debug_path_f.close()

        send_email(smtp_server=CONSTS.SMTP_SERVER,
                   sender=CONSTS.ADMIN_EMAIL,
                   receiver=f'{CONSTS.OWNER_EMAIL}',
                   subject=f'{CONSTS.WEBSERVER_NAME.upper()} job {run_number} by {email} has been failed!',
                   content=f"{email}\n\n{os.path.join(CONSTS.WEBSERVER_URL, 'results', run_number, CONSTS.RESULT_WEBPAGE_NAME)}\n"
                   f"\n{os.path.join(CONSTS.WEBSERVER_URL, 'results', run_number, 'cgi_debug.txt')}")

    except Exception as e:

        msg = 'CGI crashed before the job was submitted :('
        with open(html_path) as f:
            html_content = f.read()
        html_content = html_content.replace('QUEUED', 'FAILED')
        html_content += f'<br><br><br><center><h2><font color="red">{msg}</font><br><br>Please try to re-run your job or <a href="mailto:{CONSTS.ADMIN_EMAIL}?subject={CONSTS.WEBSERVER_NAME.upper()}%20Run%20Number%20{run_number}">contact us</a> for further information</h2></center><br><br>\n</body>\n</html>\n'
        with open(html_path, 'w') as f:
            f.write(html_content)

        exc_type, exc_obj, exc_tb = sys.exc_info()
        fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
        with open(cgi_debug_path, 'w') as cgi_debug_path_f:
            write_to_debug_file(cgi_debug_path_f,
                                f'\n{"$" * 100}\n\n{msg}\n\n{fname}: {exc_type}, at line: {exc_tb.tb_lineno}\n\ne.args[0]: {e.args[0]}\n\n{"$" * 100}')

        # Send me a notification email every time there's a failure
        try:
            email = form['email'].value.strip() if form['email'].value.strip() else 'NO_EMAIL'
        except:
            email = 'NO_EMAIL'
        send_email(smtp_server=CONSTS.SMTP_SERVER,
                   sender=CONSTS.ADMIN_EMAIL,
                   receiver=f'{CONSTS.OWNER_EMAIL}',
                   subject=f'{CONSTS.WEBSERVER_NAME.upper()} job {run_number} by {email} has been failed!',
                   content=f"{email}\n\n{os.path.join(CONSTS.WEBSERVER_URL, 'results', run_number, CONSTS.RESULT_WEBPAGE_NAME)}\n"
                   f"\n{os.path.join(CONSTS.WEBSERVER_URL, 'results', run_number, 'cgi_debug.txt')}")

        # logger.info(f'Waiting {2*CONSTS.RELOAD_INTERVAL} seconds to remove html refreshing headers...')
        # Must be after flushing all previous data. Otherwise it might refresh during the writing.. :(
        from time import sleep

        # sleep(2 * CONSTS.RELOAD_INTERVAL)
        with open(html_path) as f:
            html_content = f.read()
        html_content = html_content.replace(CONSTS.RELOAD_TAGS, f'<!--{CONSTS.RELOAD_TAGS}-->')
        with open(html_path, 'w') as f:
            f.write(html_content)

    # logging submission
    with open(CONSTS.SUBMISSIONS_LOG, 'a') as f:
        f.write(f'{email}\t{run_number}\t{ctime()}\n')

    with open(cgi_debug_path, 'a') as f:  # for cgi debugging
        f.write(f'{ctime()}: Submission was documented in \n')


def notify_user(run_number, email, job_title, msa_name, tree_name, minIR, maxIR):
    job_name = f'{job_title}\n' if job_title else ''
    notification_content = f'Your submission details are:\n\n{job_name}'

    notification_content += f'Multiple sequence alignment: {msa_name}\n'
    notification_content += f'Phylogenetic tree: {tree_name}\n'
    notification_content += f'{minIR} <= indel rate posterior <= {maxIR}\n'

    notification_content += f'Once the analysis will be ready, we will let you know! Meanwhile, you can track the ' \
        f'progress of your job at:\n{CONSTS.WEBSERVER_URL}/results/{run_number}/{CONSTS.RESULT_WEBPAGE_NAME}\n\n'

    send_email(smtp_server=CONSTS.SMTP_SERVER,
               sender=CONSTS.ADMIN_EMAIL,
               receiver=f'{email}',
               subject=f'{CONSTS.WEBSERVER_NAME.upper()} - your job has been submitted! (Run number: {run_number})',
               content=notification_content)


if __name__ == '__main__':
    run_cgi()
