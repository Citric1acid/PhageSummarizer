# Server of Phage Summarizer
# using flask server
# made by Zigang Song

import os
from flask import Flask, render_template, request, send_file

from log_handle import logging
from sequence_analysis import change_suffix, find_item
from combine_results import PhageResult

app = Flask(__name__)


def save_temp_files(file_storages):
    for f in file_storages:
        f.save("temp/" + change_suffix(f.filename.replace(' ', '_')))  # replace space with _


def check_file(name: str):
    accepted = ('.fasta', '.fa', '.fna')
    ok = any(name.endswith(a) for a in accepted)
    if not ok:
        logging.warning(f'Not accepted format: {name}')
    return ok


@app.route('/', methods=['GET'])
def index():
    return render_template('index.html')


@app.route('/upload', methods=['POST'])
def upload():
    if 'fileInput' not in request.files:
        return 'Bad Request', 400

    file_storages = request.files.getlist('fileInput')
    # check files
    if not all(check_file(f.filename) for f in file_storages):
        return 'Wrong file format', 400

    files = [change_suffix(f.filename.replace(' ', '_')) for f in file_storages]  # replace space with _
    methods = request.values.getlist('method')
    logging.info(f'Submitting files: {files}')
    logging.info(f'Using methods: {methods}')

    save_temp_files(file_storages)

    ids = []
    err = []
    if "phaster" in methods:
        try:
            ids.extend(PhageResult.submit_jobs("phaster", files))
        except Exception as Err:
            logging.error(Err, stack_info=True)
            err.append("Error happened when using PHASTER")
    if "phigaro" in methods:
        try:
            ids.extend(PhageResult.submit_jobs("phigaro", files))
        except Exception as Err:
            logging.error(Err, stack_info=True)
            err.append("Error happened when using Phigaro")
    if "VIBRANT" in methods:
        try:
            ids.extend(PhageResult.submit_jobs("VIBRANT", files))
        except Exception as Err:
            logging.error(Err, stack_info=True)
            err.append("Error happened when using VIBRANT")
    if err:
        return '\n'.join(err), 400
    else:
        return str(ids)


@app.route('/refresh', methods=['GET'])
def refresh():
    if 'ids' in request.values:
        request_ids = request.values['ids']
        if not request_ids:  # empty
            return ''

        logging.info(f'Refreshing {request_ids}')
        ids = list(map(int, request_ids.split(',')))
    else:
        ids = None
    results = PhageResult.refresh(ids)
    return ''.join(p.to_html() for p in results)


@app.route('/fasta', methods=['GET'])
def get_fasta():
    obj_id = int(request.values['id'])
    logging.info(f'Getting fasta of {obj_id}')
    phage: PhageResult = find_item(PhageResult.results, lambda x: id(x) == obj_id) or \
        find_item(PhageResult.combine_used, lambda x: id(x) == obj_id)
    if phage is None:
        logging.error(f'Id {obj_id} does not exist')
        return 'Invalid', 403
    fasta = f'temp/{phage.method}_output/{phage.name}.fasta'
    if not os.path.isfile(fasta):  # write fasta file
        phage.sequences(fasta)
    return send_file(fasta)


@app.route('/blast', methods=['GET'])
def blast():
    obj_id = int(request.values['id'])
    logging.info(f'BLASTing {obj_id}')
    phage: PhageResult = find_item(PhageResult.results, lambda x: id(x) == obj_id) or \
        find_item(PhageResult.combine_used, lambda x: id(x) == obj_id)
    if phage is None:
        logging.error(f'Id {obj_id} does not exist')
        return 'Invalid', 403

    return phage.blast()


@app.route('/<path:path>', methods=['GET'])
def get_file(path):
    if os.path.isfile(path):
        if path.endswith(".html") or path.endswith(".fasta"):  # only allow html and fasta
            return send_file(path)
        else:
            return 'Invalid', 403
    else:
        return 'Not Found', 404


if __name__ == '__main__':
    app.run(host='0.0.0.0', port=8080)
