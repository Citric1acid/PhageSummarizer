# 用 PHASTER 在基因组中找噬菌体序列，使用PHASTER的api
# http://phaster.ca/instructions#urlapi
import os.path
import logging
import requests
import json
from time import sleep
from collections import defaultdict
from sequence_analysis import PhageSequence, readfile, read_sequences, find_item

phaster_url = "http://phaster.ca/phaster_api"
phaster_url_contigs = "http://phaster.ca/phaster_api?contigs=1"
phaster_get_url = "http://phaster.ca/phaster_api?acc={}"


def post_file(file, contigs=None):
    if contigs is None:
        # decide if it contains multiple contigs automatically
        contigs = len(read_sequences(file)) > 1
    url = phaster_url_contigs if contigs else phaster_url
    ans = requests.post(url, open(file))
    # Response: {"job_id":"ZZ_023a167bf8","status":"You're next!..."}
    return json.loads(ans.content)


def get_result(url=None, job_id=None):
    if url is None:
        if job_id is None:
            raise ValueError("No url or job_id provided")
        url = phaster_get_url.format(job_id)
    ans = requests.get(url)
    # Response: {"job_id": "NC_000913.3", "status": "Complete", "url": "phaster.ca/submissions/NC_000913.3",
    #            "zip": "phaster.ca/submissions/NC_000913.3.zip", "summary": "Criteria for scoring prophage regions..."}
    data = json.loads(ans.content)
    if 'summary' in data:
        data['done'] = True
        logging.info(f'{data["job_id"]} done')
    else:
        data['done'] = False
    return data


def parse_summary(data, file=None):
    temp = f'temp/phaster_output/{file}.summary.txt'
    if isinstance(data, dict):
        if 'summary' not in data:
            return
        summary: str = data['summary']
        print(summary, file=open(temp, 'w'))  # store in temp
    else:
        # read from temp directly
        if not os.path.isfile(temp):  # no result yet
            return
        summary = readfile(temp)
    lines = summary.split('\n')
    result = defaultdict(list)  # creates a list for each key
    default_contig_name = None
    skip = True

    for l in lines:
        # phage 的信息在 ----- 行以后，前面的行忽略
        if l.strip().startswith('-----'):
            skip = False
            continue
        if skip:
            continue

        terms = l.split()
        if not terms:  # empty line
            continue
        # each line:
        # ['REGION', 'REGION_LENGTH', 'COMPLETENESS(score)', 'SPECIFIC_KEYWORD', 'REGION_POSITION', 'TRNA_NUM',
        # 'TOTAL_PROTEIN_NUM', 'PHAGE_HIT_PROTEIN_NUM', 'HYPOTHETICAL_PROTEIN_NUM', 'PHAGE+HYPO_PROTEIN_PERCENTAGE',
        # 'BACTERIAL_PROTEIN_NUM', 'ATT_SITE_SHOWUP', 'PHAGE_SPECIES_NUM', 'MOST_COMMON_PHAGE_NAME(hit_genes_count)',
        # 'FIRST_MOST_COMMON_PHAGE_NUM', 'FIRST_MOST_COMMON_PHAGE_PERCENTAGE', 'GC_PERCENTAGE']
        # REGION_POSITION: NODE_2_length_652681_cov_91.473123:299387-343631
        score = terms[2]
        position = terms[4]
        if ':' in position:
            # contig_name:start-end
            contig_name, position = position.split(':')
            contig_name = contig_name.split(',')[0]
            start, end = map(int, position.split('-'))
        else:
            # start-end
            if default_contig_name is None:  # get the id of original sequence
                default_contig_name = read_sequences(f'temp/{file}.fasta')[0].id \
                    if os.path.isfile(f'temp/{file}.fasta') else "genome"
            contig_name = default_contig_name
            start, end = map(int, position.split('-'))
        # store result as a PhageSequence
        result[contig_name].append(PhageSequence(contig_name, start, end, "score=" + score))

    return result


def upload(files: list) -> list:
    """return a list of job_ids"""
    job_ids = []
    for file in files:
        response = post_file("temp/" + file)
        job_ids.append(response['job_id'])
    return job_ids


def wait_for_results(job_ids: list):
    result = []
    try:
        while job_ids:
            update = []  # need waiting another round
            for j in job_ids:
                data = get_result(phaster_get_url.format(j))
                if data['done']:
                    result.append(parse_summary(data))
                else:
                    # not finished
                    update.append(j)
            job_ids = update
            sleep(1)
    finally:
        return result  # always return


if __name__ == '__main__':
    # file = 'C213_spades.fasta'
    # print(upload([file]))
    print(wait_for_results(['ZZ_285989e558']))
    # result:
    # [defaultdict(<class 'list'>, {
    #   'NODE_2_length_652681_cov_91.473123': [
    #       {'position': [299387, 343631], 'length': '44.2Kb', 'score': 'intact(140)'}],
    #   ...
    # })]
