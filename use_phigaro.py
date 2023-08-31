import pandas as pd
import subprocess
from collections import defaultdict
from sequence_analysis import PhageSequence, clear_name

script = "../use_phigaro.sh"
out_dir = "temp/phigaro_output/"


def run(*files) -> list:
    """return a list of output names"""
    # if the first argument is list, split this
    if isinstance(files[0], (list, tuple)):
        files = files[0]
    command = ' '.join([script, *files])
    # invoke command line
    print(subprocess.run(command, shell=True, cwd="temp/").stdout)
    output = [out_dir + clear_name(f) for f in files]
    return output


def parse_result(name) -> defaultdict:
    file = name + '.phigaro.tsv'
    data = pd.read_table(file, sep='\t')
    result = defaultdict(list)
    for i, x in data.iterrows():
        # scaffold	begin	end	transposable	taxonomy
        contig_name = x['scaffold']
        start, end = x['begin'], x['end']
        result[contig_name].append(PhageSequence(contig_name, start, end, "taxonomy=" + x['taxonomy']))
    return result


if __name__ == '__main__':
    # print(run("C213_spades.fasta"))
    print(parse_result(out_dir + "C213_spades"))
    # result:
    # defaultdict(<class 'list'>, {
    #   'NODE_2_length_652681_cov_91.473123': [
    #       {'position': [299447, 339321], 'length': 39875, 'taxonomy': 'Myoviridae'}],
    #   ...
    # })
