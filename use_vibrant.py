import pandas as pd
import subprocess
from collections import defaultdict
from sequence_analysis import PhageSequence, clear_name

script = '../use_vibrant.sh'
output_path = 'temp/VIBRANT_output/VIBRANT_{0}/VIBRANT_results_{0}/VIBRANT_integrated_prophage_coordinates_{0}.tsv'


def run(*files) -> list:
    """return a list of output names"""
    # if the first argument is list, split this
    if isinstance(files[0], (list, tuple)):
        files = files[0]
    command = ' '.join([script, *files])
    # invoke command line
    print(subprocess.run(command, shell=True, cwd="temp/").stdout)
    output = [output_path.format(clear_name(f)) for f in files]
    return output


def parse_result(file) -> defaultdict:
    data = pd.read_table(file, sep='\t')
    result = defaultdict(list)
    for i, x in data.iterrows():
        # scaffold	fragment	protein start	protein stop	protein length	nucleotide start	nucleotide stop	nucleotide length
        contig_name = x['scaffold'].split()[0]
        start, end = x['nucleotide start'], x['nucleotide stop']
        result[contig_name].append(PhageSequence(contig_name, start, end, ''))
    return result


if __name__ == '__main__':
    print(parse_result(output_path.format('Pseudomonas_aeruginosa_PAO1')))
    # defaultdict(<class 'list'>, {
    #   'NC_002516.2': [PhageSequence(contig=NC_002516.2, start=673191, end=702831, length=29641, info=)]
    # })
