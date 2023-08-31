# invoke blast
import subprocess
import pandas as pd


def blast(query, output):
    command = ' '.join(['./use_blast.sh', query, output])
    print(subprocess.run(command, shell=True).stdout)


def parse_result(output) -> pd.DataFrame:
    # qseqid sseqid stitle pident evalue length sstart send
    table = pd.read_table(output, sep='\t', header=None)
    table.columns = ['qseqid', 'sseqid', 'stitle', 'pident', 'evalue', 'length', 'sstart', 'send']
    return table
