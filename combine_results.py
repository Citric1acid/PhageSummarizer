# Phage Summarizer 的主要部分
# 汇总不同方法查找噬菌体的结果，并可视化
# 现可以运行的方法: PHASTER (API), Phigaro, VIBRANT

import os
import logging
import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from collections import defaultdict
from itertools import chain
from concurrent.futures import ThreadPoolExecutor

from sequence_analysis import PhageSequence, read_sequences, contig_length, clear_name, union_of_intervals
from draw_figs import draw_fig

# methods
import use_phaster
import use_phigaro
import use_vibrant
import use_blast

executor = ThreadPoolExecutor(max_workers=8)


class PhageResult:
    # input: temp/{name}.fasta
    # output: temp/{method}_output/{name}.fasta
    name: str
    method: str
    link: str
    data: dict
    blast_hit: dict
    # store results
    waiting = []
    results = []
    combine_used = []  # store this to prevent recombining
    temp_results = []  # store not finished jobs

    def __init__(self, name='', method='', link='', data=None, _id=0, blast_hit=None, done=1):
        self.name = name
        self.method = method
        self.link = link
        self.data = data if data is not None else {}
        self._id = _id
        self.blast_hit = blast_hit if blast_hit is not None else {}
        self.done = done  # 0 -> search -> 1 -> blast -> 2

    def __str__(self):
        return f'name: {self.name}, method: {self.method}, link: {self.link} data:\n' + str(self.data)

    def to_file(self, file):
        with open(file) as f:
            f.write(str(self))

    def to_list(self):
        if hasattr(self, '__list'):  # cached
            return self.__list

        # get the length of each contig first
        lengths = contig_length(read_sequences(f'temp/{self.name}.fasta'))

        def contig_to_list(contig, c: list):
            current = 1
            result = []
            c.sort(key=lambda x: x.start)

            for phage in c:
                start, end = phage.start, phage.end
                if start > current:
                    result.append({"contig": contig, "start": current, "end": start - 1, "type": "Bacterial"})
                result.append({"contig": contig, "start": start, "end": end, "type": "Phage"})
                current = end + 1
            # the final bacterial part
            if lengths[contig] > current:
                result.append({"contig": contig, "start": current, "end": lengths[contig], "type": "Bacterial"})
            return result

        self.__list = list(chain.from_iterable(contig_to_list(x, y) for x, y in self.data.items()))
        return self.__list

    def sequences(self, output=None):
        if self.done == 0:
            return
        seq_dict = SeqIO.to_dict(read_sequences(f'temp/{self.name}.fasta'))
        result = []
        for x, y in self.data.items():
            for phage in y:
                start, end = phage.start, phage.end
                seq = seq_dict[x].seq[start:end+1]
                seqid = f'{x}:{start}-{end}'
                result.append(SeqRecord(seq, id=seqid))

        if output is not None:
            SeqIO.write(result, output, 'fasta')
        return result

    def blast(self):
        # prepare its sequences first
        for x, y in self.data.items():
            for phage in y:
                seqid = f'{x}:{phage.start}-{phage.end}'
                # add to store blast hit
                self.blast_hit[seqid] = 'Waiting...'

        query = f'temp/{self.method}_output/{self.name}.fasta'
        output = f'temp/{self.method}_output/{self.name}.out'
        if self.method == 'Combined':
            output = f'temp/{self.method}_output/{self.name}_{self._id}.out'

        if not os.path.isfile(query):  # write fasta file
            self.sequences(query)

        if not self.data:  # empty result
            return 'No phage sequence.'

        if os.path.isfile(output):
            if not os.path.getsize(output):  # empty file
                return 'Please wait ...'

            # done
            result = use_blast.parse_result(output)
            self.done = 2
            if len(result) == 0:  # no alignment
                return 'No alignment'
            else:
                # find first alignment
                for i, row in result.iterrows():
                    seqid = row['qseqid']
                    if self.blast_hit.get(seqid) in (None, 'Waiting...'):
                        self.blast_hit[seqid] = row['stitle'].split(',')[0]
                for seqid in self.blast_hit:
                    # no alignment
                    if self.blast_hit[seqid] == 'Waiting...':
                        self.blast_hit[seqid] = '<span style="color: red">No alignment</span>'

                return result.to_html()

        else:
            logging.info(f'BLASTing {self.name}')
            executor.submit(use_blast.blast, query, output)
            return 'Please wait ...'

    def to_html(self):
        try:
            title = f'<div id={self._id} class="result"><h4>{self.name}  ({self.method})</h4>'
            link = f'<p><a href="{self.link}">Details</a></p>' if self.link else ''
            if self.done == 0:
                return title + link + str(self.data) + '</div>'

            link = f'''<p><a href="{self.link}">Details</a>  <a href="/fasta?id={id(self)}">Download Fasta</a>  <a 
                href="/blast?id={id(self)}">BLAST result</a></p>'''
            if not self.data:
                return title + link + '<p style="color: red">No phage sequence found.</p>' + '</div>'

            if isinstance(self.data, pd.DataFrame):
                table = self.data.to_html()
                fig = draw_fig(self.to_list(), height=150 + 50 * len(self.data))
                return title + link + table + fig + '</div>'

            self.blast()
            if self.blast_hit:
                # include blast result
                row = '<tr><td>{}</td> <td>{}</td> <td>{}</td> <td>{}</td> <td>{}</td>  <td>{}</td></tr>'
                table = '''<table border=1 cellspacing=0 bordercolor="gray"><thead><tr>
                        <th>Contig</th> <th>Start</th> <th>End</th> <th>Length</th> <th>Info</th> <th>BLAST hit</th>
                        </tr></thead><tbody>'''
                # write each row for the table
                blast_hit = iter(self.blast_hit.values())  # add this to each row
                for contig, locs in self.data.items():
                    if not locs:
                        continue
                    loc: PhageSequence = locs[0]
                    table += row.format(contig, loc.start, loc.end, loc.length, loc.info, next(blast_hit))
                    for loc in locs[1:]:
                        table += row.format('', loc.start, loc.end, loc.length, loc.info, next(blast_hit))

            else:
                row = '<tr><td>{}</td> <td>{}</td> <td>{}</td> <td>{}</td> <td>{}</td></tr>'
                table = '''<table border=1 cellspacing=0 bordercolor="gray"><thead><tr>
                        <th>Contig</th> <th>Start</th> <th>End</th> <th>Length</th> <th>Info</th></tr></thead>
                        <tbody>'''
                # write each row for the table
                for contig, locs in self.data.items():
                    if not locs:
                        continue
                    loc: PhageSequence = locs[0]
                    table += row.format(contig, loc.start, loc.end, loc.length, loc.info)
                    for loc in locs[1:]:
                        table += row.format('', loc.start, loc.end, loc.length, loc.info)
            table += '</tbody></table>'
            fig = draw_fig(self.to_list(), height=150 + 50 * len(self.data))
            return title + link + table + fig + '</div>'
        except Exception as Err:
            logging.error(repr(Err))
            return ''

    @classmethod
    def merge(cls, *instances):
        if len(instances) <= 1:  # no need to merge
            return
        if any(not o.done for o in instances):  # any o is not done
            return
        name = instances[0].name
        logging.info(f'merging results of {name}')
        new_data = defaultdict(list)

        for o in instances:
            # remove from result and put in combine_used
            cls.results.remove(o)
            if o.method != 'Combined':
                cls.combine_used.append(o)

            for x, y in o.data.items():
                new_data[x].extend(y)
        for x in new_data:
            new_intervals = union_of_intervals([phage.start, phage.end] for phage in new_data[x])
            new_data[x] = [PhageSequence(x, loc[0], loc[1]) for loc in new_intervals]
        # add combined result
        new_id = '+'.join(str(o._id) for o in instances)
        combined = cls(name, 'Combined', '', new_data, new_id)
        combined.blast()
        cls.results.append(combined)

    @classmethod
    def submit_jobs(cls, method, files: list) -> list:
        """submit files using method
        return the job ids or names to look up for result"""
        # make sure all files exist
        # files = list(filter(os.path.isfile, files))
        files = [f for f in files if os.path.isfile("temp/" + f)]
        if not files:  # empty input
            raise ValueError("No file submitted!")

        logging.info(f'uploading {files} using {method}')

        if method == "phaster":
            jobs = use_phaster.upload(files)
            new_jobs = [{"name": clear_name(f), "method": "phaster", "job": j} for f, j in zip(files, jobs)]

        elif method == "phigaro":
            executor.submit(use_phigaro.run, files)  # run it in background
            jobs = [use_phigaro.out_dir + clear_name(f) for f in files]
            new_jobs = [{"name": clear_name(f), "method": "phigaro", "job": j} for f, j in zip(files, jobs)]

        elif method == "VIBRANT":
            executor.submit(use_vibrant.run, files)  # run it in background
            jobs = [use_vibrant.output_path.format(clear_name(f)) for f in files]
            new_jobs = [{"name": clear_name(f), "method": "VIBRANT", "job": j} for f, j in zip(files, jobs)]

        else:
            logging.error(f'Wrong method: {method}')
            return []

        JobHandler.add_jobs(new_jobs)
        cls.waiting.extend(new_jobs)
        return [x['_id'] for x in new_jobs]

    @classmethod
    def get_result(cls, _id, name, method, job):
        if name is None and job is None:
            raise ValueError("No file or job provided")

        name = clear_name(name) if isinstance(name, str) else 'unknown'
        logging.info(f'Getting result for {name} using {method}')

        if method == "phaster":
            # read from temp if stored
            if os.path.isfile(f'temp/phaster_output/{name}.summary.txt'):
                summary = use_phaster.parse_summary(None, name)
                link = "http://phaster.ca/submissions/" + job
                return cls(name, "PHASTER", link, summary, _id)

            # get result from website
            link = use_phaster.phaster_get_url.format(job)
            data = use_phaster.get_result(link)
            link = "http://phaster.ca/submissions/" + job
            if data['done']:
                summary = use_phaster.parse_summary(data, name)
                return cls(name, "PHASTER", link, summary, _id)
            else:  # not done
                return cls(name, "PHASTER", link, data, _id, done=0)

        elif method == "phigaro":
            if not os.path.isfile(job + ".phigaro.tsv"):
                # not finished
                return cls(name, "Phigaro", '', 'Waiting', _id, done=0)
            link = job + ".phigaro.html"
            data = use_phigaro.parse_result(job)
            return cls(name, "Phigaro", link, data, _id)

        elif method == "VIBRANT":
            if not os.path.isfile(job):
                # not finished
                return cls(name, "VIBRANT", '', 'Waiting', _id, done=0)
            data = use_vibrant.parse_result(job)
            return cls(name, "VIBRANT", '', data, _id)

    @classmethod
    def make_combined_results(cls, ids):
        # classify results by name
        results_dict = defaultdict(list)
        for r in cls.results:
            if ids is not None and r._id not in ids:
                continue
            results_dict[r.name].append(r)

        for name in results_dict:
            cls.merge(*results_dict[name])

    @classmethod
    def refresh(cls, ids: list | None) -> list:
        cls.temp_results = []  # empty this first
        done = []
        for w in cls.waiting:  # {"_id": 0, "method": ..., "name": f, "job": j}
            if ids is not None and w['_id'] not in ids:
                continue
            try:
                result = cls.get_result(**w)
                if isinstance(result, cls):
                    if result.done:
                        # do BLAST automatically
                        result.blast()

                        # done, store it and remove from waiting
                        cls.results.append(result)
                        done.append(w)

                    else:  # not finished
                        cls.temp_results.append(result)
            except Exception as Err:
                logging.error(repr(Err))

        # remove these jobs from waiting
        cls.waiting = [w for w in cls.waiting if w not in done]
        # add combined results
        cls.make_combined_results(ids)
        all_results = cls.results + cls.combine_used + cls.temp_results
        return list(filter(lambda x: x._id in ids, all_results)) if ids is not None else all_results


class JobHandler:
    # this is actually a closure
    count = 0
    all_jobs = []

    def __init__(self):  # don't need its instance
        if not os.path.isfile('jobs.dat'):
            with open('jobs.dat', 'w') as f:
                f.write('_id\tname\tmethod\tjob\n')
        self.read_jobs()

    @classmethod
    def read_jobs(cls) -> list:
        # _id  name  method  job
        jobs = pd.read_table('jobs.dat', sep='\t')
        cls.all_jobs = jobs.to_dict('records')
        if len(cls.all_jobs) > 1:
            cls.count = cls.all_jobs[-1]['_id'] + 1
        return cls.all_jobs

    @classmethod
    def add_jobs(cls, jobs: dict | list[dict]):
        # this adds _id to each job
        if isinstance(jobs, dict):
            jobs = [jobs]
        with open('jobs.dat', 'a') as f:
            for j in jobs:
                j['_id'] = cls.count
                cls.all_jobs.append(j)
                f.write('{_id}\t{name}\t{method}\t{job}\n'.format(**j))
                cls.count += 1


# initialize this
JobHandler()
# shallow copy
PhageResult.waiting = JobHandler.all_jobs.copy()
