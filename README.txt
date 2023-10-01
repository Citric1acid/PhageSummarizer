Phage Summarizer

This program is used to find phage sequences in bacterial genomes.
This work summarizes and combines results from several existing phage-finding programs, in order to make better prediction on phage sequences. It combines the result from different programs, and uses BLAST to verify the results. It combines the advantage of different programs, and has higher accuracy. PHASTER, Phigaro, and VIBRANT are currently supported. This work will be posted online.

Author: Zigang Song
Contact info: songzigang@stu.pku.edu.cn

Dependencies:
  - python=3.11
  - biopython=1.81
  - Flask=2.3.3
  - numpy=1.25.0
  - pandas=2.0.3
  - plotly=5.16.1
  - requests=2.31.0


VIBRANT and Phigaro need python3.7, however
Their dependencies:
  - python=3.7
  - hmmer=3.3.2
  - phigaro=2.3.0
  - prodigal=2.6.3
  - scikit-learn=0.21.3
  - vibrant=1.2.1

Usage:
Run `python server.py` to start the server. 
The program accepts FASTA files as input. Users can upload multiple files at once for convenience. The program will run the methods as specified by the user, summarize their result and visualize. 
