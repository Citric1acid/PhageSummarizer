#!/bin/bash
# usage:
# ./use_vibrant.sh folder
# ./use_vibrant.sh file1 file2 ...
# ./use_vibrant.sh *

conda activate py37
output_dir="VIBRANT_output"
threads=8
echo "使用VIBRANT处理fasta文件"
echo "$*"

if test -d "$*" ; then
  # 文件夹
  for name in `ls "$*"/*.fasta`
  do
    echo "Analyzing $name"
    VIBRANT_run.py -i $name -folder $output_dir -t $threads
    echo "VIBRANT: $name done!"
  done
else
  # 多个文件
  for name in "$@"
  do
    if [[ $name == *.fasta ]]; then
      echo "Analyzing $name"
      VIBRANT_run.py -i $name -folder $output_dir -t $threads
      echo "VIBRANT: $name done!"
    fi
  done
fi
