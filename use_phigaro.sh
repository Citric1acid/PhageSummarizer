#!/bin/bash
# 批量处理细菌基因组寻找噬菌体片段，使用phigaro
# usage:
# ./use_phigaro.sh folder
# ./use_phigaro.sh file1 file2 ...
# ./use_phigaro.sh *

output_dir="phigaro_output"
threads=8
echo "使用phigaro处理fasta文件"
echo "$*"

if test -d "$*" ; then
  # 文件夹
  for name in `ls "$*"/*.fasta`
  do
    echo "Analyzing $name"
    phigaro -f $name -e html tsv stdout -o $output_dir -t $threads -d
    echo "Phigaro: $name done!"
  done
else
  # 多个文件
  for name in "$@"
  do
    if [[ $name == *.fasta ]]; then
      echo "Analyzing $name"
      phigaro -f $name -e html tsv stdout -o $output_dir -t $threads -d
      echo "Phigaro: $name done!"
    fi
  done
fi
