#!/bin/bash

#dataset="rb04-4096"
dataset="cw12b-1M"
#dataset="ccnews-1M"
#dataset="cw12b-4096"
#result_dir="/mnt/hgfs/Ddatabase/Result"
#data_dir="/mnt/hgfs/Ddatabase/Data"
#code_dir="/mnt/hgfs/Ddatabase/Build"

#result_dir="/media/sirius/BA167427F361AA8B1/Data/Result"
#data_dir="/media/sirius/BA167427F361AA8B1/Data/CIFF/Output"
#code_dir="/home/sirius/codes/cpp/Build"

result_dir="/mnt/home/xyzhu/result"
data_dir="/mnt/home/xyzhu/datasets"
code_dir="/mnt/home/xyzhu/codes/Build"

read_only="t"
source_dir="../pgm_index"
#epsilon=1024
index_type="pgm"
#decode_type="all"

mkdir -p "$result_dir/index/$dataset/$index_type"
mkdir -p "$result_dir/log/$dataset/$index_type"
#
## build index
#rm -r $code_dir
cmake -B $code_dir -S $source_dir
sleep 1
cd $code_dir
make -j
#16 32 62 63 64 126 127 128 254 255 256 512 1024 2048 4096 8192 16384 32768 65536 131072 262144
#64 128 256 512 1024 2048 4096 8192 16384 32768 65536 131072
#16 32 62 63 64 126 127 128 254 255 256 512 1024 2048 4096 8192 16384 32768 65536 131072 262144
#"simd" "normal" "simd_simple"
for dataset in "cw12b-1M" "ccnews-1M"
do
  for epsilon in 1
  do
    echo "————————————dataset : $dataset epsilon: $epsilon————————————"
    for decode_type in "simd"
    do
      for repeat in 1 2
      do
        mkdir -p $result_dir/index/$dataset/$index_type/$dataset-$index_type-$epsilon/
        $code_dir/pgm_decode $index_type $data_dir/$dataset/$dataset $result_dir/index/$dataset/$index_type/$dataset-$index_type-$epsilon/ $epsilon $read_only $decode_type $result_dir/log/$dataset/$index_type/$dataset-$index_type-$epsilon
      done
    done
  done
done
# build and check index

# pgm /media/sirius/BA167427F361AA8B/Data/CIFF/Output/rb04-4096/rb04-4096 /media/sirius/BA167427F361AA8B/Data/Result/index/rb04-4096/rb04-4096-pgm-64.idx 64
