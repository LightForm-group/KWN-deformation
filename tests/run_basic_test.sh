#!/bin/bash

run_analysis ()
{
# prepare scenario
rm -r ${scenario}/results
mkdir ${scenario}/results
cp ${scenario}/namelist.input .
# run model
~/bin/KWN-deform 2>&1 > ${scenario}_log.txt 
# run checks on model outputs
echo "checking output file checksums"
cd ${scenario}/results || output="no result folder" && output=$(diff ../reference_result_checksums.txt <(md5 *))
if [ -z "${output}" ];then
    echo "${scenario} outputs match reference outputs"
else
    echo ${output}
fi
cd -
}

scenario="./test_1a"
run_analysis

scenario="./test_2a"
run_analysis

echo "tests finished"



