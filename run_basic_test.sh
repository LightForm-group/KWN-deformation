#!/bin/bash

run_analysis ()
{
./run_kwn.sh ${scenario} 2>&1 > ${scenario}_log.txt 
echo "checking output file checksums"
cd ${scenario}/results || output="no result folder" && output=$(diff ../reference_result_checksums.txt <(md5 *))
if [ -z "${output}" ];then
    echo "${scenario} outputs match reference outputs"
else
    echo ${output}
fi
cd -
}

scenario="tests/test_1a"
run_analysis

scenario="tests/test_2a"
run_analysis

echo "tests finished"



