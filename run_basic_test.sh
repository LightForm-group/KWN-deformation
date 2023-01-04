scenario="tests/test_1a"

./run_kwn.sh ${scenario}

echo "checking output file checksums"

cd ${scenario}/results && output=$(diff ../reference_result_checksums.txt <(md5 *)) || output="no result folder"

if [ -z "${output}" ];then
    echo "${scenario} outputs match reference outputs"
else
    echo ${output}
fi

echo "tests finished"



