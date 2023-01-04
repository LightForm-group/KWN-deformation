scenario="tests/test_1a"

./run_kwn.sh ${scenario}

echo "checking output file checksums"

cd ${scenario}/results || output="no result folder" && output=$(diff ../reference_result_checksums.txt <(md5 *))

if [ -z "${output}" ];then
    echo "${scenario} outputs match reference outputs"
else
    echo ${output}
fi

echo "tests finished"



