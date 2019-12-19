#!/bin/bash
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null && pwd )"

echo "Installing GRASP to ${DIR} grasp"
cd ${DIR} || exit

if [ -e "${DIR}/grasp" ]; then
    echo "${DIR}/grasp already exists. Aborting."
    exit 1
fi

git clone https://github.com/mortenpi/grasp.git grasp || exit
cd grasp/ && git checkout mp/libgrasp2 || exit
./configure.sh && cd build/ || exit
make && make install || exit
