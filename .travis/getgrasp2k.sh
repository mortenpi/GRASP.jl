#!/bin/bash
TARBALL="2018-01-18-grasp2k-qed.tar.gz"

if [ -z ${TRAVIS_BUILD_DIR+x} ]; then
    echo "ERROR: \$TRAVIS_BUILD_DIR unset."
    exit 1
fi
INSTALL_DIR=`dirname ${TRAVIS_BUILD_DIR}`
echo "Installing GRASP2K to ${INSTALL_DIR}"

cd ${INSTALL_DIR}
wget https://s3-eu-west-1.amazonaws.com/massey-mortenpi/${TARBALL} || exit 2

mkdir -p grasp2k && cd grasp2k/ || exit 3
tar -xf "${INSTALL_DIR}/${TARBALL}" || exit 4

./configure.sh && cd build/ || exit 5
make && make install || exit 6
