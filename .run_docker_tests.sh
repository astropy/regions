#!/bin/bash

git submodule init
git submodule update

docker info

cat << EOF | docker run -i \
                        -v ${PWD}:/regions_src \
                        -a stdin -a stdout -a stderr \
                        astropy/affiliated-32bit-test-env:1.6 \
                        bash || exit $?

cd /regions_src

echo "Output of uname -m:"
uname -m

echo "Output of sys.maxsize in Python:"
python -c 'import sys; print(sys.maxsize)'

pip install https://github.com/astrofrog/pytest-arraydiff/archive/master.zip

python setup.py test -V -a "-s"

EOF
