#!/bin/bash
set -ev
# get coverage information
lcov --capture --directory /home/nuto/build --output-file ../coverage.info # 2> /dev/null
# filter out system stuff that we don't control
lcov -r ../coverage.info '/usr/include/*' -o ../coverage.info
# filter out external libraries
lcov -r ../coverage.info '/home/nuto/source/external/*' -o ../coverage.info
# upload to codecov
bash <(curl -s https://codecov.io/bash) -f ../coverage.info # > /dev/null 2> /dev/null
