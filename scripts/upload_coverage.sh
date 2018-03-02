#!/bin/bash
set -ev
# get coverage information
lcov --capture --directory /home/nuto/build --output-file ../coverage.info
# filter out system stuff that we don't control
lcov -r ../coverage.info '/usr/include/*' -o ../coverage.info
# filter out external libraries
lcov -r ../coverage.info '/home/nuto/source/external/*' -o ../coverage.info
# upload to codecov, -F is the flag (unit or integration test)
bash <(curl -s https://codecov.io/bash) -f ../coverage.info -F $1
# delete gcda files so that they don't show up in next report
find /home/nuto/build -name '*.gcda' -delete
