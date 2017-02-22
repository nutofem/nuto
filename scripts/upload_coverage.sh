#!/bin/bash
set -ev
if [[ "$COVERAGE" == "TRUE" ]]; then
    # get coverage information
    lcov --capture --directory build --output-file coverage.info > /dev/null
    # filter out system stuff that we don't control
    lcov -r coverage.info '/usr/include/*' -o coverage.info > /dev/null
    # upload to codecov
    bash <(curl -s https://codecov.io/bash)
fi
