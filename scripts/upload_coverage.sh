#!/bin/bash
set -ev
if [[ "$COVERAGE" == "--coverage" ]]; then
    # get coverage information
    lcov --capture --directory build --output-file coverage.info 2> /dev/null
    # filter out system stuff that we don't control
    lcov -r coverage.info '/usr/include/*' -o coverage.info
    # filter out external libraries
    EXTERNALFILES="$(pwd)/external/*"
    lcov -r coverage.info "$EXTERNALFILES" -o coverage.info
    # upload to codecov
    bash <(curl -s https://codecov.io/bash) -f coverage.info > /dev/null 2> /dev/null
fi
