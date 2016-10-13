@page Coverage Getting the coverage information for the test suite

First of all, you need to have `lcov` installed, i.e.

    sudo apt-get install lcov

Then, you need to set `ENABLE_COVERAGE` to `TRUE` in your CMake options.

At last, you'll need to run `make coverage` (or `ninja coverage` or ...).
