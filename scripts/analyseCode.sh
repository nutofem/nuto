#!/bin/bash
find . \( -name '*.cpp' -o -name '*.h' \) -exec clang-format -i +
sudo docker run -ti -v "$(pwd)":/app --workdir=/app coala/base coala
