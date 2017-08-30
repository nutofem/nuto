find -name '*.cpp' -o -name '*.h' | xargs clang-format -i

sudo docker run -ti -v $(pwd):/app --workdir=/app coala/base coala


