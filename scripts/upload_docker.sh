docker build -t nuto_docker .
docker tag nuto_docker nuto/nuto_docker:dev
docker push nuto/nuto_docker:dev
