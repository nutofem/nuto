sudo docker build -t nuto_docker .
sudo docker tag nuto_docker nuto/nuto_docker:xenial
sudo docker push nuto/nuto_docker:xenial
