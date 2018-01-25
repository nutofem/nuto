@page UploadDockerImage Uploading a new Docker image

1. Create an account on https://hub.docker.com/
2. Remember your account name and ask us to be made a "collaborator"
3. "Login" to your local Docker service using
~~~
docker login
~~~
4. Make the wanted changes to our `Dockerfile`
5. Create a new Docker image
~~~
docker build -t nuto_docker <nuto_source_directory>
~~~
6. Create a new tag (but give it a better name)
~~~
docker tag nuto_docker nuto/nuto_docker:my_fance_new_tag
~~~
7. Push it to Docker Hub
~~~
docker push nuto/nuto_docker:my_fancy_new_tag
~~~
