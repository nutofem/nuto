@page CustomRepoGuide GUIDE to a full NuTo Custom repository with Travis CI and code coverage

# GUIDE to a full NuTo Custom repository with Travis CI and code coverage

## Introduction

This guide shows you how to setup your own git repository for your custom NuTo applications on github. It will allow you to work on NuTo and your own applications without the need of reinstalling NuTo after every change in the NuTo source code.


## Basic Setup

First create your own github repository. You need to log in on www.github.com and click on the "+" in the upper right corner. Then select new repository. Follow the steps until you have created a repository. Now check out your empty repository to a location of your choice. Switch to the "code" tab of your repository and click on the https button under "quick setup". Copy the link and type the following command in a bash:  

~~~bash
git clone ${HTTPS_LINK} ${FOLDER_WHERE_YOU_WANT_YOUR_REPO}
~~~

switch to the folder and create a CMakeLists.txt

~~~bash
cd myNuToAppsFolder
touch CMakeLists.txt
~~~


A submodule is something like a link to another repository. You don't store the code inside your repository, but you always get the code, when you pull your repo. Use the following commands to create a submodule:

~~~bash
git submodule add https://github.com/nutofem/nuto.git nuto
git submodule init
git submodule update
~~~

The first command adds the repository at https://github.com/nutofem/nuto.git (NuTo) as a submodule into the folder nuto. The other two commands are for initialization. You can also just pull NuTo without using submodules, but using submodules have some advantages. For example, the current nuto branch you are working with is stored when you do a commit. When you pull your repository, the branch is also pulled (Not sure here. You might need to swicth to the nuto directory and pull it. Travis does it automatically). This is especially helpful when using Travis CI to test your code that needs NuTo (more on this later). However, a submodule is not automatically pulling the latest version on the repository. If you want to get the latest version, you have to go into the submodule folder and pull it

~~~bash
git pull
~~~

You can read more about submodules here:
https://chrisjean.com/git-submodules-adding-using-removing-and-updating/

Now you can do your first commit and push to your repository. Use:

~~~bash
git status
~~~

to get an overview about all untracked files and folders. Add them all with

~~~bash
git add ${UNTRACKED_FILE_OR_FOLDER}
~~~

Then do 

~~~bash
git commit
~~~

write a meaningful message and push to your repository:

~~~bash
git push
~~~



## CMake Setup

Open the CMakeLists with an editor of your choice. Add the header:

~~~cmake
cmake_minimum_required(VERSION 3.9)
project(NuToApps)
~~~

Since you are going to build NuTo from a directory that is not the original root directory, you need to to some adjustments, otherwise the compiler might not find some includes:

~~~cmake
find_path(EIGEN_INCLUDE_DIR NAMES Eigen/Core PATH_SUFFIXES eigen3)
include_directories(${EIGEN_INCLUDE_DIR})
include_directories("nuto")   # <---!!! might be nuto/src in older versions
~~~ 
`IMPORTANT:` Since I have not included every NuTo header so far, there might be other adjustments that need to be added in order to compile your programs

Now we want to add NuTo in our Build tree, so that it gets compiled together with our applications. Therefore we need to include a cmake module first:
~~~cmake
include(ExternalProject)
~~~

Now we can use the script to add NuTo to our build tree:
~~~cmake
ExternalProject_Add(nuto
    CMAKE_ARGS "-DCMAKE_INSTALL_PREFIX=${CMAKE_CURRENT_BINARY_DIR}"
               "-DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}"
    SOURCE_DIR "${PROJECT_SOURCE_DIR}/RELATIVE_PATH_TO_NUTO"
    BUILD_ALWAYS 1
    )
~~~
The provided CMAKE_ARGS tell cmake to compile NuTo with the same build type as our project and install the shared libs into our build directory instead of the default folder. The SOURCE_DIR command tells cmake where to find nuto. BUILD_ALWAYS makes sure, that NuTo always gets compiled before your apps.

The libs are added to a subfolder into your build directory. Therefore you need to tell cmake where to find them:
~~~cmake
link_directories("${CMAKE_CURRENT_BINARY_DIR}/lib")
~~~


Now you can build an executable using NuTo:
~~~cmake
add_executable(Test test.cpp)
add_dependencies(Test nuto)
target_link_libraries(Test
    NuToBase
    NuToMath
    NuToMechanics
    )
~~~
The `add_dependencies` makes sure, that the libs are build first, before cmake is looking for them.

This is basically everything you need to build, but we might not be done yet (see next topic). The full CMakeLists looks like this:

~~~cmake
cmake_minimum_required(VERSION 3.7)
project(NuToApps)
set(CMAKE_VERBOSE_MAKEFILE OFF)

find_path(EIGEN_INCLUDE_DIR NAMES Eigen/Core PATH_SUFFIXES eigen3)
include_directories(${EIGEN_INCLUDE_DIR})
include_directories("nuto")


include(ExternalProject)
ExternalProject_Add(NuTo
    CMAKE_ARGS "-DCMAKE_INSTALL_PREFIX=${CMAKE_CURRENT_BINARY_DIR}"
               "-DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}"
    SOURCE_DIR "${PROJECT_SOURCE_DIR}/nuto"
    BUILD_ALWAYS 1
    )


link_directories("${CMAKE_CURRENT_BINARY_DIR}/lib")
~~~




## Setting the shared library search path

In some environments it might occur that you can build your applications, but they do not find the (correct) NuTo libs. There are two possible problems:

1. The program will crash and tell you it didn't find the NuTolibs. 
2. You can build your applications and run them, but if you change something in NuTo, it won't affect your program.

Both problems are related to the fact, that NuTo is build as a shared library by default. Since we installed the libs into the build directory and not into the default system library directory, you might have to tell your executables where to look for them. Problem 1. is obviously related to this. The second problem occurs, if you have installed a NuTo version before, so the application is linked against the wrong shared libs. Both can be fixed by setting the environmental variable `LD_LIBRARY_PATH` to your build directory.

You can do this in the terminal:
~~~bash
export LD_LIBRARY_PATH=PATH_TO_YOUR_BUILD_DIRECTORY/lib
~~~

But this will only work for your current shell. To make it permant, you have to edit your .bashrc. If you are using QTCreator, you can also set an environmental variable in the "Projects"-Menu (Ctrl+5). Add the `LD_LIBRARY_PATH` in the "Build environment" section with the correct folder as value.


## Travis CI

If you want to be sure that your stuff is also running on different machines and setups, you should include travis ci into your project. This is the case for NuTo itself and there is no reason to not include it into your project, except spending some time to get it running correctly. If you follow the instructions below, this should not be an issue.

`IMPORTANT:` Since everything might change over time, the documentation might get outdated. You can always have a look at NuTo's `.travis.yml` to get an idea what might have changed and what you need to adjust. If thats the case, don't forget to update this documentation!

First thing you have to to, is to register on https://travis-ci.org . If it is not part of the registration process, you now have to link your github and your travis account. Google it, if you get stuck here, but as far as I remember, it was not hard to accomplish. If this is done, open your profile on travis-ci.org . You should now see all your repositories listed and marked with a gray "x". Just click the gray field to activate travis for your repository. All you need to do now, is to provide a file named `.travis.yml` (don't miss the dot at the front) to your projects root repository, which configures your test environment.

I will show you one possible configuration first, before we go through each command step by step:

~~~yml
language: cpp
sudo: false
services: docker
dist: trusty

matrix:
    include:
        - compiler: clang
          env: BUILD_TYPE=Release
        - compiler: g++
          env: BUILD_TYPE=Release
        - compiler: g++
          env: BUILD_TYPE=Debug COVERAGE=--coverage
        - compiler: clang
          env: BUILD_TYPE=Release WERROR=-Werror
        - compiler: g++
          env: BUILD_TYPE=Release WERROR=-Werror
    allow_failures:
        - env: BUILD_TYPE=Release WERROR=-Werror

before_install:
    - travis_retry timeout 120 docker pull nuto/nuto_docker:dev
    - docker run -itd --name dock -u 0 -v $(pwd):/home/nuto/source nuto/nuto_docker:dev

script:
    - docker exec dock cmake -DCMAKE_CXX_COMPILER=$CXX -DCMAKE_BUILD_TYPE=$BUILD_TYPE -DCMAKE_CXX_FLAGS="$COVERAGE $WERROR" ../source
    - docker exec dock make -j2
    - docker exec dock ctest
~~~

So the header is quiet self explanatory:

~~~yml
language: cpp
sudo: false
services: docker
dist: trusty
~~~

We want to test something written in cpp, do not need sudo permissions on the server and request the additional docker (virtual machine) service on a ubuntu trusty (14.04) system.

The next part sets your different test configurations

~~~yml
matrix:
    include:
        - compiler: clang
          env: BUILD_TYPE=Release
        - compiler: g++
          env: BUILD_TYPE=Release
        - compiler: g++
          env: BUILD_TYPE=Debug COVERAGE=--coverage
        - compiler: clang
          env: BUILD_TYPE=Release WERROR=-Werror
        - compiler: g++
          env: BUILD_TYPE=Release WERROR=-Werror
    allow_failures:
        - env: BUILD_TYPE=Release WERROR=-Werror
~~~

The `include:` section covers the different build types you want to use. The 'compiler:' flag is obvious. The `env:` flag sets environmental variables that you use later inside instructions and if statements. The `allow_failures:` section specifies build setups, that are allowed to fail. So the rating if your test passed or not does not depend on this setup. In the example above, we want two builds to fail, if the compiler creates any warnings, but we don't want those builds to cause a failure of the whole test series. 

Now we define what we want to do before the tests are executed using the flag `before_install:`. (Notice that all the following sections basically just take bash commands that you can run on your own system.)

~~~yml
before_install:
    - travis_retry timeout 120 docker pull nuto/nuto_docker:dev
    - docker run -itd --name dock -u 0 -v $(pwd):/home/nuto/source nuto/nuto_docker:dev
~~~

First we try to pull the current NuTo docker image (more infos regarding docker can be found in the following section):

~~~bash
docker pull nuto/nuto_docker:dev
~~~

The preceding `travis_retry timeout 120` just tells travis to retry it after 120 second in case the server is not responding. Then we start the image with 

~~~bash
docker run -itd --name dock -u 0 -v $(pwd):/home/nuto/source nuto/nuto_docker:dev
~~~

This will start the image (`docker run  ... nuto/nuto_docker:dev`) in detached mode (`-itd`). In detached mode, the image won't "take over" the bash on travis. To send commands to that image, we give it the name dock (`--name dock`). The command `-v $(pwd):/home/nuto/source` mounts our current working directory (`$(pwd)`) on travis to the folder `/home/nuto/source` inside our image. This way we are capable of sharing data between the host (Travis) and our virtual machine. In this case we share the NuTo source code, which Travis automatically pulls, with our image. This way we avoid pulling it twice. The flag `-u 0` tells the image, that we want to log in as root. At this point this is not necessary, but if you want write access to the mounted directory (see codecov section) you need active root permissions.

As stated before, you are just executing bash commands, so the following section should be somehow self explanatory:

~~~yml
script:
    - docker exec dock cmake -DCMAKE_CXX_COMPILER=$CXX -DCMAKE_BUILD_TYPE=$BUILD_TYPE -DCMAKE_CXX_FLAGS="$COVERAGE $WERROR" ../source
    - docker exec dock make -j2
    - docker exec dock ctest
~~~
 
The command `docker exec dock` just sends a bash command to the docker image that we named "dock". The only noteworthy thing is the cmake call, where we pass our environmental variables which we set in the `matrix:` section to cmake. Additionally you should have a look at the end of the cmake call: `../source`. Here we specify the path to our source files, which we mounted to "/home/nuto/source". Since we are using ".." our current location inside the image is obviously a folder located inside "/home/nuto". This folder is "build". This is the predefined starting directory of the current nuto image. You have to be careful which command you call from which directory. Calling a command from the wrong directory might cause your tests to fail. You can always check what you are doing inside the docker container by pulling the image to your computer and copying the commands from your .travis.yml. This will be discussed in the following section.



**Additional notes**

The travis .travis.yml can contain additional section that are not needed here. You can find an overview here: 

https://docs.travis-ci.com/user/customizing-the-build/

If you make any changes to your `.travis.yml` and want to know if it is valid before pushing to your repository, you can copy and paste the code into the corresponding field on 

http://lint.travis-ci.org/

Or follow the instructions here:

https://docs.travis-ci.com/user/travis-lint

# Docker

Docker is a virtual machine that we use to run an image on travis. This way we avoid installing all the necessary components for NuTo over and over again during each test run on travis. This saves a lot of time but adds a little bit more complexity to the whole test environment. However, if your test fails and the reason is not obvious due to the output, you can easily reproduce the travis results on your own computer in order to find the exact issue. First you need to install docker

~~~bash
sudo apt-get install docker
~~~

Now you can simply copy the commands from the .travis.yml. Watch out that you might need to replace the `$(pwd)` with the location of the NuTo sources on your computer when starting the image. The results should be the same as on Travis. You can start modifying/adding/removing commands and see how you get your tests to pass. If you need to see what is actually going on inside the image, just replace the `itd` with `-it` and drop the `--name dock` when starting the image. This will send you to the bash inside of the image so that you can navigate through the folders etc. a little bit more comfortable. You can also drob `docker exec dock` in front of every command. To get back to your system just type

~~~bash
exit
~~~ 

For more information about docker, just google what you want to know. If you have configured your test environment you should not need that much more to know about it.


# Codecov

Codecov is a website that tells you how much of your code is actually used or tested. Especially if you are planning to provide some new features to NuTo that you first test in your private repository (for example integrands), you should make sure that you test everything well. Using codecov is quiet simple, but there are some mean pitfalls which can drive you crazy.
First thing you need to do is to register on 

https://codecov.io 

You can do this with your github account automatically. Everything else that you need to do is now uploading a report generated by a travis build. Therefore you add the following lines to your .travis.yml at the end of the script section:

~~~yml
    - if [ [ "$COVERAGE" == "--coverage" ] ]; then docker exec dock lcov --directory . --capture --output-file ../source/coverage.info; fi
    - if [ [ "$COVERAGE" == "--coverage" ] ]; then bash <(curl -s https://codecov.io/bash); fi
~~~

The `if [ [ "$COVERAGE" == "--coverage" ] ]; then ... ;fi` part of the code tells travis to only execute these lines if the environmental variable `$COVERAGE` is set to "--coverage". This has nothing to do with codecov itself. It just prevents your travis build from failing when no codecov data is generated. If you have a look at the CMake call in the `.travis.yml`, you will find this flag: `-DCMAKE_CXX_FLAGS="$COVERAGE $WERROR"`. So the variable is also passed to Cmake. CMake itself will pass it to the compiler. The flag "--coverage" tells the compiler to generate coverage data. Without it, there will be no data to process and the codecov commands will fail (and consequently your tests fail too). Hence we use the if block. 
Let's have a look at the actual command passed to the docker image: 

~~~bash
lcov --directory . --capture --output-file ../source/coverage.info
~~~

The program lcov just processes the coverage data generated by the programs. It is important to know, that you actually need to run your programs to generate the data. Compiling them is not enough! `--directory . --capture` tells lcov to gather the data in the current directory. '--output-file ../source/coverage.info' tells lcov to write the result file into your source directory (inside the docker image we are in "/home/nuto/build"). Now here are two possible pitfalls:
First we need root access inside the image as stated in the travis section above. This is done by starting the image with `-u 0`. The reason is, that our NuTo sources are just mounted to the directory "/home/nuto/source" and in our current image you don't have write permission to that directory when starting as the default user.
The second proble is, that you need to upload the coverage report from the source directory and not the build or any other directory. Otherwise the output from the codecov website might not give you the expected results.

The second command

~~~bash
bash <(curl -s https://codecov.io/bash)
~~~

does the actual upload to the codecov website. Notice that we don't use the docker image here! Since we have written the data to the shared directory it is not necessary! And thats basically it!

**Additional notes**


If you don't want certain folders to be included into your coverage report you can just add the following line between the other two commands:

~~~yml
- if [ [ "$COVERAGE" == "--coverage" ] ]; then docker exec dock lcov --remove ../source/coverage.info '${FOLDER_TO_EXCLUDE}' --output-file ../source/coverage.info; fi
~~~

Typical folders might be "usr/*" to exclude system libraries or "*/test*" to exclude your testfiles. The latter one won't remove the code from everything that is tested, just the code from your actual testfiles itself. Usually you are not interested of how much files in your testfile are run, just how many lines of your tested objects are hit.

If you want to configure the codecov report, you can add a `.codecov.yml` to your projects root directory. This way you can specify wich coverage (total and difference) is treated as a pass and when it is considered a fail. If it fails, github will mark the current commit with a red "X" even though the travis build was successful. Since the default setup produces fails from time to time based on to restrictive parameters I recommand using a `.codecov.yml`.

Here is some example code:

~~~yml
codecov:
    branch: master
    notify:
        require_ci_to_pass: yes

coverage:
    precision: 2
    round: down
    range: "70...100"

    status:
        project:
            default:
                enabled: yes
                target: 90%
                threshold: 5%


        patch:
            default:
                enabled: yes
                target: 90%
~~~

To understand the meaning of the different sections, just have a look here:

https://github.com/codecov/support/wiki/Codecov-Yaml

Notice that the code must be written like python code. So the indentation must be correct, otherwise the file is not valid. A quick check if your file is valid can be done by running 

~~~bash
curl --data-binary @.codecov.yml https://codecov.io/validate
~~~

inside the directory where your `.codecov.yml` is located.