#!/bin/bash
set -ev
# Automatically build and upload doxygen documentation to Github pages
# inspired by https://gist.github.com/vidavidorra/548ffbcdae99d752da02
REPO=git@github.com:nutofem/nuto.git
# only upload the master documentation to github pages
if [[ "$TRAVIS_BRANCH" == "master" ]]
then
    echo Generating the github documentation

    openssl aes-256-cbc -K $encrypted_7c7dcd9a49b0_key -iv $encrypted_7c7dcd9a49b0_iv -in scripts/traviskey.enc -out traviskey -d
    chmod 600 traviskey
    eval `ssh-agent -s`
    ssh-add traviskey

    cd build
    make doc > /dev/null

    # Get the current gh-pages branch
    git clone --depth 1 -b gh-pages $REPO
    cd nuto

    # Configure git
    git config --global push.default simple
    git config user.name "Travis CI"
    git config user.email "nuto@nuto.fem"

    # Remove all the old files
    rm -rf *

    # Tell Github it should not use Jekyll on this branch
    touch .nojekyll

    # Copy the generated documentation to here
    cp -r ../doc/nuto/html/* .

    echo 'Uploading documentation to the gh-pages branch...'
    git add --all

    git commit -m "Deploy code docs to GitHub Pages Travis build: ${TRAVIS_BUILD_NUMBER}" -m "Commit: ${TRAVIS_COMMIT}" > /dev/null

    git push $REPO gh-pages
fi
