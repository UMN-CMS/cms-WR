#!/bin/bash
DOCCONF=fulldoc

if [ ! -d "doc/doxygen/${DOCCONF}/" ];then
    mkdir -p doc/doxygen/${DOCCONF}/
fi
if [ ! -d "doc/doxygen/${DOCCONF}/html" ];then
    cd doc/doxygen/${DOCCONF}/
    git clone -b gh-pages  git@github.com:UMN-CMS/cms-WR.git html
    cd -
fi

doxygen ${DOCCONF}

cd doc/doxygen/${DOCCONF}/html
ls
git remote -v 
git pull
git add *.html
git add *.css
git add *.gif
git commit -m "updated documentation" -a
git commit -m "updated documentation" -a
git push origin gh-pages:gh-pages
cd -

