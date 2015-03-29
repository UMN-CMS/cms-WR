#!/bin/bash
DOCCONF=fulldoc

if [ ! -d "doc/doxygen/${DOCCONF}/" ];then
    mkdir -p doc/doxygen/${DOCCONF}/
fi

if [ ! -d "doc/doxygen/${DOCCONF}/html" ];then
    cd doc/doxygen/${DOCCONF}/
    git clone -b gh-pages  git@github.com:UMN-CMS/cms-WR.git html
    cd -
else
    cd doc/doxygen/${DOCCONF}/html/
    if [ "`git branch | grep -c gh-pages`" == "0" ];then
	cd -
	rm doc/doxygen/${DOCCONF}/html/ -Rf
	cd doc/doxygen/${DOCCONF}/
	git clone -b gh-pages  git@github.com:UMN-CMS/cms-WR.git html
	cd -
    fi
fi

doxygen ${DOCCONF}

cd doc/doxygen/${DOCCONF}/html
ls
git remote -v 
git branch 
git pull
git add *.html
git add *.css
git add *.gif
git commit -m "updated documentation" -a
git commit -m "updated documentation" -a
git push origin gh-pages:gh-pages
cd -

