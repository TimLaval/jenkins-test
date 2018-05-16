#!/bin/bash

clone_abcd () {
	ssh-keygen  -t rsa -b 4096 -C "timothee.laval@barco.com"
	eval $(ssh-agent -s)
	ssh-add ~/.ssh/id_rsa
	clip < ~/.ssh/id_rsa.pub
	repository="git@github.com:TimLaval/gtest-jenkins.git"
	localfolder=="/abcd"
	git clone "$repository" "$localfolder"
}

do_test () {
	rm -rf abcd_build
	mkdir abcd_build
	cd abcd_build
	cmake /work/abcd/gtest-jenkins
	make test
}

do_documentation () {
	echo 'DEBUG documentation'
}

do_coverage () {
	echo 'DEBUG coverage'
}


clone_abcd
do_test
do_documentation
do_coverage
