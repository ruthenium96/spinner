#!/bin/bash

function init(){
    docker build . -t july 
    docker run --name july -v $PWD:/project -d july
}

function clean(){
    docker container stop july
    docker container rm july
}

function shell(){
    docker exec --interactive --tty july /bin/bash
}

if [ $1 != "init" -a $1 != "clean" -a $1 != "shell" ]
    then
        echo "no such command" $1
        exit 1
fi

if [ $1 == "init" ]
    then
        init
fi

if [ $1 == "clean" ]
    then
        clean
fi

if [ $1 == "shell" ]
    then
        shell
fi

