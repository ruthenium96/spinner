#!/bin/bash

function init(){
    docker build . -t july 
    docker run --name july --env CXX=clang++-13 -v $PWD:/project -d july
}

function clean(){
    docker container stop july
    docker container rm july
}

function shell(){
    docker exec --interactive --tty july /bin/bash
}

function shell_restart(){
    docker restart july
    docker exec --interactive --tty july /bin/bash
}

if [ $1 != "init" -a $1 != "clean" -a $1 != "shell" -a $1 != "shell_restart" ]
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

if [ $1 == "shell_restart" ]
    then
        shell_restart
fi
