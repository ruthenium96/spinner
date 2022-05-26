#!/bin/bash

function init(){
    docker build . -t spinner
    docker run --name spinner --env CXX=clang++-13 -v $PWD:/project -d spinner
}

function clean(){
    docker container stop spinner
    docker container rm spinner
}

function shell(){
    docker exec --interactive --tty spinner /bin/bash
}

function shell_restart(){
    docker restart spinner
    docker exec --interactive --tty spinner /bin/bash
}

if [ $1 != "init" ] && [ $1 != "clean" ] && [ $1 != "shell" ] && [ $1 != "shell_restart" ]
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
