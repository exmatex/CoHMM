#!/bin/bash
#dir=${PWD##*/}
dir=${PWD}
cd $dir
rm -r cn*
sleep 1
master=
[[ ! -f $i ]] && echo "No hostfile given!" >&2 && exit 1
for i in $(cat "$1"); do
    mkdir "$i"
    cd "$i"
    if [[ -z $master ]]; then
        nohup redis-server &
        #nohup redis-cli flushdb &
        master=$i
    else
        #echo "port 6789" > redis.conf
        {
            echo "port 6379"
            echo "slaveof $master 6379"
            echo "dbfilename dump$i.rdb"
            echo "slave-read-only no"
        } > redis.conf
        nohup ssh -n $i "cd $dir/$i; $(type -p redis-server) ./redis.conf " &
        echo "started slave on $i"
    fi
    sleep 2
    cd -
done

