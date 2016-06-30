#!/bin/bash
channel=$1
mode=$2
finalDir=$3
workingDir=$4
args=$5

mkdir -p $finalDir
cd $workingDir
eval `scramv1 runtime -sh`
eval "./bin/analysis -m $mode -c $channel  -d $finalDir $args"
