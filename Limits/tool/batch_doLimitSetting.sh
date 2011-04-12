#!/bin/sh

installLoc=$1
shift

source ${installLoc}/setuproot.sh

exec ${installLoc}/doLimitSetting.exe $*

