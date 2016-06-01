#!/bin/bash

#start in cmsWR dir
eval "cd test"
eval "sed 's@eleData@"

eval "root -l -b -q runCheckEleScalesSmearings.C"

