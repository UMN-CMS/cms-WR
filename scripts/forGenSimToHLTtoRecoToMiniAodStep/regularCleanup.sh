#!/bin/bash
#execute this script periodically to prevent core dumps and failed jobs creating large RAWSIM or GEN-SIM root files
#in local job submission directory from negatively affecting other running jobs

eval "rm core.*"
eval "rm WR_M-*.root"

