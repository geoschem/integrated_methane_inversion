#!/bin/bash

du -h -s J_emis*/sat_obs.gosat.00.m | sort -n > rundir_check.log

exit 0