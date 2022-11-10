#!/bin/bash

set -eu
#e as soon as the one single command fails the whole script stops 
#u the script breaks when using undefined variable 
set -o pipefail 
# the whole pipe fails if there is erro in in piping
