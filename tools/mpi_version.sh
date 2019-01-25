#!/bin/bash
[[ "$1" == -no* ]] && exit    # no version wanted
MpiVersion="$(mpirun --version 2>&1 | grep -w MPI | sed -e 's/.* //')"
printf "%s" _${MpiVersion:-NoNe}
