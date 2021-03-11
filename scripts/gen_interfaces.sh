#! /usr/bin/env bash
# set -x

#SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

for target in ${@}; do 
  if grep -q '!InTfX!' ${target}; then
    base=${target##*/}
    echo "#if ! defined(IN_${base%.*})"
    cat ${target} | grep '!InTfX!' | sed 's/^!![ ]*/    /'
    echo "#endif"
  fi
done

echo
echo "  interface"
for target in ${@}; do 
  if grep -q '!InTf!' ${target}; then
    base=${target##*/}
    echo "#if ! defined(IN_${base%.*})"
    # grep '!InTf!' ${target} | sed -e 's/^\t//' -e 's/^      //' -e 's/^!!//' -e 's/ [ ]*!.*//' | sed 's/[ ]*$//' | sed 's/^/    /'
    grep '!InTf!' ${target} | sed -e 's/^[\t ]*!!//' -e 's/[ ]*!.*//' -e 's/[ ]*$//' -e 's/^[ \t]*/    /'
    echo "#endif"
    echo
  fi
done
echo "  end interface"
