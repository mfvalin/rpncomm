#!/bin/bash
mkdir -p ../tools
if [[ ! -x ../tools/wrap_code.exe ]] ; then
cat <<EOT >../tools/wrap_code.f90
program wrap
  implicit none
  character (len=256) :: line, temp
  character (len=6) :: pre
  write(6,100)'      interface'
1 read(5,100,end=2) line
    if(line(1:1)=='!') then
      write(6,100)line
    else
      pre = '      '
      do while(len(trim(line))>66)
        write(6,100)pre//line(1:66)//'&'
        temp=line(67:) ; line=temp ; pre = '     &'
      enddo
      write(6,100)pre//trim(line)
    endif
  goto 1
2 continue
  write(6,100)'      end interface'
  stop
100 format(A)
  end
EOT
#s.f90 -o ../tools/wrap_code.exe ../tools/wrap_code.f90 2>/dev/null 1>/dev/null
gfortran -o ../tools/wrap_code.exe ../tools/wrap_code.f90 2>/dev/null 1>/dev/null
rm -f ../tools/wrap_code.f90
fi
if [[ -z $1 ]] ; then
  grep '!InTf!' | sed -e 's/^\t//' -e 's/^      //' -e 's/^!!//' -e 's/ [ ]*!.*//' | ../tools/wrap_code.exe | sed 's/[ ]*$//'
else
  [[ -f $1 ]] || exit 0
  if grep -q '!InTf!' $1 ; then
    echo "#if ! defined(IN_${1%.*})"
    grep '!InTf!' $1 | sed -e 's/^\t//' -e 's/^      //' -e 's/^!!//' -e 's/ [ ]*!.*//' | ../tools/wrap_code.exe | sed 's/[ ]*$//'
    echo "#endif"
  fi
fi
