#!/bin/bash
if [[ ! -x ./wrap_code.exe ]] ; then
cat <<EOT >./wrap_code.f90
program wrap
  implicit none
  character (len=256) :: line, temp
  character (len=6) :: spaces='      '
  character (len=6) :: cont  ='     &'
  character (len=1) :: cont2 ='&'
  write(6,100)spaces//'interface'
1 read(5,100,end=2) line
    if(len(trim(line))<=66) then
      write(6,101)spaces,line(1:66)
    else
      write(6,101)spaces,line(1:66),cont2
      temp=line(67:)
      line=temp
      do while(len(trim(line))>66)
        write(6,101)cont,line(1:66),cont2
        temp=line(67:)
        line=temp
      enddo
      write(6,101)cont,line(1:66)
    endif
  goto 1
2 continue
  write(6,100)spaces//'end interface'
  stop
100 format(A)
101 format(A6,A66,A1)
  end
EOT
s.f90 -o ./wrap_code.exe ./wrap_code.f90 2>/dev/null 1>/dev/null
rm -f ./wrap_code.f90
fi
grep '!InTf!' | sed -e 's/^\t//' -e 's/^      //' -e 's/^!!//' -e 's/!.*//' | ./wrap_code.exe
