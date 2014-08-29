#!/bin/bash
typeset -a array=('' '(1)' '(1,1)' '(1,1,1)' '(1,1,1,1)')
typeset -a arrayp=('' '(:)' '(:,:)' '(:,:,:)' '(:,:,:,:)')
typeset -a arrayy=('' '(lbound(what,1))' '(lbound(what,1),lbound(what,2))' '(lbound(what,1),lbound(what,2),lbound(what,3))' '(lbound(what,1),lbound(what,2),lbound(what,3),lbound(what,4))')
type=integer
kind=4
cat <<EOT
#ifndef CODE
interface RPN_COMM_ptr
#endif
EOT
for typex in integer real
do
for kind in 4 8
do

for dim in 1 2 3 4
do
  [[ $typex == i* ]] && type='i'
  [[ $typex == r* ]] && type='r'
  cat <<EOT
function RPN_COMM_ptr_${type%??????}${kind}_${dim}d(what) result(ptr)
use ISO_C_BINDING
implicit none
$typex(KIND=$kind), dimension${array[$dim]}, target :: what
type(C_PTR) :: ptr
#ifdef CODE
ptr = c_loc(what${array[$dim]})
return
#endif
end function RPN_COMM_ptr_${type%??????}${kind}_${dim}d

EOT
if [[ "1" == "2" ]] 
then
  cat <<EOT
function RPN_COMM_ptr_${type%??????}${kind}_${dim}dp(what) result(ptr)
use ISO_C_BINDING
implicit none
$typex(KIND=$kind), dimension${arrayp[$dim]}, pointer :: what
type(C_PTR) :: ptr
ptr = c_loc(what${arrayy[$dim]})
return
end function RPN_COMM_ptr_${type%??????}${kind}_${dim}dp

EOT
fi
done


done
done
cat <<EOT
#ifndef CODE
end interface RPN_COMM_ptr
#endif
EOT
