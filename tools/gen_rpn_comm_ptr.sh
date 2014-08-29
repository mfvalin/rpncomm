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
subroutine RPN_COMM_ptr_${type%??????}${kind}_${dim}d(what,ptr)
use ISO_C_BINDING
implicit none
include 'RPN_COMM_types.inc'
$typex(KIND=$kind), dimension${array[$dim]},intent(IN), target :: what
type(rpncomm_ptr), intent(OUT) :: ptr
#ifdef CODE
ptr%p = c_loc(what${array[$dim]})
return
#endif
end subroutine RPN_COMM_ptr_${type%??????}${kind}_${dim}d

EOT
if [[ "1" == "2" ]] 
then
  cat <<EOT
subroutine RPN_COMM_ptr_${type%??????}${kind}_${dim}dp(what,ptr)
use ISO_C_BINDING
implicit none
$typex(KIND=$kind), dimension${arrayp[$dim]}, pointer :: what
type(rpncomm_ptr), intent(OUT) :: ptr
ptr%p = c_loc(what${arrayy[$dim]})
return
end subroutine RPN_COMM_ptr_${type%??????}${kind}_${dim}dp

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
