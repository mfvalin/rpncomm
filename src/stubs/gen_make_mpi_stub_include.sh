
cat <<EOT
program make_mpif_includes
include 'mpif.h'
EOT
for i in MPI_2DOUBLE_PRECISION MPI_2INTEGER MPI_2REAL MPI_DOUBLE_PRECISION MPI_INTEGER \
         MPI_INTEGER4 MPI_INTEGER8 MPI_LOGICAL MPI_REAL MPI_REAL4 MPI_REAL8 MPI_SUCCESS
do
  echo "print 100,'#define $i ',$i"
done
cat <<EOT
100 format(A,I10)
end program
EOT

