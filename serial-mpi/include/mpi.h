#define MPI_2DOUBLE_PRECISION         24
#define MPI_2INTEGER         25
#define MPI_2REAL         23
#define MPI_DOUBLE_PRECISION         17
#define MPI_INTEGER          7
#define MPI_INTEGER4         10
#define MPI_INTEGER8         11
#define MPI_LOGICAL          6
#define MPI_REAL         13
#define MPI_REAL4         14
#define MPI_REAL8         15
#define MPI_ANY_SOURCE       1000
#define MPI_ANY_TAG       1001
#define MPI_BYTE       1002
#define MPI_INTEGER1       1002
#define MPI_MAX       1003
#define MPI_MIN       1004
#define MPI_PACKED       1005
#define MPI_COMM_NULL       1006
#define MPI_GROUP_NULL       1007
#define MPI_UB       1008
#define MPI_LB       1009
#define MPI_COMM_WORLD       1010
#define MPI_CHARACTER       1011
#define MPI_SUM       1012
#define MPI_PROD       1013
#define MPI_DATATYPE_NULL       1014
#define MPI_LAND       1015
#define MPI_BAND       1016
#define MPI_REAL16       1017
#define MPI_LOR       1018
#define MPI_COMPLEX       1019
#define MPI_COMPLEX8       1020
#define MPI_COMPLEX16       1021
#define MPI_BOR       1022
#define MPI_COMPLEX32       1023
#define MPI_DOUBLE_COMPLEX       1024
#define MPI_LXOR       1025
#define MPI_LBOR       1026
#define MPI_BXOR       1027
#define MPI_MAXLOC       1028
#define MPI_MINLOC       1029
#define MPI_OP_NULL       1030
#define MPI_REPLACE       1031
#define MPI_WIN_NULL       1032
#define MPI_LOCK_EXCLUSIVE       1033
#define MPI_LOCK_SHARED       1034
#define MPI_INFO_NULL       1035
#define MPI_UNDEFINED       1036
#define MPI_IN_PLACE       1037


#define MPI_SUCCESS          0

typedef int MPI_Fint;

typedef int MPI_Datatype ;
typedef int MPI_Group ;
typedef int MPI_File ;
typedef int MPI_Info ;
typedef int MPI_Op ;
typedef int MPI_Request ;
typedef int MPI_Win ;
typedef int MPI_Comm ;

MPI_Comm MPI_Comm_f2c(MPI_Fint comm);
MPI_File MPI_File_f2c(MPI_Fint file);
MPI_Group MPI_Group_f2c(MPI_Fint group);
MPI_Info MPI_Info_f2c(MPI_Fint info);
MPI_Op MPI_Op_f2c(MPI_Fint op);
MPI_Request MPI_Request_f2c(MPI_Fint request);
MPI_Datatype MPI_Type_f2c(MPI_Fint datatype);
MPI_Win MPI_Win_f2c(MPI_Fint win);

MPI_Fint MPI_Comm_c2f(MPI_Comm comm);
MPI_Fint MPI_File_c2f(MPI_File file);
MPI_Fint MPI_Group_c2f(MPI_Group group);
MPI_Fint MPI_Info_c2f(MPI_Info info);
MPI_Fint MPI_Op_c2f(MPI_Op op);
MPI_Fint MPI_Request_c2f(MPI_Request request);
MPI_Fint MPI_Type_c2f(MPI_Datatype datatype);
MPI_Fint MPI_Win_c2f(MPI_Win win);

int MPI_Comm_rank(MPI_Comm comm, int *rank);
int MPI_Comm_size(MPI_Comm comm, int *size);
int MPI_Allgather(const void *sendbuf, int  sendcount, MPI_Datatype sendtype, const void *recvbuf, int recvcount, MPI_Datatype recvtype, MPI_Comm comm);
int MPI_Allreduce(const void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm);
int MPI_Barrier(MPI_Comm comm);
int MPI_Bcast(void *buffer, int count, MPI_Datatype datatype, int root, MPI_Comm comm);
