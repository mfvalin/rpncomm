#include <mpi.h>
typedef union {
  void *ftn;
  MPI_Datatype type;
  MPI_Group group;
  MPI_File file;
  MPI_Info info;
  MPI_Op op;
  MPI_Request request;
  MPI_Win win;
  MPI_Comm comm;
} F_C_union;

typedef struct {
  F_C_union all;
  int t1;
  int t2;
} F_C_univ ;

MPI_Fint RPN_COMM_comm_c2f( F_C_univ *in) {
  return MPI_Comm_c2f(in->all.comm);
}

F_C_union RPN_COMM_comm_f2c( F_C_univ *in) {
  F_C_univ temp;
  MPI_Fint comm = in->t2;
  temp.all.comm = MPI_Comm_f2c(comm);
  return temp.all;
}

MPI_Fint RPN_COMM_type_c2f( F_C_univ *in) {
  return MPI_Type_c2f(in->all.type);
}

F_C_union RPN_COMM_type_f2c( F_C_univ *in) {
  F_C_univ temp;
  MPI_Fint type = in->t2;
  temp.all.type = MPI_Type_f2c(type);
  return temp.all;
}

MPI_Fint RPN_COMM_group_c2f( F_C_univ *in) {
  return MPI_Group_c2f(in->all.group);
}

F_C_union RPN_COMM_group_f2c( F_C_univ *in) {
  F_C_univ temp;
  MPI_Fint group = in->t2;
  temp.all.group = MPI_Group_f2c(group);
  return temp.all;
}

MPI_Fint RPN_COMM_file_c2f( F_C_univ *in) {
  return MPI_File_c2f(in->all.file);
}

F_C_union RPN_COMM_file_f2c( F_C_univ *in) {
  F_C_univ temp;
  MPI_Fint file = in->t2;
  temp.all.file = MPI_File_f2c(file);
  return temp.all;
}

MPI_Fint RPN_COMM_info_c2f( F_C_univ *in) {
  return MPI_Info_c2f(in->all.info);
}

F_C_union RPN_COMM_info_f2c( F_C_univ *in) {
  F_C_univ temp;
  MPI_Fint info = in->t2;
  temp.all.info = MPI_Info_f2c(info);
  return temp.all;
}

MPI_Fint RPN_COMM_op_c2f( F_C_univ *in) {
  return MPI_Op_c2f(in->all.op);
}

F_C_union RPN_COMM_op_f2c( F_C_univ *in) {
  F_C_univ temp;
  MPI_Fint op = in->t2;
  temp.all.op = MPI_Op_f2c(op);
  return temp.all;
}

MPI_Fint RPN_COMM_request_c2f( F_C_univ *in) {
  return MPI_Request_c2f(in->all.request);
}

F_C_union RPN_COMM_request_f2c( F_C_univ *in) {
  F_C_univ temp;
  MPI_Fint request = in->t2;
  temp.all.request = MPI_Request_f2c(request);
  return temp.all;
}

MPI_Fint RPN_COMM_window_c2f( F_C_univ *in) {
  return MPI_Win_c2f(in->all.win);
}

F_C_union RPN_COMM_window_f2c( F_C_univ *in) {
  F_C_univ temp;
  MPI_Fint win = in->t2;
  temp.all.win = MPI_Win_f2c(win);
  return temp.all;
}

MPI_Fint RPN_COMM_c2f( F_C_univ *in, int what) {
  switch(what) {
    case 0:
      return MPI_Type_c2f(in->all.type);
      break;
    case 1:
      return MPI_File_c2f(in->all.file);
      break;
    case 2:
      return MPI_Group_c2f(in->all.group);
      break;
    case 3:
      return MPI_Info_c2f(in->all.info);
      break;
    case 4:
      return MPI_Op_c2f(in->all.op);
      break;
    case 5:
      return MPI_Request_c2f(in->all.request);
      break;
    case 6:
      return MPI_Win_c2f(in->all.win);
      break;
    case 7:
      return MPI_Comm_c2f(in->all.comm);
      break;
  }
}

void *RPN_COMM_f2c( MPI_Fint in, int what) {
  F_C_univ temp;
  switch(what) {
    case 0:
      temp.all.type = MPI_Type_f2c(in);
      break;
    case 1:
      temp.all.file = MPI_File_f2c(in);
      break;
    case 2:
      temp.all.group = MPI_Group_f2c(in);
      break;
    case 3:
      temp.all.info = MPI_Info_f2c(in);
      break;
    case 4:
      temp.all.op = MPI_Op_f2c(in);
      break;
    case 5:
      temp.all.request = MPI_Request_f2c(in);
      break;
    case 6:
      temp.all.win = MPI_Win_f2c(in);
      break;
    case 7:
      temp.all.comm = MPI_Comm_f2c(in);
      break;
  }
  return temp.all.ftn;
}
