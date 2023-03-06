#include "Bastien/string_utils.i"
#include "Pierre/POP/sfit.i"

func list_inc(file,n=){
  /* DOCUMENT
     works only for paths defined like #include "Pierre/POP/sfit.i"
  */

  if(is_void(n)) n=10000;
  a=open(file);
  u=rdline(a,n);
  iu=where((strmatch(u,"#include")==1)&(strmatch(u,"//")!=1));
  niu=numberof(iu);
  _list=[];
  //for(i=1;i<=niu;i++) grow,_list,strreplace(split2words(u(i),sep="/")(0),"\";","")(1);
  for(i=1;i<=niu;i++) grow,_list,strreplace(split2words(u(i),sep="\"")(-1),"\";","")(1);
  //  for(i=1;i<=niu;i++) grow,_list,strreplace(split2words(u(i),sep="Yorick")(-1),"\";","")(1);

  for(i=1;i<=numberof(_list);i++){grow,_list,list_inc(_list(i));};
  
  for(i=1;i<=niu;i++) write,split2words(u(i),sep="/")(0);
  
  return _list;
};
