#include "spline.i"

func spline2(y,x,xp){
  /*
    DOCUMENT
    spline interp of an array in its first direction only
    

  */

  local nl,na,nx,nxp,res;

  nl=dimsof(y)(2);
  nx=dimsof(x)(2);
  na=dimsof(y)(3);
  nm=dimsof(y)(4);
  nxp=numberof(xp);
  
  if(nx!=nl) {write,"GOTO HELL";return [];};
  res=array(0.,nxp,na,nm);

  
  for(i=1;i<=na;i++){
    for(j=1;j<=nm;j++){
      
      res(,i,j)=spline(y(,i,j),x,xp);
      
    };

  };
  return res;
};

func genVC(nab,w,ofs=){
  /* DOCUMENT
     generates a variance-covariance matrix
     a square matrix with a gaussian diagonal of width w (can be a vector => variable width) and ofset ofs (can also be avector)
  */

  if (numberof(w)==1) w=w(-:1:nab);
  if (is_void(ofs)) ofs=0.;
  if (numberof(ofs)==1) ofs=ofs(-:1:nab);
  res=array(0.,nab,nab);
  for(i=1;i<=nab;i++){res(,i)=makebump(nab,i+ofs(i),w(i),N=[]);};
  return res;
};


func fourop(n){

  /* DOCUMENT
     returns a Fourier operator of size n as defined in numerical recipes
  */

  res=array(0.,n,n);
  res=exp((2*1i*pi/n)*((indgen(n)-1)(,-:1:n))*((indgen(n)-1)(-:1:n,)));
  return res;
};
  
func trace(A){
  /* DOCUMENT
     returns the trace of square matrix A
  */

  return A(*)(::numberof(A(,1)));
};
