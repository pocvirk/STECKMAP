// Tools for popcin
// Vectorized version of craptest.i





func pad(x,pad1,&pad2){

  local nr,pad2;
  nr=numberof(x(,1,1));
  if (mod(nr,2)!=0) pad2=pad1;
  if (mod(nr,2)==0) {
  pad2=pad1+1;
  write,"WARNING: spectrum size even => possible shift by fft convolution";
  };

  
  img2=[];
  img2=grow(transpose(indgen(pad2)(,-:1:na,)(,,-:1:nm)*0.),transpose(x));
  img2=grow(img2,transpose(indgen(pad2)(,-:1:na,)(,,-:1:nm)*0.));
  img2=transpose(img2);

  return img2;
};


func vdtreat(llos,xp0t,&xp01){

  local los;
  // Normalize losvd
  los=llos;
  
  los=los/(los(sum,,)(-:1:nr,,));


  // Resample losvd over spectrum range
 
   
  u3=array(0.,nr)(,-:1:na,)(,,-:1:nm);
  deltal=xd0(numberof(xd0))-xd0(1);
  xp01=span(-deltal/2.,deltal/2.,nr)(,-:1:na,)(,,-:1:nm);
  
  for(i=1;i<=na;i++){
    for(j=1;j<=nm;j++){
      
    u3(,i,j)=interp(los(,i,j),xp0t(,i,j),xp01(,i,j));
    };
  };
  

  // Renormalize losvd to have sum(u3)=1 when sum(los)=1
  // how to do this properly ?
  
  u3/=((xp0(nr,)-xp0(1,))/numberof(xp0(,1))*nr/deltal)(-:1:nr,);
  
  los=u3;

  return los;
};





func sptreat(u1,xd0,&xd01,&nr,opt=){
  /* DOCUMENT
     opt=1 => Normalization of all spectras
  */

  local nr,u2,u3,xd01;
  nr=int((xd0(numberof(xd0))-xd0(1))/min((x0(2)-x0(1))/x0));
  xd01=span(xd0(1),xd0(numberof(xd0)),nr);
  u2=interp(u1,xd0,xd01);
  if(!is_void(opt)) u2=u2/((u2(avg,,))(-:1:nr,,));
  return u2;
};

func xtreat(x,nx){
  local w;
  w=array(0.,nx)(,-:1:na)(,,-:1:nm);
  ind=indgen(nj-ni+1)+pad2+ni-1;
  w(ind,,)=x;
  return w;
};

  
  
