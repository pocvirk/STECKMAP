// Tools for popcin
// Vectorized version of craptest.i


func _pad(x,pad1,&pad2){

  local nr,pad2,na;
  nr=dimsof(x)(2);
  //if (!is_void(dimsof(x)(3))) na=dimsof(x)(3);
  if (mod(nr,2)!=0) pad2=pad1;
  if (mod(nr,2)==0) {
  pad2=pad1+0;
  //write,"WARNING: spectrum size even => possible half pixel shift by fft convolution";
  };

  
  img2=[];
  img2=grow(transpose(indgen(pad1)*0.),transpose(x));
  img2=grow(img2,transpose(indgen(pad2)*0.));
  img2=transpose(img2);

  return img2;
};

func pad(x,pad1,pad2){

  local nr,pad2;
  nr=numberof(x(,1));
  //if (mod(nr,2)!=0) pad2=pad1;
  //if (mod(nr,2)==0) {
  //pad2=pad1+1;
  //write,"WARNING: spectrum size even => possible shift by fft convolution";
  //};

  
  img2=[];
  grow,img2,array(0.,pad1),x,array(0.,pad2);
  return img2;
};



func _pad(x,pad1,pad2){

  /* DOCUMENT same as pad but no care taken about nr be even or not
   */

  local nr,na;
  nr=dimsof(x)(2);
  //write,"WARNING: spectrum size even => possible half pixel shift by fft convolution";

  
  img2=[];
  img2=grow(transpose(indgen(pad1)*0.),transpose(x));
  img2=grow(img2,transpose(indgen(pad2)*0.));
  img2=transpose(img2);

  return img2;
};




func padn1(x,pad1,pad2){

  local nr,na;
  nr=dimsof(x)(2);
  //if (!is_void(dimsof(x)(3))) na=dimsof(x)(3);
  
  
  img2=[];
  img2=grow(transpose(indgen(pad1)*0.),transpose(x));
  img2=grow(img2,transpose(indgen(pad2)*0.));
  img2=transpose(img2);

  return img2;
};


func padn2(x,pad1,pad2){

  local nr,pad2,na,img2;
  nr=dimsof(x)(2);
  na=dimsof(x)(3);
  

  
  img2=[];
  img2=grow(transpose(indgen(pad1)(,-:1:na)*0.),transpose(x));
  img2=grow(img2,transpose(indgen(pad2)(,-:1:na)*0.));
  img2=transpose(img2);

  return img2;
};


func padn3(x,pad1,pad2){

  local nr,na,img2;
  nr=dimsof(x)(2);
  na=dimsof(x)(3);
  nm=dimsof(x)(4);
  
  
  img2=[];
  img2=grow(transpose(indgen(pad1)(,-:1:na,)(,,-:1:nm)*0.),transpose(x));
  img2=grow(img2,transpose(indgen(pad2)(,-:1:na,)(,,-:1:nm)*0.));
  img2=transpose(img2);

  return img2;
};




func vdtreat(llos,xp0t,&xp01,norm=){

  local los;
  // Normalize losvd
  
  los=llos;
  if (!(norm==0)){
  los/=(los(sum,)(-:1:nr,));
  }

  // Resample losvd over spectrum range
 
   
  u3=array(0.,nr)(,-:1:na);
  deltal=xd0(0)-xd0(1);
  xp01=span(-deltal/2.,deltal/2.,nr)(,-:1:na);
  
  for(i=1;i<=na;i++){
   
      
    u3(,i)=interp(los(,i),xp0t(,i),xp01(,i));
   
  };
  

  // Renormalize losvd to have sum(u3)=1 when sum(los)=1
  // how to do this properly ?
  
  u3/=((xp0(nr,)-xp0(1,))/numberof(xp0(,1))*nr/deltal)(-:1:nr,);
  
  los=u3;

  return los;
};



func sptreat(u1,xd0,&xd01,&nr,opt=){
  /* DOCUMENT
     opt=1 => Normalization of all spectra
  */

  local nr,u2,u3,xd01;
  nr=int((xd0(0)-xd0(1))/min((x0(2)-x0(1))/x0));
  xd01=span(xd0(1),xd0(0),nr);
  u2=interp(u1,xd0,xd01);
  if(!is_void(opt)) u2=u2/((u2(avg,,))(-:1:nr,,));
  return u2;
};

func sptreat2(u1,xd0,x0,&xd01,&nr,opt=){
  /* DOCUMENT
     opt=1 => Normalization of all spectras
     created for sfit.i
  */

  local nr,u2,u3,xd01;
  nr=int((xd0(0)-xd0(1))/min((x0(2)-x0(1))/x0));
  xd01=span(xd0(1),xd0(0),nr);
  u2=interp(u1,xd0,xd01);
  if(!is_void(opt)) u2=u2/((u2(avg,,))(-:1:nr,,));
  return u2;
};

func sptreatn1(u1,xd0,&xd01,&nr,opt=){
  /* DOCUMENT
     opt=1 => Normalization of all spectra
  */

  local nr,u2,u3,xd01;
  nr=int((xd0(0)-xd0(1))/min((x0(2)-x0(1))/x0));
  xd01=span(xd0(1),xd0(0),nr);
  u2=interp(u1,xd0,xd01);
  if(!is_void(opt)) u2=u2/((u2(avg,,))(-:1:nr,,));
  return u2;
};



func xtreat1(x,nx){
  local w;
  w=array(0.,nx);
  ind=indgen(nj-ni+1)+pad2+ni-1;
  w(ind)=x;
  return w;
};

func xtreat2(x,nx,na){
  local w;
  w=array(0.,nx)(,-:1:na);
  ind=indgen(nj-ni+1)+pad2+ni-1;
  w(ind,)=x;
  return w;
};

func xtreat(x,nx){
  /* DOCUMENT 
     originally taken from vpopcin2ls.i to comply with 2d age-kni inversion
  */
  local w;
  w=array(0.,nx)(,-:1:nab)(,,-:1:nm);
  ind=indgen(nj-ni+1)+pad2+ni-1;
  w(ind,,)=x;
  return w;
};

  
  
