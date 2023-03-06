// MODULE specific to 2d age-kinematics inversion
// separated from the Qfunctions.i to avoir overcrowding
// BE CAREFUL THAT IT DOES NOT REDEFINE SOME FUNCTIONS (IT PROBABLY DOES THOUGH, via vpopcin2ls2.i instead of vpopcin2ls4.i)

#include "vpopcin2ls2.i"

func Pa(x,rra,&ag){
  /* DOCUMENT
     needs nv
  */
  ag=x*0.;
  res=0.;
  for(i=1;i<=nlos;i++){res+=P(x(i,),rra,g);ag(i,)=g;};
  return res;
};

func Pv(x,rrv,&vg){
  /* DOCUMENT
     needs nab
  */
  vg=x*0.;
  res=0.;
  for(i=1;i<=nab;i++){res+=P(x(,i),rrv,g);vg(,i)=g;};
  return res;
};

func Pav(x,rra,rrv,&g){
  /* DOCUMENT
     penalizing function for age-kinematics inversion
  */

  res=Pa(x,rra,ag)+Pv(x,rrv,vg);
  g=ag+vg;
  return res;
};
  

func AKNPEC(x,&g){
  /* DOCUMENT
     positive attitude version 2d age-kin inversion
     mono-metallic, but with NPEC
     derivative seems ok 
     NOTE: here Im applying the dust screen after the convolution, contrary to STECKMAP
     is it really what I want to do ?
     needs mtf
  */

  extern model1,rgl1,g1,dx1,wr1,pt,sp01,model;
  
//if (is_void(weight)) weight = 1.0;
  x=x;
  de=x(nab*nlos+1:);
  t=npe(de,nr,gt);
  
  x=x(:nab*nlos);
  x=reform(x,nlos,nab);
  rgl=deconv_mux*Pa(x^2,rr,ag)+deconv_muv*Pv(x^2,rr1,vg);
  dx=deconv_mux*ag+deconv_muv*vg;
  dx*=2.*x;
  //nnx=xtreat(x^2,nr+pad2+pad1);
  nnx=xtreat2(x^2,nr+pad2+pad1,nab);
  model = (1.0/numberof(mtf(,1)))*double(fft(mtf*fft(nnx,[-1,0]),[+1,0]));
  model=model(,sum);

  pt=pad(t,pad1,pad2);
  ipt=pt;
  //pt=pt(::-1);
  pt=roll(pt);
  sp01=pt*model;
  model1=sp01;
  residual = sp01 - bdata;
  //wr = weight*residual;
  wr = roll(pW(::-1))*residual;  // CHECK pW or roll(pW) ???

  //gx = pt(,-:1:nab)*(2.0/numberof(mtf(,1)))*double(fft(conj(mtf)*fft(wr,-1),[+1,0]));
  // keep in mind: in fourier space multiplication term to term is conlosolution in real space
  gx = roll(pt)*(2.0/numberof(mtf(,1)))*double(fft(conj(mtf)*fft(wr,-1),[+1,0]));
  gx=gx(pad2+ni:pad2+nj,);
  gx*=2.*x;
  gx+=dx;
  gx=gx(*);

  ge=(wr*model)(+)*(roll(bpad(gt,pad1,pad2)))(+,); ge=roll(ge);
  //ge=(wr*model)(+)*((bpad(gt,pad1,pad2)))(+,); ge=roll(ge);
  ge*=2.;
  
  grow,gx,ge;
  g=gx;

  // for debugging and plotting in main:
  model=sp01;
  
  if (debug==1) error;
  
  return sum(wr*residual) + rgl;
};

func AKNPEC2(x,&g){
  /* DOCUMENT
     positive attitude version 2d age-kin inlosersion
     argument passes as vector (nab*nlos + nde)
     mono-metallic, but with NPEC
     added normalization of NPEL (less expensive than normalization of v)
     NOTE: here Im applying the dust screen after the convolution, contrary to STECKMAP
     this makes sense as taking the NPEC as instrumental distortion
     needs cov and mucov;     
  */
  
  extern model1,rgl1,g1,dx1,wr1,pt;
 
 res=AKNPEC(x,g)+mucov*Pn(x(-nde+1:),cov,gcov);
 g(-nde+1:)+=mucov*gcov;
 return res;
};

func AKNPEC3(x,&g){
  /* DOCUMENT
     positive attitude version 2d age-kin inversion
     argument passes as vector (nab*nlos + nde)
     mono-metallic, but with NPEC
     added normalization of NPEL (less expensive than normalization of v)
     added penalization for e, controlled by mue and Le
     derivative is ok. 
     NOTE: here Im applying the dust screen after the convolution, contrary to STECKMAP
     this makes sense as taking the NPEC as instrumental distortion
     needs cov and mucov;     
  */
  
  extern model1,rgl1,g1,dx1,wr1,pt;
 
  res=AKNPEC2(x,g);
  res+=deconv_mue*P(x(-nde+1:),rr4,ge);
  g(-nde+1:)+=deconv_mue*ge;
 return res;
};



