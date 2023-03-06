// Several objective functions (wiht/without kinematics, various treatment of extinction: parametric (Calzetti ?) non-parametric, age-dependent... all kinds of stuff.

// ADDING a penalization to allow for user-supplied AMR. Best way to do it would be in the objective functions, and creating a new P.

#include "random.i"
#include "regul.i"
#include "optim_pack.i"

// needs bloc,x0, and much more : see Wqcnfaz

//func Q0(x,&g){// unclipped
//  return Wfaz1cc(x,g); };
//func Q1(x,&g){// clipped
//  return Wfaz1cc(x,g,pos=2);  };

func Zrescale(z){
  // Rescales z (to match other parameters scales before minimization)
  return (-1.e-1*log10(z/0.1));
};

func Zrescalem1(y){
  // inverse of Zrescale (back to physical z)
  return (10^(-10.*y-1));
};

func LWA(nhop,a1,a2,b=){
  /* DOCUMENT
     computes luminosity weighted age of components between a1 and a2
     needs ta and ab OR specification of a basis b ** PREFERRED **
     note: SAD is not squared inside LWA. Do it prior to applying LWA
     ages are logarithmically averaged:
     LWA=10^sum (nhop(i)*log10(b.ages))
  */
  write,dimsof(b.flux)
  if(!is_void(b)) {ta=b.ages;ab=indgen(numberof(nhop));};
  if (a1==a2) return ta(a1);
  nhop=max(nhop,0.);
  return 10^(((log10(ta(ab(a1:a2)))*nhop(a1:a2,))(sum,))/(nhop(a1:a2,)(sum,)));
};

func LWM(nhop,a1,a2){
  /* DOCUMENT
    same as LWA for metallicity
    note: SAD part in nhop is not squared
    NOTE: It should make a big difference wether you have rescaled or absolute or log metallicities in nhop!!!!
    seems to be vectorized, didnt test it though
  */
  //  nhop(1:nab)=max(nhop(1:nab),0.);  removed 11/08/06
  //  return ((nhop(nab+a1:nab+a2,)*nhop(a1:a2,))(sum,))/(nhop(a1:a2,)(sum,));
  return ((nhop(nab+a1:nab+a2)*nhop(a1:a2))(sum))/(nhop(a1:a2)(sum));
};

func bpad(sp,pad1,pad2){
  /* DOCUMENT
     pad a 2d basis
  */
  local na;
  na=dimsof(sp)(3);
  res=array(0.,numberof(sp(,1))+pad2+pad1,na);
  for(i=1;i<=na;i++){
    res(,i)=pad(sp(,i),pad1,pad2);
  };
  return res;
};



func chi2r(d,y,snr,w=){
  /* DOCUMENT
     computes reduced chi^2 of d fitted by y given snr, w is the diagonal of the inverse variance-covariance matrix
  */
  local nl;
  
  nl=numberof(d);
  if(is_void(w)) w=((snr/d)^2)/nl;
  return (((d-y)^2)*w)(sum);
};

func sqd(d,y,snr,w=){
  /* DOCUMENT
     computes the square of the euclidian distanc between d and y
  */
  local nl;
  
  return ((d-y)^2)(sum);
};


func SAD2MASS(SAD,AMR){
  /* DOCUMENT
     converts flux fractions into mass fractions
  */
  na=numberof(SAD);
  MF=[];
  //  for(i=1;i<=na;i++){
  //    grow,MF,SAD(i)*MsLinterp(AMR);
  //};
  return SAD*MsLinterp(AMR);
};

func SAD2SFR2(SAD,AMR,ba,&bab,lo=,hi=){
  /* DOCUMENT
     need to provide timebins extremities through lo and hi
     defaut is lo=10Myr and hi=15 Gyr
     check bab to make sure the intervals arent crappy
  */

  if(is_void(lo)) lo=10.;
  if(is_void(hi)) hi=18000.;
  bab=lo;
  grow,bab,ba(zcen);
  grow,bab,hi;
  baba=bab(dif);
  return SAD2MASS(SAD,AMR)/baba;
};

  


func SAD2SFR(hop,z,bab,N=){
  /* DOCUMENT
     converts SAD 2 SFR
  */
  na=numberof(bab)-1;
  nrea=(numberof(dimsof(hop))==3)?dimsof(hop)(3):1;
  rhop=hop;
  SF=rhop*0.;
  nrhop=rhop/rhop(sum,)(-:1:na,);
  for(i=1;i<=nrea;i++){SF(,i)=(MsLinterp(z(,i))*nrhop(1:na,i))/bab(dif);}
  return SF;
};


func SFR2SAD(SF,z,bab,N=){
  // converts SFR 2 SAD
  na=numberof(bab)-1;
  hop=SF;
  rhop=hop;
  nrea=(numberof(dimsof(hop))==3)?dimsof(hop)(3):1;
  SAD=rhop*0.;
  nrhop=rhop/rhop(sum,)(-:1:na,);
  for(i=1;i<=nrea;i++){SAD(,i)=(1./((MsLinterp(z(,i))/bab(dif))))*nrhop(1:na,i);}
  return SAD;
};


  

func _MsLinterp(a,_z){
  // INTERPOLATES MsL values. requires MsL array
  
  zmax=1.e1*max(_m);
  zmin=1.e-1*min(_m);
  _z=max(_z,zmin);
  _z=min(_z,zmax);

  iz1=is_void(where(_z>=_m))?where(_z>=_m)(1):1;
  iz2=is_void(where(_m>=_z))?where(_m>=_z)(1):numberof(_m);
  
  A=(_m(iz1)-_z)/(_m(iz1)-_m(iz2));
  B=(-_m(iz2)+_z)/(_m(iz1)-_m(iz2));
  res=A*MsL(a,iz2)+B*MsL(a,iz1);
  return res;
};

func MsLinterp(z){
  b=array(0.,nab);
  for(i=1;i<=nab;i++){
    b(i)=_MsLinterp(i,z(i));
  };
  return b;
};




//func genREGUL(ni) {
  // generates regularization matrix D3
//      res = array(0.,ni-3,ni);
//for(i=1;i<=ni-3;i++){
//res(i,i)=1;res(i,i+1)=-3;res(i,i+2)=3;res(i,i+3)=-1;}
//return res;
//};

func genREGUL(ni,s=){
  if (is_void(s)) s="D2";
  
  if (s=="I") {
    res = diag(array(1.,ni));};
  
  if (s=="D1") {
    res = array(0.,ni-1,ni);
    for(i=1;i<=ni-1;i++){
      res(i,i)=-1.;res(i,i+1)=1.;};
  };
  
  if (s=="D2") {
    res = array(0.,ni-2,ni);
    for(i=1;i<=ni-2;i++){
      res(i,i)=-1.;res(i,i+1)=2.;res(i,i+2)=-1.;};
  };
  
  if (s=="D3") {
    res = array(0.,ni-3,ni);
    for(i=1;i<=ni-3;i++){
      res(i,i)=1;res(i,i+1)=-3;res(i,i+2)=3;res(i,i+3)=-1;};
  };
  

  return res;
};




func P(x,rr,&g){
  /* DOCUMENT
     classical tikhonov regul by kernel rr so that P(x)=xtRx and gradent.
  */
  
  local g;
  g=2.*(rr(,+)*x(+));
  return (x(+)*(rr(,+)*x(+))(+));
};

func logP(x,rr,&g){
  /* DOCUMENT
     classical regul by kernel rr and gradient, except
  */
}


func __P(x,rr,&g,rrp=,mup=,prior=){
  /* DOCUMENT
     classical tikhonov regul by kernel rr so that P(x)=xtRx and gradent.
     optional arguments prior is a prior, same size as x, rrp is the kernel definig the scalar product between x and the prior. Identity is default, mup is the weight of the prior distribution.
  */
  
  local g;

//if(!is_void(mup)&is_void(rrp)) rrp=diag(array(1.,numberof(x)));
  gp=0.;
  res=0.;
  if (is_void(mup)) mup=0.;
  if(!is_void(prior)){
    res=mup*(((x-prior)(+)*rrp(,+))(+)*(x-prior)(+));
    gp=mup*2.*(rrp(,+)*(x-prior)(+));
  };
    
  g=2.*(rr(,+)*x(+))+gp;
  return (x(+)*(rr(,+)*x(+))(+))+res;
};


func _P(x,L,&g){
  /* DOCUMENT
     classical tikhonov regul by kernel rr so that P(x)=xtRtRx and gradent.
  */
  
  local g;
  g=2.*((L(+,)*L(+,))(,+)*x(+));
  return (x(+)*((L(+,)*L(+,))(,+)*x(+))(+));
};

func Pn(x,c,&g){
  /* DOCUMENT
     minimizing the objective function + Pn will provide a solution with the constraint sumi xi =c
  */

  g=(2.*(x(sum)-c))(-:1:numberof(x));
  return (x(sum)-c)^2;

};
  



func _spinterp(a,_z,&g){
  /* DOCUMENT
     WARNING! should be better than linear interpolation
     derivative is not continousou
     old version
  */


  //_z=max(_z,0.0004);
  //_z=min(_z,0.05);
  zmax=1.e1*max(_m);
  zmin=1.e-1*min(_m);

  //added 23may05
  {
    //zmax=_m(2);
    //zmin=_m(-1);
  };
  
  _z=max(_z,zmin);
  _z=min(_z,zmax);
  //iz1=(where(_z>=_m))(1);
  //iz2=iz1-1;

  if(backup=="yes"){
  iz1=is_void(where(_z>=_m))?where(_z>=_m)(1):1;  // AARRGGHH
  iz2=is_void(where(_m>=_z))?where(_m>=_z)(1):numberof(_m); // AAARRGGGHH too
  };

  // that was added 23may05
  {
  iz1=is_void(dimsof(where(_z>=_m)))?1:where(_z>=_m)(1);   
  iz2=is_void(dimsof(where(_m>_z)))?numberof(_m):where(_m>_z)(0);
  };


  //write,_z;
  //write,iz1,iz2;
  
  A=(_m(iz1)-_z)/(_m(iz1)-_m(iz2));
  B=(-_m(iz2)+_z)/(_m(iz1)-_m(iz2));
  res=A*bloc(,a,iz2)+B*bloc(,a,iz1);
  g=(1./(_m(iz2)-_m(iz1)))*(bloc(,a,iz2)-bloc(,a,iz1));
  return res;
};

func spinterp(a,_z,&g){
  /* DOCUMENT
     WARNING! could be better than linear interpolation
     derivative is not continous
     CHECK THIS FOR THE NEXT TIME I CAN CODE
  */


  //_z=max(_z,0.0004);
  //_z=min(_z,0.05);
  zmax=1.e0*max(_m)-1.e-5;
  zmin=1.e0*min(_m)+1.e-5;

  //added 23may05
  {
    //zmax=_m(2);
    //zmin=_m(-1);
  };
  
  _z=max(_z,zmin);
  _z=min(_z,zmax);
  iz1=is_void(dimsof(where(_z>=_m)))?1:where(_z>=_m)(0);   
  iz2=is_void(dimsof(where(_m>_z)))?numberof(_m):where(_m>_z)(1);

  //  write,_m;
  //write,is_void(dimsof(where(_z>=_m)));
  //write,where(_z>=_m);
  //write,_z;
  //write,iz1,iz2;
  
  A=(_m(iz1)-_z)/(_m(iz1)-_m(iz2));
  B=(-_m(iz2)+_z)/(_m(iz1)-_m(iz2));
  res=A*bloc(,a,iz2)+B*bloc(,a,iz1);
  g=(1./(_m(iz2)-_m(iz1)))*(bloc(,a,iz2)-bloc(,a,iz1));
  return res;
};





func spinterpa(_a,_z,&ga,&gz){
  /* DOCUMENT
     WARNING!
     same warnings as spinterp
     derivative seems ok, as long as not on the sample points
     DERIVATIVE IS OK EXCEPT ON SAMPLE POINTS

     
  */
  
  amax=ta(-1);
  amin=ta(2);
  _a=max(_a,amin);
  _a=min(_a,amax);

  ia1=is_void(dimsof(where(_a<=ta)))?numberof(ta):where(_a<=ta)(1);
  ia2=is_void(dimsof(where(ta<_a)))?1:where(ta<_a)(0);

  //write,where(_a<=ta);
  
  A=(ta(ia1)-_a)/(ta(ia1)-ta(ia2));
  B=(-ta(ia2)+_a)/(ta(ia1)-ta(ia2));
  ias1=spinterp(ia1,_z,gz1);
  ias2=spinterp(ia2,_z,gz2);
  res=A*ias2+B*ias1;
  //  write,A,B;
  ga=(1./(ta(ia2)-ta(ia1)))*(ias2-ias1);
  gz=A*gz2+B*gz1;
  return res;
};

func spinterpatest(_a,_z,&ga,&gz){
  /* DOCUMENT
     WARNING!
     same warnings as spinterp
     derivative seems ok, as long as not on the sample points
     DERIVATIVE IS OK EXCEPT ON SAMPLE POINTS
     requires ta,
     
  */
  
  amax=ta(-1);
  amin=ta(2);
  _a=max(_a,amin);
  _a=min(_a,amax);

  ia1=is_void(dimsof(where(_a<=ta)))?numberof(ta):where(_a<=ta)(1);
  ia2=is_void(dimsof(where(ta<_a)))?1:where(ta<_a)(0);

  //write,where(_a<=ta);
  
  A=(ta(ia1)-_a)/(ta(ia1)-ta(ia2));
  B=(-ta(ia2)+_a)/(ta(ia1)-ta(ia2));
  ias1=spinterp(ia1,_z,gz1);
  ias2=spinterp(ia2,_z,gz2);
  res=A*ias2+B*ias1;
  //  write,A,B;
  ga=(1./(ta(ia2)-ta(ia1)))*(ias2-ias1);
  gz=A*gz2+B*gz1;
  return res;
};





func buildb2(z,&g){
  b=array(0.,dd,nab);
  //b=rspb*0.;
  g=b;
  for(i=1;i<=nab;i++){
    b(,i)=(spinterp(i,z(i),spg));
    g(,i)=(spg);
  };
  return b;
};

func bound2(z,&g){
  // quadratic out of bounds, flat inside
  local g;
  g=0.*z;
  res=0.;
  for(i=1;i<=numberof(z);i++){
    res+=_bound2(z(i),b2g);
    g(i)=b2g;
  };
  return res;
};

func _bound2(z,&g){
  //z1=0.0008;
  //z2=0.01;
  // log Z values  (_m=-log10(_m/0.1))
  //z1=0.1e-1;
  //z2=2.9e-1;
  // FIX boundaries from _m

  local g;
  
  z1=is_void(zlim)?min(base.met):zlim(1);  
  z2=is_void(zlim)?max(base.met):zlim(2);
  
    b=(z<=z1)?((z1-z))^2:(z>=z2)?(z-z2)^2:0;
  g=(z<=z1)?(-2.*(z1-z)):(z>=z2)?(2.*(z-z2)):0;
  return b;};

func _abound2(a,&g){
  //z1=0.0008;
  //z2=0.01;
  // log Z values  (_m=-log10(_m/0.1))
  //z1=0.1e-1;
  //z2=2.9e-1;
  // FIX boundaries from _m

  local g;
  
  a1=is_void(alim)?min(base.met):alim(1);  
  a2=is_void(alim)?max(base.met):alim(2);
  
    b=(a<=a1)?((a1-a))^2:(a>=a2)?(a-a2)^2:0;
  g=(a<=a1)?(-2.*(a1-a)):(a>=a2)?(2.*(a-a2)):0;
  return b;};



func ke(l){
  // ke for a starburst galaxy from calzetti et al. l: lambda in microns
  return 1.17*(-1.857+1.04/l)+1.78;};

func keUV(l){
  // ke in the UV for a starburst galaxy from calzetti et al. l: lambda in microns
  return 1.17*(-2.156+1.509/l-0.198/l^2+0.011/l^3)+1.78;};


func ds(x,x0,&g){
  // derivative ok
  // WARNING!! x is E(B-V)gas
  // and E(B-V)star=0.44E(B-V)gas in Calzetti 2001
  // uses ke
  local nl;
  x=max(x,-145.);   //   to avoid 
  x=min(x,5.e1);    //       pow errors
  nl=numberof(x0);
  t=10^(-0.4*x*ke(x0/1.e4));
  N=t(sum)/nl;
  res=t/N;
  dt_de=-0.4*ke(x0/1.e4)*t*log(10);
  dN_de=dt_de(sum)/nl;
  g=-(1./N^2)*(dN_de*t-N*dt_de);
  return res;};

func dsUV(x,x0,&g){
  // same as ds for UV extinction law
  // derivative ok
  local nl;
  x=max(x,-145.);   //   to avoid 
  x=min(x,5.e1);    //       pow errors
  nl=numberof(x0);
  t=10^(-0.4*x*keUV(x0/1.e4));
  N=t(sum)/nl;
  res=t/N;
  dt_de=-0.4*keUV(x0/1.e4)*t*log(10);
  dN_de=dt_de(sum)/nl;
  g=-(1./N^2)*(dN_de*t-N*dt_de);
  return res;};


func fds(x0,x,&grad,deriv=){
  /* DOCUMENT
     used to find E(B-V) from non-parametric exintction law
  */
    res=ds(x,x0,g);
  if(deriv==1) grad=[g];
  return res;
};

func fdsUV(x0,x,&grad,deriv=){
  /* DOCUMENT
     same as fds in the UV
   */
    res=dsUV(x,x0,g);
  if(deriv==1) grad=[g];
  return res;
};

func vds(x,x0,&g){
  // derivative ok
  local na,nl;
  na=numberof(x);
  nl=numberof(x0);
  x=max(x,-145.);   //   to avoid 
  x=min(x,5.e1);    //       pow errors
  t=10^(-0.4*x(-:1:nl,)*ke(x0/1.e4)(,-:1:na));
  N=t(sum,)/nl;
  res=t/N(-:1:nl,);
  dt_de=-0.4*(ke(x0/1.e4)(,-:1:na))*t*log(10);
  dN_de=dt_de(sum,)/nl;
  g=-(1./(N^2)(-:1:nl,))*(dN_de(-:1:nl,)*t-N(-:1:nl,)*dt_de);
  return res;
};

func GENdecdde(nde,nx0){
  /* DOCUMENT
     computes an array necessary for the computation of the derivative of npe
     note: would be nice to do it wiht cubic splines rahter than linear interpolation
     needs external paramter bnpec
  */
  local q;
  decdde=array(0.,nx0,nde);
  for(i=1;i<=nde;i++){
    q=array(0.,nde);
    q(i)=1.
    if(bnpec=="tri")  decdde(,i)=interp(q,indgen(nde),span(1,nde,nx0));   // "triangle" functions
    if(bnpec=="spl") decdde(,i)=spline(q,indgen(nde),span(1,nde,nx0));   // proper cubic spline
    
  };
  return decdde;
};

func GENdx0dx1(x0,x1){
  /* DOCUMENT
     generates dx0dx1 such that plh,sp0,x0, and plh,sp0(+)*dx0dx1(,+),x1 are the same.
     REMARK: we also have:dx0dx1(,+)*x1(+)=x0 (nearly: only along the smallest domain)
     needs external paramter bnpec


  */
  local q;
  nx0=numberof(x0);
  dx0dx0=array(0.,nx0,nx0);
  for(i=1;i<=nx0;i++){
    q=array(0.,nx0);
    q(i)=1.;
    if (bnpec=="tri") dx0dx0(,i)=interp(q,x0,x1);
    if (bnpec=="spl") dx0dx0(,i)=spline(q,x0,x1);
  };
  return dx0dx0;
};


func npe(de,nx0,&g){
  /* DOCUMENT
     npe for Non Parametric Extinction
     output is a non-parametric extinction curve with numberof(de) anchor points regularly spaced along the wavalength axis x0
     NOTE: creates decdde if it doesnot exist with GENdecdde
     nx0=numberof(x0)
     derivative checked and seems ok
     needs external paramter bnpec
  */

  local nde;
  nde=numberof(de);
  if (is_void(decdde)) decdde=GENdecdde(nde,nx0);
  if ((dimsof(decdde)([2,3])==[nx0,nde])(sum)!=2)
    {
    write,"dimensions changed, recreate decdde";
    error;
  };
  g=decdde;
  if(bnpec=="tri") return interp(de,double(indgen(nde)),span(1.,double(nde),nx0));
  if(bnpec=="spl") return spline(de,double(indgen(nde)),span(1.,double(nde),nx0));
  
};





func polf(c,x0,&g){
  local d,res;
  d=numberof(c);
  res=x0*0.;
  p=  g=res(,-:1:d); 
  for(i=1;i<=d;i++){res+=c(i)*x0^(i-1);g(,i)=x0^(i-1);};
  return res;
};

func _polf(c,x0,&g){
  local d,res;
  d=numberof(c);
  res=x0*0.;
  p=  g=res(,-:1:d); 
  for(i=1;i<=d;i++){res+=c(i)*x0^i;g(,i)=x0^(i-1);};
  return res;
};


func npolf(c,x0,&g,&nc){
  local d,res;
  d=numberof(c);
  res=x0*0.;
  g=res(,-:1:d);
  N=0.;
  for(i=1;i<=d;i++){ N+=(c(i)/i)*(((x0(0))^(i))-x0(1)^(i)) ;};
  // c(i)/i and not c(i)/i+1 because c(1) is a0
  //N=N/(x0(0)-x0(1));
  nc=(x0(0)-x0(1))*c/N;
  res=polf(nc,x0);
  for(i=1;i<=d;i++){g(,i)=(x0(0)-x0(1))*(x0^(i-1))*(N-c(i)*((x0(0)^i-x0(1)^i)/i))/N^2;};
  return res;
};

func npolf2(c,x0,&g,&nc){
  /* DOCUMENT
     returns the value of the normalized polynomial a0 + sum (c>=1) ci x0^i with a0 defined so that this polynomial is effectively normalized. c begins with coefficient a1.
  */
  local d,res,nc;
  d=numberof(c);
  res=x0*0.;
  g=0.*res(,-:1:d);
  dc0dc=c*0.;
  c0=0.;
  for(i=1;i<=d;i++){ c0+=(c(i)/(i+1))*(x0(0)^(i+1)-x0(1)^(i+1));};
  c0*=-1./(x0(0)-x0(1));
  for(i=1;i<=d;i++){dc0dc(i)=-(x0(0)-x0(1))*(c(i)/(i+1))*(x0(0)^(i+1)-x0(1)^(i+1));};
  nc=c0;
  grow,nc,c;
  res=polf(nc,x0,gp);

  for(i=1;i<=d;i++){g(,i)=gp(,1)*dc0dc(i) + gp(,i+1);};

  return res;
};

  func npolf3(c,x0,&g,&nc){

  // DERIVATIVE OK // but possible degeneracy: different coefficients can give same output
  
  local d,res;
  d=numberof(c);
  res=x0*0.;
  g=res(,-:1:d);
  N=0.;
  for(i=1;i<=d;i++){ N+=(c(i)/i)*(((x0(0))^(i))-x0(1)^(i)) ;};
  //N=N/(x0(0)-x0(1));
  nc=(x0(0)-x0(1))*c/N;
  res=polf(nc,x0);
  for(j=1;j<=d;j++){
    g(,j)=(x0(0)-x0(1))*(x0^(j-1))/N;
    for(i=1;i<=d;i++){
      g(,j)+=-(x0(0)-x0(1))*((c(i)/(j*N^2))*(x0(0)^j-x0(1)^j)*x0^(i-1));}
  };
  return res;
};


func npenalty1(x,&g){
  /* DOCUMENT
     penalty1 revisited for the needs of LF inversion
  */
  //  res=0.5*((spd(,+)*(x^2)(+)-d)^2)(sum);
  model=spd(,+)*(x^2)(+);
  resi=model-d;
  wresi=resi*W;
  res=0.5*(wresi*resi)(sum);
  g=(spd*(2.*x)(-,))(+,)*wresi(+);
  return res;
};
  



func penalty1(x,positive,para,&g){
  /* DOCUMENT
     objective function for the linear regularized problem. needs spd,d,mu,rr=LtL
     para= 1 => no reparameterization, linear
           2 => log reparameterization
           3 => quad reparameterization

     positive=2: gradient clipping    
           
  */
  
  if(positive) x=max(x,0.);

  // linear space
if (para==1){
  res=0.5*((spd(,+)*x(+)-d)^2)(sum)+0.5*mu*(x(+)*(rr(,+)*x(+))(+));
  g=spd(+,)*(spd(,+)*x(+)-d)(+)+mu*(rr(,+)*x(+));
  };
  

  // logarithmic space

if (para==2){
  res=0.5*((spd(,+)*(exp(x))(+)-d)^2)(sum)+0.5*mu*((exp(x))(+)*(rr(,+)*(exp(x))(+))(+));
  g=(spd*(exp(x))(-,))(+,)*(spd(,+)*(exp(x))(+)-d)(+)+0.5*2.*mu*((rr*(exp(x))(-,))(,+)*(exp(x))(+));
};

  

  // quadratic space

if (para==3){ 
  res=0.5*((spd(,+)*(x^2)(+)-d)^2)(sum)+0.5*mu*((x^2)(+)*(rr(,+)*(x^2)(+))(+));
  g=(spd*(2.*x)(-,))(+,)*(spd(,+)*(x^2)(+)-d)(+)+0.5*mu*((rr*(4.*x)(-,))(,+)*(x^2)(+));
  };


  
  if (positive >= 2 && is_array((i = where((x==0.)*(g >0.0))))) g(i) = 0.0;
  
  return res;
};

func linpenalty0(x,&g){
  // plain linear case
  return penalty1(x,0,1,g);
};

func linpenalty2(x,&g){
  // positivity through  gradient clipping
  return penalty1(x,2,1,g);
};

func logpenalty(x,&g){
  //log reparamterisation
  return penalty1(x,0,2,g);
};

func quadpenalty(x,&g){
  // quadratic reparameterisation
  return penalty1(x,0,3,g);
};

func quadpenaltyPlog(x,&g){
  /* DOCUMENT
  quadratic reparameterization, penalization is in log space, e.e. log(..) is t be smooth
  uses mup
  derivative seems ok as long as Ds is small enough and mup=0
  */

  res=quadpenalty(x,gx)+mup*P(log(x^2),rr,gpx);
  //g=gx+(1./x)*2.*x*gpx;
  g=gx+2.*mup*gpx;
  return res;
};

func quadpenaltyPlog2(x,&g){
  /* DOCUMENT
  quadratic reparameterization, penalization is like in 2d age -kin stuff
  derivative seems ok as long as Ds is small enough and mup=0
  */

  //  info,x;
  //res=quadpenaltypenalty(x,gx);
  res=npenalty1(x,gx);
  
  //+mup*P(log(x^2),rr,gpx);
  rx=reform(x,nz,nl);
  res+=mua*Pa(rx^2,rra,ag)+mub*Pv(rx^2,rrb,bg);
  //g=gx+(1./x)*2.*x*gpx;
  //  g=gx+2.*mup*gpx;
  g=gx+2.*mua*((rx*ag)(*))+2.*mub*((rx*bg)(*));
  return res;
};



func para1(az,&g){
  /* DOCUMENT
     straight parametric reconstruction in age, z
     find solutions brutally wiht map
     tam(1,,)=ta
     tam(2,,)=(_m(-:1:68,))
     tam=array(0.,2,68,7)
  */
  
  g=0.*az;
  model=spinterpatest(az(1),az(2),ga,gz);
  resi=model-d;
  wresi=W*resi;
  ki2=(wresi*resi)(sum);
  g(1)=2.*(ga(+)*wresi(+));
  g(2)=2.*(gz(+)*wresi(+));

  bou=muz1*_bound2(az(2),bgz);
  abou=muab1*_abound2(az(1),bga);
  g(1)+=muab1*bga;
  g(2)+=muz1*bgz;
    
  // CLIP Z
  //if ((g(2)>0)&(az(2)<1.e-6)) g(2)=0.;
  
  return ki2+bou;
};

func npecpara(x,&g){

  /* DOCUMENT
     requires _x0,W,d
  */
     
  extern model;
  
  g=0.*x;
  de=x(3:);
  mo1=spinterpatest(x(1),x(2),ga,gz);
  model=mo1;
  t=npe(de,numberof(_x0),gt);
  model=t*model;
  resi=model-d;
  wresi=W*resi;
  ki2=(wresi*resi)(sum);
  g(1)=2.*((t*ga)(+)*wresi(+));
  g(2)=2.*((t*gz)(+)*wresi(+));

  gde=2.*(gt(+,)*(mo1*wresi)(+));
  g(3:)=gde;
  PL,x(2),x(1),msize=0.1;
  pause,1;
  return ki2;
};

func npecparak(x,&g){
  
  /* DOCUMENT
     parametric in age,z, non parametric in losvd, NPEC
     NOTE: NPEL is applied before convolution, like the rest in 1d stuff and unlike the age-kin 2d inversion
     requires _x0,W,and bdata, and ta=b.ages
     derivative looks ok
     
  */
  
  extern model;

  x=x;
  g=0.*x;
  xv=x(3+nde:);
  de=x(3:2+nde);
  t=npe(de,numberof(_x0),gt);
  nemodel=spinterpatest(x(1),x(2),gage,gz);
  rmodel=t*nemodel;
  prmodel=pad(rmodel,pad1,pad2);
  mtf=fft(prmodel);
  nnx=xtreat1(xv^2,nr+pad2+pad1);
  
  model = (1.0/numberof(mtf))*double(fft(mtf*fft(nnx,[-1]),[+1]));
  residual = model - bdata;
  wresi = roll(pW)*residual;
  ki2=(wresi*residual)(sum);
  
  ga=(2.0/numberof(mtf))*double(fft(conj(fft(nnx,[-1,0]))*(fft(wresi,[-1,0])),[+1,0]));
  ga=ga(pad2+2:nr+pad2+1)(::-1);  // HOW STRANGE ??!!

  gv = (2.0/numberof(mtf))*double(fft(conj(mtf)*fft(wresi,[-1,0]),[+1,0]));
  gv=gv(pad2+ni:pad2+nj,);

  //wresi=wresi(pad2+2:nr+pad2+1)(::-1);
  wresi=wresi(pad2+2:nr+pad2+1);

  
  gx1=((t*gage)(+)*(ga)(+));
  gx2=((t*gz)(+)*(ga)(+));
  gde=((gt)(+,)*(nemodel*ga)(+));
  
  gv*=2.*xv;

  g=[];
  grow,g,gx1,gx2,gde,gv;
  
  PL,x(2),x(1),msize=0.1;
  pause,1;
  return ki2;
};


func npecparak2(x,&g){
  
  /* DOCUMENT
     parametric in age,z, NPEC.
     tricky: losvd is specified from outside, through los. Hence gradient is not computed.
     NOTE: NPEL is applied before convolution, like the rest in 1d stuff and unlike the age-kin 2d inversion
     requires _x0,W,and bdata, and ta=b.ages
     derivative looks ok
     to be used in fnpecpara2, copy of fnpecpara
     
  */
  
  extern model;

  x=x;
  g=0.*x;
  xv=los;
  de=x(3:2+nde);
  t=npe(de,numberof(_x0),gt);
  nemodel=spinterpatest(x(1),x(2),gage,gz);
  rmodel=t*nemodel;
  prmodel=pad(rmodel,pad1,pad2);
  mtf=fft(prmodel);
  nnx=xtreat1(xv^2,nr+pad2+pad1);
  
  model = (1.0/numberof(mtf))*double(fft(mtf*fft(nnx,[-1]),[+1]));
  residual = model - bdata;
  wresi = roll(pW)*residual;
  ki2=(wresi*residual)(sum);
  
  ga=(2.0/numberof(mtf))*double(fft(conj(fft(nnx,[-1,0]))*(fft(wresi,[-1,0])),[+1,0]));
  ga=ga(pad2+2:nr+pad2+1)(::-1);  // HOW STRANGE ??!!

  //  gv = (2.0/numberof(mtf))*double(fft(conj(mtf)*fft(wresi,[-1,0]),[+1,0]));
  //gv=gv(pad2+ni:pad2+nj,);

  //wresi=wresi(pad2+2:nr+pad2+1)(::-1);
  wresi=wresi(pad2+2:nr+pad2+1);

  
  gx1=((t*gage)(+)*(ga)(+));
  gx2=((t*gz)(+)*(ga)(+));
  gde=((gt)(+,)*(nemodel*ga)(+));
  
  //gv*=2.*xv;

  g=[];
  //grow,g,gx1,gx2,gde,gv;
  grow,g,gx1,gx2,gde;

  
  PL,x(2),x(1),msize=0.1;
  pause,1;
  return ki2;
};
  
  


func fpara1(d,R=,u0=){
  /* DOCUMENT
     performs a brute force exploration of parameter space, takes the absolute minimum and refines with conjugate gradient
  */

  dic=1.e2;
  tam=array(0.,2,40,5);
  tam(1,,)=ta;
  tam(2,,)=(_m(-:1:40,));
  ftam=tam(1,,)*0.;
  for(i=1;i<=numberof(ta);i++){for(j=1;j<=numberof(_m);j++){ftam(i,j)=para1(tam(,i,j));};};
  u0=(tam(,*))(,(ftam(*)(mnx)));
  u=optim_driver(para1,u0,verb=100,frtol=1.e-20,fatol=1.e-30,fmin=1.e-10,maxeval=500,ndirs=2);
  r=para1(u);
  
  if(!is_void(R)){
    ws;
    plk,ftam,_m,ta;
    //plcir,-log10(para1(u))/dic,u(2),u(1),color="green";
    myPL,u(2),u(1),msize=-log10(para1(u));
    for(i=1;i<=R;i++){
    u0=tam(,((int(abs(random(1))*numberof(tam(1,*))))(1)));
    u1=optim_driver(para1,u0,verb=100,frtol=1.e-20,fatol=1.e-30,fmin=1.e-10,maxeval=500,ndirs=2);
    grow,r,para1(u1);
    myPL,u1(2),u1(1),msize=-log10(para1(u1));
    //plcir,-log10(para1(u1))/dic,u1(2),u1(1),color="green";
    pause,1;
    grow,u,u1;
    };
  };

  if(!is_void(R)) u=reform(u,2,1+R);
  return u(,r(mnx));
};
  

func fnpecpara(d,RG=,u0=,s=,meval=){
  /* DOCUMENT
     performs a brute force exploration of parameter space, takes the absolute minimum and refines with conjugate gradient and npec, R times random guess
     requires _m,ta,nde,
     uses npecpara, which requires _x0,W,d
     uses interpspatest,which uses spinterp, which requires _m, ta, bloc
  */

  extern tam,ftam,r;

  if(is_void(s)) s=1.e-20;
  if(is_void(RG)) RG=-1;
  r=[];
  dic=1.e2;
  tam=array(0.,2,nb,numberof(_m));
  tam(1,,)=ta(,-:1:numberof(_m));
  //tam(1,,)=ta;
  tam(2,,)=(_m(-:1:nb,));
  ftam=tam(1,,)*0.;
  for(i=1;i<=numberof(ta);i++){
    for(j=1;j<=numberof(_m);j++){
      ftam(i,j)=para1(tam(,i,j));
    };
  };
  u0=(tam(,*))(,(ftam(*)(mnx)));
  grow,u0,array(1.,nde);
  
  ws;plk,ftam,_m,ta;
  myPL,u0(2),u0(1),msize=-log10(npecpara(u0));
  pause,1;
  
  if(RG==0) {write,"ok";return u0;};
  
  u=optim_driver(npecpara,u0,verb=100,frtol=1.e-20,fatol=1.e-20,fmin=[],maxeval=meval,ndirs=2);
  r=npecpara(u);
  myPL,u(2),u(1),msize=-log10(npecpara(u));
  
  if(r<=s) {write,"ok";return u;};
  
  if(RG>0){
    for(i=1;i<=RG;i++){
      u0=tam(,((int(abs(random(1))*numberof(tam(1,*))))(1)));
      grow,u0,array(1.,nde);
      u1=optim_driver(npecpara,u0,verb=100,frtol=1.e-20,fatol=1.e-20,fmin=[],maxeval=meval,ndirs=[]);
      //write,u1;
      r1=npecpara(u1);
      //write,r1;
      grow,r,r1;
      myPL,u1(2),u1(1),msize=-log10(npecpara(u1)),color="blue";
      limits,min(ta),max(ta);
      range,min(_m),max(_m);
      //plcir,-log10(para1(u1))/dic,u1(2),u1(1),color="green";
      pause,1;
      grow,u,u1;
      if (r1<=s) break;
    };
  };
  
  if(i>=RG) i=RG;
  if(RG>0) {u=reform(u,2+nde,i+1);write,"ok";return u(,r(mnx));};
  write,"ok";
  return u;
};

func fnpecpara2(d,RG=,u0=,s=,meval=){
  /* DOCUMENT
     performs a brute force exploration of parameter space, takes the absolute minimum and refines with conjugate gradient and npec, R times random guess
     requires _m,ta,nde,
     uses npecpara, which requires _x0,W,d
     uses interpspatest,which uses spinterp, which requires _m, ta, bloc
  */

  extern tam,ftam,r;

  if(is_void(s)) s=1.e-20;
  if(is_void(RG)) RG=-1;
  r=[];
  dic=1.e2;
  tam=array(0.,2,nb,numberof(_m));
  tam(1,,)=ta(,-:1:numberof(_m));
  //tam(1,,)=ta;
  tam(2,,)=(_m(-:1:nb,));
  ftam=tam(1,,)*0.;
  for(i=1;i<=numberof(ta);i++){
    for(j=1;j<=numberof(_m);j++){
      qi=tam(,i,j);
      grow,qi,array(1.,nde);
      ftam(i,j)=npecparak2(qi);
    };
  };
  u0=(tam(,*))(,(ftam(*)(mnx)));
  grow,u0,array(1.,nde);
  
  ws;plk,ftam,_m,ta;
  myPL,u0(2),u0(1),msize=-log10(npecparak2(u0));
  pause,1;
  
  if(RG==0) {write,"ok";return u0;};
  
  u=optim_driver(npecparak2,u0,verb=100,frtol=1.e-20,fatol=1.e-20,fmin=[],maxeval=meval,ndirs=2);
  r=npecparak2(u);
  myPL,u(2),u(1),msize=-log10(npecparak2(u));
  
  if(r<=s) {write,"ok";return u;};
  
  if(RG>0){
    for(i=1;i<=RG;i++){
      u0=tam(,((int(abs(random(1))*numberof(tam(1,*))))(1)));
      grow,u0,array(1.,nde);
      u1=optim_driver(npecparak2,u0,verb=100,frtol=1.e-20,fatol=1.e-20,fmin=[],maxeval=meval,ndirs=[]);
      //write,u1;
      r1=npecparak2(u1);
      //write,r1;
      grow,r,r1;
      myPL,u1(2),u1(1),msize=-log10(npecparak2(u1)),color="blue";
      limits,min(ta),max(ta);
      range,min(_m),max(_m);
      //plcir,-log10(para1(u1))/dic,u1(2),u1(1),color="green";
      pause,1;
      grow,u,u1;
      if (r1<=s) break;
    };
  };
  
  if(i>=RG) i=RG;
  if(RG>0) {u=reform(u,2+nde,i+1);write,"ok";return u(,r(mnx));};
  write,"ok";
  return u;
};


func fnpecparak(d,RG=,u0=,s=,meval=,w0=){
  /* DOCUMENT
     same as fnpecpara with kinematics
     requires _m,ta,nde,
     uses npecparak, which requires _x0,W,bdata and ta.
     uses interpspatest,which uses spinterp, which requires _m, ta, bloc
     w0 is the first PSF FWHM tried
     doesnt converge well.
  */

  extern tam,ftam,r;
  
  if(is_void(s)) s=1.e-20;
  if(is_void(RG)) RG=-1;
  r=[];
  dic=1.e2;
  tam=array(0.,2,nb,numberof(_m));
  tam(1,,)=ta(,-:1:numberof(_m));
  //tam(1,,)=ta;
  tam(2,,)=(_m(-:1:nb,));
  ftam=tam(1,,)*0.;
  for(i=1;i<=numberof(ta);i++){
    for(j=1;j<=numberof(_m);j++){
      q=tam(,i,j);
      grow,q,array(1.,nde),makebump(nlos,int(nlos/2),w0,N=1);
      ftam(i,j)=npecparak(q);
    };
  };
  u0=(tam(,*))(,(ftam(*)(mnx)));
  grow,u0,array(1.,nde),makebump(nlos,int(nlos/2),w0,N=1);
  
  ws;plk,ftam,_m,ta;
  myPL,u0(2),u0(1),msize=-log10(npecparak(u0));
  pause,1;
  
  if(RG==0) {write,"ok";return u0;};
  
  u=optim_driver(npecparak,u0,verb=100,frtol=1.e-20,fatol=1.e-20,fmin=[],maxeval=meval,ndirs=2);
  r=npecparak(u);
  myPL,u(2),u(1),msize=-log10(npecparak(u));
  
  if(r<=s) {write,"ok";return u;};
  
  if(RG>0){
    for(i=1;i<=RG;i++){
      u0=tam(,((int(abs(random(1))*numberof(tam(1,*))))(1)));
      grow,u0,array(1.,nde),makebump(nlos,int(nlos/2),w0,N=1);
      u1=optim_driver(npecparak,u0,verb=100,frtol=1.e-20,fatol=1.e-20,fmin=[],maxeval=meval,ndirs=[]);
      //write,u1;
      r1=npecparak(u1);
      //write,r1;
      grow,r,r1;
      myPL,u1(2),u1(1),msize=-log10(npecparak(u1)),color="blue";
      limits,min(ta),max(ta);
      range,min(_m),max(_m);
      //plcir,-log10(para1(u1))/dic,u1(2),u1(1),color="green";
      pause,1;
      grow,u,u1;
      if (r1<=s) break;
    };
  };
  
    if(i>=RG) i=RG;
    info,u;
    if(RG>0) {u=reform(u,2+nde+nlos,int(numberof(u)/(2+nde+nlos)));write,"ok";return u(,r(mnx));};
  write,"ok";
  return u;
};




func fpop(X,&g){
  /*DOCUMENT
    mono-metallic pop, dust screen from Calzetti 1999*/
  extern ki2,rx;
  g=0.*X;
  x=X(:nab);
  ebv=X(0);
  t=ds(ebv,x0,tg);
  fd=d*t;
  ki2=0.5*((spb(,+)*x(+)-fd)^2)(sum);
  rx=0.5*mu*(x(+)*(rr(,+)*x(+))(+));
  res=ki2+rx;
  g(:nab)=spb(+,)*(spb(,+)*x(+)-fd)(+)+mu*(rr(,+)*x(+));
  g(0)=(tg*d*(fd-spb(,+)*x(+)))(sum);
  return res;
};





func faz(X,&g){
  /*DOCUMENT
    1dAMR and dust screen from Calzetti 1999*/

  extern ki2,rx,rz;
  g=X*0.;
  x=X(:nab)
  z=X(nab+1:2*nab);
  ebv=X(0);
  t=ds(ebv,x0,gt);
  

  rspb=buildb2(z,bg);
  ki2=0.5*((rspb(,+)*x(+)-fd)^2)(sum);
  rx=0.5*mu*(x(+)*(rr(,+)*x(+))(+));
  rz=0.5*muZ*(z(+)*(rr(,+)*z(+))(+));
  res=ki2+rx+rz+mub*bound2(z,gbz);
  gkiz=1.*(bg*((rspb(,+)*x(+)-fd)(,-:1:numberof(x))*x(-:1:numberof(d),))); //right
  g(nab+1:2*nab)=gkiz(sum,)+muZ*(rr(,+)*z(+))+mub*gbz;  // z part 
  g(:nab)=rspb(+,)*(rspb(,+)*x(+)-fd)(+)+mu*(rr(,+)*x(+)); // x part
  g(0)=(gt*d*(fd-rspb(,+)*x(+)))(sum);  // extinction part of gradient
  //g=mus*g;
  return res;
};

func qfaz(X,&g){

  extern ki2,rx,rz;
  g=X*0.;
  x=X(:nab)
  z=X(nab+1:2*nab);
  ebv=X(0);
  t=ds(ebv,x0,gt);
  
  // FORSAKEN

  
  rspb=buildb2(z,bg);
  model=t*((rspb)(,+)*(x^2)(+));
  ki2=0.5*((model-d)^2)(sum);
  rx=0.5*mu*(x(+)*(rr(,+)*x(+))(+));
  rz=0.5*muZ*(z(+)*(rr(,+)*z(+))(+));
  res=ki2+rx+rz+mub*bound2(z,gbz);
  gkiz=1.*(bg*((model-d)(,-:1:numberof(x))*x(-:1:numberof(d),))); //right
  g(nab+1:2*nab)=gkiz(sum,)+muZ*(rr(,+)*z(+))+mub*gbz;  // z part 
  g(:nab)=rspb(+,)*(model-d)(+)+mu*(rr(,+)*x(+)); // x part
  g(0)=(gt*d*(fd-rspb(,+)*x(+)))(sum);  // extinction part of gradient
  //g=mus*g;
  return res;
};



func qnfaz(X,&g){
  
 extern ki2,rx,rz,rebv,model;
  g=X*0.;
  x=X(:nab);
  z=X(nab+1:2*nab);
  ebv=X(2*nab+1:3*nab);
  t=vds(ebv,_x0,gt);
  

  
  rspb=buildb2(z,bg);
  model=(t*rspb)(,+)*(x^2)(+);
  ki2=0.5*(((t*rspb)(,+)*(x^2)(+)-d)^2)(sum);
  rx=0.5*mu*((x^2)(+)*(rr(,+)*(x^2)(+))(+));
  rz=0.5*muZ*(z(+)*(rr(,+)*z(+))(+));
  rebv=0.5*muebv*(ebv(+)*(rr(,+)*ebv(+))(+));
  res=ki2+rx+rz+rebv+mub*bound2(z,gbz);
  gkiz=1.*(bg*t*(((t*rspb)(,+)*(x^2)(+)-d)(,-:1:numberof(x))*(x^2)(-:1:numberof(d),))); //right
  g(nab+1:2*nab)=gkiz(sum,)+muZ*(rr(,+)*z(+))+mub*gbz;  // z part 
  g(:nab)=((t*rspb)*(2.*x)(-,))(+,)*((t*rspb)(,+)*(x^2)(+)-d)(+)+0.5*mu*((rr*(4.*x)(-,))(,+)*(x^2)(+)); // x part
  gebv=1.*(gt*rspb*(((t*rspb)(,+)*(x^2)(+)-d)(,-:1:numberof(x))*(x^2)(-:1:numberof(d),)));
  g(2*nab+1:3*nab)=gebv(sum,)+muebv*(rr(,+)*ebv(+));
  return res;
};



func Wqnfaz(X,&g){
  
 extern ki2,rx,rz,rebv,model;
  g=X*0.;
  x=X(:nab);
  z=X(nab+1:2*nab);
  ebv=X(2*nab+1:3*nab);
  t=vds(ebv,_x0,gt);
  

  
  rspb=buildb2(z,bg);
  model=(t*rspb)(,+)*(x^2)(+);
  ki2=0.5*((((t*rspb)(,+)*(x^2)(+)-d)^2)*W)(sum);
  rx=0.5*mu*((x^2)(+)*(rr(,+)*(x^2)(+))(+));
  rz=0.5*muZ*(z(+)*(rr(,+)*z(+))(+));
  rebv=0.5*muebv*(ebv(+)*(rr(,+)*ebv(+))(+));
  res=ki2+rx+rz+rebv+mub*bound2(z,gbz);
  gkiz=1.*(bg*t*((((t*rspb)(,+)*(x^2)(+)-d)*W)(,-:1:numberof(x))*(x^2)(-:1:numberof(d),))); //right
  g(nab+1:2*nab)=gkiz(sum,)+muZ*(rr(,+)*z(+))+mub*gbz;  // z part 
  g(:nab)=((t*rspb)*(2.*x)(-,))(+,)*(((t*rspb)(,+)*(x^2)(+)-d)*W)(+)+0.5*mu*((rr*(4.*x)(-,))(,+)*(x^2)(+)); // x part
  gebv=1.*(gt*rspb*((((t*rspb)(,+)*(x^2)(+)-d)*W)(,-:1:numberof(x))*(x^2)(-:1:numberof(d),)));
  g(2*nab+1:3*nab)=gebv(sum,)+muebv*(rr(,+)*ebv(+));
  return res;
};



func Wqcnfaz(X,pos,&g){

  // quadratic in age, positivity of E(B-v) through clipping
  // needs bloc,_x0,W,mu,muZ,mue,mub,nab,

  

  extern ki2,rx,rz,rebv,model;
  
  g=X*0.;
  x=X(:nab);
  z=X(nab+1:2*nab);
  ebv=X(2*nab+1:3*nab);
  if (pos) ebv=max(ebv,0.);   //E(B-v) clipping
  t=vds(ebv,_x0,gt);
  

  
  rspb=buildb2(z,bg);
  model=(t*rspb)(,+)*(x^2)(+);
  ki2=0.5*((((t*rspb)(,+)*(x^2)(+)-d)^2)*W)(sum);
  rx=0.5*mu*((x^2)(+)*(rr(,+)*(x^2)(+))(+));
  rz=0.5*muz*(z(+)*(rr(,+)*z(+))(+));
  rebv=0.5*mue*(ebv(+)*(rr(,+)*ebv(+))(+));
  res=ki2+rx+rz+rebv+mub*bound2(z,gbz);
  gkiz=1.*(bg*t*((((t*rspb)(,+)*(x^2)(+)-d)*W)(,-:1:numberof(x))*(x^2)(-:1:numberof(d),))); //right
  g(nab+1:2*nab)=gkiz(sum,)+muz*(rr(,+)*z(+))+mub*gbz;  // z part 
  g(:nab)=((t*rspb)*(2.*x)(-,))(+,)*(((t*rspb)(,+)*(x^2)(+)-d)*W)(+)+0.5*mu*((rr*(4.*x)(-,))(,+)*(x^2)(+)); // x part
  gebv=1.*(gt*rspb*((((t*rspb)(,+)*(x^2)(+)-d)*W)(,-:1:numberof(x))*(x^2)(-:1:numberof(d),)));
  if (!is_void(pos)&&(is_array((i = where((!ebv)*(gebv(sum,) >= 0.0)))))) gebv(i) = 0.0; // clipping E(B-V)
  g(2*nab+1:3*nab)=gebv(sum,)+mue*(rr(,+)*ebv(+));
  
  return res;
};

func Wqcfaz(X,pos,&g){

  // quadratic in age, positivity of E(B-v) through clipping. Only 1 E(B-V)
  // needs bloc,_x0,W,mu,muZ,mue,mub,nab,

  // CURRENTLY IN WORK

  extern ki2,rx,rz,rebv,model;
  
  g=X*0.;
  x=X(:nab);
  z=X(nab+1:2*nab);
  ebv=X(0);
  if (pos) ebv=max(ebv,0.);   //E(B-v) clipping
  t=vds(ebv,_x0,gt);
  

  
  rspb=buildb2(z,bg);
  model=(t*rspb)(,+)*(x^2)(+);
  ki2=0.5*((((model-d)^2)*W))(sum);
  rx=0.5*deconv_mux*((x^2)(+)*(rr(,+)*(x^2)(+))(+));
  rz=0.5*deconv_muz*(z(+)*(rr(,+)*z(+))(+));
  //rebv=0.5*deconv_mue*(ebv(+)*(rr(,+)*ebv(+))(+));
  rebv=0.;
  res=ki2+rx+rz+rebv+mub*bound2(z,gbz);
  gkiz=1.*(bg*t*((((t*rspb)(,+)*(x^2)(+)-d)*W)(,-:1:numberof(x))*(x^2)(-:1:numberof(d),))); //right
  g(nab+1:2*nab)=gkiz(sum,)+deconv_muz*(rr(,+)*z(+))+mub*gbz;  // z part 
  g(:nab)=((t*rspb)*(2.*x)(-,))(+,)*(((t*rspb)(,+)*(x^2)(+)-d)*W)(+)+0.5*deconv_mux*((rr*(4.*x)(-,))(,+)*(x^2)(+)); // x part
  gebv=1.*(gt*rspb*((((t*rspb)(,+)*(x^2)(+)-d)*W)(,-:1:numberof(x))*(x^2)(-:1:numberof(d),)));   //**********************************
  gebv=gebv(sum);
  
  if (!is_void(pos)&&(is_array((i = where((!ebv)*(gebv >= 0.0)))))) gebv(i) = 0.0; // clipping E(B-V)
  g(2*nab+1)=gebv
    
  return res;
}; 


func Wqcfaz2(X,pos,&g){

  // quadratic in age, positivity of E(B-v) through clipping. Only 1 E(B-V)
  // needs bloc,_x0,W,mu,muZ,mue,mub,nab,
  // changing treatment of dust screen. Should speed up by several
  // derivative seems ok.
  
  extern ki2,rx,rz,rebv,model;
  
  g=X*0.;
  x=X(:nab);
  z=X(nab+1:2*nab);
  ebv=X(0);
  if (pos) ebv=max(ebv,0.);   //E(B-v) clipping
  t=ds(ebv,_x0,gt);
  

  
  rspb=buildb2(z,bg);
  model=(rspb(,+)*(x^2)(+))*t;
  resi=model-d;
  wresi=resi*W;
  ki2=(wresi*resi)(sum)
  rx=deconv_mux*P(x^2,rr,gPx);
  rz=deconv_muz*P(z,rr3,gPz);
  //rebv=0.5*deconv_mue*(ebv(+)*(rr(,+)*ebv(+))(+));
  rebv=0.;
  res=ki2+rx+rz+rebv+mub*bound2(z,gbz);
  gkiz=2.*((bg*t(,-:1:numberof(x))*(x^2)(-:1:numberof(d),))(+,)*(wresi(,-))(+,)); //right
  g(nab+1:2*nab)=gkiz+deconv_muz*gPz+mub*gbz;  // z part 
  g(:nab)=2.*(((rspb*(t(,-)))*(2.*x(-,)))(+,)*(wresi)(+))+deconv_mux*2.*x*gPx; // x part
  //gebv=1.*((gt(,-)*rspb)*(wresi(,-:1:numberof(x))*(x^2)(-:1:numberof(d),)));
  //**********************************
  gebv=2.*(gt(+)*((rspb(,+)*(x^2)(+))*wresi)(+));
  //gebv=gebv(sum);
  
  if ((pos)&&((!ebv)*gebv >= 0.0)) gebv = 0.0; // clipping E(B-V)
  g(2*nab+1)=gebv
    
  return res;
}; 


func Wqcfaz3(X,pos,&g){

  // extinction and spectroiphotometric errors modelized by a poynomial npolf3
  // quadratic in age
  // needs bloc,_x0,W,mu,muZ,mue,mub,nab,
  // derivative seems ok
  // but it doesnt converge very well

  
  extern ki2,rx,rz,rebv,model;
  
  g=X*0.;
  x=X(:nab);
  z=X(nab+1:2*nab);
  ebv=X(2*nab+1:);
  if (pos) ebv=max(ebv,0.);   //E(B-v) clipping
  t=npolf3(ebv,_x0,gt);
  

  
  rspb=buildb2(z,bg);
  model=(rspb(,+)*(x^2)(+))*t;
  resi=model-d;
  wresi=((((model-d))*W));
  ki2=(wresi*resi)(sum)
  rx=deconv_mux*P(x^2,rr,gPx);
  rz=deconv_muz*P(z,rr3,gPz);
  //rebv=0.5*deconv_mue*(ebv(+)*(rr(,+)*ebv(+))(+));
  rebv=0.;
  res=ki2+rx+rz+rebv+mub*bound2(z,gbz);
  gkiz=2.*((bg*t(,-:1:numberof(x))*(x^2)(-:1:numberof(d),))(+,)*(wresi(,-))(+,)); //right
  g(nab+1:2*nab)=gkiz+deconv_muz*gPz+mub*gbz;  // z part 
  g(:nab)=2.*(((rspb*(t(,-)))*(2.*x(-,)))(+,)*(wresi)(+))+deconv_mux*2.*x*gPx; // x part
  //gebv=1.*((gt(,-)*rspb)*(wresi(,-:1:numberof(x))*(x^2)(-:1:numberof(d),)));
  //**********************************
  gebv=2.*(gt(+,)*((rspb(,+)*(x^2)(+))*wresi)(+));
  //gebv=gebv(sum);
  
  if ((pos)&&(is_array((!ebv)*gebv >= 0.0))) gebv = 0.0; // clipping E(B-V)
  g(2*nab+1:)=gebv
    
  return res;
}; 

func Wqcfaz4(X,pos,&g){

  // quadratic in age, and spectrophotometric error +extinction is nonparametric
  // through npe
  // needs bloc,_x0,W,mu,muZ,mub,nab,mue,nde
  // 
  // derivative seems ok.
  
  extern ki2,rx,rz,rebv,model;
  
  g=X*0.;
  x=X(:nab);
  z=X(nab+1:2*nab);
  de=X(2*nab+1:2*nab+nde);
  if (pos) de=max(de,0.);   //E(B-v) clipping
  t=npe(de,numberof(_x0),gt);

  
  rspb=buildb2(z,bg);
  model=(rspb(,+)*(x^2)(+))*t;
  resi=model-d;
  wresi=resi*W;
  ki2=(wresi*resi)(sum)
  rx=deconv_mux*P(x^2,rr,gPx);
  rz=deconv_muz*P(z,rr3,gPz);
  rde=0.;
  res=ki2+rx+rz+rde+mub*bound2(z,gbz);
  gkiz=2.*((bg*t(,-:1:numberof(x))*(x^2)(-:1:numberof(d),))(+,)*(wresi(,-))(+,)); //right
  g(nab+1:2*nab)=gkiz+deconv_muz*gPz+mub*gbz;  // z part 
  g(:nab)=2.*(((rspb*(t(,-)))*(2.*x(-,)))(+,)*(wresi)(+))+deconv_mux*2.*x*gPx; // x part
  //gebv=1.*((gt(,-)*rspb)*(wresi(,-:1:numberof(x))*(x^2)(-:1:numberof(d),)));
  //**********************************
  gde=2.*(gt(+,)*((rspb(,+)*(x^2)(+))*wresi)(+));
  //gebv=gebv(sum);
  
  if ((pos)&&((!ebv)*gebv >= 0.0)) gebv = 0.0; // clipping E(B-V)
  g(2*nab+1:2*nab+nde)=gde
    
  return res;
}; 

func Wqcfaz5(X,pos,&g){

  // quadratic in age, and spectrophotometric error +extinction is nonparametric
  // through npe
  // needs bloc,_x0,W,mu,muZ,mub,nab,mue,nde,muc,co
  // adds a penalization on x so that x is normalized to c
  // that is x(:nab)(sum)=c

  extern ki2,rx,rz,rebv,model;

  res=Wqcfaz4(X,pos,g)+muc*Pn((X(:nab))^2,co,gc);
  gc=gc*2.*X(:nab);
  g(:nab)+=muc*gc;
  return res;
};
  
  


  func Wfaz1cc(X,positive, &g)
{

  //needs _x0,deconv_muv,deconv_mux,deconv_mue,deconv_muz,mub,
  // _m,nab,pad1,pad2,d
  
  
  extern model,bdata,pW;

  g=X*0.;
  if (positive) X = max(0.0, x);  // clips all variables
  x=X(:nab);
  z=X(nab+1:2*nab);
  ebv=X(2*nab+1:3*nab);
  xv=X(3*nab+1:);
  t=vds(ebv,_x0,gt);
  //bdata=roll(pad(d,pad1,pad2))(::-1);
  //pW=roll(pad(W,pad1,pad2))(::-1);
  

  rglx=deconv_mux*ip_roughness(x,1,dx,which=1,order=2);
  rglv=deconv_muv*ip_roughness(xv,1,dxv,which=1,order=2);
  rgle=deconv_mue*ip_roughness(ebv,1,de,which=1,order=2);
  rglz=deconv_muz*ip_roughness(z,1,dz,which=1,order=2)+mub*bound2(z,gbz);
  nnx=xtreat1(xv,nr+pad2+pad1);

  sp=buildb2(z,bg);
  sp0=sp;
  sp=t*sp;
  ma=sp(,+)*x(+);
  
  pma=pad(ma,pad1,pad2);
  mtf=fft(pma);
  
  model = (1.0/numberof(mtf))*double(fft(mtf*fft(nnx,-1),+1));
  

  residual = model - bdata;
  wr = roll(pW)*residual;

  ga=(2.0/numberof(mtf))*double(fft(conj(fft(nnx,[-1,0]))*(fft(wr,[-1,0])),[+1,0]));
  ga=ga(pad2+2:nr+pad2+1)(::-1);  // HOW STRANGE ??!!

  gx=ga(+)*sp(+,);
  gx+=  deconv_mux*dx;

  gz=ga(+)*(bg*t)(+,);
  gz*=x;
  gz+=deconv_muz*dz+mub*gbz;

  ge=ga(+)*(sp0*gt)(+,);
  ge*=x;
  ge+=deconv_mue*de;
  
  
  gv = (2.0/numberof(mtf))*double(fft(conj(mtf)*fft(wr,[-1,0]),[+1,0]));
  gv=gv(pad2+ni:pad2+nj,);
  gv +=deconv_muv*dxv;

  grow,gx,gz,ge,gv;
  g=gx;
  if (positive >= 2 && is_array((i = where((!X)*(g > 0.0))))) g(i) = 0.0;
  return sum(wr*residual) +rglv+rglx+rglz;
};

  func Wfaz1cq(X,positive, &g)
{

  // quadratic in sfh only
  //needs _x0,deconv_muv,deconv_mux,deconv_mue,deconv_muz,mub,
  // _m,nab,pad1,pad2,bloc,
  
  extern wr,model,rglx,rglv,rgle,rglz;

  g=X*0.;
  if (positive) X = max(0.0, x);
  x=X(:nab);
  z=X(nab+1:2*nab);
  ebv=X(2*nab+1:3*nab);
  xv=X(3*nab+1:);
  t=vds(ebv,_x0,gt);
  
  

  rglx=deconv_mux*ip_roughness(x^2,1,dx,which=1,order=2);
  dx*=2.*x;
  rglv=deconv_muv*ip_roughness(xv,1,dxv,which=1,order=2);
  rgle=deconv_mue*ip_roughness(ebv,1,de,which=1,order=2);
  rglz=deconv_muz*ip_roughness(z,1,dz,which=1,order=2)+mub*bound2(z,gbz);
  nnx=xtreat1(xv,nr+pad2+pad1);

  sp=buildb2(z,bg);
  sp0=sp;
  sp=t*sp;
  ma=sp(,+)*(x^2)(+);
  
  pma=pad(ma,pad1,pad2);
  mtf=fft(pma);
  
  model = (1.0/numberof(mtf))*double(fft(mtf*fft(nnx,-1),+1));
  

  residual = model - bdata;
  wr = roll(pW)*residual;

  ga=(2.0/numberof(mtf))*double(fft(conj(fft(nnx,[-1,0]))*(fft(wr,[-1,0])),[+1,0]));
  ga=ga(pad2+2:nr+pad2+1)(::-1);  // HOW STRANGE ??!!

  gx=2.*x*(ga(+)*sp(+,));
  gx+=  deconv_mux*dx;

  gz=ga(+)*(bg*t)(+,);
  gz*=(x^2);
  gz+=deconv_muz*dz+mub*gbz;

  ge=ga(+)*(sp0*gt)(+,);
  ge*=(x^2);
  ge+=deconv_mue*de;
  
  
  gv = (2.0/numberof(mtf))*double(fft(conj(mtf)*fft(wr,[-1,0]),[+1,0]));
  gv=gv(pad2+ni:pad2+nj,);
  gv +=deconv_muv*dxv;

  grow,gx,gz,ge,gv;
  g=gx;
  if (positive >= 2 && is_array((i = where((!X)*(g > 0.0))))) g(i) = 0.0;
  return sum(wr*residual) +rglv+rglx+rglz+rgle;
 

  
};

  func Wfaz1cqq(X,positive, &g)
{

  // quadratic in sfh and losvd,clipping in E(B-V)
  //needs _x0,deconv_muv,deconv_mux,deconv_mue,deconv_muz,mub,
  // _m,nab,pad1,pad2,
  
  extern wr,model,rglx,rglv,rgle,rglz;

  g=X*0.;
  x=X(:nab);
  z=X(nab+1:2*nab);
  ebv=X(2*nab+1:3*nab);
  if (!is_void(pos)) ebv=max(ebv,0.);   //E(B-v) clipping
  xv=X(3*nab+1:);
  t=vds(ebv,_x0,gt);
  
  

  rglx=deconv_mux*P(x^2,rr,dx);
  dx*=deconv_mux*2.*x;
  //rglv=deconv_muv*ip_roughness(xv^2,1,dxv,which=1,order=2);
  rglv=deconv_muv*P(xv^2,rr1,dxv);
  dxv*=deconv_muv*2.*xv;
  rgle=deconv_mue*P(ebv,rr,de);
  de*=deconv_mue;
  rglz=deconv_muz*P(z,rr3,dz)+mub*bound2(z,gbz);
  dz*=deconv_muz;
  dz+=mub*gbz;
  nnx=xtreat1(xv^2,nr+pad2+pad1);
  nxv=xtreat1(xv,nr+pad2+pad1);

  sp=buildb2(z,bg);
  sp0=sp;
  sp=t*sp;
  ma=sp(,+)*(x^2)(+);
  
  pma=pad(ma,pad1,pad2);
  mtf=fft(pma);
  
  model = (1.0/numberof(mtf))*double(fft(mtf*fft(nnx,-1),+1));
  

  residual = model - bdata;
  wr = roll(pW)*residual;

  ga=(2.0/numberof(mtf))*double(fft(conj(fft(nnx,[-1,0]))*(fft(wr,[-1,0])),[+1,0]));
  ga=ga(pad2+2:nr+pad2+1)(::-1);  // HOW STRANGE ??!!

  gx=2.*x*(ga(+)*sp(+,));
  gx+=dx;

  gz=ga(+)*(bg*t)(+,);
  gz*=(x^2);
  gz+=dz;

  ge=ga(+)*(sp0*gt)(+,);
  ge*=(x^2);
  ge+=de;
  if (!is_void(pos)&&(is_array((i = where((!ebv)*(ge >= 0.0)))))) gebv(i) = 0.0; // clipping E(B-V)
  
  
  gv = (2.0/numberof(mtf))*double(fft(conj(mtf)*fft(wr,[-1,0]),[+1,0]));
  gv=gv(pad2+ni:pad2+nj,);
  gv*=2.*xv;
  gv+=dxv;

  grow,gx,gz,ge,gv;
  g=gx;
  
  //if (positive >= 2 && is_array((i = where((!X)*(g > 0.0))))) g(i) = 0.0;
  return sum(wr*residual) +rglv+rglx+rglz+rgle;
 
};

  func B2R(X,positive, &g)
{

  //needs _x0,deconv_muv,deconv_mux,deconv_mue,deconv_muz,mub,
  // _m,nab,pad1,pad2,d
  // NO POSITIVITY NO REPARAMETERIZAZTION NO BOUNDARIES FOR Z NO CLIPPING
  // supposed to work with the new optimizatuin  routines by Eric
  // where the bounds are given to the optimization driver.
  
  
  extern model,bdata,pW,rglv,rglx,rglz,rgle,resid;

  g=X*0.; 
  x=X(:nab);
  z=X(nab+1:2*nab);
  ebv=X(2*nab+1:3*nab);
  xv=X(3*nab+1:);
  t=vds(ebv,_x0,gt);
  //bdata=roll(pad(d,pad1,pad2))(::-1);
  //pW=roll(pad(W,pad1,pad2))(::-1);
  

  rglx=deconv_mux*ip_roughness(x,1,dx,which=1,order=2);
  rglv=deconv_muv*ip_roughness(xv,1,dxv,which=1,order=2);
  rgle=deconv_mue*ip_roughness(ebv,1,de,which=1,order=2);
  rglz=deconv_muz*ip_roughness(z,1,dz,which=1,order=2);
  nnx=xtreat1(xv,nr+pad2+pad1);

  sp=buildb2(z,bg);
  sp0=sp;
  sp=t*sp;
  ma=sp(,+)*x(+);
  
  pma=pad(ma,pad1,pad2);
  mtf=fft(pma);
  
  model = (1.0/numberof(mtf))*double(fft(mtf*fft(nnx,-1),+1));
  

  residual = model - bdata;
  wr = roll(pW)*residual;

  ga=(2.0/numberof(mtf))*double(fft(conj(fft(nnx,[-1,0]))*(fft(wr,[-1,0])),[+1,0]));
  ga=ga(pad2+2:nr+pad2+1)(::-1);  // HOW STRANGE ??!!

  gx=ga(+)*sp(+,);
  gx+=  deconv_mux*dx;

  gz=ga(+)*(bg*t)(+,);
  gz*=x;
  gz+=deconv_muz*dz;

  ge=ga(+)*(sp0*gt)(+,);
  ge*=x;
  ge+=deconv_mue*de;
  
  
  gv = (2.0/numberof(mtf))*double(fft(conj(mtf)*fft(wr,[-1,0]),[+1,0]));
  gv=gv(pad2+ni:pad2+nj,);
  gv+=deconv_muv*dxv;

  grow,gx,gz,ge,gv;
  g=gx;

  resid=sum(wr*residual);
  
  return resid +rglv+rglx+rglz+rgle;
};


  func W1faz1cqq(X,positive, &g)
{

  // quadratic in sfh and losvd,clipping in E(B-V)
  //needs _x0,deconv_muv,deconv_mux,deconv_mue,deconv_muz,mub,
  // _m,nab,pad1,pad2,
  // Here there is only one dust screen for the whole population like
  // everyone does with SDSS.
  
  extern wr,model,rglx,rglv,rgle,rglz;

  g=X*0.;
  x=X(:nab);
  z=X(nab+1:2*nab);
  ebv=X(2*nab+1);
  if (!is_void(positive)) ebv=max(ebv,0.);   //E(B-v) clipping
  xv=X(2*nab+2:);
  t=vds(ebv(-:1:nab),_x0,gt);
  
  
  rglx=deconv_mux*P(x^2,rr,dx);
  dx*=deconv_mux*2.*x;
  //rglv=deconv_muv*ip_roughness(xv^2,1,dxv,which=1,order=2);
  rglv=deconv_muv*P(xv^2,rr1,dxv);
  dxv*=deconv_muv*2.*xv;
  rgle=0.;
  de=0.;
  rglz=deconv_muz*P(z,rr3,dz)+mub*bound2(z,gbz);
  dz*=deconv_muz;
  dz+=mub*gbz;
  nnx=xtreat1(xv^2,nr+pad2+pad1);
  nxv=xtreat1(xv,nr+pad2+pad1);

  sp=buildb2(z,bg);
  sp0=sp;
  sp=t*sp;
  ma=sp(,+)*(x^2)(+);
  
  pma=pad(ma,pad1,pad2);
  mtf=fft(pma);
  
  model = (1.0/numberof(mtf))*double(fft(mtf*fft(nnx,-1),+1));
  

  residual = model - bdata;
  wr = roll(pW)*residual;

  ga=(2.0/numberof(mtf))*double(fft(conj(fft(nnx,[-1,0]))*(fft(wr,[-1,0])),[+1,0]));
  ga=ga(pad2+2:nr+pad2+1)(::-1);  // HOW STRANGE ??!!

  gx=2.*x*(ga(+)*sp(+,));
  gx+=  deconv_mux*dx;

  gz=ga(+)*(bg*t)(+,);
  gz*=(x^2);
  gz+=dz;

  ge=ga(+)*(sp0*gt)(+,);
  ge*=(x^2);
  ge=ge(sum);
  if (!is_void(positive)&&(is_array((i = where((!ebv)*(ge >= 0.0)))))) ge(i) = 0.0; // clipping E(B-V)
  
  
  gv = (2.0/numberof(mtf))*double(fft(conj(mtf)*fft(wr,[-1,0]),[+1,0]));
  gv=gv(pad2+ni:pad2+nj,);
  gv*=2.*xv;
  gv+=deconv_muv*dxv;

  grow,gx,gz,ge,gv;
  g=gx;
  
  //if (positive >= 2 && is_array((i = where((!X)*(g > 0.0))))) g(i) = 0.0;
  return sum(wr*residual) +rglv+rglx+rglz+rgle;
 
};

 func W1faz1cqq4(X,positive, &g)
{

  // Non Parametric Extinction Law. Like wqcfaz4 but with losvd search
  // quadratic in sfh and losvd,clipping in E(B-V)
  //needs _x0,deconv_muv,deconv_mux,deconv_mue,deconv_muz,mub,
  // _m,nab,pad1,pad2,
  // Here there is only one dust screen for the whole population like
  // everyone does with SDSS.
  // derivative seems ok
  
  extern wr,model,rglx,rglv,rgle,rglz;

  g=X*0.;
  x=X(:nab);
  z=X(nab+1:2*nab);
  de=X(2*nab+1:2*nab+nde);
  if (positive) de=max(de,0.);   //E(B-v) clipping
  xv=X(2*nab+nde+1:);
  t=npe(de,numberof(_x0),gt);
  
  
  rglx=deconv_mux*P(x^2,rr,dx);
  dx*=deconv_mux*2.*x;
  //rglv=deconv_muv*ip_roughness(xv^2,1,dxv,which=1,order=2);
  rglv=deconv_muv*P(xv^2,rr1,dxv);
  dxv*=deconv_muv*2.*xv;
  rgle=0.;
  rglz=deconv_muz*P(z,rr3,dz)+mub*bound2(z,gbz);
  dz*=deconv_muz;
  dz+=mub*gbz;
  nnx=xtreat1(xv^2,nr+pad2+pad1);
  nxv=xtreat1(xv,nr+pad2+pad1);  // WHAT's THAT ??

  sp=buildb2(z,bg);
  sp0=sp;
  sp=t*sp;
  sp01=sp0(,+)*(x^2)(+);
  ma=sp(,+)*(x^2)(+);
  
  pma=pad(ma,pad1,pad2);
  mtf=fft(pma);
  
  model = (1.0/numberof(mtf))*double(fft(mtf*fft(nnx,[-1]),[+1]));
  

  residual = model - bdata;
  wr = roll(pW)*residual;

  ga=(2.0/numberof(mtf))*double(fft(conj(fft(nnx,[-1,0]))*(fft(wr,[-1,0])),[+1,0]));
  ga=ga(pad2+2:nr+pad2+1)(::-1);  // HOW STRANGE ??!!

  gx=2.*x*(ga(+)*sp(+,));
  gx+=  deconv_mux*dx;

  gz=ga(+)*(bg*t)(+,);
  gz*=(x^2);
  gz+=dz;

  // ge=ga(+)*(sp0*gt)(+,);
  //ge*=(x^2);
  ge= (sp01*ga)(+)*gt(+,);
  if (!is_void(pos)&&(is_array((i = where((!ebv)*(ge >= 0.0)))))) gebv(i) = 0.0; // clipping E(B-V)
  
  
  gv = (2.0/numberof(mtf))*double(fft(conj(mtf)*fft(wr,[-1,0]),[+1,0]));
  gv=gv(pad2+ni:pad2+nj,);
  gv*=2.*xv;
  gv+=deconv_muv*dxv;
  //gv+=dxv;

  grow,gx,gz,ge,gv;
  g=gx;

  //if (dbg==1) error;
  
  //if (positive >= 2 && is_array((i = where((!X)*(g > 0.0))))) g(i) = 0.0;
  return sum(wr*residual) +rglv+rglx+rglz+rgle;
 
};

func W1faz1cqq5(X,positive, &g)
{

  // Non Parametric Extinction Law. Like wqcfaz4 but with losvd search
  // quadratic in sfh and losvd,clipping in E(B-V)
  //needs _x0,deconv_muv,deconv_mux,deconv_mue,deconv_muz,mub,
  // _m,nab,pad1,pad2,muc,co
  // Here there is only one dust screen for the whole population like
  // everyone does with SDSS.
  // there is a penalization on x Pn(x,c) so that x satisfies x(sum)=c.
  // so to lift the degeneracy in norms between x and losvd and non parametric extinction law.

  
  extern wr,model,rglx,rglv,rgle,rglz;

  res=W1faz1cqq4(X,positive,g)+muc*Pn((X(:nab))^2,co,gc)+mucov*Pn((X(2*nab+nde+1:))^2,cov,gcov);
  
  gc=gc*2.*X(:nab);
  g(:nab)+=muc*gc;

  gcov=gcov*2.*X(2*nab+nde+1:);
  g(2*nab+nde+1:)+=mucov*gcov;
  
  return res;
};

func W1faz1cqq7(X,positive, &g)
{
  // like W1faz1cqq5 but with a prior in AMR
  // Non Parametric Extinction Law. Like wqcfaz4 but with losvd search
  // quadratic in sfh and losvd,clipping in E(B-V)
  //needs _x0,deconv_muv,deconv_mux,deconv_mue,deconv_muz,mub,
  // _m,nab,pad1,pad2,muc,co
  // Here there is only one dust screen for the whole population like
  // everyone does with SDSS.
  // there is a penalization on x Pn(x,c) so that x satisfies x(sum)=c.
  // so to lift the degeneracy in norms between x and losvd and non parametric extinction law.

  
  extern wr,model,rglx,rglv,rgle,rglz;
  
  res=W1faz1cqq5(X,positive,g);
  if (AMRp!="none") {
    res=res+muAMRp*P(X(nab+1:2*nab)-AMRp,rrAMRp,gp);
    g(nab+1:2*nab)=g(nab+1:2*nab)+muAMRp*gp;
  };
  return res;
};


func W1faz1cqq6(X,positive, &g)
{

  // Parametric Extinction Law. Like W1faz1cqq but with losvd search
  // quadratic in sfh and losvd,clipping in E(B-V)
  //needs _x0,deconv_muv,deconv_mux,deconv_mue,deconv_muz,mub,
  // _m,nab,pad1,pad2,muc,co
  // Here there is only one dust screen for the whole population like
  // everyone does with SDSS.
  // there is a penalization on the losvd Pn(x,c) so that x satisfies x(sum)=c.
  // so to lift the degeneracy in norms between x and losvd and non parametric extinction law.

  
  extern wr,model,rglx,rglv,rgle,rglz;

  res=W1faz1cqq(X,positive,g)+mucov*Pn((X(2*nab+2:))^2,cov,gcov);
  
  gcov=gcov*2.*X(2*nab+2:);
  g(2*nab+2:)+=mucov*gcov;
  
  return res;
};


func W1faz1cqq8(X,positive, &g)
{

  // like W1faz1cqq6 but with prior in AMR
  // Parametric Extinction Law. Like W1faz1cqq but with losvd search
  // quadratic in sfh and losvd,clipping in E(B-V)
  //needs _x0,deconv_muv,deconv_mux,deconv_mue,deconv_muz,mub,
  // _m,nab,pad1,pad2,muc,co
  // Here there is only one dust screen for the whole population like
  // everyone does with SDSS.
  // there is a penalization on the losvd Pn(x,c) so that x satisfies x(sum)=c.
  // so to lift the degeneracy in norms between x and losvd and non parametric extinction law.

  extern wr,model,rglx,rglv,rgle,rglz;
  
  res=W1faz1cqq6(X,positive,g);
  if (AMRp!="none") {
    res=res+muAMRp*P(X(nab+1:2*nab)-AMRp,rrAMRp,gp);
    g(nab+1:2*nab)=g(nab+1:2*nab)+muAMRp*gp;
  };
  return res;
};




  func W1faz1cqq2(X,positive, &g)
{

  // quadratic in sfh and losvd,clipping in E(B-V)
  //needs _x0,deconv_muv,deconv_mux,deconv_mue,deconv_muz,mub,
  // _m,nab,pad1,pad2,
  // Here there is only one dust screen for the whole population like
  // everyone does with SDSS.
  // ext is  polynomial to account for spectrophotometric error
  // derivative seems ok except for first polynomial coefficient ???!?
  
  extern wr,model,rglx,rglv,rgle,rglz;

  g=X*0.;
  x=X(:nab);
  z=X(nab+1:2*nab);
  ebv=X(2*nab+1:2*nab+3);
  if (positive) ebv=max(ebv,0.);   //E(B-v) clipping
  xv=X(2*nab+4:);
  t=npolf3(ebv,_x0,gt);
  
  
  rglx=deconv_mux*P(x^2,rr,dx);
  dx*=deconv_mux*2.*x;
  //rglv=deconv_muv*ip_roughness(xv^2,1,dxv,which=1,order=2);
  rglv=deconv_muv*P(xv^2,rr1,dxv);
  dxv*=deconv_muv*2.*xv;
  rgle=0.;
  de=0.;
  rglz=deconv_muz*P(z,rr3,dz)+mub*bound2(z,gbz);
  dz*=deconv_muz;
  dz+=mub*gbz;
  nnx=xtreat1(xv^2,nr+pad2+pad1);
  nxv=xtreat1(xv,nr+pad2+pad1);

  sp=buildb2(z,bg);
  sp0=sp;
  sp=t*sp;
  ma=sp(,+)*(x^2)(+);
  
  pma=pad(ma,pad1,pad2);
  mtf=fft(pma);
  
  model = (1.0/numberof(mtf))*double(fft(mtf*fft(nnx,-1),+1));
  

  residual = model - bdata;
  wr = roll(pW)*residual;

  ga=(2.0/numberof(mtf))*double(fft(conj(fft(nnx,[-1,0]))*(fft(wr,[-1,0])),[+1,0]));
  ga=ga(pad2+2:nr+pad2+1)(::-1);  // HOW STRANGE ??!!

  gx=2.*x*(ga(+)*sp(+,));
  gx+=  deconv_mux*dx;

  gz=ga(+)*(bg*t)(+,);
  gz*=(x^2);
  gz+=dz;

  //ge=ga(+)*(sp0*gt)(+,);
  //ge*=(x^2);

  ge=(((sp0(,+)*(x^2)(+))*ga)(+)*gt(+,));
  //gebv=2.*(gt(+,)*((rspb(,+)*(x^2)(+))*wresi)(+));

  
  //ge=ge(sum);
  if (!is_void(pos)&&(is_array((i = where((!ebv)*(ge >= 0.0)))))) gebv(i) = 0.0; // clipping E(B-V)
  
  
  gv = (2.0/numberof(mtf))*double(fft(conj(mtf)*fft(wr,[-1,0]),[+1,0]));
  gv=gv(pad2+ni:pad2+nj,);
  gv*=2.*xv;
  gv+=deconv_muv*dxv;

  grow,gx,gz,ge,gv;
  g=gx;
  
  //if (positive >= 2 && is_array((i = where((!X)*(g > 0.0))))) g(i) = 0.0;
  return sum(wr*residual) +rglv+rglx+rglz+rgle;
 
};


  func Rkin(X,positive, &g)
{

  // Tentative to get kinematics of observfations of higher resolutions than model
  // ie sigmavpop^2=sigmavmodel^2-sigmavlos^2
  // where the model is PEGASE at rest and we minimize data*los-model
  // quadratic in sfh and losvd,clipping in E(B-V)
  //needs _x0,deconv_muv,deconv_mux,deconv_mue,deconv_muz,mub,
  // _m,nab,pad1,pad2,
  // Here there is only one dust screen for the whole population like
  // everyone does with SDSS.
  
  extern wr,model,rglx,rglv,rgle,rglz;

  g=X*0.;
  x=X(:nab);
  z=X(nab+1:2*nab);
  ebv=X(2*nab+1);
  if (!is_void(pos)) ebv=max(ebv,0.);   //E(B-v) clipping
  xv=X(2*nab+2:);
  t=ds(ebv,_x0,gt);
    
  nnx=xtreat1(xv^2,nr+pad2+pad1);
  

  sp=buildb2(z,bg);
  model=(sp(,+)*(x^2)(+))*t;
  sp0=sp;
    
  pmodel=pad(model,pad1,pad2);
  pd=pad(d,pad1,pad2);
  mtf=fft(pd);
  
  bda = (1.0/numberof(mtf))*double(fft(mtf*fft(nnx,-1),+1));
  fresi=roll(pmodel(::-1)) -bda;
  wfresi=fresi*roll(pW(::-1));
  resi = pmodel - roll(bda(::-1));
  wr=pW*resi;
  wresi = wr(pad2+2:nr+pad2+1); 
  
  ki2=(wr*resi)(sum);
  rx=deconv_mux*P(x^2,rr,gPx);
  gPx*=deconv_mux*2.*x;
  rz=deconv_muz*P(z,rr3,gPz);
  rv=deconv_muv*P(xv^2,rr1,gPv);
  gPv*=deconv_muv*2.*xv;
  res=ki2+rx+rz+rv+mub*bound2(z,gbz);

  gkiz=2.*((bg*t(,-:1:numberof(x))*(x^2)(-:1:numberof(d),))(+,)*(wresi(,-))(+,)); //right
  gz=gkiz+deconv_muz*gPz+mub*gbz;
  gz=gz(,1);// z part 
  gx=2.*(((sp*(t(,-)))*(2.*x(-,)))(+,)*(wresi)(+))+gPx; // x part
  ge=2.*(gt(+)*((sp(,+)*(x^2)(+))*wresi)(+));
if (!is_void(positive)&&(!ebv*ge >= 0.0)) ge = 0.0; // clipping E(B-V)
  
  gv = -(2.0/numberof(mtf))*double(fft(conj(mtf)*fft(wfresi,[-1,0]),[+1,0]));
  gv=gv(pad2+ni:pad2+nj,); // probably have to roll(gv(::-1)) or st like that
  gv=roll(gv(::-1));
  gv*=2.*xv;
  gv+=gPv;


  grow,gx,gz,ge,gv;
  g=gx;

  //grow,1,2;
  return res;

  // Derivative seems to be ok
};

func tburp(ma,xv,positive, &g){
  // just a test
  extern wr,model;

  
  //if (positive) x = max(0.0, x);
  //xa=x(:nab);
  //xv=x(nab+1:);
  //rgla=deconv_mua*ip_roughness(xa,1,dxa,which=1,order=2);
  //rglv=deconv_muv*ip_roughness(xv,1,dxv,which=1,order=2);
  nnx=xtreat1(xv,nr+pad2+pad1);

  //ma=sp(,+)*xa(+);
  pma=pad(ma,pad1,pad2);
  mtf=fft(pma);
  
  model = (1.0/numberof(mtf))*double(fft(mtf*fft(nnx,-1),+1));
  

  residual = model - bdata;
  wr = roll(pW)*residual;

  ga=(2.0/numberof(mtf))*double(fft(conj(fft(nnx,[-1,0]))*(fft(wr,[-1,0])),[+1,0]));
  //ga=(2.0/numberof(mtf))*double(fft((fft(nnx,[-1,0]))*(fft(wr,[-1,0])),[+1,0]));
  //psp=padn2(sp,pad1,pad2);
  //rpsp=roll(psp,[dimsof(psp)(2)/2,0]);
  //for(i=1;i<=na;i++){rpsp(,i)=roll(psp(,i));};
  //rpsp=psp;
  //ga=ga(+)*rpsp(+,);
  //ga=ga(sum,) + deconv_mua*dxa;
  //  ga=ga + deconv_mua*dxa;

  //gv = (2.0/numberof(mtf))*double(fft(conj(mtf)*fft(wr,[-1,0]),[+1,0]));
  //gv=gv(pad2+ni:pad2+nj,);
  //gv +=deconv_muv*dxv;

  //grow,ga,gv;
  
  // if (positive >= 2 && is_array((i = where((!x)*(g > 0.0))))) g(i) = 0.0;
   g=ga(pad2+2:nr+pad2+1)(::-1);
   return sum(wr*residual); //+rglv+rgla;
 

  
};

//func f(x,arg,&g){g=arg*(x^(arg-1.));return x^arg;};
//func g(x){return x^3;};

func checkder(f,n,sol=,Ds=,arg=,noplot=,i1=,i2=) {
  /* DOCUMENT
     derivative check tool for scalar functions of several variables
     first comes analytic gradient and then numerical gradient
     the check is performed for the arguments i1 to i2
     will probably only work for vector argumnents, not matrices
     nder is red
  */
  

  local u1,Ds,nder,ader,u2,g,res;
  
  if(is_void(sol)) sol=1.e-1*abs(sin(indgen(n)/10.));
  if(is_void(Ds)) Ds=1.e-6;
  if(is_void(i1)) i1=1;
  if(is_void(i2)) i2=n;
  n=i2-i1+1;
  
  u1=sol;
  nder=u1*0.;
  ader=nder;
  for(i=i1;i<=i2;i++){
    u2=u1;
    u2(i)+=Ds;
    if(is_void(arg)) nder(i)=(f(u2)-f(u1,g))/Ds;
    if(!is_void(arg)) nder(i)=(f(u2,arg)-f(u1,arg,g))/Ds;
    ader(i)=g(i);
    //error;
  };

  //res=array(0.,n,2);
  //res(,1)=ader;
  //res(,2)=nder;

  res=ader;
  grow,res,nder;
  
  if(is_void(noplot)){
    ws;
    plh,ader;
    plh,nder,color="red";
    limits,i1,i2;
  };
  
  return res;

};
  


// check derivative

#if 0
u1=0.*sol+indgen(20)^3;
//u1(11:20)=sol(11:20);
Ds=1.e-4;
nder=u1*0.;
ader=nder;
for(i=1;i<=20;i++){
  u2=u1;
  u2(i)+=Ds;
  nder(i)=(P(u2,rr)-P(u1,rr,g))/Ds;
  ader(i)=g(i);
};
#endif

#if 0
u1=sol;
//u1(11:20)=sol(11:20);
ds=1.e-5;
nder=x0(,-:1:3);
ader=nder;
for(i=1;i<=3;i++){
  u2=u1;
  u2(i)+=ds;
  nder(,i)=(npolf3(u2,x0)-npolf3(u1,x0,g))/ds;
  ader(,i)=g(,i);
};
#endif


