// read synthetic spectra from my own runs of PEGASE-HR
#include "STECKMAP/Pierre/POP/sfit.i"

//bak=convert_all("/Users/pedro/work/STECKMAP_SAV/Pat_coelho/coelho+07.ms_agb.p02p04.chab.9gyr.psb.fits");

//models=["agb","rgb","tpagb"];
allalphas=[0.,0.4];

for(ialpha=1;ialpha<=2;ialpha++){

alpha=allalphas(ialpha);

 if(alpha==0.)  ll=exec("ls coelho*p00.chab*.fits");
 if(alpha==0.4)  ll=exec("ls coelho*p04.chab*.fits");

 if(strmatch(ll(1),"ms_agb")) models="ms_agb";
 if(strmatch(ll(1),"ms_rgb")) models="ms_rgb";
 if(strmatch(ll(1),"ms_tpagb")) models="ms_tpagb";



//ll=ll(:2);
nll=numberof(ll);
_m=[];
_ms=[];
bloc=[];

zs1=split2words(ll,sep=".")(,3);
for(i=1;i<=nll;i++){
  zs=strcut(zs1(i),3)(1);
  sign=strmatch(zs,"p");
  if(sign==0) sign=-1.;
  zs=strreplace(zs,"m","");
  zs=strreplace(zs,"p","");
  zs=str2float(zs);
  zs=zs*sign;
  grow,_ms,zs;

  if(0){
  alphas=strcut(zs1(i),3)(2);
  sign=strmatch(alphas,"p");
  if(sign==0) sign=-1.;
  alphas=strreplace(alphas,"m","");
  alphas=strreplace(alphas,"p","");
  alphas=str2float(alphas);
  alphas=alphas*sign;
  grow,_alphas,alphas;
  };
  
 };
Fes=_ms*0.1;
_ms=_ms*0.1+(alpha==0?0.:0.3);

//ages=str2float((split2words(strreplace(ll,".fits",""),sep="t")(,2)));
ages=split2words(ll,sep=".g")(,6);
ages=str2float(ages); // in Gyr
//ages;

//error

nm=3;
nages=10;

bloc=[];
for(i=1;i<=nll;i++){
  a=fits_read(ll(i),h);
  grow,bloc,a;
  //b=fits_read(ll(i),hb,hdu=3);
  //nages=numberof(*b(1));
  nw=fits_get(h,"NAXIS1");
  wave=fits_get(h,"CRVAL1")+(float(indgen(nw))-float(fits_get(h,"CRPIX1")))*fits_get(h,"CDELT1");
  //  COM=fits_get(h,"COMMENT");
  //  write,COM(24);
  //  grow,_m,str2float(split2words(COM(24))(3));
  //  INFO,*b(1);
  //  INFO,wave;
  
};

nw=fits_get(h,"NAXIS1");
wave=fits_get(h,"CRVAL1")+(float(indgen(nw))-float(fits_get(h,"CRPIX1")))*fits_get(h,"CDELT1");


//rbloc=reform(bloc,nw,nages,nm);
 rbloc=reform(bloc,nw,nages*nm);
 ind=sort(ages-_ms*100.);
 rbloc=rbloc(,ind);
 rbloc=reform(rbloc,nw,nages,nm);
 rages=ages(ind);
 rages=reform(rages,nages,nm);

 rems=_ms(ind);
 rems=reform(rems,nages,nm);
 
 
//error

_m=0.02*10^(_ms);
_m=_m(1::nages);
rbloc=rbloc(,,sort(_m)(::-1));
_m=_m(sort(_m))(::-1);
bloc=rbloc;

//error

if(1){

  _a=indgen(nages);
  ta=rages(:nages)*1.e3; // in Myr
  //Res=(_x0(0)-_x0(1))/(_x0(dif)(avg));
  //  Res=1./(0.9/wave(avg));
  Res=1./wave(avg); // 1 angstrom resolution
  Rdlambda=wave(dif)(avg); // sampling in angstroms
  _x0=wave;
 
  f=createb("coelho07."+models+".alpha="+pr1(alpha)+".yor");
 save,f,bloc,_x0,_a,_m,ta,Res,Rdlambda;
 close,f;
};

 };



