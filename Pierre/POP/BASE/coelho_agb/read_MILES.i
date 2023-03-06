// read synthetic spectra from my own runs of PEGASE-HR
#include "STECKMAP/Pierre/POP/sfit.i"

ll=exec("ls *.fits");

//ll=ll(:2);
nll=numberof(ll);
_m=[];
bloc=[];

zs=split2words(ll,sep="ZT")(,2);
sign=strmatch(zs,"p");
sign(where(sign==0))=-1.;
//zs=split2words(zs,sep="mp");
zs=strreplace(zs,"m","");
zs=strreplace(zs,"p","");
zs=str2float(zs);
zs=zs*sign;
//zs;
//error

//ages=str2float((split2words(strreplace(ll,".fits",""),sep="t")(,2)));
ages=split2words(ll,sep="T");
ages=strreplace(ages,".fits","")(,0);
ages=str2float(ages); // in Gyr
//ages;

//error

nm=7;
nages=50;

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

rbloc=reform(bloc,nw,nages,nm);
_m=0.02*10^(zs);
_m=_m(1::nages);
rbloc=rbloc(,,sort(_m)(::-1));
_m=_m(sort(_m))(::-1);
bloc=rbloc;



if(1){

  _a=indgen(nages);
  ta=ages(:nages)*1.e3;
  //Res=(_x0(0)-_x0(1))/(_x0(dif)(avg));
  Res=1./(0.9/wave(avg));
  Rdlambda=wave(dif)(avg); // sampling in angstroms
  _x0=wave;
 
 f=createb("MILES_SSPs_kroupa_universal_v8.0.yor");
 save,f,bloc,_x0,_a,_m,ta,Res,Rdlambda;
 close,f;
};





