// read synthetic spectra from my own runs of PEGASE-HR
#include "STECKMAP/Pierre/POP/sfit.i"

ll=exec("ls MILES/*.fits");

//ll=ll(:2);
nll=numberof(ll);
_m=[];
bloc=[];

zs=str2float(split2words(ll,sep="zt")(,2));
ages=str2float((split2words(strreplace(ll,".fits",""),sep="t")(,2)));

nm=6;
nages=46;

bloc=[];
for(i=1;i<=nll;i++){
  a=fits_read(ll(i),h);
  grow,bloc,a;
  nw=fits_get(h,"NAXIS1");
  wave=fits_get(h,"CRVAL1")+(float(indgen(nw))-float(fits_get(h,"CRPIX1")))*fits_get(h,"CDELT1");
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

  _a=indgen(numberof(*b(1)));
  ta=ages(:nages)*1.e3;
Res=(_x0(0)-_x0(1))/(_x0(dif)(avg));
 Rdlambda=_x0(dif)(avg) // smapling in angstroms
   _x0=wave;
 
 f=createb("MILES_SSPs.yor");
 save,f,bloc,_x0,_a,_m,ta,Res,Rdlambda;
 close,f;
};





