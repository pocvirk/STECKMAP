#include "ascii2ls.i"

func tfilter(name,x0,s=){
  /* DOCUMENT
     return the transmission of known filter sampled on wavelengths x0 (interpolated) from the file filters.dat in PEGASE.2
     if name=="list" or void, returns the list of filters in the file
     BE CAREFUL TO FINISH THE FILTERS.DAT FILE WITH A CARRIAGE RETURN
  */
  if(name=="list") return(readfilter("list"));
  filt=readfilter(name);
  if(is_void(s)){
  if(x0(0)<filt(,1)(0)) write,"WARNING! domains not matched";
  if(x0(1)>filt(,1)(1)) write,"WARNING! domains not matched";
  };
  return interp(filt(,2),filt(,1),x0);
};


func readfilter(name) {
/* DOCUMENT reads the filter name from filters.dat
     returns a 2 column array, wavelength and transmission
     if name=="list",returns list of filters
     requires file "./DATA/filters.dat"
     bug fixed by Ariane Lancon
*/

  fn=steckmapdir+"Pierre/POP/DATA/filters.dat";
  nlines_fn=0;            // AL
  ndummy = sread(exec("wc "+fn), nlines_fn);   // AL : cf Eric/system.i
  f=open(fn);
  // a=rdline(f,10000);
  a=rdline(f,nlines_fn);   // AL
  if(name=="list") {
    is=where(strmatch(a,"'")==1);
    return a(is);
  };
  is=where(strmatch(a,name)==1)(1);

// Ariane's fix is troublesome for filters with stoopid names such as 300 (F300W) so only sensible fix is making more senible names, such as F300W
  is=where(strmatch(a,name)==1);                     // AL : keep allmatches
  iis =1;                                            // AL
  if(0){ // we'll take only the first match and expect the user to specify accurately the filter name he's looking for
      while (split2words(a(is(iis)),sep="'")(2) != name) {   // AL
          iis++;                                           // AL
      }                                                 // AL
      if (iis>numberof(is)) error, "ERROR: readfilter - filter name missing";
  };
  //is= is(iis);
  nl=int(str2float(split2words(a(is),sep=" ")(1)));
  //nl=0;                                  // AL
  //ndummy = sread(a(is), nl);             // AL, commented by PO
  filt1=array(0.,[2,nl,2]);
  for(i=1;i<=nl;i++){
    // filt1(i,)=str2float(split2words(a(is+i),sep=" ")(1:2));
    ndummy = sread(a(is+i),filt1(i,1),filt1(i,2));   // AL
  };
  return filt1;
};

func readfilter_old(name){
  /* DOCUMENT
     reads the filter name from filters.dat
     returns a 2 column array, wavelength and transmission
     if name=="list",returns list of filters
     requires file "./DATA/filters.dat"
  */
  fn=steckmapdir+"Pierre/POP/DATA/filters.dat";
  f=open(fn);
  a=rdline(f,10000);
  if(name=="list") {
    is=where(strmatch(a,"'")==1);
    return a(is);
  };
  is=where(strmatch(a,name)==1)(1);
  nl=int(str2float(split2words(a(is),sep=" ")(1)));
  filt1=array(0.,nl,2);
  for(i=1;i<=nl;i++){
    filt1(i,)=str2float(split2words(a(is+i),sep=" ")(1:2));
  };
  return filt1;
};
  
func fflux_old(sp,x0,name,s=){
  /* DOCUMENT
     returns the inverse of the flux of sp through filter name
     be careful that the spectral domain spans the whole bandpass of the filter
  */

  filt=tfilter(name,x0,s=s);
  // THIS IS A BIT LOWER
  return 1./((integ(sp*filt,x0,x0(0))/integ(filt,x0,x0(0)))*(x0(0)-x0(1)));
  //return 1./((integ(sp*filt,x0,x0(0))/integ(filt,x0,x0(0))));
  //return 1./integ(sp*filt,x0,x0(0));  // THAT SHOULD BE OK BUT IS A BIT LARGE
};

func fflux(sp,x0,name,s=){
  /* DOCUMENT
     returns the inverse of the flux of sp through filter name
     be careful that the spectral domain spans the whole bandpass of the filter
  */

  //INFO,x0;
  filt=tfilter(name,x0,s=s);
  // THIS IS A BIT LOWER
//  return 1./((integ(sp*filt,x0,x0(0))/integ(filt,x0,x0(0)))*(x0(0)-x0(1)));
  return 1./((integ(sp*filt,x0,x0(0))/integ(filt,x0,x0(0))));
  // careful this is problematic for flat_wide, since it has throughput of 1 on a very large domain // no because b.wave is given through x0 so its cut just to the size of the basis
  //return 1./integ(sp*filt,x0,x0(0));  // THAT SHOULD BE OK BUT IS A BIT LARGE
};



func bfflux(b,name,s=){
  /* DOCUMENT
     same as fflux for a basis
     be careful that domains match and basis is in correct units and not normalized
     returns an array
     not very clever: reads filters.dat many times
  */
  local nb,nm;
  bloc=b.flux;
  x0=b.wave;
  nb=dimsof(bloc)(3);
  nm=dimsof(bloc)(4);

  res=array(0.,nb,nm);
  for(i=1;i<=nb;i++){
    for(j=1;j<=nm;j++){
      res(i,j)=fflux(bloc(,i,j),x0,name,s=s);
    };
  };

  return res;
};
  

func bflux(b,name,s=){
  /* DOCUMENT
     same as fflux for a basis
     be careful that domains match and basis is in correct units and not normalized
     returns an array
     not very clever: reads filters.dat many times
  */
  local nb,nm;
  bloc=b.flux;
  x0=b.wave;
  nb=dimsof(bloc)(3);
  nm=dimsof(bloc)(4);

  res=array(0.,nb,nm);
  for(i=1;i<=nb;i++){
    for(j=1;j<=nm;j++){
      res(i,j)=fflux(bloc(,i,j),x0,name,s=s);
    };
  };

  return res;
};
  
  

if(0){
  //Masses for solar z from BC03, required to compute the M/LI correctly
  _mm=[];
  tabmasses=array(0.,numberof(bBC03.ages),6);
  //tabmasses=array(0.,220,6);

  for(i=6;i>=1;i--){
    fn="/raid/ocvirk/BC03/bc03/models/Padova1994/chabrier/bc2003_lr_m"+pr1(i+1)+"2_chab_ssp.4color";
    f=open(fn);
    a=rdline(f,249);
    iZ=where(strmatch(a,"Z="))(1);
    grow,_mm,str2float(split2words(a(iZ),sep="Z=I ")(-1));
    is=where(strmatch(a,"#"))(0)+1;
    labs=split2words(a(is-2),sep=" ");
    bigtab=str2float(split2words(a(is:),sep=" "));
    // masses in 7th position
    masses=interp(bigtab(,7),bigtab(,1),log10(bBC03.ages)+6);
    //masses=bigtab(,7);
    tabmasses(,i)=masses;
  };

  
  // with these definition MsL <- 
}

//bBC03=bRbasis(ager,nbins=nb,basisfile=models(2),dlambda=10,R=[],inte=inter,br=br,zr=1,wavel=[6500,13000],N=0,navg=1);write,"ok";write,dimsof(bBC03.flux);
//write,"WARNING! bBC03 modified";
//ML=bfflux(bBC03.flux,bBC03.wave,"IJ");
// with these definition MsL <- ML*tabmasses(tabmasses(,:5));
//MsL= ML*(tabmasses(,:5));
//sMsL=MsL; //save it !
