// A few more tools for stellar populations study

#include "pop_paths.i"
#include "Qfunctions.i"
#include "spline.i"
#include "lmfit.i"
#include "colors.i"

func mycorrel(x,y){
  /* DOCUMENT
     normalizes x and y and then xycorrel
  */

  x-=fft_smooth(x,numberof(x)/5);
  y-=fft_smooth(y,numberof(y)/5);
    
  return xycorrel(x,y);
};


func setmus(snr,R=){
  /* DOCUMENT
     returns value for vector mus containing mu for each variable
     by now only should only be used for simulation purposes
     values putative
  */

  local snr,lum,lus,mux;
  
  if(is_void(snr)) return [2.e2, 1.e0,1.e2,1.e2,1.e4,1.e8];
  //lum=[10.,1.e-1,5.e-3,1.e-4]; //look up table for mux// used to be ok for nab=20;
  lum=[1.e1,1.e0,5.e-2,1.e-2];
  lus=[10.,25.,50.,100.];
  mux=interp(lum,lus,snr);
  return [2.e2, mux ,1.e2,1.e2,1.e4,1.e8];
};



func mkPENAL(ni,s=){
  /* DOCUMENT
     builds a penalization tikhonov(s="tikho"), gradient (s="grad") or laplacian (s="lap") defaut tikho so that P(X)=xtRtRx. consistent with function GCV
     prefer genREGUL
  */

  if (s=="tikho") return diag(array(1.,ni));
  if (s=="grad") {rr=genREGUL(ni);return rr;};
  if (s=="lap") {rr=genREGUL(ni); return rr(+,)*rr(+,);};

};

func genREGULAZ(na,nz,s=){
  /* DOCUMENT
     generates penalizations for 2d age-Z linear inversions
     na age bins
     nz metallcity bins
  */

  local L;

  l1=genREGUL(na,s=s);
  d=dimsof(l1);
  L=array(0.,d(2)*nz,d(3)*nz);
  for(i=1;i<=nz;i++){
    L(1+d(2)*(i-1):d(2)*i,1+d(3)*(i-1):d(3)*i)=l1;
  };
  return L;
};


func selage(A,s=,mit=){
  /* DOCUMENT
     returns a sequence of index ind such that for each i (A(,ind)(,dif))^2(sum,)>s, i.e. the distance between the elements of A(,ind) will be larger than s
  */

  local mit,s,nA,res;
  
  if (is_void(mit)) mit=1;
  if (is_void(s)) s=((A(,dif))^2)(sum,)(avg);

  nA=A;
  for(j=1;j<=mit;j++){
    na=dimsof(nA)(3);
    i1=1;
    i2=2;
    res=i1;
    for(i2=1;i2<=na;i2++){
      v=(((nA(,i1)-nA(,i2))^2)(sum)>=s);
      grow,res,((v==1)?i2:[]);
      i1=(v==1)?i2:i1;
    };
    nA=nA(,res);
  };
  return res;
  
}

func _GSO(A,Cm1=){

  /* DOCUMENT
     gram schmidt orthogonolization with variable metric C such that <x,y>=xtCy
     Cm1=inv(C). Can be computationnally very expensive because C can be big matrix (typically 10000*10000 in my case)
  */
  
  local na,Cm1;
  
  na=dimsof(A)(3);
  nd=dimsof(A)(2);
  if(is_void(Cm1)) Cm1=diag(array(1.,nd));
  
  B=A;  
  B(,1)=Cm1(,+)*(A(,1))(+);  // not normalized shit

  for(i=1;i<=na;i++){
    
    B(,i)=Cm1(,+)*(A(,i))(+);
    
    for(j=1;j<i;j++){
      
      B(,i)=B(,i)+(-(A(,i)(+)*B(,j)(+)))*B(,j);
    };

    B(,i)/=sqrt((B(,i)(+)*Cm1(+,))(+)*B(,i)(+));
        
  };

  return B;
};
  
func GSO(A,&v,Cm1=){

  /* DOCUMENT
     gram schmidt orthogonolization with diagonal metric given by Cm1 (metric vector instead of a metric matrix)
     additional output v contains projection of A onto B, that is passage matrix from A to B
     
  */

  local na;
  
  na=dimsof(A)(3);
  nd=dimsof(A)(2);
  if(is_void(Cm1)) Cm1=(array(1.,nd));

  
  B=A;
  B(,1)=A(,1)/sqrt((A(,1)*Cm1)(+)*A(,1)(+)); 

  for(i=1;i<=na;i++){
    
    B(,i)=Cm1*(A(,i));
    
    for(j=1;j<i;j++){
      
      B(,i)=B(,i)+(-(A(,i)(+)*B(,j)(+)))*B(,j);
    };

    B(,i)/=sqrt((B(,i)*Cm1)(+)*B(,i)(+));
        
  };

  v=B(+,)*A(+,);
  
  
  return B;
};
  

  

  
func bRbasis(ages,&FWHM,mets=,nbins=,R=,wavel=,N=,basisfile=,dlambda=,zr=,br=,inte=,navg=,list=,intages=){

  /* DOCUMENT

  bRbasis(list=1) will give a list of the available SSP models in basisdir
  
     builds a basis with ages (in yrs) between ages(1) and ages(2) (ages(1)>ages(2)), with nbins bins, averaging each ssp over the bin, as in the paper- should be better than buildb. Resolution (FWHM) is specified by R, and the spectral range by wavel (array of 2). Wanted metallicities in the basis are given by mi (array of 2, decreasing)
     Also builds an array MsL containing the M/L ratio of each age bin for each metallicity. A file with a vector ta containing the ages should be provided.
     When dlambda is given the basis is interpolated at dlambda.
     When dlambda is set to 0 the dlambda is automatically chosen as wavel(avg)/R.
     zr=1 rescales metallicities through Zrescale function defined in Qfunctions.i
     
     WARNING !! extremities of the spectrum are badly defined when R!=10000
     that's why br is for. br=1 -> automatic removal of wrong boundaries, over a range defined by deltalambda=(FWHM/2.) in final spectrum
     inte is the number of ages bins on which the original basis will be interpolated
     N=0 gives an unnormalized flux output (young is bright, old is faint)
     navg=1 inhibits the averaging process, so that output base.flux(,i,m) is just the instantaneous SSP at age i and not the average over (i-1/2,i+1/2).
     SET navg=1 and N=0 to get unnormalized basis

     intages=1 probably gives ages more reasonable than intages=0
     
     WARNING!! metallicities in decreasing order
     WARNING!! SOMETHING WRONG WITH br=1 sometimes
     CAREFUL!! Kroupa is default for bc03 ?, "PHR" is salpeter by default
  */

  extern MsL,bab,ibab,di;

  if (list==1) {exec("ls "+basisdir);return;};
  if (is_void(basisfile)) basisfile=basisdir+"ELO3.yor";
  if (basisfile=="LCB") basisfile=basisdir+"LCB_krou.yor";
  if (basisfile=="PHR") basisfile=basisdir+"ELO3_salp.yor";
  if (basisfile=="ELO3s") basisfile=basisdir+"ELO3_salp.yor";
  if (basisfile=="ELO3old") basisfile=basisdir+"ELO3.old";
  if (basisfile=="SPEED") basisfile=basisdir+"SPEED_kroupa.yor";
  if (basisfile=="bc03") basisfile=basisdir+"bc03_Pa94_Cha_raw.yor";
  if (basisfile=="BC03") basisfile=basisdir+"bc03_Pa94_Cha_raw.yor";
  if (basisfile=="GD05sp") basisfile=basisdir+"GD05_salpeter_padova.yor";
  if (basisfile=="GD05sg") basisfile=basisdir+"GD05_salpeter_geneva.yor";
  if (basisfile=="MILES") basisfile=basisdir+"MILES_SSP.yor";
  if (is_void(nbins)) nbins=10;
  if (is_void(mets)) mets=[6.e-2,4.e-4];
  if (is_void(ages)&(strmatch(basisfile,"ELO3"))) ages=[1.e7,2.e10];
  if (is_void(ages)) ages="full";
  //if (is_void(wavel)) wavel=[4000.,6700.];
  if (is_void(dlambda)) dlambda=[]; 
  if (!is_void(R)) {if (R>10000.) {write, "you dreamin', man !";error;};};
  if (is_void(N)) N=1;
  if(is_void(navg)) navg=0;
  if(is_void(zr)) zr=1;
  if(is_void(intages)) intages=0;
  
  upload,basisfile,s=1;
  dib=dimsof(bloc);
  dd=dib(2); 
  if(R==Res) R=[];

  
  
  // TEST PART what to do if wavel requested exceeds actual domain of the model

  if(is_void(wavel)) wavel=[_x0(1),_x0(0)];
  if(!is_void(wavel)){
    wavel=wavel;
  
  wavel(1)=max(wavel(1),_x0(1));
  wavel(2)=min(wavel(2),_x0(0));
  
  wai=where((_x0<=wavel(2))&(_x0>=wavel(1)));
  bloc=bloc(wai,,);
  _x0=_x0(wai);
  dib=dimsof(bloc);
  dd=dib(2);
  };

  //**********************************************************
  // interpolate basis in log age range (because ages are not log-sampled in ELO3.yor)
  //**********************************************************
 if (!is_void(inte)){
   nta=spanl(ta(2),ta(0),int(inte));
   //ibloc=array(0.,dib(2),int(inte),dib(4));
   ibloc=interp(bloc,ta,nta,2);
   ta=nta;
   bloc=ibloc;
 };
 
 dib=dimsof(bloc);
 dd=dib(2); 

 if(ages=="full") ages=[ta(2),ta(0)]*1.e6;
 
 //********************************************************
 // FIRST STEP, form the vector ab from the bins vector bab
 //********************************************************
 
 bab=10^(span(log10(ages(1)),log10(ages(2)),nbins+1)); //classic age bins
 //    bab=10^(span(log10(2.e9),log10(8.e9),nab+1));
 fibab=abs((ta*1.e6)(,-:1:numberof(bab))-bab(-:1:numberof(ta),));
 ibab=fibab(mnx,);
 // if (intages==1) ibab=nta;
 
 //*************************************************
 // find index in bloc of the required metallicities
 //*************************************************
 imet=where((_m<=mets(1))&(_m>=mets(2)));

 //************************************************
 // Compute FWHM (pix) to convolve with to obtain R
 //************************************************
if ((!is_void(R))&&(R<10000.)){
  FWHM1=(_x0(avg)/Res)/(_x0(dif)(avg));
  FWHM2=(wavel(avg)/R)/(_x0(dif)(avg));
  FWHM=sqrt((FWHM2^2-FWHM1^2));                        
};

 asp=bloc(:dd,ibab(:-1),imet); // just to form the table
 MsL=array(0.,nbins,dimsof(asp)(4));

 //***********************************************
 // compute M/L ratio for each bin and metallicity
 // NB: the mass loss of the initial population is
 // not taken into account in this computation.
 // i.e. the mass remains 1 Msolar
 // the L is the flux integrated in the whole
 // spectral domain considered.
 // it's easy to get the magnitudes through filter
 // by using a function I dont remember the name of
 //***********************************************
 
 for(i=1;i<=nbins;i++){
   MsL(i,)=1./(integ(bloc(sum,ibab(i):ibab(i+1),imet),ta( ibab(i):ibab(i+1)),ta( ibab(i):ibab(i+1))(0),1)/(ta(ibab(i+1))-ta(ibab(i))));
 };
 
 if (N==1) bloc(:dd,,)/=bloc(:dd,,)(avg,,)(-:1:dd,,);   // Normalize
 
 for(i=1;i<=nbins;i++){
   if (navg!=1) asp(,i,)=integ(bloc(:dd,ibab(i):ibab(i+1),imet),ta( ibab(i):ibab(i+1)),ta( ibab(i):ibab(i+1))(0),2);
   if (!is_void(R)) {for(j=1;j<=numberof(imet);j++)

     {asp(,i,j)=fft_smooth(asp(,i,j),FWHM);write,i,j;};};
 };

 if (N==1) asp/=asp(avg,,)(-:1:dd,,);   // Normalize

 if (!is_void(dlambda)) {
   dlambda=dlambda==0?(wavel(avg))/(2.*(is_void(R)?Res:R)):dlambda;
   ndd=int((wavel(0)-wavel(1))/dlambda);
   nwave=wavel(1)+dlambda*indgen(ndd);
   nasp=array(0.,ndd,nbins,numberof(imet));
   for(i=1;i<=nbins;i++){
     for(j=1;j<=numberof(imet);j++){
       nasp(,i,j)=interp(asp(,i,j),_x0,nwave);
     };
   };
   asp=nasp;
   _x0=nwave
 };
 
   // CUT BLUE AND RED BOUNDARIES

   if (!is_void(br)){
     if(is_void(dlambda)) dlambda=(_x0(dif))(avg);
   np=int((wavel(0)+wavel(1))/(2.*(is_void(R)?Res:R))*(1./dlambda))+1;
   asp=asp(np:-np,,);
   };
   ndd=dimsof(asp)(2);
   
   
   if (N==1) asp/=asp(avg,,)(-:1:ndd,,);   // Normalize
   dd=ndd;

 di=[dd,nbins,numberof(imet)];
 include, "base_struct.i",1;
 
 base=basStruct();
 base.filename=basisfile;
 base.flux=asp;
 base.wave=_x0(:dd);
 base.nages=int(ibab(:-1)+ibab(dif)/2);
 base.ages=float(ta)(base.nages);
 if(intages==1) base.ages=10^(span(log10(ages(1)),log10(ages(2)),nbins));
 base.met=is_void(zr)?_m(imet):Zrescale(_m(imet));
 base.R=is_void(R)?Res:R;

 return base;
};


func bRbasis2(ages,&FWHM,&bab,mets=,nbins=,R=,wavel=,N=,basisfile=,dlambda=,zr=,br=,inte=,navg=,list=,intages=,dirsfr=){

  /* DOCUMENT
  Prefer this one over bRbasis
  bRbasis(list=1) will give a list of the available SSP models in basisdir
  
     builds a basis with ages (in yrs) between ages(1) and ages(2) (ages(1)>ages(2)), with nbins bins, averaging each ssp over the bin, as in the paper- should be better than buildb. Resolution (FWHM) is specified by R, and the spectral range by wavel (array of 2). Wanted metallicities in the basis are given by mi (array of 2, decreasing)
     dirsfr=1 -> basis is in Msol/yr

     Also builds an array MsL containing the M/L ratio of each age bin for each metallicity. A file with a vector ta containing the ages should be provided.
     When dlambda is given the basis is interpolated at dlambda.
     When dlambda is set to 0 the dlambda is automatically chosen as wavel(avg)/R.
     zr=1 rescales metallicities through Zrescale function defined in Qfunctions.i
     
     WARNING !! extremities of the spectrum are badly defined when R!=10000
     that's why br is for. br=1 -> automatic removal of wrong boundaries, over a range defined by deltalambda=(FWHM/2.) in final spectrum
     inte is the number of ages bins on which the original basis will be interpolated
     N=0 gives an unnormalized flux output (young is bright, old is faint)
     navg=1 inhibits the averaging process, so that output base.flux(,i,m) is just the instantaneous SSP at age i and not the average over (i-1/2,i+1/2).
     SET navg=1 and N=0 to get unnormalized basis

     intages=1 probably gives ages more reasonable than intages=0
     also returns the bins vector bab (size nbins+1) in bab
     DO (10^bab)(dif) to get the time elapsed in each bin in years
     
     WARNING!! metallicities in decreasing order
     WARNING!! SOMETHING WRONG WITH br=1 sometimes
     CAREFUL!! Kroupa is default for bc03 ?, "PHR" is salpeter by default
     CAREFUL!! CHECK THE UNITS OF b.ages:  Myr or yr ??
  */

  extern ibab,di;

  if (list==1) {exec("ls "+basisdir);return;};
  if (is_void(basisfile)) basisfile=basisdir+"ELO3.yor";
  if (basisfile=="LCB") basisfile=basisdir+"LCB_krou.yor";
  if (basisfile=="PHR") basisfile=basisdir+"ELO3_salp.yor";
  if (basisfile=="ELO3s") basisfile=basisdir+"ELO3_salp.yor";
  if (basisfile=="ELO3old") basisfile=basisdir+"ELO3.old";
  if (basisfile=="SPEED") basisfile=basisdir+"SPEED_kroupa.yor";
  if (basisfile=="bc03") basisfile=basisdir+"bc03_Pa94_Cha_raw.yor";
  if (basisfile=="BC03") basisfile=basisdir+"bc03_Pa94_Cha_raw.yor";
  if (basisfile=="GD05sp") basisfile=basisdir+"GD05_salpeter_padova.yor";
  if (basisfile=="GD05sg") basisfile=basisdir+"GD05_salpeter_geneva.yor";
  if (basisfile=="MILES") basisfile=basisdir+"MILES_SSP.yor";
  if (is_void(nbins)) nbins=10;
  if (is_void(mets)) mets=[6.e-2,4.e-4];
  if (is_void(ages)&(strmatch(basisfile,"ELO3"))) ages=[1.e7,2.e10];
  if (is_void(ages)) ages="full";
  //if (is_void(wavel)) wavel=[4000.,6700.];
  if (is_void(dlambda)) dlambda=[]; 
  if (!is_void(R)) {if (R>10000.) {write, "you dreamin', man !";error;};};
  if (is_void(N)) N=1;
  if(is_void(navg)) navg=0;
  if(is_void(zr)) zr=1;
  if(is_void(intages)) intages=0;
  if(is_void(ninter)) ninter=10;
  
  upload,basisfile,s=1;
  dib=dimsof(bloc);
  dd=dib(2); 
  if(R==Res) R=[];
  
  // TEST PART what to do if wavel requested exceeds actual domain of the model

  if(is_void(wavel)) wavel=[_x0(1),_x0(0)];
  if(!is_void(wavel)){
    wavel=wavel;
  
  wavel(1)=max(wavel(1),_x0(1));
  wavel(2)=min(wavel(2),_x0(0));
  
  wai=where((_x0<=wavel(2))&(_x0>=wavel(1)));
  bloc=bloc(wai,,);
  _x0=_x0(wai);
  dib=dimsof(bloc);
  dd=dib(2);
  };

  //**********************************************************
  // interpolate basis in log age range (because ages are not log-sampled in ELO3.yor)
  //**********************************************************
   if (!is_void(inte)){
   nta=spanl(ta(2),ta(0),int(inte));
   ibloc=array(0.,dib(2),int(inte),dib(4));
   ibloc=interp(bloc,ta,nta,2);
   ta=nta;
   bloc=ibloc;
   };
 
 dib=dimsof(bloc);
 dd=dib(2); 

 if(ages=="full") ages=[ta(2),ta(0)]*1.e6;
 
 //********************************************************
 // FIRST STEP, form the vector ab from the bins vector bab
 //********************************************************
 
 //bab=10^(span(log10(ages(1)),log10(ages(2)),nbins+1)); //classic age bins
 delta=log10(ages(2)/ages(1))/(nbins-1);
 bab=-delta/2+delta*(indgen(nbins+1)-1)+log10(ages(1));
 
 //fibab=abs((ta*1.e6)(,-:1:numberof(bab))-bab(-:1:numberof(ta),));
 //bab=fibab(mnx,);
 // if (intages==1) ibab=nta;
 
 //*************************************************
 // find index in bloc of the required metallicities
 //*************************************************
 imet=where((_m<=mets(1))&(_m>=mets(2)));

 //************************************************
 // Compute FWHM (pix) to convolve with to obtain R
 //************************************************
if ((!is_void(R))&&(R<10000.)){
  FWHM1=(_x0(avg)/Res)/(_x0(dif)(avg));
  FWHM2=(wavel(avg)/R)/(_x0(dif)(avg));
  FWHM=sqrt((FWHM2^2-FWHM1^2));                        
};

 asp=bloc(:dd,indgen(nbins),imet); // just to form the table
 // MsL=array(0.,nbins,dimsof(asp)(4));

 //***********************************************
 // compute M/L ratio for each bin and metallicity
 // NB: the mass loss of the initial population is
 // not taken into account in this computation.
 // i.e. the mass remains 1 Msolar
 // the L is the flux integrated in the whole
 // spectral domain considered.
 // it's easy to get the magnitudes through filter
 // by using a function I dont remember the name of
 //***********************************************
 
 for(i=1;i<=nbins;i++){
   //   MsL(i,)=1./(integ(bloc(sum,ibab(i):ibab(i+1),imet),ta( ibab(i):ibab(i+1)),ta( ibab(i):ibab(i+1))(0),1)/(ta(ibab(i+1))-ta(ibab(i))));
   
 };
 
 if (N==1) bloc(:dd,,)/=bloc(:dd,,)(avg,,)(-:1:dd,,);   // Normalize
 
 for(i=1;i<=nbins;i++){
   nsup=span(bab(i),bab(i+1),ninter);
   //   write,nsup;
   ibloc=interp(bloc(:dd,,imet),log10(ta+1.e-6)+6,nsup,2);
   asp(,i,)=integ(ibloc,nsup,nsup(0),2)/(nsup(0)-nsup(1));
   //if (navg!=1) asp(,i,)=integ(bloc(:dd,ibab(i):ibab(i+1),imet),ta( ibab(i):ibab(i+1)),ta( ibab(i):ibab(i+1))(0),2);
   if (!is_void(R)) {
     for(j=1;j<=numberof(imet);j++){
       asp(,i,j)=fft_smooth(asp(,i,j),FWHM);
       write,i,j;
     };
   };
 };

 if(db==1) error;
 
 if (N==1) asp/=asp(avg,,)(-:1:dd,,);   // Normalize

 if (!is_void(dlambda)) {
   dlambda=dlambda==0?(wavel(avg))/(2.*(is_void(R)?Res:R)):dlambda;
   ndd=int((wavel(0)-wavel(1))/dlambda);
   nwave=wavel(1)+dlambda*indgen(ndd);
   nasp=array(0.,ndd,nbins,numberof(imet));
   for(i=1;i<=nbins;i++){
     for(j=1;j<=numberof(imet);j++){
       nasp(,i,j)=interp(asp(,i,j),_x0,nwave);
     };
   };
   asp=nasp;
   _x0=nwave
 };
 
   // CUT BLUE AND RED BOUNDARIES

   if (!is_void(br)){
     if(is_void(dlambda)) dlambda=(_x0(dif))(avg);
   np=int((wavel(0)+wavel(1))/(2.*(is_void(R)?Res:R))*(1./dlambda))+1;
   asp=asp(np:-np,,);
   };
   ndd=dimsof(asp)(2);
   
   
   if (N==1) asp/=asp(avg,,)(-:1:ndd,,);   // Normalize
   dd=ndd;

 di=[dd,nbins,numberof(imet)];
 include, "base_struct.i",1;
 
 base=basStruct();
 base.filename=basisfile;
 base.flux=asp;
 if(dirsfr==1) base.flux=asp*((((10^bab)(dif))(-:1:numberof(_x0),,-:1:numberof(_met))));
 base.wave=_x0(:dd);
 base.nages=indgen(nbins);
   //int(ibab(:-1)+ibab(dif)/2);
 //base.ages=float(ta)(base.nages);
 //if(intages==1) base.ages=10^(span(log10(ages(1)),log10(ages(2)),nbins));
 base.ages=(10^(bab(zcen)))/1.e6;
 base.met=is_void(zr)?_m(imet):Zrescale(_m(imet));
 base.R=is_void(R)?Res:R;
 base.age_unit="yr";
 base.basistype="normalized, for flux fractions";
 if (dirsfr==1) base.basistype="flux/(Msol/yr)";
 if ((N==0)&(dirsfr==0)) base.basistype="flux/Msol";

 return base;
};


func _bRbasis(ages,mets=,nbins=,R=,wavel=,N=,basisfile=,dlambda=,zr=,br=,inte=,navg=){

  /* DOCUMENT
     builds a basis with ages (in yrs) between ages(1) and ages(2) (ages(1)>ages(2)), with nbins bins, averaging each ssp over the bin, as in the paper- should be better than buildb. Resolution (FWHM) is specified by R, and the spectral range by wavel (array of 2). Wanted metallicities in the basis are given by mi (array of 2, decreasing)
     Also builds an array MsL containing the M/L ratio of each age bin for each metallicity. A file with a vector ta containing the ages should be provided.
     When dlambda is given the basis is interpolated at dlambda.
     When dlambda is set to 0 the dlambda is automatically chosen as wavel(avg)/R.
     zr=1 rescales metallicities through Zrescale function defined in Qfunctions.i
     
     WARNING !! extremities of the spectrum are badly defined when R!=10000
     that's why br is for. br=1 -> automatic removal of wrong boundaries, over a range defined by deltalambda=(FWHM/2.) in final spectrum
     inte is the number of ages bins on which the original basis will be interpolated
     N=0 gives an unnormalized flux output (young is bright, old is faint)
     navg=1 inhibits the averaging process, so that output base.flux(,i,m) is just the instantaneous SSP at age i and not the average over (i-1/2,i+1/2).
     SET navg=1 and N=0 to get unnormalized basis

     WARNING!! metallicities decreasing
     

     
  */

  extern MsL,bab,ibab,di;

if (is_void(basisfile)) basisfile="~/perso/modeles/galaxie/CIN+POP2/ELO3.yor";
if (is_void(nbins)) nbins=10;
if (is_void(mets)) mets=[6.e-2,4.e-4];
if (is_void(ages)) ages=[50.e6,2.e10];
if (is_void(wavel)) wavel=[4000.,6690.];
 if (is_void(dlambda)) dlambda=0.; 
 if (!is_void(R)) {if (R>10000.) {write, "you dreamin', man !";error;};};
 if (is_void(N)) N=1;
 if(is_void(navg)) navg=0;
  
upload,basisfile,s=1;
dib=dimsof(bloc);
dd=dib(2); 

 if(R==Res) R=[];
 
// interpolate basis in log age range (because ages are not log-sampled in ELO3.yor)
 if (!is_void(inte)){
   nta=spanl(ta(2),ta(0),int(inte));
   ibloc=array(0.,dib(2),int(inte),dib(4));
   ibloc=interp(bloc,ta,nta,2);
   ta=nta;
   bloc=ibloc;
 };

 wai=where((_x0<=wavel(2))&(_x0>=wavel(1)));
 bloc=bloc(wai,,);
 _x0=_x0(wai);
 dib=dimsof(bloc);
 dd=dib(2);
 
// FIRST STEP, form the vector ab from the bins vector bab
 bab=10^(span(log10(ages(1)),log10(ages(2)),nbins+1)); //classic age bins
 //    bab=10^(span(log10(2.e9),log10(8.e9),nab+1));
 fibab=abs((ta*1.e6)(,-:1:numberof(bab))-bab(-:1:numberof(ta),));
 ibab=fibab(mnx,);

// find index in bloc of the required metallicities
imet=where((_m<=mets(1))&(_m>=mets(2)));

// Compute FWHM (pix) to convolve with to obtain R
if ((!is_void(R))&&(R<10000.)){
  FWHM1=(((_x0(0)+_x0(1))/2.)/Res)/(_x0(2)-_x0(1));
  FWHM2=(((wavel(0)+wavel(1))/2.)/R)/(_x0(2)-_x0(1));
  FWHM=sqrt(FWHM2^2-FWHM1^2);
};

 asp=bloc(:dib(2),ibab(:-1),imet); // just to form the table
 MsL=array(0.,nbins,dimsof(asp)(4));

 // compute M/L ratio for each bin and metallicity
 for(i=1;i<=nbins;i++){
   MsL(i,)=1./(integ(bloc(sum,ibab(i):ibab(i+1),imet),ta( ibab(i):ibab(i+1)),ta( ibab(i):ibab(i+1))(0),1)/(ta(ibab(i+1))-ta(ibab(i))));
 };
 
 if (N==1) bloc(:dd,,)/=bloc(:dd,,)(avg,,)(-:1:dd,,);   // Normalize
 
 for(i=1;i<=nbins;i++){
   if (navg!=1) asp(,i,)=integ(bloc(:dd,ibab(i):ibab(i+1),imet),ta( ibab(i):ibab(i+1)),ta( ibab(i):ibab(i+1))(0),2);
   if (!is_void(R)) {for(j=1;j<=numberof(imet);j++)

     {asp(,i,j)=fft_smooth(asp(,i,j),FWHM);};};
 };

 if (N==1) asp/=asp(avg,,)(-:1:dd,,);   // Normalize

 if (!is_void(dlambda)) {
   dlambda=dlambda==0?(wavel(0)+wavel(1))/(2.*(is_void(R)?Res:R)):dlambda;
   ndd=int((wavel(0)-wavel(1))/dlambda);
   nwave=wavel(1)+dlambda*indgen(ndd);
   nasp=array(0.,ndd,nbins,numberof(imet));
   for(i=1;i<=nbins;i++){
     for(j=1;j<=numberof(imet);j++){
       nasp(,i,j)=interp(asp(,i,j),_x0,nwave);
     };
   };
   asp=nasp;
   
   // CUT BLUE AND RED BOUNDARIES

   if (!is_void(br)){
   np=int((wavel(0)+wavel(1))/(2.*(is_void(R)?Res:R))*(1./dlambda))+1;
   asp=asp(np:-np,,);
   };
   ndd=dimsof(asp)(2);
   
   
   if (N==1) asp/=asp(avg,,)(-:1:ndd,,);   // Normalize
   dd=ndd;
   _x0=nwave
 };
 

 di=[dd,nbins,numberof(imet)];
 include, "teststruct.i",1;
 
 base=basStruct();
 base.filename=basisfile;
 base.flux=asp;
 base.wave=_x0(:dd);
 base.nages=int(ibab(:-1)+ibab(dif)/2);
 base.ages=ta(base.nages);
 base.met=is_void(zr)?_m(imet):Zrescale(_m(imet));
 base.R=is_void(R)?Res:R;

 return base;
};


func MSL(base,filter){
  /*DOCUMENT
    returns MsL for a given basis by integrating it through the specified filter.
    default is Bessell_V
  */
  if(filter=="list") return readfilter("list");
  if(is_void(filter)) filter="V_B90";
  return bfflux(base.flux,base.wave,filter);
};



func countlines(u,&wid,&w,&isc,s=){
  /* DOCUMENT
     result is a vector containing absolute value of all ups and downs of u
     takes value -100 if only one sign change is detectd
     do a histo1d to see the distributions of the depth of the lines in a spectrum
     wid is the width where the sign change takes place, takes value -100 when only one sign change is detected.
     w is 1 where sign change happens and 0 elsewhere
     isc is a list of indexes where sign change is detected
     s is threshold for sign change detection
  */

  if (is_void(s)) s=-1.e-70;
  
  du=u(dif);
  ddu=du(2:)*du(:-1); //derivative sign change detection
  isc=where(ddu<=s);
  if(is_void(dimsof(isc))) {w=[]; return -100;};
  isc=isc+1;
  wid=numberof(isc)>1?isc(dif):-100.;
  uisc=u(isc);
  res=numberof(isc)>1?abs(uisc(dif)):-100;
  w=array(0.,numberof(u));
  w(isc)=1.;
  return res;

};


func makebump(n,in,w,h=,N=){
  /* DOCUMENT
     makes a vector of 0s of length n with bumps of width w at indices in
of height h  */

  local v,nb,iv,i;
  
  nb=numberof(in);
  if (numberof(w)==1) w=w(1)(-:1:nb);
  if (is_void(h)) h=1.(-:1:nb);
  
  v=array(0.,n);
  iv=indgen(n);
  for(i=1;i<=nb;i++){v+=h(i)*exp(-((iv-in(i))^2/(2.*w(i)^2)));};
  if(N==1) v/=v(sum);
  return v;
};


func findlmax(q,&w,n=){
  /* document
     finds all local maxima including boarders */
  
  local t,r,s,i,ma,w,n;
  
  w=q*0.;
  // ok for interior
  t=countlines(q,r,s,i);
  if(is_void(s)) {w=array(0.,numberof(q)); ma=[];};

  if (!is_void(s)){
  ma=where((q(dif)>=0)&(s(2:)==1));
  ma=is_void(dimsof(ma))?[]:ma+1;
  //ma=where((q(dif)>=0)&(s(:-1)==1));
  };
  
    // check boarders
    if(q(dif)(1)<0.) grow,ma,1;
  if(q(dif)(0)>0.) grow,ma,numberof(q);

  w(ma)=1.;
  if(is_void(n)) n=numberof(ma);
  n=min(n,numberof(ma));
  
  return ma(sort(ma))(:n);
};

func _bumpweights(q,&lb){
  /* DOCUMENT
     divides a SAD in bumps and returns their weights and location.
     buggy
  */

  local lmax,t,r,s,i,u,wb,lb;
  
  //if(is_void(n)) n=2;
  lmax=findlmax(q,w);
  //a12=findlmax(q)(sort(q(findlmax(q)))(:-n+1:-1)); this selects only the 2 highest bumps we dont want that
  t=countlines(q,r,s,i);
  u=(w+!s>=1.);
  //if(q(dif)(1)<=0.) grow,ma,1;
  //if(q(dif)(0)>=0.) grow,ma,numberof(q);
  u(1)=0;
  u(0)=0;
  wb=array(0.,numberof(lmax));
  lb=array(0.,2,numberof(lmax));
  nu=numberof(where(u==0));
  for(i=1;i<=numberof(lmax);i++){
    lb(1,i)=where(u==0)(i);
    lb(2,i)=min(numberof(q),where(u==0)(min(nu,i+1)));
    wb(i)=max(0.,q(lb(1,i):lb(2,i)))(sum);
  };
  
  return wb;
};

func bumpweights(q,&lb){
  /* DOCUMENT
     divides a SAD in bumps and returns their weights and location.
     seems to work better than _bumpweights
  */

  local lmax,t,r,s,i,u,wb,lb;
  
  //if(is_void(n)) n=2;
  lmax=int(findlmax(q,w));
  smax=lmax;

  if(lmax(1)!=1) {smax=1;grow,smax,lmax;};
  if(lmax(0)!=numberof(q)) {grow,smax,numberof(q);}

  lmin=[];
  for(i=1;i<=numberof(smax)-1;i++){
    grow,lmin,(q(smax(i):smax(i+1))(mnx)-1+smax(i));
  };

  smin=lmin;
  if(lmin(1)!=1) {smin=1;grow,smin,lmin;};
  if(lmin(0)!=numberof(q)) {grow,smin,numberof(q);}
  

  
  wb=array(0.,numberof(lmax));
  lb=array(0,2,numberof(lmax));
  for(i=1;i<=numberof(lmax);i++){
    lb(1,i)=smin(i);
    lb(2,i)=smin(i+1);
    wb(i)=max(0.,q(lb(1,i):lb(2,i)))(sum);
  };
  
  return wb;
};




func findsub(q,b,n=,a=){
  /* DOCUMENT
     returns the LWA of the n most luminous sub-populations found in the solution q
     requires specification of a basis b to compute the LWA
     LWA specified in Qfunctions.i
     by default result sorted in increasing luminous weight
     if a=1 sorted in increasing age
     
  */

  local n,a,res,wb,lb,n,i;

  if(is_void(n)) n=2;
  if(is_void(a)) a=0;
  res=array(-100.,n);
  wb=bumpweights(q,lb);
  lb=lb(,where(wb>0.));
  wb=wb(where(wb>0.));
  n=min(n,numberof(wb));
  lb=lb(,sort(wb)(-(n-1):));
   wb=wb(sort(wb))(-(n-1):);
  
  res=array(0.,n);
  for(i=1;i<=n;i++){res(i)=LWA(q,lb(1,i),lb(2,i),b=b);};
  if (a==1) res=res(sort(res)(::-1));
  return res;
};


func crosschi2(sp,snr,&v,nodis=,N=,vages=){

  /* DOCUMENT
     crosschi2(sp,snr,&v,N=)
  */


  local ndl,bloc;
  bloc=sp;
  ndl=dimsof(bloc)(3);
  nlambda=dimsof(bloc)(2);
v=array(0.,ndl,ndl);
 al=v; 
if (!is_void(N)) blocn=bloc/bloc(sum,,)(-:1:dimsof(bloc)(2),,); //normalize
else  blocn=bloc; 

 for(i=1;i<=dimsof(blocn)(3);i++){
  for(j=1;j<=dimsof(blocn)(3);j++){
    alpha=(1./(((blocn(,j)^2)/blocn(,i)^2)(sum)))*((blocn(,j)/blocn(,i))(sum));
    v(i,j)=(((blocn(,i)-alpha*blocn(,j))^2)/(blocn(,i))^2)(sum);
  };
};

 v=v*snr^2;
 k2_90=0.5*(1.2816+sqrt(2*nlambda-1))^2;

 if (is_void(nodis)){
 ws;
 //plk,log(max(v,0.1*k2_90))/log(10);
 //plk,v;
 x=(is_void(vages)?indgen(ndl):vages)(,-:1:ndl);
 y=(is_void(vages)?indgen(ndl):vages)(-:1:ndl,);
 plf,log10(max(v,0.1*k2_90)),x,y;
 plc,v,y,x,levs=k2_90*[1.],width=3;
 //if(!is_void(vages)) logxy,1,1; 
 };
 return k2_90;
};

func Picard(b,spb){
/* DOCUMENT
     Picard(b,spb) plots discrete picard condition. b:data.
     also returns spb conditioning number.
  */
  s=SVdec(spb,u,v);
  co1=u(+,)*b(+);
  co2=co1/s;
  plh,log10(s),color="blue";
  plh,log(abs(co1)),color="green";
  plh,log10(abs(co2)),color="red";
  s=SVdec(spb(+,)*spb(+,));
  return log10(max(s)/min(s));
};

func degaz(sp,blocn,snr,&v,nodis=,N=){
/* DOCUMENT
   Representation of age-metallicity degeneracy
*/

  nda=dimsof(blocn)(3);
  ndz=dimsof(blocn)(4);
  nlambda=dimsof(blocn)(2);
  v=array(0.,nda,ndz);


 for(i=1;i<=nda;i++){
  for(j=1;j<=ndz;j++){
    alpha=(1./(((blocn(,i,j)^2)/sp^2)(sum)))*((blocn(,i,j)/sp)(sum));
    v(i,j)=(((sp-alpha*blocn(,i,j))^2)/(sp)^2)(sum);
  };
};

 v*=snr^2;
 k2_90=0.5*(1.2816+sqrt(2*nlambda-1))^2;

 if (is_void(nodis)){
   ws;
   plk,log(max(v,0.1*k2_90))/log(10);
   x=indgen(nda)(,-:1:ndz);
   y=indgen(ndz)(-:1:nda,);
   plc,v,y,x,levs=k2_90*[0.7,1.,1.3],width=3;
 };
 return k2_90;
};


func invm(sp,mu,&aAtil,rr=,pl=,of=,ix0=,iy=){
  /* DOCUMENT
     returns inverse model for linear problem represented by sp, regularized by kernel rr, with smoothing parameter mu.
     try _invm if it doesnt work
  WARNING: result is transposed
     
  */
  local pp,Atil,aAtil,sp,rr,ix0;

  if(is_void(ix0)) ix0=indgen(dimsof(sp)(2));
  if(is_void(iy)) iy=indgen(dimsof(sp)(3));
  if(is_void(of)) of=[2.,4.,1.33];
  if(dimsof(of)(1)==0) of=[of(1),4.,1.33];
  if(dimsof(of)(2)==2) of=[of(1),of(2),1.33];
  
  
  if(is_void(pl)) pl=[0,0,0];
  
  if(is_void(rr)) {
    pp=dimsof(sp)(3);
    rr=mkPENAL(pp,s="grad");
  };
  pp=dimsof(sp)(3);
  Atil=(LUsolve(sp(+,)*sp(+,)+mu*(rr(+,)*rr(+,))))(,+)*sp(,+);
  Atil=transpose(Atil);
  adyrange=(iy(0)-iy(1))/pp;
  
  if (pl(1)!=0) {
    ws;
    plk,Atil,iy(,-:1:numberof(ix0)),ix0(-:1:pp,);
    plh,iy(int(pp/of(1)))+pl(1)*adyrange*sp(,int(pp/2.)),ix0,color="green",width=3;
    //pltitle,"influence matrix mu="+pr1(mu);
  };

  if(pl(2)!=0){
    plh,iy(int(pp/of(2)))+pl(2)*adyrange*sp(,int(pp/4.)),ix0,color="green",width=3;
    plh,iy(int(pp/of(3)))+pl(2)*adyrange*sp(,int(pp/1.33)),ix0,color="green",width=3;
  };

  aAtil=abs(Atil(,ptp));
  aAtil/=aAtil(avg);
  if(pl(3)!=0){
    plh,iy(1)+pl(3)*adyrange*aAtil,ix0,color="blue",width=3;
  };
  
  return (Atil);
};


func _invm(sp,mu,&aAtil,rr=,pl=,ix0=,N=){
  /* DOCUMENT
     simpler version of invm
     returns influence matrix for linear problem represented by sp, regularized by kernel rr, with smoothing parameter mu.
     N=1 norms lines of Atil (is it useful ? right or wrong?)
     WARNING: SAME WARNINGS AS invm
  */
  local pp,Atil,aAtil,sp,rr,ix0;

  if(is_void(ix0)) ix0=indgen(dimsof(sp)(2));
  if(is_void(iy)) iy=indgen(dimsof(sp)(3));
  if(is_void(of)) of=[2.,4.,1.33];
  if(dimsof(of)(1)==0) of=[of(1),4.,1.33];
  if(dimsof(of)(2)==2) of=[of(1),of(2),1.33];
  
  
  if(is_void(pl)) pl=[0,0,0];
  
  if(is_void(rr)) {
    pp=dimsof(sp)(3);
    rr=mkPENAL(pp,s="grad");
  };
  Atil=(LUsolve(sp(+,)*sp(+,)+mu*(rr(+,)*rr(+,))))(,+)*sp(,+);
  Atil=transpose(Atil);
  //adyrange=(iy(0)-iy(1))/pp;
  sAtil=Atil; // for return;
  if(is_complex(Atil)) Atil=abs(Atil);
  if(!is_void(N)) Atil/=Atil(sum,)(-,);

  
  
  if (pl(1)==1) {
    ws;
    plk,Atil,iy(,-:1:numberof(ix0)),ix0(-:1:pp,);
  };

  if(pl(1)==2){
  plh,(int(pp/of(1)))*sp(,int(pp/2.)),ix0,color="green",width=3;
    //pltitle,"influence matrix mu="+pr1(mu);
  };

  if(pl(2)!=0){
    plh,(int(pp/of(2)))*sp(,int(pp/4.)),ix0,color="green",width=3;
    plh,(int(pp/of(3)))*sp(,int(pp/1.33)),ix0,color="green",width=3;
  };

  aAtil=abs(Atil(,ptp));
  aAtil/=aAtil(avg);
  if(pl(3)!=0){
    plh,iy(1)+aAtil,ix0,color="blue",width=3;
  };
  
  return (sAtil);
};



func rGCVsims(file,&snrs,&sngres,n=){

  /* DOCUMENT
     just to read and plot the GCV simus */

  //extern gres,sgres,rgres;
  if(is_void(n)) n=1000;
  n+=10;  // cause the header is size 10
  f=open(file,"r");
  t=rdline(f,n);
  close,f;
  nsnrs=strmatch(t,"snrs")(sum);
  snrs=array(0.,nsnrs);
  snrs=str2float(split2words(t(:nsnrs),sep="=")(,2));

  nf=str2int(split2words(t(nsnrs+1),sep="  ")(2));
  nmc=str2int(split2words(t(nsnrs+2),sep="=")(2));
  nad=str2int(split2words(t(nsnrs+3),sep="  ")(2));

  //gres=array(0.,min(nmc*nsnrs,n),2);

  gres=str2float(split2words(t(nsnrs+4:min(nsnrs+3+nmc*nsnrs,n)),sep="  ")(,1));
  sngres=str2float(split2words(t(nsnrs+4:min(nsnrs+3+nmc*nsnrs,n)),sep="  ")(,2));
  ng=is_void(dimsof(where(sngres==0)))?[]:(where(sngres==0))(1)-1;

  if (!is_void(ng)) {gres=gres(:ng);sngres=sngres(:ng);};
    
  sgres=gres(sort(sngres));
  sngres=sngres(sort(sngres));
  rgres=!is_void(ng)?((nmc*nsnrs>ng)?sgres:reform(sgres,nmc,nsnrs)):reform(sgres,nmc,nsnrs);

  return rgres;
};

  
func GCV(y,A,CC,amu,&err1,nodis=){
  /* DOCUMENT
     y data
     A model matrix such that y=Ax
     CC penalization (when searching min((y-Ax)t(y-Ax)+mu xtCCtCCx)
     err1: optional, contains the whole GCV curve
     amu is in powers of 10
     a GCV that runs ok
     
  */
  err1=amu*0.;
  //if (is_void(amu)) amu=span(-10.,10.,10);
  //if(is_void(CC)) CC=mkPENAL(dimsof(A)(3),s="tikho");
  mu1=10^(amu);
  
  AtA= transpose(A)(,+)*A(+,);
  CtC=(CC(+,)*CC(+,));
  
  
  for(imu=1; imu<=dimsof(amu)(2); imu++){
    
    inv= AtA + mu1(imu)*CtC;
    inv =LUsolve(inv);
    Atil=inv(,+)*A(,+);
    xtil=Atil(,+)*y(+);

    num=A(,+)*xtil(+)-y;
    num=num(+)*num(+);

    den=0.;
    for(i=1;i<=dimsof(y)(2);i++){
      den+=1-(A(i,)*Atil(,i))(sum);
    };
    den*=den;
    
    err1(imu) = (num/den)*dimsof(y)(2);  // why /dimsof(y)(2) ??
    
  };

  res=mu1(err1(mnx));
  
  if(is_void(nodis)) plg, log(err1), amu;
  return log10(res);
};

func mGCV(A,SNR,&_mu,CC=,amu=,vec=,nodisp=,nres=){

  /*DOCUMENT
    returns the median GCVmu for kernel A, SNR, penal CC defined as in GCV
  */

  local vec,nA,d,_m,amu,nres,CC,d0;
  nA=dimsof(A)(3);
  if(is_void(vec)) vec=makebump(nA,int(nA/3.),int(nA/15.),N=1);
  if(is_void(nres)) nres=10;
  if(is_void(CC)) CC=mkPENAL(dimsof(A)(3),s="tikho");
  if (is_void(amu)) amu=span(-10.,10.,10);
  d0=A(,+)*vec(+);
  d=d0*(1.+(1./SNR)*random_normal(dimsof(A)(2)));
  _mu=[];
  for(i=1;i<=nres;i++){
    grow,_mu,GCV(d,A,CC,amu);};
  plh,histo1d(_mu,amu),amu(:-1);
  return median(_mu);

};



func fitpol(u,x0,n,&pa,w=){
  /* DOCUMENT
     returns the polynomial of degree n fitting best u over x0 with weight w
     doesnot work properly. try
     plh,fitpol(random_normal(2655),span(0.,1.,2655),1,w=w),color="black"
     plh,fitpol(random_normal(2655),span(1.,2.,2655),1,w=w),color="red"
     WARNING:  x0 should not be integer but float

     prefer fitpol2
     
  */
  local x,p,v,s;
  if(is_void(w)) w=array(1.,numberof(x0));
  x=array(0.,n);
  Cm=array(0.,n,n);
  Y=0.*x
    for(i=1;i<=n;i++){
      for(j=1;j<=n;j++){
        Cm(i,j)=((x0^i)*w*x0^j)(sum);};
      Y(i)=(u*w*x0^i)(sum);
    };
  pa=SVsolve(Cm,Y,1.e-100);
  res=u*0.;
  for(i=1;i<=n;i++){ res+=pa(i)*x0^i;};
  return res;
};

func fitpol2(u,x0,n,&pa,w=){
  /* DOCUMENT
     returns the polynomial of degree n-1 fitting best u over x0 with weight w
     doesnot work properly. try
     plh,fitpol(random_normal(2655),span(0.,1.,2655),1,w=w),color="black"
     plh,fitpol(random_normal(2655),span(1.,2.,2655),1,w=w),color="red"
     WARNING:  x0 should not be integer but float and should contain no exact zero.
     
  */
  local x,p,v,s;
  if(is_void(w)) w=array(1.,numberof(x0));
  x=array(0.,n);
  Cm=array(0.,n,n);
  Y=0.*x
    for(i=1;i<=n;i++){
      for(j=1;j<=n;j++){
        Cm(i,j)=((x0^(i-1))*w*x0^(j-1))(sum);};
      Y(i)=(u*w*x0^(i-1))(sum);
    };
  pa=SVsolve(Cm,Y,1.e-100);
  res=u*0.;
  for(i=1;i<=n;i++){ res+=pa(i)*x0^(i-1);};
  return res;
};





func wnormb(b,n,&c){
  /* DOCUMENT
     normalizes a basis (or a spectra with wnorms ) using wavelet-transforms, the continuum being computed as the nth order of the wavelet transform
     The continuum is also returned in c, whcih has a basis-like structure
     check the borders and mask if necessary
  */
  local b1,c,n;
  if(is_void(n)) n=7;
  wt=yeti_wavelet(b.flux,n,which=1,border=5);
  b1=b;
  b1.flux=b.flux/wt(,,,n+1);
  c=b;
  c.flux=wt(,,,n+1);
  wt=[];
  return b1;
};

func wnorms(s,n,&c){
  /* DOCUMENT
     same as wnormb for a single spectrum

  */
local s1,c,n;
  if(is_void(n)) n=7;
  wt=yeti_wavelet(s,n,border=5);
  c=wt(,n+1);
  s1=s/c;
  return s1;
};

func fftnorm(s,n){
  local s1,c,n;
  if(is_void(n)) n=500;
  ss=fft_smooth(s,n);
  s1=s/ss;
  return s1;
};
  

func fftnormb(b,n,&c){
  /* DOCUMENT
     normalizes a basis (or a spectra with fftnorms ) using fft_smooth, 
     The continuum is also returned in c
     check the borders and mask if necessary
  */
  local b1,c,n;
  if(is_void(n)) n=500;
  c=b;
  for(i=1;i<=dimsof(b)(3);i++){
    for(j=1;j<=dimsof(b)(4);j++){
      c(,i,j)=fft_smooth(b(,i,j),n);
    };
  };
  b1=b;
  b1=b/c;
  return b1;
};



func _wnormb(s,n,&c){
   /* DOCUMENT
     same as wnormb for a basis that is not in an array struct

  */
local s1,c,n;
  if(is_void(n)) n=7;
  wt=yeti_wavelet(s,n,which=1,border=5);
  c=wt(,,,n+1);
  s1=s/c;
  return s1;
};


func cchi2(x,y,snr=,w=){
  /* DOCUMENT
     computes chi2 as sum (((xi-yi)^2)/(w^2))/numberof(x). If w is not given it is taken as (x/snr) */
    
    if (is_void(w)) w= x/snr;
  res=((((x-y)/w)^2)/numberof(x))(sum);
  return sqrt(res);
};

func resmx(A,L,mu){

  /* DOCUMENT
     returns the resolution matrix (not sure its the correct name) Achapeau such that
     xmu = Achapeau * x
  */
  local i1,Atil,res;
  i1=LUsolve(A(+,)*A(+,)+mu*(L(+,)*L(+,)));
  Atil=i1(,+)*A(,+);
  res=Atil(,+)*A(+,);
  return res;
};

func transmx(A,B,L,mu,W=){
  /* DOCUMENT
     returns the transfer matrix in model space from SSP base matrixA to SSP base matrix B
     as invmodel(A,L,mu) . B
  */
  return invmodel(A,L,mu,W=W)(,+)*B(+,);
};

     




func invmodel(A,L,mu,W=){
  /* DOCUMENT
     returns the inverse model matric Atilde so that
     xmu=Atilde . y
     W is a weight vector
  */
  local i1,Atil,res;
  if(is_void(W)){return invmodel(A,L,mu,W=array(1.,dimsof(A)(2)));};
      i1=LUsolve((A*W(,-:1:dimsof(A)(3)))(+,)*A(+,)+mu*(L(+,)*L(+,)));
    Atil=i1(,+)*A(,+);
    Atil=Atil*W(-:1:dimsof(A)(3),);
    return Atil;
};
  

func invmx(A,L,mu){
  /* DOCUMENT
     returns inverse model of A regularized by L, so that
     xmu = Atild*y
  */
  return transpose(invm(A,mu,rr=L));
};


func infmx(A,L,mu){
  /* DOCUMENT
     return the influence matrix that is Asharp so that
     ymu = Asharp 
  */
};





func TSVD(A,y,k,&At){
  /* DOCUMENT
    TSVD solution of the problem y=Ax with SVD truncated at rank k
   */

  s=SVdec(A,u,v);
  tsm1=0.*s;
  tsm1(:k)=1./s(:k);
  tx=u(+,)*y(+);
  tx=tx*tsm1;
  tx=v(+,)*tx(+,);

  ts=s*0.;ts(:k)=s(:k);
  At=diag(ts)(,+)*v(+,);
  At=u(,+)*At(+,);

  return tx;
};

func MTSVD(A,y,k,L,mu,&At){
  /* DOCUMENT
     computes the MTSVD (Modified TSVD) solution to y=Ax i.e. the min of
     (trA-y)t(trA-y) + mu xtLtLx     where trA is A truncated at rank k (done by function TSVD)
  */

  rr=L(+,)*L(+,);
  qtsvd=TSVD(A,y,k,At);
  qmtsvd=SVsolve(At(+,)*At(+,)+mu*rr,At(+,)*y(+),1.e-40);
  return qmtsvd;
};



func Lcurve(A,y,L,&gres,amu=,noplot=){
  /* DOCUMENT
    plots the L-curve of problem y=Ax regularized by penalization
    P(x)=XtLtLx
    optionnally returns an array [amu,res,reg]
    returns the mu at the corner. (maximum of second derivative)
  */
  
  if(is_void(amu)) amu=10^(float(indgen(21)-10));
  rr=L(+,)*L(+,);
  res=array(0.,numberof(amu));
  reg=res;
  for(i=1;i<=numberof(amu);i++){
    mu=amu(i);
    x=SVsolve(A(+,)*A(+,)+mu*rr,A(+,)*y(+),1.e-40);
    res(i)=((A(,+)*x(+)-y)^2)(sum);
    reg(i)=(rr(,+)*x(+))(+)*x(+);
  };

  gres=array(0.,3,numberof(amu));
  gres(1,)=amu;gres(2,)=res;gres(3,)=reg;

  if(is_void(noplot)){
    ws;
    PL,reg,res,msize=.5,color="black";
    plg,reg,res,color="black";
    logxy,1,1;
    lab=swrite(amu,format="%02E");
    for(i=1;i<=numberof(amu);i++){plt,lab(i),res(i),reg(i),tosys=1;};
    xyleg,"reg norm","residuals norm";
  };

  // make max curvature detection

  a=spline(log10(abs(res)),log10(abs(reg)));
  da=spline(a,log10(abs(reg)));
  mmu=amu(da(mxx));
  

  
  return mmu;
};

  
func isolve(A,y,L,mu,stop=,guess=,frtol=,verb=,maxeval=){

  /* DOCUMENT iterative solver of A.x=y regularized through mu L
     positivity through quad reparameterization
     uses quadpenalty in defined Qfunction.i
     stop is a stop criterion
     it uses the function
  */

  local nab;
  spd=A;
  nab=dimsof(spd)(3);
  d=y;
  rr=L(+,)*L(+,);
  if(is_void(guess))  guess=array(0.1,nab);
  pq=(optim_driver(quadpenalty,guess,verb=verb,frtol=frtol,maxeval=maxeval))^2;
  return pq;
};

func isolveP(A,y,mup,L,stop=,guess=,frtol=,verb=,maxeval=){

  /* DOCUMENT iterative solver of A.x=y regularized through a function
     Pisolve(x,&g)
     positivity through quad reparameterization
     uses quadpenalty in defined Qfunction.i
     stop is a stop criterion
     it uses the function
  */

  local nab;
  spd=A;
  nab=dimsof(spd)(3);
  d=y;
    rr=L(+,)*L(+,);
    nlos=nz;
    nab=nl;
    La=genREGUL(nl,s="D2");rra=La(+,)*La(+,);mua=1.e5;
    Lb=genREGUL(nz,s="D2");rrb=Lb(+,)*Lb(+,);mub=1.e5;
  //rr=L;
  mu=0.;
  if(is_void(guess))  guess=array(1.e0,dimsof(spd)(3));
  pq=(optim_driver(quadpenaltyPlog2,guess,verb=verb,frtol=frtol,maxeval=maxeval))^2;
  return pq;
};


 
func lcor(d,di,ncl,&_m,&im,&nind,ninhib=,noplot=,l=,nw=,nosub=){
  /* DOCUMENT

     prefer lcor2
     
     small scale wavelength correction attempt between d and di
     the spectrum is divided into segments of length ncl which are individually recalibrated in wavelength to match di
     
     nosub=1: subpixel scale shifts forbidden. more robust
     nosub=void => allows subpixel scale shifts. Careful, determination of the max of the correlation peak is buggy: offset by a few. use very small nw to get small shift errors. But nw must be larger than the real ofset

     Advice: first run lcor with nosub=1 to check for large shifts, and then when the shift has been brought back to subpixel scale, run with nosub=[]
     
plots in 0 and 1
     
  */

  if(is_void(ninhib)) ninhib=0;
  if(is_void(noplot)) p=1;
  if(is_void(nw)) nw=3;
  if(!is_void(nosub)) nosub=1;
  if(is_void(somepower)) somepower=4;
  w=array(0.,ncl);  // weights for the gaussian fit
  w(ncl/2-nw:ncl/2+nw)=1.;
  nl=numberof(d);
  nic=int(nl/ncl);
  
  
  _m=[];
  if(p) ws,1;
  for(i=1;i<=nic;i++){
    //coi=xycorrel(d((i-1)*ncl+1:i*ncl),di((i-1)*ncl+1:i*ncl));
    coi=mycorrel(d((i-1)*ncl+1:i*ncl),di((i-1)*ncl+1:i*ncl));
    if(p) plh,-10*i+roll(coi);
    if (is_void(nosub)) {
      pars=gaussfit(roll(coi)^somepower,indgen(ncl),w);
      y=pars(1)*exp(-[(indgen(numberof(a))-pars(2))/pars(3)]^2);
      y2=(y(,1))^(1./somepower);

      grow,_m,pars(2);
      if(p) plh,y2-10*i,color="red",width=3;
    };
        
    if (nosub==1) grow,_m,roll(coi)(mxx); // only pixel scale shifting
  };
  ws,0;
  _m-=ncl/2+1;
  write,_m;

  // follow various implementations of the correction, not equally efficient

  if(1){ // _m=> fonction constante par morceaux sur tout le spectre.
    // not bad but oscillates
    inm=_m(-:1:ncl,);
    inm=inm(*);
    if ((nl-nic*ncl)!=0) grow,inm,array(0.,nl-nic*ncl);  
  };

  
  if(0){
    inm=interp(_m,indgen(nic),span(.5,nic+.5,nic*ncl)); // brutal interpolation
    grow,inm,array(0.,nl-nic*ncl);
  };

  if(0){  //one way  works bof
    trm=0.;
    grow,trm,_m,0.;write,trm;
    _m=trm;
    ind1=0.5;grow,ind1,(indgen(nic)+.5),nic+1;
    inm=spline(_m,ind1,span(.5,nic+1,nl)); 
  };
  
  if(0){  // another way trhat works quite ok
    ind1=indgen(nic);
    inm=spline(_m,ind1,span(.5,nic+.5,nic*ncl)); 
    if ((nl-nic*ncl)!=0) grow,inm,array(0.,nl-nic*ncl);  
  };
    
    
  if(ninhib!=0) {// lock boundaries
    inm(:ninhib)=0;
    inm(nic*ncl+1:)=0;
  };
  //inm(1880:1920)=0.;   
  if(p) plh,inm+1,color="blue",type=4;
  pause,0;
  nind=indgen(nl)+inm;
  return spline(d,indgen(nl),nind);
  
};





func lcor2(d,di,ncl,&_m,&inm,&nind,ninhib=,disp=,l=,nw=,nosub=,wS=,nwS=,d1=,tr=,wN=){
  /* DOCUMENT

     returns d piecewise shifted to match di
     
     small scale wavelength correction attempt between d and di
     the spectrum is divided into segments of length ncl which are individually recalibrated in wavelength to match di
     
     nosub=1: subpixel scale shifts forbidden. more robust
     nosub=void => allows subpixel scale shifts. Careful, determination of the max of the correlation peak is buggy: offset by a few. use very small nw to get small shift errors. But nw must be larger than the real ofset

     Advice: first run lcor with nosub=1 to check for large shifts, and then when the shift has been brought back to subpixel scale, run with nosub=[]
     
plots in 0 and 1

nwS is the number of fourier coefficients used to estimate the signal model

disp=0 -> no display
disp=1 -> fft, noise and model
disp=2 -> several quotients and fits

     
  */

  if(is_void(ninhib)) ninhib=0;
  if(is_void(noplot)) p=1;
  if(!is_void(nosub)) nosub=1;
  nic=int(nl/ncl);
  nl=numberof(d);
  
  _m=[];
  if(p) ws,1;
  for(i=1;i<=nic;i++){
    wS=array(0.,numberof(d((i-1)*ncl+1:i*ncl)));
    wS(1:nwS)=1.;//INFO,wS;
    window,2;
    //if(i==2) error;
    grow,_m,findshift(d((i-1)*ncl+1:i*ncl),di((i-1)*ncl+1:i*ncl),wf=1,disp=disp,wS=wS,wN=wN,pix=nosub,nw=nw,d1=d1,tr=tr,pars,wief,sig,quo,y2);
    //INFO,quo;
    if(disp==2) {window,1;plh,-10*i+quo;plh,-10*i+y2,color="red";};
  };
  ws,0;
  //_m-=int(ncl/2.)+1;
  write,_m;
  
  // follow various implementations of the correction, not equally efficient
  // more are given in lcor
  
  if(1){ // _m=> fonction constante par morceaux sur tout le spectre.
    // not bad but oscillates
    inm=_m(-:1:ncl,);
    inm=inm(*);
    if ((nl-nic*ncl)!=0) grow,inm,array(0.,nl-nic*ncl);  
  };
  
  if(ninhib!=0) {// lock boundaries
    inm(:ninhib)=0;
    inm(nic*ncl+1:)=0;
  };
  
  if(disp==2) plh,inm+1,color="blue",type=4;
  pause,0;
  nind=indgen(nl)+inm;
  return spline(d,indgen(nl),nind);
  
};



func _wlcor(u,ncl,&_m,&inm,&nind,ninhib=,disp=,l=,nw=,nosub=,wS=,nwS=,wN=,d1=,tr=){
  /* DOCUMENT
     small scale wavelength correction attempt between the data and its model pointed by u (galStruct)
     uses lcor2 and inherits all options from it.
     returns a new structure pointing to a new filename, resfile etc...
     writes the file named filename containing the piece-shifted spectrum.
     WARNING! erases u
  */

  
  local ba;

  su=u; // save extern variable;
  upload,u.resfile(1),s=1;
  wave=x0;
  ba=lcor2(d,pmodel1,ncl,_m,inm,nind,ninhib=ninhib,disp=disp,l=l,nw=nw,nosub=nosub,wS=wS,nwS=nwS,d1=d1,tr=tr,wN=wN);

  flux=ba;
  //wave=x0;
  
  //wave+=inm;
    
  u.resfile=split2words(u.resfile,sep=".")(1)+"_wlc."+split2words(u.resfile,sep=".")(2);
  u.filename=split2words(u.filename,sep=".")(1)+"_wlc."+split2words(u.filename,sep=".")(2);
  u.name=u.name+"_wlc";
  gal=u;
  f=createb(u.filename(1));
  save,f,gal,flux,wave,sigm,mask,R;

  u(1:0)=su(1:0); // as Dave said, it's a feature, not a bug
  return gal;
};
  
func wlcor(u,ncl,&_m,&inm,&nind,ninhib=,disp=,l=,nw=,nosub=,wS=,nwS=,wN=,d1=,tr=){
  /* DOCUMENT
     small scale wavelength correction attempt between the data and its model pointed by u (galStruct)
     uses lcor2 and inherits all options from it.
     returns a new structure pointing to a new filename, resfile etc...
     writes the file named filename containing the piece-shifted spectrum.
     WARNING! erases u
  */

  
  local ba;

  su=u; // save extern variable;
  upload,u.resfile(1),s=1;
  wave=x0;
  ba=lcor(d,pmodel1,ncl)
    
  flux=ba;
  //wave=x0;
  
  //wave+=inm;
    
  u.resfile=split2words(u.resfile,sep=".")(1)+"_wlc."+split2words(u.resfile,sep=".")(2);
  u.filename=split2words(u.filename,sep=".")(1)+"_wlc."+split2words(u.filename,sep=".")(2);
  u.name=u.name+"_wlc";
  gal=u;
  f=createb(u.filename(1));
  save,f,gal,flux,wave,sigm,mask,R;

  u(1:0)=su(1:0); // as Dave said, it's a feature, not a bug
  return gal;
};
  




func fquotient(a,b){
  /* DOCUMENT
     returns crosscorrelation between a and b by fourier quotient
     a and b need to be similar and have the same mean. (i.e. no template mismatch)
  */

  return roll(fft(fft(a)/fft(b),[-1]));
};
  
func findshift(a,b,&pars,&wief,&sig,&quo,&y2,nw=,tr=,disp=,wf=,nodispwf=,wS=,wN=,d1=,d2=,pix=,somepower=,nwS=){
  /* DOCUMENT
     new attempt to determine shift between a and b using fourier quotient
     nw is the number of pixels around the middle for fitting a gaussian on the crosscorrelation function.
     tr is the truncature rank for the fft of noisy a if one does not wish to wiener filter
     it's a bit violent
     Options for wiener filtering
     wf is off by default
     wf=1 => wiener filtering on
     automatic by default
     other wise user must specify
     wS=
     wN=
     d1=
     d2=
     if pix=1 the shift is given only in an integer number of pixels. This is more robust
     since the fourier quotient quo does not have a gaussian shape, we actually fit quo^somepower (it looks much more like a gaussian and thus the measurement of the shift is much better)
     disp=0 -> no display
     disp=1 -> wienerf display   CHECK THAT NOISE AND SIGNAL MODELS ARE OK
     disp=2 -> fourier quo display CHECK THAT THE BUMP IS WELL FITTED

     optional output sig is the width of the main bump of the quotient, i.e. an estimate of the width of the PSF that gave a from b.
     
  */

  extern w;
  local na;
  
  if(is_void(disp)) disp=1;
  na=numberof(a);
  if(is_void(nw)) nw=int(na/2.)-1;
  if(is_void(somepower)) somepower=8;
  
  //if (wf==1) {
  //  if(is_void(wS)) {wS=array(0.,na);wS(1:na/30)=1.;}
  //  if(is_void(wN)) {wN=array(0.,na);wN(na/5:-na/5)=1.;};
  //  if(is_void(d1)) d1=4;
  //  if(is_void(d2)) d2=2;

  //  wf=wienerf(a,wS=wS,wN=wN,nodisp=(disp==1?[]:1),d1=d1,d2=d2);
  //  wief=wf;
  //};

  if (wf==1) {
    wf=wienerf(a,wS=wS,wN=wN,nwS=nwS,nodisp=(disp==1?[]:1),d1=d1,d2=d2);
    wief=wf;
  };

  
  w=array(0.,numberof(a));
  w(int(numberof(a)/2.)-nw:int(numberof(a)/2.)+nw)=1.;
  fa=fft(a);
  if(!is_void(tr)) fa(tr:)=0.;
  if(!is_void(wf)) fa=fa*wf;  // wiener filter applied here
  fb=fft(b);  
  quo=abs(roll(fft(fa/fb,[-1])));
  
  pars=gaussfit(quo^somepower,indgen(numberof(a)),w);
  y=pars(1)*exp(-[(indgen(numberof(a))-pars(2))/pars(3)]^2);
  pars2=gaussfit(y(,1)^(1./somepower),indgen(numberof(a)),w);
  sig=pars2(3);
  y2=(y(,1))^(1./somepower);
  
  if(disp==2) {
    plh,quo,width=3;
    plh,y2,color="red",width=3;
    pltitle,"fourier quotient and log-normal fit";
  };
  //pars=gaussfit(y^(1./somepower),indgen(numberof(a)),w);
  return (pix==1?quo(mxx):pars(2))-int(numberof(a)/2.)-1.;
};




func wienerf(a,&s1,&n1,&ni,wS=,wN=,nwS=,nodisp=,d1=,d2=,nguess=){
  /* DOCUMENT
     returns the wiener filter 
     derived from the fft of a.
     ws and wn are same size as a.
     ws is 1 where the signal is thought to dominate,0 elsewhere
     wn is 1 where the noise is thought to dominate.
     the signal power spectrum is fitted by a log-polynomial of degree d1-1 over wS.
     the noise power spectrum is fitted by a log-polynomial of degree d2-1 over wN.
     check the display and tweak until it looks good
     default: is automatic mode
     default is d1=4 and d2=2
     is wS is not given an automatic search for the transition between signa dominated and noise dominated regions is done
     a symbol where the transition is detected.
     nwS is a guess for where to the signal dominates
  */

  local fa,na;
  fa=abs(fft(a));
  na=numberof(a);
  
  
  if(is_void(d1)) d1=4;
  if(is_void(d2)) d2=2;
  if(is_void(nguess)) nguess=20;
  if(is_void(wN)) {wN=array(0.,na);wN(na/5:-na/5)=1.;};


  // automatic detection of signal-noise transition (attempt)
  
  if(is_void(wS)&is_void(nwS)) {
    wS=array(1.,na);
    wi=wienerf(a,s1,n1,wS=wS,wN=wN,nodisp=1,d1=100,d2=2);
    wS*=0.;
    wS(1:nguess)=1.; // 
    s2=10^min(fitpol2(log10(s1),span(0.1,1.,numberof(a)),2,pa,w=wS),100.);
    // the beginning of the spectrum (up to nguess) is fitted by a fonction affine 
    
    ni=where(s2<=n1)(1); // detection of signal-noise transition
    //ni*=.6; // have to guess that
    ni=int(ni);
    wS=array(0.,na);
    wS(:int(ni))=1.;
  };

  if(!is_void(nwS)) {wS=array(0.,na);wS(:nwS)=1.;};
    
  s1=10^min(fitpol2(log10(fa),span(0.1,1.,numberof(a)),d1,pa,w=wS),100.);
  n1=10^fitpol2(log10(fa),span(0.1,1.,numberof(a)),d2,pa,w=wN);

  if(is_void(nodisp)) {
    ws;
    plh,fa;
    plh,s1,color="blue",width=3;
    plh,n1,color="red",width=3;
    logxy,0,1;
    pltitle,"fft, signal and noise models";
    if(!is_void(ni)) PL,s1(ni),ni,marker=1,msize=1,color="red",incolor="red";
  };
  return (s1^2/(s1^2+n1^2));
};

func fndate{
  /* DOCUMENT
     returns date usable for a file name
     WARNING: english convention mm-dd-yy and hour minute
  */
  return strreplace(exec("date +%D%H%M"),"/","-")(1);
};
  
func intersectS(wave0,wave1){
  /* DOCUMENT
     returns the intersection of the supports
     finest max sampling wins
     wave0 and wave 1 MUST be increasing arrays
  */

  wave0=wave0; // to avoid external crap (never know)
  wave1=wave1; // itou 

  if (wave0(0)<wave1(1)) return (wave0(0)+wave1(1))/2.;
  if (wave1(0)<wave0(1)) return (wave1(0)+wave0(1))/2.;
  
  if(wave0(dif)(max)>wave1(dif)(max)) return intersectS(wave1,wave0);
  
  wmin=(wave0(1)<wave1(1))?(where(wave0<=wave1(1)))(0)+1:1; 

  wmax=(wave0(0)>wave1(0))?(where(wave0>=wave1(0)))(1)-1:0;

  wave0=wave0(wmin:(wmax==0?numberof(wave0):wmax));
  
  return wave0;

};






