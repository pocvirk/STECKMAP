// A FEW EXAMPLES HOW TO CONVERT the .fits to .pdb file AND fill the galStruct structure
// 20/07/2006 SOME EXTENSIONS OF convert_all


#include "fits.i"
#include "spec2ls.i"


// Structure definition of galStruct

struct galStruct {
  string   name;
  string   filename;
  string   result_dir;
  string   resfile;
  string   wavesampling;   //  possible values are "LIN" or "LOG"
  double   redshift;
  double   SNR;
};            

struct sdssgalStruct {
  // similar to galStruct but we store more parameters, such as sigmav if available
  string   name;
  string   filename;
  string   result_dir;
  string   resfile;
  string   wavesampling;   //  possible values are "LIN" or "LOG"
  double   redshift;
  double   sigmav;
  double   SNR;
};    


// 3 Examples of conversion to pdb routines

func convertSDSS_old(filelist,cut=,noplot=,s=)
/* DOCUMENT 
   EXAMPLE  ll=exec("find /home5/ocvirk/spectro -name '*.fit'")
   convertSDSS(filelist);
   SEE ALSO:
   
   // first row is the spectrum,
   // second row is the wavelengths   // third the error, and
   // fourth a bitmask 
   */
  
{
  if(is_void(s)) s=0;
  if(is_void(cut)) cut=5;
  gals=[];
  for(i=1;i<=min(cut,numberof(filelist));i++)
    {
      u=fitsRead(filelist(i),h);
      wave=10^(h.crval(1)+(indgen(h.axis(1))-h.crpix(1))*h.CD(1,1));

      if (is_void(noplot)) {
      fma;
      plh,u(,1),wave,color=__rgb(,(i%20)+10);
      plg,u(,1)-u(,2),wave,color=-10,type=2,width=4;
      plh,u(,4)/u(avg,4),wave,color="blue";
      limits;//  mouse();
      };
      z= getRedshift(filelist(i));
      sigm=u(,3); flux=u(,1); mask=(u(,4)>1e7); ww=where(mask); sigm(ww) *=100;
      gal=galStruct();
      gal.redshift=z;
      gal.SNR=getSNR(filelist(i));
      gal.name= split2words(filelist(i),sep="./")(-1);
      gal.filename=strtok(filelist(i),".")(1)+".pdb";
      gal.filename;
      gal.resfile=strtok(filelist(i),".")(1)+".res1";
      gal.result_dir="/"+strjoin(split2words(gal.filename,sep="/")(1:-1),"/")+"/";
      gal.wavesampling="LOG";
      save,createb(gal.filename),gal,flux,wave,sigm,mask;
      gal.redshift;
      grow,gals,gal;
    }

  return gals;
}

func convertSDSS(filelist,cut=,noplot=,s=,putatrest=)
/* DOCUMENT 
   EXAMPLE  ll=exec("find /home5/ocvirk/spectro -name '*.fit'")
   convertSDSS(filelist);
   SEE ALSO:
   
   // first row is the spectrum,
   // second row is the wavelengths   // third the error, and
   // fourth a bitmask 
   */
  
{
  if(is_void(putatrest)) putatrest=0;
  if(is_void(s)) s=0;
  if(is_void(cut)) cut=5;
  gals=[];
  for(i=1;i<=min(cut,numberof(filelist));i++)
    {
      u=fits_read(filelist(i),h);
      //      wave=10^(h.crval(1)+(indgen(h.axis(1))-h.crpix(1))*h.CD(1,1));

      wave=10^(fits_get(h,"CRVAL1")+(indgen(fits_get(h,"NAXIS1"))-fits_get(h,"CRPIX1"))*fits_get(h,"CD1_1"));

      z=fits_get(h,"Z   ");
      if(putatrest==1) {
        write,"putting galaxy at rest using header redshift z=",z;
                        wave/=(1.+0.*z+0.00018); // is this due to a vacuum/air wavelength definition ?  0.00018 estimated from the spectrum // that seems to work better but does it depend on the galaxy ?
        //wave/=(1.+0.*z+0.0003); // is this due to a vacuum/air wavelength definition ?
        // refraction index of air = 1.0003. thats why
        // and we dont really put the gals at rest cause seems they are already at rest...
      };
      
      if (is_void(noplot)) {
      fma;
      plh,u(,1),wave,color=__rgb(,(i%20)+10);
      plg,u(,1)-u(,2),wave,color=-10,type=2,width=4;
      plh,u(,4)/u(avg,4),wave,color="blue";
      limits;//  mouse();
      };
      sigm=u(,3); flux=u(,1); mask=(u(,4)>1e7); ww=where(mask); sigm(ww) *=100;
      gal=sdssgalStruct();
      gal.redshift=z;
      gal.sigmav=fits_get(h,"VEL_DIS");
      //gal.SNR=100.; //getSNR(filelist(i));
      gal.SNR=avg([fits_get(h,"SN_G"),fits_get(h,"SN_R"),fits_get(h,"SN_I")]);
      gal.name= split2words(filelist(i),sep="./")(-1);
      //      gal.filename=strtok(filelist(i),".")(1)+".pdb";
      gal.filename=strreplace(filelist(i),".fit",".pdb");
      //      gal.filename;
      //      error;
      //      gal.resfile=strtok(filelist(i),".")(1)+".res1";
      gal.resfile=strreplace(filelist(i),".fit",".res");
      gal.result_dir="/"+strjoin(split2words(gal.filename,sep="/")(1:-1),"/")+"/";
      gal.wavesampling="LOG";
      save,createb(gal.filename),gal,flux,wave,sigm,mask;
      gal.SNR;
      grow,gals,gal;
    }

  return gals;
}



func convertVAKU(filelist,cut=,noplot=) 
/* DOCUMENT 
SEE ALSO: convertSDSS
*/

{
  if(is_void(cut)) cut=10;
  gals=[];
  for(i=1;i<=min(cut,numberof(filelist));i++)
    {
      u=fitsRead(filelist(i),h);
      upload,examplesdir+"/VAKUzlist.yor",s=1;
      //wave=10^(h.crval(1)+(indgen(h.axis(1))-h.crpix(1))*h.CD(1,1));
      wave=h.crval(2)+h.cdelt(2)*(indgen(h.axis(2)) -1.);  // CHANGED RECENTLY ON LEDA ??
      flux=u(75:300,)(sum,); // sum lines of CCD array
      
      if (is_void(noplot)) {
      fma;
      plh,flux,wave,color=__rgb(,(i%20)+10);
      limits;//  mouse();
      };
      z=zl(where((h.object)==namel)(1))
        //flux=u(,75:300)(,sum); // sum lines of CCD array
      mask=flux*0.+0.;
      sigm=flux*(1./SNRl(where((h.object)==namel)(1)));
      gal=galStruct();
      gal.redshift=z;
      gal.SNR= SNRl(where((h.object)==namel)(1)); // getSNR(filelist(i));
      gal.name= h.object;
      gal.filename=strtok(filelist(i),".")(1)+".pdb";
      gal.resfile=strtok(filelist(i),".")(1)+".res";
      if (strmatch(filelist(i),"/")==1) gal.result_dir=strjoin(split2words(gal.filename,sep="/")(1:-1),"/")+"/";
      if (strmatch(filelist(i),"/")==0) {
        di=array(string,[1,2]);
        di(1)=(exec("pwd"))(1);
        di(2)="/";
        gal.result_dir=strjoin(di);
      };
      gal.wavesampling="LIN";
      save,createb(gal.filename),gal,flux,wave,sigm,mask;
      gal.redshift;
      grow,gals,gal;
    }

  return gals;
}



func convertJ(filelist,cut=,noplot=,sm=,R=,nosav=) 
/* DOCUMENT
   EXAMPLE  ll=exec("find /home5/evan/spectro -name '*.fit'")
   R= -> the resolution you want the data in
   uses NSCtab1.yor where the published data is stored
   necessary since JAKOB's data is in a strange format
   WARNING! does weird stuff with the variance-covariance matrix

   this one is more complicated because of the handling of the velocity dispersions and possibility of smoothing before writing the .pdb.
   
   SEE ALSO: convertSDSS, convertVAKU
*/

{
  if(is_void(cut)) cut=15;
  if(is_void(nosav)) nosav=0;
  gals=[];
  for(i=1;i<=min(cut,numberof(filelist));i++)
    {
      u=fitsRead(filelist(i),h);
      flux=u;
      wave=10^(h.crval(1)+h.cdelt(1)*(indgen(h.axis(1)) -1.));
      upload,"/Users/Pedro/Yorick/Pierre/POP/EXAMPLES/NSCtab1.yor",s=1;
      //upload,examplesdir+"NSCtab1.yor",s=1;
      //idn=str2double(split2words(filelist(i),sep="NGC b.")(-1));
      idn=str2double(split2words((split2words(filelist(i),sep="NGC"))(0),sep="b")(1));
      sv=sigmav(where(idn==idcard));
      R1=(c/sv)(1);
      SNR=SNRp(where(idn==idcard))(1);
      
      
      if(!is_void(R)){
        fwhm1=(wave(avg)/R1)/wave(dif)(avg);
        fwhm2=(wave(avg)/R)/wave(dif)(avg);

        if(R<R1){
          fwhm=sqrt((fwhm2^2-fwhm1^2));
          flux=fft_smooth(flux,fwhm);
          //SNR=SNR*sqrt(fwhm);  // NO!!! only if resampled accroding to the new R
        };
        if(R>R1) R=R1;
      };
        
      
      if (is_void(noplot)) {
      fma;
      plh,flux,wave,color=__rgb(,(i%20)+10);
      limits;//  mouse();
      };
      z=0.
        //sigm=u(,3);
      
      mask=flux*0.+0.;
      sigm=mask*0.;
        //  flux*(1./SNRl(where((h.object)==namel)(1)));
      gal=galStruct();
      gal.redshift=z;
      gal.SNR=SNR(1);
        //SNRl(where((h.object)==namel)(1)); // getSNR(filelist(i));
      gal.name= h.object;
      //gal.filename=strtok(filelist(i),".")(1)+".pdb";  29/08/06
      gal.filename=strreplace(filelist(i),".fits",".pdb");
      gal.resfile=strreplace(filelist(i),".fits",".res");
      if (strmatch(filelist(i),"/")==1) gal.result_dir=strjoin(split2words(gal.filename,sep="/")(1:-1),"/")+"/";
      if (strmatch(filelist(i),"/")==0) {
        di=array(string,[1,2]);
        di(1)=(exec("pwd"))(1);
        di(2)="/";
        gal.result_dir=strjoin(di);
      };
      gal.wavesampling="LOG";
      if(nosav==0){
        f=createb(gal.filename);
        save,f,gal,flux,wave,sigm,mask;
        if(!is_void(R1)) save,f,R1;
        if(!is_void(R)) save,f,R;
      };
      
      gal.redshift;
      grow,gals,gal;
    };

  return gals;
};

func convertR(filelist,cut=,noplot=){ 
  /* DOCUMENT
     For converting Reynier's data NGC4030 fits to pdb
     good example for conversion from fits to pdb. behaves well.
     SEE ALSO: convertSDSS, convertVAKU, convertJ
  */

  if(is_void(cut)) cut=10;
  gals=[];
  for(i=1;i<=min(cut,numberof(filelist));i++)
    {

      a=fits_read(filelist(i),h,hdu=1);
      //wave=fits_get(h,"CRVAL1")+indgen(dimsof(a)(2))*fits_get(h,"CDELT1");
      wave=fits_get(h,"CRVAL1")+indgen(fits_get(h,"NAXIS1"))*fits_get(h,"CDELT1");
      u=a;
      flux=u(,sum); // sum lines of CCD array
      //flux=u(,640:670)(,sum);
      
      if (is_void(noplot)) {
      fma;
      plh,flux,wave,color=__rgb(,(i%20)+10);
      limits;//  mouse();
      };
      //z=zl(where((h.object)==namel)(1))
        //flux=u(,75:300)(,sum); // sum lines of CCD array
      z=0.00500; //measured by me from spectrum 
      //z=0.0048
      SNR=1.e2;
      mask=flux*0.+0.;
      sigm=flux*(1./SNR);
      gal=galStruct();
      gal.redshift=z;
      gal.SNR= SNR; // getSNR(filelist(i));
      gal.name= fits_get(h,"OBJECT");
      gal.filename=strreplace(filelist(i),".fits",".pdb");
      gal.resfile=strreplace(filelist(i),".fits",".res");
      if (strmatch(filelist(i),"/")==1) gal.result_dir=strjoin(split2words(gal.filename,sep="/")(1:-1),"/")+"/";
      if (strmatch(filelist(i),"/")==0) {
        di=array(string,[1,2]);
        di(1)=(exec("pwd"))(1);
        di(2)="/";
        gal.result_dir=strjoin(di);
      };
      gal.wavesampling="LIN";
        save,createb(gal.filename),gal,flux,wave,sigm,mask;
      gal.redshift;
      grow,gals,gal;
    }

  return gals;
};





func convertJ2(filelist,cut=,noplot=,sm=,R=) 
/* DOCUMENT
   EXAMPLE  ll=exec("find /home5/evan/spectro -name '*.fit'")
   R=-> resamples the data on n bins so that dlambda=0.5*lambda/R
   uses rebin
   cleaner for the variance-covariance matrix handling
   SEE ALSO:
*/


{
  if(is_void(cut)) cut=15;
  gals=[];
  for(i=1;i<=min(cut,numberof(filelist));i++)
    {
      u=fitsRead(filelist(i),h);
      flux=u;
      wave=10^(h.crval(1)+h.cdelt(1)*(indgen(h.axis(1)) -1.));
      upload,"/home6/ocvirk/perso/modeles/galaxie/POP/SDSS/JAKOB/NSCtab1.yor",s=1;
      idn=str2double(split2words(filelist(i),sep="NGC b.")(-1));
      sv=sigmav(where(idn==idcard));
      R1=(c/sv)(1);
      SNR=SNRp(where(idn==idcard))(1);
      
      
      if(!is_void(R)){
        dl=0.5*wave(avg)/R;
        dl0=wave(dif)(avg);
        n=int(dl/dl0);
        flux=rebin(flux,wave,n,x1);
        wave=x1;
        SNR=SNR*sqrt(n);
      };
        
      
      if (is_void(noplot)) {
      fma;
      plh,flux,wave,color=__rgb(,(i%20)+10);
      //plg,u(,1)-u(,2),wave,color=-10,type=2,width=4;
      //plh,u(,4)/u(avg,4),wave,color="blue";
      limits;//  mouse();
      };
      z=0.
        //sigm=u(,3);
      
      mask=flux*0.+0.;
      sigm=mask*0.;
        //  flux*(1./SNRl(where((h.object)==namel)(1)));
      gal=galStruct();
      gal.redshift=z;
      gal.SNR=SNR(1);
        //SNRl(where((h.object)==namel)(1)); // getSNR(filelist(i));
      gal.name= h.object;
      gal.filename=strtok(filelist(i),".")(1)+".pdb";
      gal.resfile=strtok(filelist(i),".")(1)+".res";
      gal.result_dir="/"+strjoin(split2words(gal.filename,sep="/")(1:-1),"/")+"/";
      f=createb(gal.filename);
      save,f,gal,flux,wave,sigm,mask;
      if(!is_void(R1)) save,f,R1;
      if(!is_void(R)) save,f,R;
         
      gal.redshift;
      grow,gals,gal;
      write,"R1=",R1;
      write,"dl0=",dl0,"\n";
      write,"dl=",dl;
      
    }

  return gals;
}

func convertJ3(filelist,cut=,noplot=,sm=,R=) 
/* DOCUMENT
   EXAMPLE  ll=exec("find /home5/evan/spectro -name '*.fit'")
   R= -> the resolution you want the data in
   WARNING! does weird stuff with the variance-covariance matrix
   smoothes and then resamples
   
   SEE ALSO:
*/


{
  if(is_void(cut)) cut=15;
  gals=[];
  for(i=1;i<=min(cut,numberof(filelist));i++)
    {
      u=fitsRead(filelist(i),h);
      flux=u;
      wave=10^(h.crval(1)+h.cdelt(1)*(indgen(h.axis(1)) -1.));
      upload,"/home6/ocvirk/perso/modeles/galaxie/POP/SDSS/JAKOB/NSCtab1.yor",s=1;
      idn=str2double(split2words(filelist(i),sep="NGC b.")(-1));
      sv=sigmav(where(idn==idcard));
      R1=(c/sv)(1);
      SNR=SNRp(where(idn==idcard))(1);
      
      
      if(!is_void(R)){
        fwhm1=(wave(avg)/R1)/wave(dif)(avg);
        fwhm2=(wave(avg)/R)/wave(dif)(avg);

        if(R<R1){
          fwhm=sqrt((fwhm2^2-fwhm1^2));
          flux=fft_smooth(flux,fwhm);
          //SNR=SNR*sqrt(fwhm);  // NO!!! only if resampled according to the new R
        };
        if(R>R1) R=R1;
      };

       if(!is_void(R)){
        dl=0.5*wave(avg)/R;
        dl0=wave(dif)(avg);
        n=int(dl/dl0);
        flux=rebin(flux,wave,n,x1);
        wave=x1;
        SNR=SNR*sqrt(n);
      };
      
      if (is_void(noplot)) {
      fma;
      plh,flux,wave,color=__rgb(,(i%20)+10);
      //plg,u(,1)-u(,2),wave,color=-10,type=2,width=4;
      //plh,u(,4)/u(avg,4),wave,color="blue";
      limits;//  mouse();
      };
      z=0.
        //sigm=u(,3);
      
      mask=flux*0.+0.;
      sigm=mask*0.;
        //  flux*(1./SNRl(where((h.object)==namel)(1)));
      gal=galStruct();
      gal.redshift=z;
      gal.SNR=SNR(1);
        //SNRl(where((h.object)==namel)(1)); // getSNR(filelist(i));
      gal.name= h.object;
      gal.filename=strtok(filelist(i),".")(1)+".pdb";
      gal.resfile=strtok(filelist(i),".")(1)+".res";
      gal.result_dir="/"+strjoin(split2words(gal.filename,sep="/")(1:-1),"/")+"/";
      f=createb(gal.filename);
      save,f,gal,flux,wave,sigm,mask;
      if(!is_void(R1)) save,f,R1;
      if(!is_void(R)) save,f,R;
         
      gal.redshift;
      grow,gals,gal;
    }

  return gals;
};


//func fits2pdb(fname)
//{
//  return strtok(fname,".")(1)+".pdb";
//}


//func pdb2fits(fname)
//{
//  return strtok(fname,".")(1)+".fits";
//}

func wait4me(nh){
  // prints a counter for nh hours 
  for(i=1;i<=nh;i++){
    pause,3600000;print,i;
  };
};
    
func mergepdbs(pdblist,nameout){
  /*DOCUMENT
    BRUTAL SUMMING OF .pdb, assumes that wavelength calibration is exactly identical, redshift also etc...
    give nameout without .pdb
    returns a galstruct
    wont work if the sizes in the pdbs are not EXACTLY identical
  */

  nl=numberof(pdblist);
  upload,pdblist(1),s=1;
  rflux=flux*0.;
  for (i=1;i<=nl;i++){
    upload,pdblist(i),s=1;
    rflux+=flux;
  };
  flux=rflux;

  gal.filename=nameout+".pdb";
  gal.resfile=nameout+".res";
  if (strmatch(nameout,"/")==1) gal.result_dir=strjoin(split2words(nameout,sep="/")(1:-1),"/")+"/";
  if (strmatch(nameout,"/")==0) {
    di=array(string,[1,2]);
    di(1)=(exec("pwd"))(1);
    di(2)="/";
    gal.result_dir=strjoin(di);
  };
  save,createb(gal.filename),gal,flux,wave,sigm,mask;
  nl;
  gal.redshift;
  return gal;
};
  
  



func convert_all(filelist,&wave,cut=,noplot=,log=,z0=,SNR0=,wav=,wavaxis=,xs=,xe=,hdu=,fsigm=,errorfile=,nosav=){ 
  /* DOCUMENT
     Attempt at making a comprehensive conversion routine that looks at the proper keywords and warns if not present. is quite verbose at the moment. Will try to comply to fitscws someday. At present its ok but not in the details
     
     OPTIONS:
     if fits file contains multiple hdus, specify with hdu=
     cut: filelist is cut at cut-th file, default is 10
     noplot   default is noplot=0
     log: log=1 enforces log wave sampling in case fits header info is inaccurate
     z0: if redshift is not supplied in fits header, can be given by user as z0
     SNR0: same as z0 for SNR. Best way though would be to provide sigma along with spectrum as in SDSS
     wav: ??
     wavaxis: 1 or 2. Useful if data is provided as 2d frame, typical for long-slit spectroscopy. CHOICE IS NOT IMPLEMENTED
     if wavaxis=[] we try to guess the wavelength axis as the one with largest dimension.
     xs, xe: start and end of the stacking in the spatial direction (again useful in long slit spectroscopy).
     sigm : sigm=1 -> forces sigm=1 instead of d/(S/N)
     errorfile: specify a name for an error file, and uses it to fill the sigm vector
     will only work for 1D data right now
     
     
     SEE ALSO: convertSDSS, convertVAKU, convertJ
  */
  
  local nl;
  if(is_void(nosav)) nosav=1;
  if(is_void(SNR0)) SNR0=100.;
  if(is_void(fsigm)) fsigm=0;
  if(is_void(cut)) cut=10;
  if(is_void(z0)) z0=0.;
  if(is_void(xs)) xs=1;
  if(is_void(xe)) xe=0;
  if(is_void(hdu)) hdu=1;
  gals=[];
  for(i=1;i<=min(cut,numberof(filelist));i++)
    {
      
      a=fits_read(filelist(i),h,hdu=hdu);
      //wave=fits_get(h,"CRVAL1")+indgen(dimsof(a)(2))*fits_get(h,"CDELT1");
      if (is_void(fits_get(h,"CRVAL1")))
        {
          write,"wavelengths origin CRVAL1 not provided     EXITED";
          return 0;
          continue;
        };

      ndims=fits_get(h,"NAXIS");
      
      if (is_void(fits_get(h,"NAXIS1")))
        {
          write,"number of pixels NAXIS1 not provided, will guess from largest dimension of data array";
          nl=max(dimsof(a)); write,nl;
        };
      if (!is_void(fits_get(h,"NAXIS1"))) {nl=fits_get(h,"NAXIS1");};
      nw=nl;
      cdeltst="CDELT1";
      ctypst="CTYPE1";
      crvalst="CRVAL1";
      crpixst="CRPIX1";
      u=a;
      flux=u;
      //      ndims=1;
      
      if (fits_get(h,"NAXIS")==2) {
        nc=fits_get(h,"NAXIS2");
        //        ndims=2;
        if (((nc>nl)&(is_void(wavaxis)))|(wavaxis==2)) {
          write, "NAXIS2>NAXIS1, assuming AXIS2 wavelength axis";
          cdeltst="CDELT2";
          ctypst="CTYPE2";
          crvalst="CRVAL2";
          nw=nc;
          //flux=transpose(flux);
          if (is_void(xe)) xe=nl;
          flux=(flux(xs:xe,))(sum,);
        };
        if((nl>nc)|(wavaxis==1)){
          nw=nl;
          if (is_void(xe)) xe=nc;
          flux=(flux(,xs:xe))(,sum);
        };
      };
      
      if(log==0) wavesampling="LIN";
      if(log==1) wavesampling="LOG";

      
      if (is_void(fits_get(h,cdeltst)))
        {
          ocdeltst=cdeltst;
          ncdeltst="CD1_1";   // FIX ME FOR WHEN DATA IS TRANSPOSED
          write,"wavelength increment "+cdeltst+"  not provided.... will try "+ncdeltst;
          cdeltst=ncdeltst;
          if (is_void(fits_get(h,cdeltst))) {
            write,"wavelength increment "+cdeltst+"  not provided.... exiting";
            return 0;
          };
        };
      
      if (is_void(fits_get(h,ctypst)))
        {
          write,ctypst+" not provided, will assume linear";
          wavesampling="LIN";
        };

      wavesampling="LIN";

      crpix1=fits_get(h,crpixst);
      
      if (is_void(fits_get(h,crpixst))) {
        write,crpixst+ "not provided, assuming 1" 
          crpix1=1;
      };
          
      if (!is_void(fits_get(h,ctypst))){
        if ((strmatch(fits_get(h,ctypst),"LOG")==1)) wavesampling="LOG";
      };
      if(log==1) wavesampling="LOG";
      
      if((wavesampling=="LIN")&(log!=1)){
        wave=fits_get(h,crvalst)+(float(indgen(nw))-float(crpix1))*fits_get(h,cdeltst);
      };

      if((wavesampling=="LOG")|log==1){
        wave=10^(fits_get(h,crvalst)+fits_get(h,cdeltst)*(float(indgen(nw)) -float(crpix1)));
      };
      
      //li=fits_list_header(h,pr=0);
      li="truc";
      
      if (is_void(noplot)) {
        //fma;
        if (ndims==2) {ws,0;plia,a;};
        ws,1;plh,flux,wave,color=__rgb(,(i%20)+10);
        xyleg,"Angstroms";
        limits;//  mouse();
      };
      //z=zl(where((h.object)==namel)(1))
      //flux=u(,75:300)(,sum); // sum lines of CCD array
      
      if(is_void(li(where(strmatch(li,"REDSHIFT"))))) {
        write,"no redshift in fits header... will assume 0 or take user input";
      };
      if(!is_void(li(where(strmatch(li,"REDSHIFT"))))) {
        z0=fits_get(h,"REDSHIFT");
      };
      z=z0;
      
      if(is_void(li(where(strmatch(li,"SNR"))))) {
        write,"no S/N in fits header... will assume 100 or take user input";
      };
      if(!is_void(li(where(strmatch(li,"SNR"))))) {
        SNR0=fits_get(h,"SNR");
      };
      
      SNR=SNR0;
      
      mask=flux*0.+0.;
      if(fsigm==0) sigm=((abs(flux))^1.)*(1./SNR);
      if(fsigm==1) sigm=flux*0.+1.;
      if(!is_void(errorfile)) {
          sigm=fits_read(errorfile); // need to normalize ?
      };

       if (is_void(noplot)) {
        //fma;
        if (ndims==2) {ws,0;plia,a;};
        ws,1;plh,flux,wave,color="black";
        plh,sigm,wave,color="cyan";
        xyleg,"Angstroms";
        limits;//  mouse();
      };

      
          
      gal=galStruct();
      gal.redshift=z;
      gal.SNR= SNR; // getSNR(filelist(i));
      gal.name= is_void(fits_get(h,"OBJECT"))?"unknown":fits_get(h,"OBJECT");
      gal.filename=strreplace(filelist(i),".fits",".pdb");
      if(gal.filename==filelist(i)) gal.filename=strreplace(filelist(i),".fit",".pdb");
      gal.resfile=strreplace(filelist(i),".fits",".res");
      if(gal.resfile==filelist(i)) gal.resfile=strreplace(filelist(i),".fit",".res");
      if (strmatch(filelist(i),"/")==1) gal.result_dir=strjoin(split2words(gal.filename,sep="/")(1:-1),"/")+"/";
      if (strmatch(filelist(i),"/")==0) {
        di=array(string,[1,2]);
        di(1)=(exec("pwd"))(1);
        di(2)="/";
        gal.result_dir=strjoin(di);
      };
      gal.wavesampling=wavesampling;
      if(nosav==0){
        f=createb(gal.filename);
        save,f,gal,flux,wave,sigm,mask;
        close,f;
      };
      write,"redshift "+pr1(gal.redshift);

      
      grow,gals,gal;
    }

  return gals;
};



func convert_from_mina(filelist,cut=,noplot=,log=,z0=,SNR0=,wav=,wavaxis=,xs=,xe=){ 
  /* DOCUMENT
     Attempt at making a comprehensive conversion routine that looks at the proper keywords and warns if not present. is quite verbose at the moment. Will try to comply to fitscws someday. At present its ok but not in the details
     
     OPTIONS:
     cut: filelist is cut at cut-th file, default is 10
     noplot   default is noplot=0
     log: log=1 enforces log wave sampling in case fits header info is inaccurate
     z0: if redshift is not supplied in fits header, can be given by user as z0
     SNR0: same as z0 for SNR. Best way though would be to provide sigma along with spectrum as is SDSS
     wav: ??
     wavaxis: 1 or 2. Useful if data is provided as 2d frame, typical for long-slit spectroscopy. CHOICE IS NOT IMPLEMENTED
     if wavaxis=[] we try to guess the wavelength axis as the one with largest dimension.
     xs, xe: start and end of the stacking in the spatial direction (again useful in long slit spectroscopy).
     
     SEE ALSO: convertSDSS, convertVAKU, convertJ
  */
  
  local nl;
  
  if(is_void(SNR0)) SNR0=100.;
  if(is_void(cut)) cut=10;
  if(is_void(z0)) z0=0.;
  if(is_void(xs)) xs=1;
  if(is_void(xe)) xe=0;
  gals=[];
  for(i=1;i<=min(cut,numberof(filelist));i++)
    {
      
      a=fits_read(filelist(i),h,hdu=1);
      wave=a(,1);
      flux=a(,2);
      wavesampling="LIN";
      
      li="truc";
      
      if (is_void(noplot)) {
        //fma;
        ws,1;plh,flux,wave,color=__rgb(,(i%20)+10);
        xyleg,"Angstroms";
        limits;//  mouse();
      };
      
      if(is_void(li(where(strmatch(li,"REDSHIFT"))))) {
        write,"no redshift in fits header... will assume 0 or take user input";
      };
      if(!is_void(li(where(strmatch(li,"REDSHIFT"))))) {
        z0=fits_get(h,"REDSHIFT");
      };
      z=z0;
      
      if(is_void(li(where(strmatch(li,"SNR"))))) {
        write,"no S/N in fits header... will assume 100 or take user input";
      };
      if(!is_void(li(where(strmatch(li,"SNR"))))) {
        SNR0=fits_get(h,"SNR");
      };
      
      SNR=SNR0;
      
      mask=flux*0.+0.;
      sigm=((abs(flux))^1.)*(1./SNR);
      gal=galStruct();
      gal.redshift=z;
      gal.SNR= SNR; // getSNR(filelist(i));
      gal.filename=strreplace(filelist(i),".fits",".pdb");
      gal.name=gal.filename;
      gal.resfile=strreplace(filelist(i),".fits",".res");
      if (strmatch(filelist(i),"/")==1) gal.result_dir=strjoin(split2words(gal.filename,sep="/")(1:-1),"/")+"/";
      if (strmatch(filelist(i),"/")==0) {
        di=array(string,[1,2]);
        di(1)=(exec("pwd"))(1);
        di(2)="/";
        gal.result_dir=strjoin(di);
      };
      gal.wavesampling=wavesampling;
      save,createb(gal.filename),gal,flux,wave,sigm,mask;
      gal.redshift;
      grow,gals,gal;
    }

  return gals;
};



func convert_4388(filelist,cut=,noplot=,log=,z0=,SNR0=,wav=,wavaxis=,xs=,xe=){ 
  /* DOCUMENT
     Attempt at making a comprehensive conversion routine that looks at the proper keywords and warns if not present. is quite verbose at the moment. Will try to comply to fitscws someday. At present its ok but not in the details
     
     OPTIONS:
     cut: filelist is cut at cut-th file, default is 10
     noplot   default is noplot=0
     log: log=1 enforces log wave sampling in case fits header info is inaccurate
     z0: if redshift is not supplied in fits header, can be given by user as z0
     SNR0: same as z0 for SNR. Best way though would be to provide sigma along with spectrum as is SDSS
     wav: ??
     wavaxis: 1 or 2. Useful if data is provided as 2d frame, typical for long-slit spectroscopy. CHOICE IS NOT IMPLEMENTED
     if wavaxis=[] we try to guess the wavelength axis as the one with largest dimension.
     xs, xe: start and end of the stacking in the spatial direction (again useful in long slit spectroscopy).
     
     SEE ALSO: convertSDSS, convertVAKU, convertJ
  */
  
  local nl;
  
  if(is_void(SNR0)) SNR0=100.;
  if(is_void(cut)) cut=10;
  if(is_void(z0)) z0=0.;
  if(is_void(xs)) xs=1;
  if(is_void(xe)) xe=0;
  gals=[];
  for(i=1;i<=min(cut,numberof(filelist));i++)
    {
      
      a=fits_read(filelist(i),h,hdu=1);
      //wave=fits_get(h,"CRVAL1")+indgen(dimsof(a)(2))*fits_get(h,"CDELT1");
      if (is_void(fits_get(h,"CRVAL1")))
        {
          write,"wavelengths origin CRVAL1 not provided     EXITED";
          //return 0;
          continue;
        };

      ndims=fits_get(h,"NAXIS");
      
      if (is_void(fits_get(h,"NAXIS1")))
        {
          write,"number of pixels NAXIS1 not provided, will guess from largest dimension of data array";
          nl=max(dimsof(a)); write,nl;
        };
      if (!is_void(fits_get(h,"NAXIS1"))) {nl=fits_get(h,"NAXIS1");};
      nw=nl;
      cdeltst="CDELT1";
      ctypst="CTYPE1";
      crvalst="CRVAL1";
      crpixst="CRPIX1";
      u=a;
      flux=u;
      //      ndims=1;
      
      if (fits_get(h,"NAXIS")==2) {
        nc=fits_get(h,"NAXIS2");
        //        ndims=2;
        if (nc>nl) {
          write, "NAXIS2>NAXIS1, assuming AXIS2 wavelength axis";
          cdeltst="CDELT2";
          ctypst="CTYPE2";
          crvalst="CRVAL2";
          nw=nc;
          //flux=transpose(flux);
          if (is_void(xe)) xe=nl;
          flux=(flux(xs:xe,))(sum,);
          info,flux;
        };
        if(nl>nc){
          write,"NAXIS1>NAXIS2,  assuming AXIS1 wavelength axis";
          nw=nl;
          if (is_void(xe)) xe=nc;

          flux=flux(,1,);
          flux=(flux(,xs:xe))(,,sum);
          write,dimsof(flux);
          info,flux;
        };
      };
      
      if(log==0) wavesampling="LIN";
      if(log==1) wavesampling="LOG";

      
      if (is_void(fits_get(h,cdeltst)))
        {
          ocdeltst=cdeltst;
          ncdeltst="CD1_1";   // FIX ME FOR WHEN DATA IS TRANSPOSED
          write,"wavelength increment "+cdeltst+"  not provided.... will try "+ncdeltst;
          cdeltst=ncdeltst;
          if (is_void(fits_get(h,cdeltst))) {
            write,"wavelength increment "+cdeltst+"  not provided.... exiting";
            return 0;
          };
        };
      
      if (is_void(fits_get(h,ctypst)))
        {
          write,ctypst+" not provided, will assume linear";
          wavesampling="LIN";
        };

      wavesampling="LIN";

      crpix1=fits_get(h,crpixst);
      
      if (is_void(fits_get(h,crpixst))) {
        write,crpixst+ "not provided, assuming 1" 
          crpix1=1;
      };
          
      if (!is_void(fits_get(h,ctypst))&(strmatch(fits_get(h,ctypst),"LOG")==1)&(log!=0)) wavesampling="LOG";
      
      if((wavesampling=="LIN")&(log!=1)){
        wave=fits_get(h,crvalst)+(float(indgen(nw))-float(crpix1))*fits_get(h,cdeltst);
      };
      
      if(wavesampling=="LOG"){
        wave=10^(fits_get(h,crvalst)+fits_get(h,cdeltst)*(float(indgen(nw)) -float(crpix1)));
      };
      
      //li=fits_list_header(h,pr=0);
      li="truc";
      
      if (is_void(noplot)) {
        //fma;
        if (ndims==2) {ws,0;plk,a;};
        ws,1;plh,flux,wave,color=__rgb(,(i%20)+10);
        xyleg,"Angstroms";
        limits;//  mouse();
      };
      //z=zl(where((h.object)==namel)(1))
      //flux=u(,75:300)(,sum); // sum lines of CCD array
      
      if(is_void(li(where(strmatch(li,"REDSHIFT"))))) {
        write,"no redshift in fits header... will assume 0 or take user input";
      };
      if(!is_void(li(where(strmatch(li,"REDSHIFT"))))) {
        z0=fits_get(h,"REDSHIFT");
      };
      z=z0;
      
      if(is_void(li(where(strmatch(li,"SNR"))))) {
        write,"no S/N in fits header... will assume 100 or take user input";
      };
      if(!is_void(li(where(strmatch(li,"SNR"))))) {
        SNR0=fits_get(h,"SNR");
      };
      
      SNR=SNR0;
      
      mask=flux*0.+0.;
      sigm=((abs(flux))^1.)*(1./SNR);
      gal=galStruct();
      gal.redshift=z;
      gal.SNR= SNR; // getSNR(filelist(i));
      gal.name= is_void(fits_get(h,"OBJECT"))?"unknown":fits_get(h,"OBJECT");
      gal.filename=strreplace(filelist(i),".fits",".pdb");
      gal.resfile=strreplace(filelist(i),".fits",".res");
      if (strmatch(filelist(i),"/")==1) gal.result_dir=strjoin(split2words(gal.filename,sep="/")(1:-1),"/")+"/";
      if (strmatch(filelist(i),"/")==0) {
        di=array(string,[1,2]);
        di(1)=(exec("pwd"))(1);
        di(2)="/";
        gal.result_dir=strjoin(di);
      };
      gal.wavesampling=wavesampling;
      save,createb(gal.filename),gal,flux,wave,sigm,mask;
      gal.redshift;
      grow,gals,gal;
    }

  return gals;
};

