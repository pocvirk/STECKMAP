// TO DO:
//fix kin=0 when nMC!=0
// check the gradient clipping of E(B-V) in W1faz1cqq  
// clean redundancies in L1-L2-L3 definitions  ok
// put an automatic cov <-1 and mucov if epar=3 and kin=1  ok
// OR give an option for interpolating data on basis support, i.e. Quick and dirty
// see how to limit installed dependencies ok
// The handling of the 2d age-kin routine is a bit messy, due to differences in handling the fft and the rolls and that kind of stuff. still LOSVD seems inverted in 2d compared to 1d...
// QUESTION WHILE SPLINING NPEC: is (splineinterp(g))'=splineinterp(g') ? seems so since it works....
// Why does the fit in model sampling window always look better than the one in initial sampling window ? I might not be going back properly to the initial data sampling because of the resampling, then normalization steps and so on...
// probably rather than resampling the data towards the model support I should leave the data untouched and resample the high res model towards the data support by using a resampling matrix to preserve the ability of computing analytical gradients.
// ADDED user-supplied AMR through penalization with a prior:
// replaced W1faz1cqq5 and 6 respectively with 7 and 8

//#include "digit2.i"
#include "STECKMAP/Pierre/POP/pop_paths.i"
#include "string_utils.i"
#include "STECKMAP/Eric/plot.i"
#include "STECKMAP/Chris/utls.i"
// includes a lot of stuff: 
#include "STECKMAP/Chris/linterp.i"
#include "STECKMAP/Pierre/POP/Qfunctions_sfit.i"
#include "STECKMAP/Eric/invpb.i"
#include "STECKMAP/Eric/optim_pack.i"
#include "STECKMAP/Eric/fft_utils.i"
// Eric/fft_utils.i requires Eric/utils.i
#include "STECKMAP/Eric/random.i"
// requires gamma.i, which is ok.
#include "STECKMAP/Pierre/POP/2DAK2ls.i"
#include "STECKMAP/Pierre/POP/vpopcin2ls4.i"  
//#include "Pierre/POP/moments.i"      
#include "STECKMAP/Pierre/POP/constantes.i"
#include "STECKMAP/Pierre/POP/pop2ls.i"
#include "STECKMAP/Pierre/POP/math2ls.i"
#include "STECKMAP/Pierre/POP/conv2ls.i"
#include "STECKMAP/Pierre/POP/gal2ls.i"
#include "STECKMAP/Pierre/POP/plot2ls.i"
#include "STECKMAP/Pierre/POP/moments.i"
#include "STECKMAP/Pierre/POP/file2ls.i"


//******************* Objective functions  ***********

func Q0(x,&g,fpar){// unclipped
  // NO KIN
  if ((fpar==[0,1])(avg)==1)  return Wqcfaz2(x,0,g);
  if ((fpar==[0,0])(avg)==1)  return Wqcfaz3(x,0,g);
  if ((fpar==[0,3])(avg)==1)  return Wqcfaz5(x,0,g);
  

  // WITH KIN
  if ((fpar==[1,0])(avg)==1)  return W1faz1cqq2(x,0,g);
  if ((fpar==[1,1])(avg)==1)  return W1faz1cqq8(x,0,g);
  if ((fpar==[1,2])(avg)==1)  return Wfaz1cqq(x,0,g);
  if ((fpar==[-1,1])(avg)==1) return Rkin(x,0,g);
  //if ((fpar==[1,3])(avg)==1)  return W1faz1cqq5(x,0,g); //without priors
  if ((fpar==[1,3])(avg)==1)  return W1faz1cqq7(x,0,g);  // with AMR prior
  if ((fpar==[2,3])(avg)==1)  return AKNPEC3(x,g); 
  
};

func Q1(x,&g,fpar){// clipped in E(B-V)
  // NO KIN
  if ((fpar==[0,1])(avg)==1)  return Wqcfaz2(x,2,g);
  if ((fpar==[0,0])(avg)==1)  return Wqcfaz3(x,2,g);
  if ((fpar==[0,3])(avg)==1)  return Wqcfaz5(x,0,g); // no need to clip the NPextinction Curve

  // WITH KIN
  if ((fpar==[1,0])(avg)==1)  return W1faz1cqq2(x,2,g);
  if ((fpar==[1,1])(avg)==1)  return W1faz1cqq8(x,2,g);
  if ((fpar==[1,2])(avg)==1)  return Wfaz1cqq(x,2,g);
  if ((fpar==[-1,1])(avg)==1) return Rkin(x,2,g);
  //if ((fpar==[1,3])(avg)==1)  return W1faz1cqq5(x,2,g); //without priors
  if ((fpar==[1,3])(avg)==1)  return W1faz1cqq7(x,2,g);  // with AMR prior

  if ((fpar==[2,3])(avg)==1)  return AKNPEC3(x,g);
};

//***************************************************************************
//   Just an example file
//   ref: Vazdekis et al.2001, ApJL, 551, L127
//***************************************************************************
fV=examplesdir+"NGC4621-1.fits";


//===========================================================================
//************************* THE FITTING ENGINE ******************************
//===========================================================================

func sfit(gallist,base,&_ki,&_resi,nMC=,MC=,verb=,maxit=,maxitMC=,noskip=,meval=,mevalMC=,meval3=,kin=,epar=,MASK=,guess=,mus=,vlim=,sav=,asc=,RMASK=,zlim=,rmin=,bf=,padv=,frtol1=,N=,nor=,nde=,sigdef=,L1=,L2=,L3=,L4=,muc=,co=,pl=,mucov=,cov=,parage=,RG=,muv=,mux=,muz=,mue=,mub1=,mub2=,z2d=,bnpec=,splin=,__sigm=,AMRp=,muAMRp=,rrAMRp=,Nb=,lsp=,plavgW=,nedgemask=,dMC=){
  /* DOCUMENT

  ***************************************************************************
  =============================== DOCUMENTATION =============================
  ***************************************************************************
  
     gallist is a list of galaxies (structure described in SDSS-gal.i)
     base is the SSP basis to interpret gallist (base struct described above)
     _ki variable to store all MC residuals 
     nMC: nb of MC
     maxit: maximum number of iterations for the first solution
     maxitMC: same for the MC exploration
     meval: maximum number of evaluations fot the first solution
     mevalMC: same for MC exploration
     dMC: choice of the mock data used for MonteCarlo exploration:
         "bestmodel" -> the best model is used, noised to match gal.snr. expect chi2~1
         "data" -> the initial data is used, noised to match gal.snr or using the error spectrum if available. expect chi2~2 following dispersion(normal(mu,sigma)+normal(mu,sigma))= sqrt(2)*sigma
         default is "bestmodel"
         
         noskip: force analysis even if results file already exists
     epar: type of extinction:
         2 -> n ebv
         1 -> 1 ebv
         0 -> continuum correction by degree 2 polynomial (not working well)
         3 -> non parametric flux correction with nde anchor points (NPEC)
     kin: type of kinematic treatment.
         0 -> nope.
         1 -> 1d kinematics.
         2 -> 2d kinematics (not yet done) (2d kin are only available with epar=3 and paratmeter z2d must be given by user as the number of the metallicity to be used)
     
     guess: starting point for the iterative minimization.
     mus: Vector containing the smoothing and binding parameters
     muv, mux, muz, mue, mub1, mub2
     Default mus=[2.e2, 1.e0,1.e2,1.e2,1.e4,1.e8]
     when epar=3 mue sets the smoothing of the NPEC

     vlim: velocity limits of the losvd. array of 2 double [vinf,vsup]

     sav: if 1 saves into resfile, default is 0, so not saving
     asc: if 1 converts the results file .res into an ascii file with explicit variables
     NB: if asc=1 then sav=1

     RMASK: user-provided mask in rest frame: pixels corresponding to a wavelength listed in RMASK are masked RMASKS are of the form [[linf1,lsup1],[linf2,lsup2],...]
     NOTE: if only one region is to be masked, it is still needed to write
     RMASK=[[linf1,lsup1]]
     zlim: limits in transformed metallicities (useful if one does not want to use the full z range of the basis or probe age-z degeneracy effects) array of 2 double [zinf, zsup], by default set to the limits of the basis
     MC: type of MC to perform:
         "MC" makes regular MC (i.e. inversion of a noised model)
         "RG" makes inversion with random guesses (bf*random_normal). Useful to probe secondary extrema  
     defaut is MC="MC"

     rmin: In "RG" mode the minimization is performed until the objective function is equal or smaller than rmin. If rmin is real, then the specified value is used as is. If rmin is a character (for example rmin="1.05") then the minimization is performed down to rmin*res1, where res1 is the value of the objective function at the first solution found.
     
     pad=[pad1,pad2]: left and right padding lengths in pixels, needs to be tuned only if algo experiences unlucky slow execution. By default pad="AUTO" so that spectrum is padded to a length favorable for fft.

     N: N=1 normalizes the input spectrum so that d(avg)=1 before fitting.
     by default N=1
     but *********   SET N=0 FOR SIMULATIONS !!! ***********

     *********************************
     ===========  IMPORTANT ==========
     *********************************
     
     Nb: Nb=1 normalizes the basis same way as the data
     by default Nb=1 -> working with FLUX FRACTIONS BY DEFAULT!!
     sso set Nb=0 to work directly in masses or anything else
     or Nb="smart"
          
     L1 is the form of penalization to use for the sad and can be "D1","D2" or "D3"
     "D1"-> Gradient 
     "D2"-> Laplacian
     "D3"-> 3rd order
     L2 is the penalization for the losvd
     L3 is the penalization for the AMR
     L4 is the penalization for the non parametric flux recalibration (this penalization is only effective with kin=2 though... My belief is that it wont change dramatically what happens with the kin=1 runs)
     note that setting L3 to D1 and entering a large muz (see mus) will yield a mono-metallic solution.
     If L3=D2 and muz is large, muz actually sets the maximum slope of the AMR
     The same can be done for E(B-V)

     IF kin=2 then L1 is the penalization in the age direction and L2 is the penalization in the velocity direction, quite consistently with the previous definitions
     z2d is the assumed metallicity for the whole population.
     
     nor gives the level of smoothing used for normalization of spectrum
     This has until now to be done in internal, when data and basis have same wavelength sampling otherwise doesnot work properly
     nde is the number of anchor points for the non parametric flux recalibration method (or NPEC: Non Parametric Extinction Curve)

     bnpec gives the choice of basis function used for expanding the non-parametric extinction curve. "tri" will give "triangular" functions (corrsponding to linear interpolation) while "spl" will use cubic splines
     default is cubic splines this parameter is global and used as is by the functions GENdecdde, GENdx0dx1 and npe
     
     bf: amplitude of the noise for the random guess MC ("RG"), default is auto, i.e. scaled on best solution
     splin: type of interpolation used for resampling of the data (AND ONLY THE DATA)if needed:
     1: spline interpolation
     0: linear interpolation
     NB: the basis is always spline-interpolated

     parage: parametric fitting of age and Z. output is 1age, 1Z. No kinematics are possible in this mode. (IMPLEMENT THAT SOMEDAY)

     AMRp=,muAMRp=,rrAMRp= are used to supply respectively a prior for the AMR, the weighting parameter of the prior and the kernel for computing prior-solution distances.

     lsp: 1-> large display for the spectrum in model sampling default is 1
     plavgW: 1-> plots an W/W(avg) for better dynamic range in the window
     nedgemask -> number of masked pixels on the edge of the spectrum (default is 25, mandatory when kin=1)

     
======================================================================
     There are several steps:
     
     - cross spectral range. The highest resolution from data and models is kept
     - blueshift observed spectrum according to gallist.REDSHIFT
     - log-resample wavelengths wave0...
     - Genereate weights (W) (from variance-covariance matrix and masks)
     - Generate smoothing kernels L1,L2,L3 as needed
     - Minimize the objective function:
          1) with little metallicity constraints
          2) with stronger  metallicity constraints
          3) a third useless step ? used to be there for E(B-V) clipping

     - save
     - Repeat minimization steps for MC
     - plots the best solution
     - save
     - convert to ascii if specified so by user
     
     WARNING: the stellar age distribution part of the result has to be squared
     WARNING: will not work properly for lists of files. Loop outside sfit

     MORE COMMENTS: The log-resampling has to be done inside the sfit since another resampling is involved when intersecting the wavelengths domains of data and model
     ABOUT SAD AND LOSVD: THE SAVED VALUES ARE SQRT(SAD) and SQRT(LOSVD)
     
  */
  
  extern bloc,_x0,x0,dd,ab,_m,d,rr,rr1,rr3,W,u,mub,pW,bdata,d,model0,q0,res0,model1,q1,res1,dfile,rv1,rv2,riv,nlos,pmodel0,pmodel1,ni,nj,nr,pad1,pad2,nab;

  if(is_void(base)) {write,"No basis provided";return [];};

  
  if (is_void(maxit)) maxit=1;  // ????   NOT USED !!
  if (is_void(meval)) meval=500;
  if (is_void(meval3)) meval3=meval;
  if (is_void(maxitMC)) maxitMC=maxit;
  if (is_void(mevalMC)) mevalMC=meval;
  if (is_void(noskip)) noskip=0;
  if (is_void(kin)) kin=0;  // 0 -> no kin  ; 1 -> 1d kin ; 2 -> 2d kin (to come) ; -1 -> when data is higher resolution than model, doesnt work yet
  if (is_void(epar)) epar=1;  // 1 -> extinction calzetti-like; 0 -> polynomial continuum correction; 2 -> age-dependent extinction 3 -> NPEC
  // if (is_void(mus)) {mus=setmus(gallist.SNR(1)/sqrt(base.R(1)/10000.));}; // whatever ??? changed on 24 feb 05
  if (is_void(mus)) {mus=setmus(gallist.SNR(1));};
  if (kin==0) vlim=[];
  if (kin==2) epar=3;
  if (is_void(z2d)) z2d=1; // DEFAULT TO HIGHEST METALLICITY IN 2D Age-kin inversion
  if (is_void(vlim)) {vlim=[-900.,900.];};
  if (is_void(sav)) sav=0;
  if (is_void(asc)) asc=0;
  if (asc==1) sav=1;
  if (is_void(zlim)) zlim=[min(base.met),max(base.met)];
  if (is_void(MC)) MC="MC";
  if (is_void(dMC)) dMC="bestmodel";
  if (!((dMC=="bestmodel")|(dMC=="data"))) {
      write,"dMC has illegal value. set to \"bestmodel\" or \"data\" ";
      error;
  }; // ADDED PO 04/11/2008
  if (is_void(bf)) bf="auto";
  if (is_void(rmin)) rmin=0.;
  if (is_void(padv)) padv="auto";
  if (is_void(MASK)) MASK="nope";
  if (is_void(frtol1)) frtol1=1.e-20;
  if (is_void(N)) N=1;   //   Normalization of DATA!
  if (base.basistype=="flux/(Msol/yr)") Nb=0;
  if (base.basistype=="flux/Msol") Nb=0;
  if (is_void(Nb)) Nb=1;  // Normalization of BASIS!
  if (is_void(guessepar)) guessepar=[1.,1.,-5.e-5];  // initialization for polynomial flux calibration correction 
  if (is_void(nde)) nde=10;
  if (is_void(bnpec)) bnpec="spl";
  if (is_void(sigdef)) sigdef=1;
  if (!is_void(L1)|is_void(L2)) L2=L1;
  if (!is_void(L1)&is_void(L3)) L3=L1;
  if (!is_void(L1)&is_void(L4)) L4=L1;
  if (is_void(muc)) muc=0.; // normalize the SAD ?
  if (is_void(co)) co=0.;  // normalize the SAD ?
  if (epar==3) {co=1.;muc=1.e4;}; // NORMALIZE by default if NPEC is used
  if ((kin==1) & is_void(mucov)) mucov=1.e4;  // normalize the LOSVD by default if kin=1
  if ((kin==1) & is_void(cov)) cov=1.;  
  if (is_void(mucov)) mucov=0.;  // normalize the LOSVD ?
  if (is_void(cov)) cov=0.; // normalize LOSVd ?
  if (is_void(pl)) pl=1;
  if (is_void(parage)) parage=0;
  if (is_void(RG)) RG=[];
  if (is_void(verb)) verb=100;
  if (is_void(splin)) splin=1;
  if (is_void(__sigm)) sigm=1;
  if (is_void(lsp)) lsp=1;
  if (is_void(plavgW)) plavgW=1;
  if (!is_void(AMRp)&is_void(rrAMRp)){
    dib=dimsof(base.flux);
    nab=dib(3);
    rrAMRp=diag(array(1.,nab));
  };
  if (is_void(muAMRp)&!is_void(AMRp)) muAMRp=1.e5;
  if (is_void(AMRp)) AMRp="none";
  if (is_void(muAMRp)) muAMRp="none";
  if (is_void(rrAMRp)) rrAMRp="none";
  if (is_void(nedgemask)) nedgemask=25;
  
  
  // setting the smoothing parameters according to mus vector
  if (!is_void(muv)) muv=muv;
  else muv=mus(1);
  if(!is_void(mux)) mux=mux;
  else mux=mus(2);
  if (!is_void(muz)) muz=muz;
  else muz=mus(3);
  if (!is_void(mue)) mue=mue;
  else mue=mus(4);
  if (!is_void(mub1)) mub1=mub1;
  else mub1=mus(5);
  if (!is_void(mub2)) mub2=mub2;
  else mub2=mus(6);
  
  // setting v1 and v2 according to vlim vector
  v1=vlim(1);
  v2=vlim(2);

  //setting padding lengths according to padv vector if no auto padding
  if(padv!="auto") {pad1=padv(1);pad2=padv(2);};

  // VARIOUS INITIALIZATIONS
  dib=dimsof(base.flux);
  nab=dib(3);
  mpar=[kin,epar];
  _nq=500;  // length of gres array containing all solutions. Not critical
  //decdde=[]; // Initialization of decdde array for npel
  ns=((dimsof(gallist)(1)))?dimsof(gallist)(2):1;
  //gres=array(0.,3*nab+_nq,ns*(((is_void(nMC)?0:nMC) +1))); //array for storing the results
  gres=[];
  _ki=array(0.,ns*(((is_void(nMC)?0:nMC) +1)));
  _resi=_ki;


  // START WORKING
  for(j=1;j<=ns;j++){
    u=gallist(j);
    f=openb(u.filename(1));                
    restore,f;
    close,f; 

    decdde=[]; // Initialization of decdde array for npel  ADDED 11th 11 2005

    
    // SAVE ORIGINAL VALUES

    sflux=flux;
    swave=wave;
    ssigm=sigm;
    
    if (is_void(sigm)) sigm=1./(u.SNR/flux);
    if (!is_void(__sigm)&(__sigm==1)) sigm=flux*0.+1.;
    if (is_void(mask)) mask=0.*flux;

    
    // TEST existence of resfile and content
    if ((isgres2(u.resfile(1)))&((noskip==0))) {  // if MC exists 
      print,mp_rank," skipping ",u.filename(1);
      continue;
    };
    
    //blueshift wavelengths
    wave0=wave/(1.+u.redshift);

      
    // cross spectral ranges
    wave00=intersectS(wave0,base.wave);
    write,"spectral ranges crossed with ",(wave00(dif))(avg)," Angstrom step, (and dispersion of step ",(wave00(dif))(rms);

    //spline-interpolate data on new support
    if(splin==1) write, " spline-interpolating data on finest common wavelength range";
    if(splin!=1) write, " linear-interpolating data on finest common wavelength range";

    flux=((splin==1)?spline(flux,wave0,wave00):interp(flux,wave0,wave00));
    sigm=((splin==1)?spline(sigm,wave0,wave00):interp(sigm,wave0,wave00));
    mask=((splin==1)?spline(mask,wave0,wave00):interp(mask,wave0,wave00));
    wave0=wave00;    

    if(1){ // this is the old way of doing it.... not sure how right that is...
        //see, we build an exponemtial wavelength support instead of a log support....
        // LOG resample ALWAYS if kin>1
        // log step is the half the finest of lin step
        //if ((u.wavesampling=="LIN")&(kin>0)){
        if(kin>0) {
            write,"this is old days kin"
            ldl=log10(1.+0.5*(abs(wave0(dif))(min))/wave0(1)); // that used to be my choice.... maybe it undersamples ?
            ldl=log10(1.+0.5*(abs(wave0(dif))(min))/wave0(0)); // changed PO 07/11/2008 this way it should oversample
            nlogwave=int(log10(wave0(0)/wave0(1))/ldl);
            logwave=wave0(1)*10^(ldl*indgen(nlogwave));
            flux=((splin==1)?spline(flux,wave0,logwave):interp(flux,wave0,logwave));
            sigm=((splin==1)?spline(sigm,wave0,logwave):interp(sigm,wave0,logwave));
            mask=((splin==1)?spline(mask,wave0,logwave):interp(mask,wave0,logwave));
            wave0=logwave;  
        };
    };

    if(0){ // this is the NEW way of doing the log resample.
        // LOG resample ALWAYS if kin>1
        // log step is the half the finest of lin step
        //if ((u.wavesampling=="LIN")&(kin>0)){
        if(kin>0) {
            write,"this is new kin";
            ldl=log10(1.+0.5*(abs(wave0(dif))(min))/wave0(1)); // that used to be my choice.... maybe it undersamples ?
            ldl=log10(1.+0.5*(abs(wave0(dif))(min))/wave0(0)); // changed PO 07/11/2008 this way it should oversample or at least it leads to a smaller ldl...
            nlogwave=int(log10(wave0(0)/wave0(1))/ldl);
            logwave_log=log10(wave0(1))+ldl*indgen(nlogwave);
            logwave=10^(logwave_log);
            flux=((splin==1)?spline(flux,wave0,logwave):interp(flux,log10(wave0),logwave_log)); // well thats just the same!
            sigm=((splin==1)?spline(sigm,wave0,logwave):interp(sigm,log10(wave0),logwave_log));
            mask=((splin==1)?spline(mask,wave0,logwave):interp(mask,log10(wave0),logwave_log));
            wave0=logwave;  
        };
    };

    if(0){ // this is the NEW way of doing the log resample.
        //I dont believe this should make a difference but I decided to replace all log10 with log
        // LOG resample ALWAYS if kin>1
        // log step is the half the finest of lin step
        //if ((u.wavesampling=="LIN")&(kin>0)){
        if(kin>0) {
            write,"this is new kin with log neper instread of log10";
            ldl=log(1.+0.5*(abs(wave0(dif))(min))/wave0(1)); // that used to be my choice.... maybe it undersamples ?
            ldl=log(1.+0.5*(abs(wave0(dif))(min))/wave0(0)); // changed PO 07/11/2008 this way it should oversample or at least it leads to a smaller ldl...
            nlogwave=int(log(wave0(0)/wave0(1))/ldl);
            logwave_log=log(wave0(1))+ldl*indgen(nlogwave);
            logwave=exp(logwave_log);
            flux=((splin==1)?spline(flux,wave0,logwave):interp(flux,log(wave0),logwave_log)); // well thats just the same!
            sigm=((splin==1)?spline(sigm,wave0,logwave):interp(sigm,log(wave0),logwave_log));
            mask=((splin==1)?spline(mask,wave0,logwave):interp(mask,log(wave0),logwave_log));
            wave0=logwave;  
        };
    };

    
    
    wmin=1;
    wmax=numberof(wave0);
        
    //spline-interpolate basis on new support
    // not a good idea to linearly interpolate models anyway... i know but....spline may be adding spectral shape out of the blue...
    if(splin==1){
        write,"spline-interpolating basis in lambda";
        bloc=spline2(base.flux,base.wave,wave0);
    };
    if(splin!=1){
        write,"linear-interpolating basis in lambda";
        bloc=interp2(base.flux,base.wave,wave0);
    };
    //bloc=spline2(base.flux,base.wave,wave0); // commented 10/11/2008 ...
    
    
    // renormalize basis on new support (for working with flux fractions)
    if (Nb==1) bloc/=bloc(avg,,)(-:1:dimsof(wave0)(2),,);
    if (Nb=="smart") bloc/=bloc(avg);

    // optional: build a continuum-subtracted basis
    bloc=(is_void(nor)?bloc:_wnormb(bloc,nor,bc));

    flux=flux(wmin:wmax);

    // use the sigma given in the data file if sigdef=1 (otherwise use the SNR estimate to build variance-covariance matrix sigdef=1 is default)
    if(sigdef==1)  {
      sigm=sigm(wmin:wmax);
      //sigma(where(sigma<=1.e-5))=1.e10;
    };

    // normalize data and sigma ? by default yes
    if (N==1) {
      sigm/=flux(avg);
      flux/=flux(avg);
    };

    //    if(sigdef==0) sigm=1./(u.SNR/flux);
    
    _x0=wave0;x0=_x0;dd=numberof(_x0);nr=dd;
    _m=base.met;
    d=flux;      // DATA
    // build continuum-subtracted data ? see wnorms
    d=(is_void(nor)?d:wnorms(d,nor,dc));

    if(1){ // this is to be used with log10 and not log neper
    // Generate losvd support
    deltal=log10(x0(0)/x0(1));
    //xp01=span(-deltal/2.,deltal/2.,nr); // used to be so as of 10/11/2008
    xp01=log10(x0)-log10(x0)(avg); // changed 10/11/2008 PO // works sort of the same as line above // looks ok
    ul1=log10(1.+v1/c);
    ul2=log10(1.+v2/c);
    ni=(where(xp01<ul1))(0);
    nj=(where(xp01>ul2))(1);
    rv1=c*(10^(xp01(ni))-1.);
    rv2=c*(10^(xp01(nj))-1.);
    riv=c*(10^(xp01(ni:nj))-1.);
    nlos=numberof(riv); // size of losvd in pixels
    };

    if(0){ // this is to be used with log10 and not log neper
    // Generate losvd support
    deltal=log(x0(0)/x0(1));
    //xp01=span(-deltal/2.,deltal/2.,nr); // used to be so as of 10/11/2008
    xp01=log(x0)-log(x0)(avg); // changed 10/11/2008 PO // works sort of the same as line above // looks ok
    ul1=log(1.+v1/c);
    ul2=log(1.+v2/c);
    ni=(where(xp01<ul1))(0);
    nj=(where(xp01>ul2))(1);
    rv1=c*(exp(xp01(ni))-1.);
    rv2=c*(exp(xp01(nj))-1.);
    riv=c*(exp(xp01(ni:nj))-1.);
    nlos=numberof(riv); // size of losvd in pixels
    };





    
    // automatic padding
    if(padv=="auto"){
      targ=fft_best_dim(nr+2*nlos);
      pad1=(targ-nr)/2;
      //pad2=(targ-nr+1)/2;
      pad2=targ-nr-pad1;
    };
    
    bdata=roll(pad(d,pad1,pad2))(::-1);
    if(kin==2) bdata=roll(pad(d,pad1,pad2));
    //data=roll(bdata)(pad2+2:nr+pad2+1)(::-1);


    // normalize the mus according to the number of pixels of the losvd, the SAD, the AMR and age-dependent extinction.
    deconv_muv=muv/nlos;
    deconv_mux=mux/nab;
    deconv_muz=muz/nab;
    deconv_mue=mue/nab;

    
    
    // PLOT the resampled basis and data for a check
    if (!mp_rank){
      if(pl!=0){
        ws,0;
        plb,bloc(,,1),wave0;
        plh,flux,wave0,color="blue";
      };
    };

    // ===============================================================
    // ***************************** Build W **********************
    // ================================================================

    // CHECK FOR ZERO ELEMENTS IN sigm
    // and give them huge sigma
    Mw=max(x0(dif));
    if(!is_void(dimsof(where(sigm==0)))) sigm(where(sigm==0))=1.e10;
    W=1./sigm^2;
    if((!is_void(__sigm))&(__sigm==1)) W=sigm*0.+numberof(W);

    
    zw1=0.;  // zeroweight
    zw2=0.;  // zeroweight

    //***************************************************************
    //================ Apply User-provided mask at rest =============
    //***************************************************************
    
    if (!is_void(RMASK)){
      mwave1=[];
      nRMASK=dimsof(RMASK)(3);
      for(i=1;i<=nRMASK;i++){grow,mwave1,intersectS(RMASK(,i),wave0);};
      for(i=1;i<=numberof(mwave1);i++){
        if (!is_void(dimsof(where(abs(wave0-mwave1(i))<=Mw)))) W(where(abs(wave0-mwave1(i))<=Mw))=zw2;};
    };
    
    //******************* ALWAYS MASK EXTREMITIES *******************
    W(:nedgemask)=zw1;
    W(dd-nedgemask:)=zw1;

    //***************************************************************
    //===================== instrument mask  ========================
    //================ Useful for dead pixels and such ============== 
    //***************************************************************
    
    if (!is_void(MASK)&(!(MASK=="nope"))){
      mwave/=(1.+u.redshift);
      
      // ============================================================
      // ******** USER-provided MASK (wal or SDSS) OBSOLETE *********
      // ============================================================
      
      for(i=1;i<=numberof(mwave);i++){
        if (!is_void(dimsof(where(abs(wave0-mwave(i))<=Mw)))) W(where(abs(wave0-mwave(i))<=Mw))=zw2;};
    };

    // Divide by number of dof
    W=W/max(1.,(numberof(where(W>=1.e-5))));
    // now W is like (1./(dd*sigma^2));

    // Padding the weight vector W
    pW=pad(W,pad1,pad2);
    pW=pW(::-1); //  WARNING: ALWAYS check consistency with functions in Qfunctions.i
      
    if (kin==0) {pW=W;bdata=d;};

    //*******************************************************
    // ===================== parametric fit =================
    // ===================== age & metallicity ====== =======
    //*******************************************************
    if(parage==1){
      ta=log10(base.ages);
      nb=numberof(ta);
      muz1=0;
      muab1=0;
      ws,0;
      if(kin==0) q=fnpecpara(d,RG=RG,s=1.e-5,meval=meval);
      if(kin==1) {v=sfit(u,base,noskip=1,meval=500,kin=1,epar=3,MASK=MASK,RMASK=RMASK,sav=0);
      
      //      if(kin==1) q=fnpecparak(bdata,RG=RG,s=1.e-5,meval=meval,w0=0.5);
      los=v(-nlos+1:,1);
      q=fnpecpara2(bdata,RG=RG,s=1.e-5,meval=meval);
      };
      
      if (sav) {
        if(kin==0) {ki2=npecpara(q); los="none";losvd="none";};
        if(kin==1) {ki2=npecparak2(q);losvd=los;los=los^2;};
        LWAge=10^q(1);
        LWMet=Zrescalem1(q(2));
        ade=q(3:nde+2);
        
        if(epar==3){    // FIT POLYNOMIAL TO NPEC TO GET an average E(B-V)
          W1=W*0.+1.;
          W1(where(W==0))=0;
          W1=fft_smooth(W1,100)^4;
          //INFO,x0;

          pmodel1=model;  // TWEAKED 17/06/06
          npec=npe(q(3:3+nde-1),numberof(x0))/npe(q(3:3+nde-1),numberof(x0))(avg)*pmodel1(avg);
          ebv=[0.1];truc=lmfit(fds,x0,ebv,npec,W1,deriv=1);  // mean E(B-V)
          sigmaebvgas=((ds(ebv,x0)-npec)*W1)(rms);   // flux calibration error
          ebvgas=ebv;
          ebvstar=0.44*ebvgas;
        };
                
        basisname=base.filename;
        R=base.R;
        RG=is_void(RG)?"void":RG;
        _ki=ki2;
        info,q;
        g=createb(u.resfile(1));
        dfile=u.filename(1);
        save,g,d,model,basisname,dfile,W,x0,R,_m,r,ta,RG,parage,q,ki2,_ki,nde,kin,tam,ftam,los,losvd,LWAge,LWMet,ebvgas,ebvstar,sigmaebvgas,npec,ade,bdata,muAMRp,rrAMRp,AMRp;
       
        close,g;
        write,mp_rank,"saved",u.resfile(1);
      };
      return q;
      
    };

    
    
    //=============================================================
    //************** NON-parametric inverse method ****************
    //=============================================================
    
    // CREATE DEFAULT GUESS (and solution format)
    
    sol=0.01*array(sqrt(1./double(nab)),nab);  // SAD
    //z=array(0.15,nab); // AMR
    z=array(b.met(avg),nab); // AMR // changed 10 jul 2008 PO
    ebv=(epar==1)?0.5:array(0.5,nab); // extinction: dust screen or age-dependent ?
    if(epar==0) ebv=[1.,1.e0,-1.e-5]; // polynomial
    if(epar==3) ebv=array(1.,nde); // non-parametric extinction law
    los=(kin==0)?[]:array(0.1,nlos); // LOSVD
    grow,sol,z,ebv,los;  // put it all together

    if(kin==2) {
      // NOTE THAT THESE 2 GUESSES LEAD TO SGNIFICANTLY DIFFERENT SOLUTIONS
      sol=array(0.1,nlos*nab+nde);   // flat guess
      sol=makebump(nlos,int(nlos/2.),int(nlos/10.),N=1)(,-:1:nab);  //gaussian guess, flat in age
      sol=sol(*);
      grow,sol,array(1.,nde);
      if(numberof(dimsof(z2d))==1){  // changed PO 14/11/2008
          write,"using metallicity number ",z2d," i.e. log(Z/Zsun)=",log10(Zrescalem1(b.met)(z2d)/Zsun);
          img=bpad(bloc(,,z2d),pad1,pad2);  
      };

      if(numberof(dimsof(z2d))==2){
          write,"using user-provided AMR";
          img=bpad(buildb2(z2d),pad1,pad2);
      };

      mtf = fft(img,[-1,0]);
    };
    
    sol=is_void(guess)?sol:guess;  // USE user-provided GUESS OR NOT

    //********************************************************************
    //==================    Build penalizations   ========================
    //********************************************************************
    
    rr=genREGUL(nab);
    rr= rr(+,)*rr(+,);  // first derivative 
    rr=rr(,+)*rr(+,);  // second derivative regularization

    // Build L1 smoothing kernel for SAD
    if(!is_void(L1)) {rr=genREGUL(nab,s=L1);rr=rr(+,)*rr(+,);}; 
    
    // the same for the losvd
    rr1=genREGUL(nlos);
    rr1= rr1(+,)*rr1(+,);  // first derivative  sometimes better (??)
    rr1=rr1(,+)*rr1(+,);  // second derivative penalization

    // Build L2 smoothing kernel for LOSVD
    if(!is_void(L2)) {rr1=genREGUL(nlos,s=L2);rr1=rr1(+,)*rr1(+,);};

    // Build L3 smoothing kernel for AMR
    rr3=rr;
    if(!is_void(L3)) {rr3=genREGUL(nab,s=L3);rr3=rr3(+,)*rr3(+,);};

    // Build L4 smoothing kernel for the Non parametric extinction law
    rr4=genREGUL(nde,s=L4);rr4=rr4(+,)*rr4(+,);


    //***********************************************************************
    //*****************  START MINIMISATION STEPS ***************************
    //***********************************************************************
    mub=mub1;
    q=sol;
    for(it=1;it<=maxit;it++){
    q=optim_driver(Q0,q,verb=verb,maxeval=meval,ndirs=40,frtol=frtol1,farg=mpar);
    };

    mub=mub2;
    for(it=1;it<=maxit;it++){
    q=optim_driver(Q0,q,verb=verb,maxeval=meval,ndirs=40,frtol=frtol1,farg=mpar);
    };
    q0=q;

    for(it=1;it<=maxit;it++){
    q=optim_driver(Q1,q,verb=verb,maxeval=meval3,ndirs=40,frtol=frtol1,farg=mpar);
    };
    q1=q;


    if(kin!=2){
    sad=q1(:nab);
    amr=q1(nab+1:2*nab);
    ade=q1(2*nab+1:2*nab+numberof(ebv));
    ade=max(ade,0.); // BECAUSE q1 might not be positive
    losvd=(kin==0)?[]:q1(2*nab+numberof(ebv)+1:2*nab+numberof(ebv)+numberof(los));
    q1=sad;
    grow,q1,amr,ade,losvd;         
        
    res0=Q0(q0,g,mpar);
    model0=model;
    
    //pmodel0=(kin==0)?model0:roll(model0)(pad1+1:pad2+dd)(::-1);
    pmodel0=(kin==0)?model0:(mod(nr,2)==1)?roll(model0)(pad2+1:pad1+dd+1)(::-1):roll(model0)(pad1+1:pad2+dd)(::-1);
    };
    
    res1=Q1(q1,g,mpar);
    if(kin==2) {
      model=model(::-1);
      q0=q1;
      ade=q1(nab*nlos+1:nab*nlos+nde);
      ade=max(ade,0.);
      res0=res1;
      pmodel0=model;
//      pmodel1=model;
    };

    model1=model;

    //pmodel1=(kin==0)?model1:roll(model1)(pad1+1:pad2+dd)(::-1);
    pmodel1=(kin==0)?model1:(mod(nr,2)==1)?roll(model1)(pad2+1:pad1+dd+1)(::-1):roll(model1)(pad1+1:pad2+dd)(::-1);

    npec=0;
      if((epar==3)&(kin!=2)) {
        npec=npe(q1(2*nab+1:2*nab+nde),numberof(x0))/npe(q1(2*nab+1:2*nab+nde),numberof(x0))(avg)*pmodel1(avg);
      };

      if((epar==3)&kin==2) {
        npec=npe(q1(nab*nlos+1:nab*nlos+nde),numberof(x0))/npe(q1(nab*nlos+1:nab*nlos+nde),numberof(x0))(avg)*pmodel1(avg);
      };
        
      spmodel1=spline(pmodel1,x0*(1.+u.redshift),swave);
      //redshift spmodel1 to the galaxy's redshift
      nspmodel1=spmodel1*sflux(avg);
            
      ssd=spline(d,x0,swave);
      sW=spline(W,x0*(1.+u.redshift),swave); //redshift W as well
      nssigm=ssigm;
      nssigm(where(sW<=2.e-2))=1.e10;
      
      //************************************************************************
      //=============================== PLOTS ==================================
      //************************************************************************

      
      
      if(pl!=0){
        if (!mp_rank){
          
          ws,1;
          //plh,ssigm/ssd(avg),swave,color="blue";
          //plh,ssd/ssd(avg),swave,color="black";
          //plh,d,x0,color="black",type=2;
          //plh,spmodel1,swave,color="red";

          plh,sflux,swave;
          plh,nspmodel1,swave,color="red";
          plh,ssigm,swave,color="blue";
          pltitle,"intial data sampling";

          //window,2,style=(lsp==1?sdir+"Gist/large.gs":sdir+"Gist/bboxed.gs");
          if (lsp==1) WSL,2;
          ws,2;
          plh,(plavgW==1?W/W(avg):W),x0,color="blue";
          plh,d,x0,color="black";
          //plh,q1;
          //plh,sol,color="green";
          plh,pmodel1,x0,color="red";
          if(epar==3) plh,npec,x0,color="green";
          pltitle,"model sampling";
          pause,1;
        

          if(kin!=2){
          ws,3;
          window,3,style=sdir+"Gist/win31.gs";
          plsys,3;
          plh,sad^2,log10(base.ages)+6,width=3;
          plh,sol(:nab)^2,log10(base.ages)+6,width=3,color="green";
          //xytitle;
          plt1,"SAD",0.18,0.83,tosys=0;
          plt1,"log10(age[yr])",0.32,0.89,tosys=0;
          plt1,"flux fractions",0.66,0.84,tosys=0,orient=3;
          plsys,2;
          plh,log10(Zrescalem1(amr)/0.02),log10(base.ages)+6,width=3;
          plh,log10(Zrescalem1(sol(nab+1:2*nab))/0.02),log10(base.ages)+6,width=3,color="green";
          plt1,"AMR",0.18,0.68,tosys=0;
          plt1,"log10(Z/Zsun)",0.66,0.7,tosys=0,orient=3;
          
          plsys,1;
          if(kin==1) {
            plh,losvd^2,riv,width=3;
            plh,sol(2*nab+nde+1:2*nab+nde+nlos)^2,riv,width=3,color="green";
            plt1,"LOSVD",0.18,0.53,tosys=0;
            xyleg,"v[km/s]","";};
          
          //plh,W,x0,color="blue";
          //plh,d,x0,color="black";
          //plh,q0;
          //plh,sol,color="green";
          //plh,pmodel0,x0,color="red";
          //if(epar==3) plh,npec,x0,color="green";
          //pltitle,"no clip solution";
          pltitle,u.name;
          pause,1;
          };
          
          if(kin==2){
            ws,3;
            //            window,3,style=usdir+"/Yorick/Gist/bboxed.gs";
            // tweaked that cause of missing gist files
            window,3,style="boxed2.gs";
            avd=(reform(q1(:nlos*nab),nlos,nab))^2; // reform age-velocity distribution
            plk,avd,log10(b.ages)+6,riv;
            pltitle,"age-velocity distribution";
            xyleg,"v[km/s]","log(age[yr])";
          };
        };
      };
      
      if(is_void(gres)) gres=array(0.,numberof(q),ns*(((is_void(nMC)?0:nMC) +1))); //array for storing the results
      // STORE chi2 and residuals
      gres(,(1+(j-1)*(is_void(nMC)?0:(nMC)+1)))=q;
      _resi((1+(j-1)*(is_void(nMC)?0:(nMC)+1)))=res1;
      _ki((1+(j-1)*(is_void(nMC)?0:(nMC)+1)))=chi2r(pmodel1,d,1.,w=W);
      
      //(((bdata-model)^2)*((kin==0?(pW):roll(PW))))(sum);
      
      

      //===================================================================
      //*********** MONTE CARLO exploration of the likelihood surface *****
      //===================================================================
      sd=d;  //  store original d for later
      // Create automatic rmin for random guess exploration
      rmin=is_real(rmin)?rmin:str2float(rmin)*res1;
      
      if (!is_void(nMC)){
        
      ki2=array(0.,nMC);

      //*******************************************************
      // START MC LOOP
      for(i=1;i<=nMC;i++){
        
        write,mp_rank,MC,"ing",u.filename(1),"step",i,"of",nMC;
        
        if (MC=="MC") { // REGULAR MC
          // 9 july 2008 -- MM & PO 
          // Two options, choose either the best model or the data to 
          // estimate the errors.
          if(dMC=="bestmodel") {
              d=pmodel1; // choose the best the model to start the MC simulations
              write,"MCing using best model, expect chi2~1";
          };

          if(dMC=="data")  {
              d=sd;  // choose the data to start the MC simulatins; dont forget to 
          // in that case what should be the chi2 1 or 2 ?
          // keep in mind that the sum of 2 gaussians (mu,sigma) is a gaussian (2mu, sqrt(2)*sigma)
          // so chi2 should be 2
              write,"MCing using noised data, expect chi2~2";
          };
          
          write,"no user-provided sigma, using global SNR to noise data";
          snr=u.SNR;
          // using the error spectrum to re-sample the spectrum -- MM & PO
          // d=d+sigm*random_normal(numberof(d)); // nice, provided it has been properly provided by user
          d*=1.+(1./snr)*random_normal(numberof(d)); // otherwise just use constant snr approximation

          if((kin==2)|(kin==1)) bdata=roll(pad(d,pad1,pad2))(::-1);
          if(kin==2) bdata=bdata(::-1) // added 15/01/2008 this has to do with a slightly different treatment of fft in kin=1 and kin=2 // I dunno really why tho
          q=sol;
        };
        
        if (MC=="RG") { // RANDOM GUESS tests secondary extrema
          if(bf=="auto") bf=avg(sqrt(sad^2));
          q=abs(bf*random_normal(numberof(q)),0.);
          //q=sol;
        };
        
        // START MINIMISATION SEQUENCE
        mub=mub1;
        for(it=1;it<=maxitMC;it++){
          q=optim_driver(Q0,q,verb=100,maxeval=mevalMC,ndirs=40,frtol=frtol1,farg=mpar,fmin=rmin);
        };
        
        mub=mub2;
        for(it=1;it<=maxitMC;it++){
          q=optim_driver(Q0,q,verb=100,maxeval=mevalMC,ndirs=40,frtol=frtol1,farg=mpar,fmin=rmin);
        };
        
        for(it=1;it<=maxitMC;it++){
          q=optim_driver(Q1,q,verb=100,maxeval=mevalMC,ndirs=40,frtol=frtol1,farg=mpar,fmin=rmin);
        };
        
        gres(:numberof(q),i+(j-1)*(nMC+1)+1)=q;
        res1=Q1(q,g,mpar);
        _resi(i+(j-1)*(nMC+1)+1)=res1;
        _ki(i+(j-1)*(nMC+1)+1)=((((bdata-model)^2)/(max(model^2,1.e-5)))*roll(pW))(sum)*(((u.SNR)^2)/sum(pW));

      };
      
    };
      // END MC loop
      //****************************************************


      
      //*************************************************************
      //========== COMPUTE INTEGRATED PHYSICAL QUANTITIES ===========
      //*************************************************************

      if(kin==2){
          sad=avd(sum,);sad/=sad(sum);sad=sqrt(sad);
          if(numberof(dimsof(z2d))==1){amr=array(b.met(z2d),nab);};
          if(numberof(dimsof(z2d))==2){amr=z2d;};
      };
          
      
      LWAge=LWA(sad^2,1,nab,b=base); // luminosity weighted age
      sadamr=sad^2;grow,sadamr,amr; 
      //LWMet=log10(Zrescalem1(LWM(sadamr,1,nab))/Zsun); // luminosity weighted metallicity
      LWMet=Zrescalem1(LWM(sadamr,1,nab)); // luminosity weighted metallicity
      
      if(epar==1) {ebvdummy=q1(2*nab+1); sigmaebvgas="none";};
      if(epar==3){    // FIT POLYNOMIAL TO NPEC TO GET avg E(B-V)
      W1=W*0.+1.;
      W1(where(W==0))=0;
      W1=fft_smooth(W1,100)^4;
//      INFO,x0;
//      write,"determining E(B-V)";
      ebvdummy=[0.1];
      INFO,ds(ebvdummy,x0);
      truc=lmfit(fds,x0,ebvdummy,npec,W1,deriv=1);  // mean E(B-V)
      sigmaebvgas=((ds(ebvdummy,x0)-npec)*W1)(rms);   // flux calibration error
      };

      ebvgas=ebvdummy;
      ebvstar=0.44*ebvdummy;  // if using calzetti et al.

      write, "ebv_star = ", ebvstar;
        
      if(kin==2) {losmoments=[];for(im=1;im<=nab;im++){grow,losmoments,moments((avd(,im))^2,riv);};losmoments=reform(losmoments,4,nab);};
      if(kin==1) losmoments=moments(losvd^2,riv);  // moments of losvd
      if(kin==0) losmoments="none";
      
      //====================================================================
      //********* Integrated quantities for the MC solutions ***************
      //====================================================================
      if(!is_void(nMC)){
      sq1=q1; // save q1 to give it back its value at the end of this loop
      ssad=sad;
      samr=amr;
      sade=ade;
      
      
      if(kin==1) slosvd=losvd;
      if(epar==3) snpec=npec;
      
        for(i=1;i<=nMC;i++){
          q1=gres(,i+1);
          sad=q1(:nab);
          amr=q1(nab+1:2*nab);
          ade=q1(2*nab+1:2*nab+numberof(ebv));
          ade=max(ade,0.); // BECAUSE q1 might not be positive
          losvd=(kin==0)?[]:q1(2*nab+numberof(ade)+1:2*nab+numberof(ade)+numberof(riv));
          if(epar==3) {
            npec=npe(q1(2*nab+1:2*nab+nde),numberof(x0))/npe(q1(2*nab+1:2*nab+nde),numberof(x0))(avg)*pmodel1(avg);
          };
          
          grow,LWAge,LWA(sad^2,1,nab,b=base); // luminosity weighted age
          sadamr=sad^2;grow,sadamr,amr; 
          //LWMet=log10(Zrescalem1(LWM(sadamr,1,nab))/Zsun); // luminosity weighted metallicity
          grow,LWMet,Zrescalem1(LWM(sadamr,1,nab)); // luminosity weighted metallicity
          
          if(epar==1) {grow,ebvdummy,q1(2*nab+1); grow,sigmaebvgas,"none";};
          if(epar==3){    // FIT POLYNOMIAL TO NPEC TO GET avg E(B-V)
            W1=W*0.+1.;
            W1(where(W==0))=0;
            W1=fft_smooth(W1,100)^4;
            //INFO,x0;
            write,"determining E(B-V)";
            ebvdummy=[0.1];truc=lmfit(fds,x0,ebvdummy,npec,W1,deriv=1);  // mean E(B-V)
            grow,sigmaebvgas,((ds(ebvdummy,x0)-npec)*W1)(rms);   // flux calibration error
          };
          
          grow,ebvgas,ebvdummy;
          grow,ebvstar,0.44*ebvdummy;  // if using calzetti et al.
        
          if(kin==1) grow,losmoments,moments(losvd^2,riv);  // moments of losvd
          if(kin==0) grow,losmoments,"none";  
          
        };
        
        //***** RE-allocate original values ***************
      q1=sq1; // get back saved value for q1;
      sad=ssad;
      amr=samr;
      ade=sade;
      
      if(kin==1) losvd=slosvd;
      if(epar==3) npec=snpec;
      //************* OK ********************************
      };

      //**************************************************************
      //************ COMPUTE SFRS, MASSES, SAVE AS TXT FILES *********
      //**************************************************************

      // SAVE SAD AS TXT FILE

      sads=sfrs=amrs=masses=gres(:nab,)*0.;
      
        f=open(u.resfile(1)+"-SAD.txt","w");
        write,f,"Ages (Myr)",numberof(base.ages);
        write,f,base.ages;
        write,f,"SAD (flux fractions)",numberof(sad);
        write,f,sad^2;
        sads(,1)=sad^2;
        if(!is_void(nMC)){
          for(iMC=1;iMC<=nMC;iMC++){
            write,f,"MC ",iMC," of ",nMC;
            sads(,iMC)=(gres(:nab,iMC+1))^2;
            write,f,(gres(:nab,iMC+1))^2;
          };
        };

        // SAVE AMR AS TXT FILE

        f=open(u.resfile(1)+"-AMR.txt","w");
        write,f,"Ages (Myr)",numberof(base.ages);
        write,f,base.ages;
        write,f,"AMR (0.02 is solar)",numberof(sad);
        write,f,Zrescalem1(amr);
        amrs(,1)=Zrescalem1(amr);
        if(!is_void(nMC)){
          for(iMC=1;iMC<=nMC;iMC++){
            write,f,"MC ",iMC," of ",nMC;
            write,f,Zrescalem1((gres(nab+1:2*nab,iMC+1)));
            amrs(,iMC)=Zrescalem1((gres(nab+1:2*nab,iMC+1)));
          };
        };

        // SAVE MASSES AS TEXT FILE
        MsL=base.MsLratio;
        f=open(u.resfile(1)+"-MASS.txt","w");
        write,f,"Ages (Myr)",numberof(base.ages);
        write,f,base.ages;
        write,f,"Masses in each time bin",numberof(sad);
        write,f,SAD2MASS(sad^2,amr);
        masses(,1)=SAD2MASS(sad^2,amr);
        if(!is_void(nMC)){
          for(iMC=1;iMC<=nMC;iMC++){
            write,f,"MC ",iMC," of ",nMC;
            write,f,SAD2MASS((gres(:nab,iMC+1))^2,gres(nab+1:2*nab));
            masses(,iMC)=SAD2MASS((gres(:nab,iMC+1))^2,gres(nab+1:2*nab));
          };
        };
        
         // SAVE SFR AS TEXT FILE
        MsL=base.MsLratio;
        f=open(u.resfile(1)+"-SFR.txt","w");
        write,f,"Ages (Myr)",numberof(base.ages);
        write,f,base.ages;
        write,f,"SFR (Unnormalized Msol/yr)",numberof(sad);
        write,f,SAD2SFR(sad^2,amr,10^base.bab,N=1);
        sfrs(,1)=SAD2SFR(sad^2,amr,10^base.bab,N=1);
        if(!is_void(nMC)){
          for(iMC=1;iMC<=nMC;iMC++){
            write,f,"MC ",iMC," of ",nMC;
            write,f,SAD2SFR((gres(:nab,iMC+1))^2,gres(nab+1:2*nab),10^base.bab,N=1);
            sfrs(,iMC)=SAD2SFR((gres(:nab,iMC+1))^2,gres(nab+1:2*nab),10^base.bab,N=1);
          };
        };

              
          
      //====================================================
      //******************* SAVE results *******************
      //====================================================
      
      if (sav) {
        g=createb(u.resfile(1));
        dfile=u.filename(1);
        
        if (0) {  // FOR DEBUG ONLY
          write,mp_rank,"saving",u.resfile(1);
          info,g;
          info,d;
          info,pmodel0;
          info,q0;
          info,res0;
          info,pmodel1;
          info,res1;
          info,dfile;
          info,W;
          info,x0;
          info,riv;
          info,ab;
          info,_m;
          info,mux;
          info,muv;
          info,muz;
          info,mue;
          info,mub1;
          info,mub2;
        };
        d=sd; // reattribute d its correct (non-MC) value
        basisname=base.filename;
        ab=base.nages;
        ages=base.ages;
        R=base.R;
        if (is_void(losvd)) losvd="none";
        if (is_void(avd)) avd="none";

        save,g,d,sigm,pmodel0,q0,res0,pmodel1,q1,sad,amr,ade,losvd,res1,dfile,W,pW,x0,riv,ab,ages,R,_m,mus,mux,muv,mue,muz,mub1,mub2,gres,meval,zlim,MC,_ki,MASK,rmin,_resi,basisname,kin,epar,nde,maxitMC,mevalMC,bf,muc,co,cov,mucov,wave,npec,wave0,sflux,swave,ssigm,spmodel1,nspmodel1,wmin,wmax,sW,nssigm,LWAge,LWMet,ebvgas,ebvstar,losmoments,avd,sads,masses,amrs,sfrs;
        close,g;
        write,mp_rank,"saved",u.resfile(1);
        if (asc==1) plouc=pdb2asc(u.resfile(1));

      };
      
      
  };
  return gres; 
  
};





// for plot puropses
#if 0 
hop=convertVAKU(gall,noplot=1);
for(i=1;i<=6;i++){ws;upload,hop.resfile(i),s=1;plh,q1,width=3;pltitle,hop.name(i);hcp_file,"/home4/ocvirk/perso/modeles/galaxie/POP/VAKU/VAKU221003-3-"+pr1(i)+".ps";hcp;hcp_finish;};
#endif


// for debugging purposes
if(0){
a=convertVAKU(fV);
//b=bRbasis(basisfile="PHR",wavel=[4000.,7000.],nbins=20);
b=bRbasis(nbins=20);
// b=bRbasis([1.e8,2.e10],nbins=20);
   
 v=sfit(a,b,kin=2,nde=40,epar=3,parage=[],sav=1,RG=[],noskip=1,verb=10,meval=300,RMASK=[[5100.,5200.]],z2d=1,bnpec=[],muv=1.e9,mux=1.e7);
 // v=sfit(a,b,kin=1,nde=40,epar=3,parage=[],sav=1,RG=[],noskip=1,verb=10,meval=300,RMASK=[],z2d=1,bnpec=[],muv=1.e5,mux=[]);
};

//mx=[3.,0.1];grow,mx,de,los;npecparak(mx);
if(0){
q=[3.,0.1];grow,q,array(1.,nde),makebump(nlos,int(nlos/2),3.,N=1);
 npecparak(q);
 bdata=model;
 nq=[2.,0.15];grow,nq,array(1.,nde),makebump(nlos,int(nlos/2),1.,N=1);
 npecparak(q);
 npecparak(nq);
 u1=optim_driver(npecparak,nq,verb=100,maxeval=500,frtol=1.e-20);
 u1q=u1;
 uq=q;
 unq=nq;
 };

if(0){
// derivative checking
  truc=checkder(AKNPEC,nab*nlos+nde,sol=sol,Ds=1.e-8,i1=30,i2=50);
};
