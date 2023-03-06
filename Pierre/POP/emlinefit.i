#include "STECKMAP/Pierre/POP/sfit.i"
#include "Bastien/mathbast.i"

// comments: tend to get lots of SIGFPE error (division by 0 or the likes) so the only way to process the whole sample is to rerun it with firsguesssiamg=4 and 5.  oh well... eventually you get thenm alll...


func fitgaussian(line,wave,firstguess,&linefit,auto=){

  /* DOCUMENT
     tool for fitting emission lines with a gaussian profile
     to get best results one should select a domain centered on the line of interest to avoid confusion
     The firstguess is: [amplitude, mean, sigma]
     if auto==1 amplitude first guess is taken as line(avg)
  */

  sspar=array(0.,3);
  sspar(3)=1./(2.*firstguess(3)^2); // sigma
  sspar(1)=firstguess(2); // mean
  sspar(2)=firstguess(1); // amplitude

  if(auto==1) sspar(2)=line(avg);

  linefitpars=gausslmfit(line,wave,par=sspar);
  linefit=fgaussienne(wave,linefitpars);
  
  sigma=sqrt(1./(2.*abs(linefitpars(3))));
  mean=linefitpars(1);
  amplitude=linefitpars(2);

  return [amplitude,mean,sigma];
};

func fitandstoreemlines(file,sav=,skip=){
  /* DOCUMENT
     needs an existing .res file
     loads the corresponindg .res file, fits the emission lines, then saves d, pmodel1, x0 and the emission lines model in a new .em.pdb file including the emission line. This way
     skip=1 allows to skip the files already processed
     dont forget to play with sigma (firstguess(3)) if the code crashes
  */

  resemfile=strreplace(file,".pdb",".em.pdb");
  if ((exec("ls "+resemfile)(1)==resemfile)) {
    write,"skipping";
    return [];
  };
  resfile=strreplace(file,".pdb",".res");
  write,"opening "+resfile;
  
  upload,resfile,s=1;

nwidth=6; // half width of the wave support
emlinesnames=["OII","Halpha","NII","SII_1","SII_2","NII_6548.","Hbeta","Hdelta","Hgamma","OIII_4959","OIII_5007","NeIII_3869","HeI_5876","OI_6300"];
emlineswave=[3727.,6563.,6582.25,6716.,6731.,6547.,4861.,4101.,4340.4,4959.,5007.,3869.,5876.,6300.];
nem=numberof(emlineswave);
emlinemodel=x0*0.;

//i=13;
// fit OII 3727
for(i=1;i<=nem;i++){
nOII=(abs(x0-emlineswave(i)))(mnx);
wave=x0(max(nOII-nwidth,1):min(nOII+nwidth,numberof(x0)));
line=(d-pmodel1)(max(nOII-nwidth,1):min(nOII+nwidth,numberof(x0)));

firstguess=[line(avg),emlineswave(i),4.];
linefitpars=fitgaussian(line,wave,firstguess,linefit);

 emlinemodel(max(nOII-nwidth,1):min(nOII+nwidth,numberof(x0)))+=max(linefit,0.);

 };
 
 ws,0;
 plh,line,wave;
 plh,linefit,wave,color="red";
 
 ws,1;
 plh,d-pmodel1,x0;
 plh,emlinemodel,x0,color="red";

 ws,2;
 plh,W,x0,color="blue";
 plh,d,x0,color="red";
 plh,d-emlinemodel,x0;
 plh,emlinemodel,x0,color="green";
 
 ssigm=sigm;
 
 if(sav==1){
   upload,file;
   f=createb(resemfile);
   mask=interp(mask,wave,x0);
   wave=x0;
   sigm=ssigm;
   flux=d-emlinemodel;
   save,f,gal,flux,wave,sigm,mask,emlinemodel;
   close,f;
   write,"wrote ",resemfile;
 };
 return "done";
 };




if(1){
  ll=exec("ls /Users/pedro/work/Es_SFH/sdss_gals/fullsample/*-???.pdb"); // fancy regexp needed in order not to confuse the files em.pdb with the .pdb
  //  ill=37;
  //  ll(ill);
  nll=numberof(ll);
  for(ill=1;ill<=nll;ill++){
    ill;
    fitandstoreemlines(ll(ill),sav=1,skip=1);
  };
 };





if(0){

ll=exec("ls /Users/pedro/work/Es_SFH/sdss_gals/fullsample/*.res");
ill=10;// 9 is a good example with nice em lines
for(ill=1;ill<=numberof(ll);ill++){
  ill;
upload,ll(ill),s=1;
ws,1;
//plh,d,x0;
//plh,pmodel1,x0,color="red";
plg,d-pmodel1,x0;

ws,2;
plh,d,x0;
plh,pmodel1,x0,color="red";

nwidth=6; // half width of the wave support
emlinesnames=["OII","Halpha","NII","SII_1","SII_2","NII_6548.","Hbeta","Hdelta","Hgamma","OIII_4959","OIII_5007","NeIII_3869","HeI_5876","OI_6300"];
emlineswave=[3727.,6563.,6582.25,6716.,6731.,6547.,4861.,4101.,4340.4,4959.,5007.,3869.,5876.,6300.];
nem=numberof(emlineswave);
 emlinemodel=x0*0.;

//i=13;
// fit OII 3727
for(i=1;i<=nem;i++){
nOII=(abs(x0-emlineswave(i)))(mnx);
wave=x0(max(nOII-nwidth,1):min(nOII+nwidth,numberof(x0)));
line=(d-pmodel1)(max(nOII-nwidth,1):min(nOII+nwidth,numberof(x0)));

firstguess=[line(avg),emlineswave(i),3.];
linefitpars=fitgaussian(line,wave,firstguess,linefit);

 emlinemodel(max(nOII-nwidth,1):min(nOII+nwidth,numberof(x0)))+=max(linefit,0.);
 
ws,0;
plh,line,wave;
plh,linefit,wave,color="red";

 };

window,1;
plh,emlinemodel,x0,color="red";

};
 
 };

if(0){

n=100;

line=makebump(n,n/2,1.,N=[])*(1.+0.1*random_normal(n))+0.2*random_normal(n);
wave=span(4000.,5000.,n);

sigmaguess=100.; // it seems the 1st guess in sigma should be large. it leads to less failures than a small sigmaguess IF THE MEANGUESS IS GOOD AND THERE IS NO NOISE!!!
meanguess=wave(avg)+10.*random_normal(1);
write,meanguess;
ampguess=1.;
 firstguess=[ampguess,meanguess,sigmaguess];

 //linefitpars=gausslmfit(line,wave,par=sspar);
 linefitpars=fitgaussian(line,wave,firstguess,linefit);
   //linefit=fgaussienne(wave,linefitpars);

 sigma=linefitpars(3)
mean=linefitpars(2);
amplitude=linefitpars(1);

ws;
plg,line,wave;
plg,linefit,wave,color="red";

limits,4400.,4600.;
  
 };





if(0){

n=100;

line=makebump(n,n/2,1.,N=[])*(1.+0.1*random_normal(n))+0.2*random_normal(n);
wave=span(4000.,5000.,n);

sspar=[0.,0.7,2.e-5]*(1.+0.*random_normal(3));
sigmaguess=100.; // it seems the 1st guess in sigma should be large. it leads to less failures than a small sigmaguess IF THE MEANGUESS IS GOOD
meanguess=wave(avg)+10.*random_normal(1);
write,meanguess;
ampguess=1.;
sspar(3)=1./(2.*sigmaguess^2);
sspar(1)=meanguess;
sspar(2)=ampguess;


linefitpars=gausslmfit(line,wave,par=sspar);
linefit=fgaussienne(wave,linefitpars);

sigma=sqrt(1./(2.*abs(linefitpars(3))));
mean=linefitpars(1);
amplitude=linefitpars(2);

ws;
plg,line,wave;
plg,linefit,wave,color="red";

limits,4400.,4600.;
  
 };
