func getki2(file,change=,s=){
/* DOCUMENT 
     returns ki2 of fitted spectra
     WARNING : restores a score of variables
   SEE ALSO:
 */

  extern gres,pmodel0,d;
  
  pmodel0=[];
  gres=[];
  
  if (numberof(file)>1)
    {
      tt=array(0.,numberof(file),2);
      for(i=1;i<=numberof(file);i++) tt(i,)=is_void(change)?getki2(file(i)(1),s=1):getki2(strtok(file(i),".")(1)+".res1",s=1);
      return tt;
    };

  // TEST EXISTENCE OF RESFILE

  resfile=is_void(change)?file:strtok(file,".")(1)+".res1";
  craps="ls "+resfile(1);
  crap=exec(craps(1));
  exist=(crap(1)==resfile(1)?1:0);

  
  if (exist){
    if (is_void(s)) upload,(resfile);
    if (!is_void(s)) upload,(resfile),s=1;
  };
  if (is_void(pmodel0)) {res0=-1.;res1=-1.;};
  return [res0,res1];
};


func isgres(file){
  /* DOCUMENT 
     returns numberof MCs for datafiles with existing resfiles, 0 if resfile constains no MC (no gres) doesnt exist
     WARNINGWARNINGWARNING : restores a score of variables
     WHICH ACTUALLY MASK CRITICAL VARIABLES FOR THE FIT
     DONT USE IT
     USE isgres2
   SEE ALSO:
 */
  
if (numberof(file)>1)
    {
      tt=array(0.,numberof(file),3);
      for(i=1;i<=numberof(file);i++) tt(i,)=isgres(file(i));
      return tt;
    };

 gres=[];

 u=getki2(file,s=1);
 if(u(1)){return is_void(gres)?[0.,0.,0.]:dimsof(gres);};  // resfile exists, but does not/does contain MC

 return [-1.,-1.,-1.];  // resfile does not exist
 
};


func isgres2(file){
  /* DOCUMENT 
     returns 0 if resfile constains no MC (no gres)
     returns -1 if resfile doesnt exist
     returns 1 if resfile contains MC whatever the number of MC
     GOOD THING: doesnt restore any variable
     less dangerous than isgres
     doesnt use getki2
   SEE ALSO:
 */
  
  if (numberof(file)>1)
    {
      tt=array(0.,numberof(file),3);
      for(i=1;i<=numberof(file);i++) tt(i,)=isgres2(file(i));
      return tt;
    };

// TEST existence 
  craps="ls "+file(1);
  crap=exec(craps(1));
  if (!is_void(crap)) exist=(crap(1)==file(1)?1:0);
  if (is_void(crap)) exist=0;
  if (!exist) return 0;

 
  u=openb(file);
  pvn=get_vars(u);
  vn=(*pvn(1));

  return !is_void(dimsof(where(vn=="gres")));
};


func HugRes(file,&hx0,n=){
  // returns results from hughes ascii files.
  // first row 2IDs and one metallicity, then mass fractions, and then sfrs
  // IDs are PLATEID and FIBID 
  //  PlateIDArr(ip1),ind0(iadd)+1 :          367         543
  //TypeArr(ip1,ind0(iadd)),Prob : GALAXY                    0
  //SN_MEDIAN,D4,ZFITS,ValueMass      38.7819      1.50199    0.0291693
  //    9.89120
  //ibnpMet :        0
  //f1,f2,f3,f4,f5,f6 :      0.606997      0.00000      0.00000     0.367277
  //    0.00000    0.0257266
  //sfr1,sfr2,sfr3,sfr4,sfr5,sfr6 :      0.300490      0.00000      0.00000
  //    1.00000      0.00000     0.700468


  //   Ce sont :

  //ibnpMet=best fit Z (0=solaire, 1=sous-solaire (40 %), 2=supra solaire (1.5 
  //Zsol)) 

  //f1,f2 ... f6 = fractions de la totale masse stellaire formee sur les 6 
  //populations que je resous (13.5-8, 8-4, 4-1.5, 1.5-0.5, 0.5-0.1, 0.1-1 
  //Gyr en loockback time) 

  //sfr1 ... sfr6 =taux de formation stellaires correspondants (normalisation 
  //abritraire dependant de chaque galaxie, de facon a ce que max(sfr_i)=1 
  //pour chaque galaxie)  

  
  local f,u,v,prob,iprob,res,hx0;

  hx0=[11000.,6000.,2500.,1000.,250.,50.];  // Hughes age support (from above)
  
  if(is_void(n)) n=10000;
  f=open(file,"r");
  u=rdline(f,n);
  close,f;
  v=strmatch(u,"TypeArr");
  indv=where(v==1);
  //iprob=str2int(strpart(u(indv),strlen(u(indv(1)))-5:strlen(u(indv(1)))));
  iprob=str2int(strpart(u(indv),0:0));
  ires=where(iprob==0);
  ni=numberof(ires);
  res=array(0.,6,3,ni);

  for(i=1;i<=ni;i++){
    res(1:2,1,i)=str2int((split2words(strtok((u(indv(ires(i))-1)),":")(2))));
    res(3,1,i)=str2float(strtok((u(indv(ires(i))+3)),":")(2));
    res(1:4,2,i)=str2float(split2words(strtok((u(indv(ires(i))+4)),":")(2)));
    res(5:6,2,i)=str2float(split2words(u(indv(ires(i))+5)));
    ns1=numberof(str2float((split2words(strtok((u(indv(ires(i))+6)),":")(2)))));
    res(1:ns1,3,i)=str2float((split2words(strtok((u(indv(ires(i))+6)),":")(2))));
    ns2=numberof(str2float(split2words(u(indv(ires(i))+7))));
    res(ns1+1:ns1+ns2,3,i)=str2float(split2words(u(indv(ires(i))+7)));
  };

  return res;

  // to plot and compare hughes results to mine, after selection of thecommon galaxies.
  #if 0
  // should build a basis b to have a bab before doing that and a correct MsL array
  b=bRbasis(nbins=10,zr=1,mets=[0.05,0.004]);
  nab=numberof(b.ages);
  _m=b.met;

  // in Hres(4,1,i) we want to store the galaxy rank in sll

  for(i=1;i<=numberof(sll);i++){
    _FIBID=[];
    _PLAID=[];
    u=convertSDSS(sll(i),noplot=1);grow,_FIBID,h.FIBERID;grow,_PLAID,h.PLATEID;
  };

  for(i=1;i<=dimsof(Hres)(4);i++){
    Hres(4,1,i)=where((_FIBID==Hres(2,1,i))&(_PLAID==Hres(1,1,i)))(1);
  };
        
  for(i=1;i<=51;i++){ws,1;u=convertSDSS(sll(int(Hres(4,1,i))),noplot=1);upload,u.resfile(1);SF=(q1(1:10)^2)*MsLinterp(q1(11:20))/bab(dif)*0.3e9;SF/=max(SF);plh,SF,(ta(ab)),width=3;plh,Hres(,3,i),(hx0),color="red",width=3;pltitle,u.filename(1);logxy,1,0;hcp_file,"hcomSF"+pr1(i)+".ps";hcp;hcp_finish;ws,0;plh,d,x0;plh,pmodel1,x0,color="red";plh,W,x0,color="blue";pltitle,u.filename(1);range,0.5,1.5;if (1){ window,1;hcp_file,"/raid/ocvirk/SDSS/Hughes/plots/hcomSP"+pr1(Hres(1,1,i))+"-"+pr1(Hres(2,1,i))+".ps";hcp;hcp_finish;};};

  #endif

};

func wait4me(nh){
  // prints a counter for nh hours 
  for(i=1;i<=nh;i++){
    pause,3600000;print,i;
  };
};
    


