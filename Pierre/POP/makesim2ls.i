// TOOLS FOR GENERATING LARGE PSEUDO-DATASETS

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



func MakMok2(sq,v,kin,epar,snr,nR,patt,fseed,&pmodel1,N=,nabvec=,pl=,nde=){
  /* DOCUMENT 
     creates mock data from model parameters q, using basis v, with signal to noise ratio snr. A number nR of mock spectra is created, and recorded in generic files beginning with fseed.
     nab=number of age bins.
     Careful! the SFH part of q is squared before computing the corresponding spectrum is nab is not given
     if nab is given this is prevenmted
     this is avoided 
     N=1 switches to normalized outputs.
     *********   SET N=0 FOR SIMULATIONS !!! ***********
     default is N=1
     WARNING fseed should finish with / and should be an ABSOLUTE PATH otherwise crac boum hue
     
  */

  local urgl,R,q;
  q=sq;
  
  extern spmodel1;

  //  nabvec=[];
  if(is_void(nabvec)) nabvec=numberof(v.ages);
  
  if(is_void(patt)) patt="0patt";
  if(is_void(N)) N=1;
  if(is_void(pl)) pl=1;

  if (patt=="SDSS") {
    write,"IM a big dork";
    pathna="/raid/ocvirk/SDSS/spectro/1d_20/0266/1d/spSpec-51630-0266-174.fit";
    u=convertSDSS(pathna,noplot=1);};

  if(patt=="0patt") u=_mok(v,snr);

  hop=sfit(u,v,meval=5,noskip=1,kin=kin,epar=epar,pl=0,nde=nde);//mus=[1.,1.,1.,1.,1.,1.]);
  if(numberof(q)!=numberof(where(hop!=0.))) {
    error, "bad size for q, debug & see hop";};
  deconv_mux=deconv_muv=deconv_muz=deconv_mue=mub=0.;
  base=v;

  //write,vec(:40)(sum) ////////////////// DEBUG /////////////////////

  write,"INFO nabvec:";
  info(nabvec);
  if(!is_void(nabvec)) {
    write,"WARNING: nabvec not void! Will compute Q0(sqrt(SAD))! (i.e. sad must be given physical, not square rooted)";
    q(:nabvec)=sqrt(q(:nabvec));
  };

  //write,vec(:40)(sum) ////////////////// DEBUG /////////////////////

  
  urgl=Q0(q,g,[kin,epar]);
  
  //error;
  pmodel1=(kin==0)?model:roll(model)(pad1+1:pad2+dd)(::-1);
  gals=[];
  spmodel1=pmodel1; // FOR DEBUG
  if (N==1) pmodel1/=pmodel1(avg);
  spmodel2=pmodel1;  // FOR DEBUG
  //error
  
  mok_dir="/"+strjoin(split2words(fseed,sep="/")(1:0),"/")+"/";
  blurk=exec("mkdir -p "+mok_dir);
  
  
  for (i=1;i<=nR;i++){
    
    flux=pmodel1*(1.+(1./snr)*random_normal(numberof(pmodel1)));
    //error;
    wave=x0;
    //sigm=pmodel1/snr;   works fine
      sigm=array(pmodel1(avg)/snr,numberof(pmodel1));
      mask=0.*flux;
      
      gal=galStruct();
      gal.redshift=0.;
      gal.SNR=snr;
      gal.name=fseed+pr1(i);
      gal.filename=fseed+pr1(i)+".pdb";
      gal.resfile=fseed+pr1(i)+".res13";
      gal.result_dir=mok_dir;
      R=v.R; // resolution of the basis
      f=createb(gal.filename);
      save,f,gal,flux,wave,sigm,mask,R;
      close,f;
      write,"saved"+gal.filename;
      grow,gals,gal;
  };

  if(pl!=0){
    ws;
    plh,flux;
  };
  
  f=createb((fseed+"q.yor"));
  save,f,q,gals;
  close,f;
  write,"saved"+fseed+"q.yor";
  
  return gals;
};

func _mok(v,snr) {
/* DOCUMENT 
   sub-routine for mock data used by MakMok2

*/

  local srange,gal,dl;
  flux=array(0.1,numberof(v.wave));  
  wave=v.wave;
  z=0.;
  mask=flux*0.+0.;
  sigm=flux/snr;
  //  flux*(1./SNRl(where((h.object)==namel)(1)));
  gal=galStruct();
  gal.redshift=0.;
  gal.SNR=100.;
  //SNRl(where((h.object)==namel)(1)); // getSNR(filelist(i));
  gal.name= "crapmok";
  gal.filename="crapmok.pdb";
  gal.resfile="crapmok.res";
  gal.result_dir=exec("pwd")(1);
  save,createb(gal.filename),gal,flux,wave,sigm,mask;
  
  return gal;
};


func dburstmok(base,__nr=,ns=,zbs=,snr=,patt=,fseed=,simdir=,w=){
  /* DOCUMENT
     produces a bunch of double burst mock data in with fseed as seed names.
     the bursts are chosen between ages(1) and ages(2), same conventions as bRbasis, ns is the number of steps between ages(1) and ages(2). nr, number of realizations of the same expriment.
     snr is snr per pixel
     files and subdirectories are dumped into simdir
     the metallicities are given by vector zbs
     the widths of the bursts given by w
     simdir should be absolute name otherwise crash
  */

  local __nr;
  
  if (is_void(ns)) ns=10;
  if(is_void(zbs)) zbs=[0.02,0.02];
  if(is_void(snr)) snr=100.;
  if (is_void(__nr)) __nr=1;
  if (is_void(fseed)) fseed="truc";
  if (is_void(simdir)) {
    curdir=exec("pwd");
    simdir=curdir+"/sims/";
    simdir=simdir(1);
  };
  if (is_void(w)) w=0.1;

  nab=numberof(base.ages);
  age1=base.ages(min);
  age2=base.ages(max);
  bages=spanl(age1,age2,ns);
  ibages=[];
  for(i=1;i<=ns;i++){
    grow,ibages,(abs(bages(i)-base.ages))(mnx);
  };

  a1=ibages;
  na1=numberof(a1);
  nb=nab;


  _data=[];

  if(1){
  for(i=2;i<=na1;i++){
    for(j=1;j<i;j++){
      write,i,j;
      ca1=a1(i);
      ca2=a1(j);
      fseed=simdir+"/a1="+pr1(ca1)+"/a2="+pr1(ca2)+"/";
      //q=sqrt(makebump(nb,[ca1,ca2],0.6,N=1));
        q=sqrt(makebump(nb,[ca1,ca2],w,N=1));
        
        //m1=max(abs(random(1))*_m(max),_m(min));
        //m2=max(abs(random(1))*_m(max),_m(min));
        m1=Zrescale(zbs(1));
        m2=Zrescale(zbs(2));
        grow,q,span(m1(1),m2(1),nb);
        grow,q,0.;
        data=MakMok2(q,base,0,1,snr,__nr,patt,fseed,N=0);
        grow,_data,data;
        write,fseed;
    };
  };
  };
        

  write,ibages;
  return _data;
  
};



func mburstmok(base,__nr=,ns=,zbs=,snr=,patt=,fseed=,simdir=,w=){
  /* DOCUMENT
     produces a bunch of mono burst mock data in with fseed as seed names.
     the bursts are chosen between the extreme ages of the basis, same conventions as bRbasis, ns is the number of steps between ages(1) and ages(0). nr, number of realizations of the same expriment.
     snr is snr per pixel
     files and subdirectories are dumped into simdir
     the metallicities are given by vector zbs
     the widths of the bursts given by w
     simdir should be absolute name otherwise crash
  */
  
  local __nr;
  
  if (is_void(ns)) ns=10;
  if(is_void(zbs)) zbs=[0.02,0.02];
  if(is_void(snr)) snr=100.;
  if (is_void(__nr)) __nr=1;
  if (is_void(fseed)) fseed="truc";
  if (is_void(simdir)) {
    curdir=exec("pwd");
    simdir=curdir+"/sims/";
    simdir=simdir(1);
  };
  if (is_void(w)) w=0.1;
  
  nab=numberof(base.ages);
  age1=base.ages(min);
  age2=base.ages(max);
  bages=spanl(age1,age2,ns);
  ibages=[];
  for(i=1;i<=ns;i++){
    grow,ibages,(abs(bages(i)-base.ages))(mnx);
  };
  
  a1=ibages;
  na1=numberof(a1);
  nb=nab;
  
  
  _data=[];
  
  if(1){
    for(i=1;i<=na1;i++){
      write,i;
      ca1=a1(i);
      fseed=simdir+"/a1="+pr1(ca1)+"/";
      q=sqrt(makebump(nb,ca1,w,N=1));
      
      //m1=max(abs(random(1))*_m(max),_m(min));
      //m2=max(abs(random(1))*_m(max),_m(min));
      m1=Zrescale(zbs(1));
      m2=Zrescale(zbs(2));
      grow,q,span(m1(1),m2(1),nb);
      grow,q,0.;
      data=MakMok2(q,base,0,1,snr,__nr,patt,fseed,N=0);
      grow,_data,data;
      write,fseed;
    };
  };
  
  
  
  write,ibages;
  return _data;
  
};

  

