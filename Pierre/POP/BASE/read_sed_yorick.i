// To read BC2003 BC03 ASCII SED files and make a big DATA cube out of it!!
// ctrl c c left click to comment stuff
// Age is in Myr like in PEGASE


nab=221
bsp=array(0.,6954,221,6);
_m=[0.0001, 0.0004, 0.004, 0.008, 0.02, 0.05](::-1);
for(j=1;j<=6;j++){

  fn="bc2003_hr_m"+pr1(1+j)+"2_chab_ssp.ised_ASCII";
  f=open(fn,"r");
  a=rdline(f,1000000);
  close,f;
  
  stra=a(1:37);
  ages=str2float(split2words(strjoin(stra,"  ")))(2:)/1.e6;
  comment=a(42);
  z=str2float(split2words(comment,sep="Z=")(0));
  strw=a(43:679);
  wave=str2float(split2words(strjoin(strw,"  ")))(2:);
  
  strf=a(680:338214);
  //fl=str2float(split2words(strjoin(strf,"  "))); // too big, have to cut it in pieces
  //rstrf=reform(strf,6137,55);
  //jstrf=[];ff=[];
  //for(i=1;i<=6137;i++){grow,jstrf,strjoin(rstrf(1,),"  ");write,i};
  //for(i=1;i<=6137;i++){
  //  for(j=1;j<=55;j++){grow,ff,str2float(split2words(rstrf(i,j)));write,i;};
  //};
  // LOUSY
  
  ff=[];ff1=[];
  for(i=1;i<=numberof(strf);i++){grow,ff,str2float(split2words(strf(i),sep="  "));write,((mod(i,1000)==0)?i:[]);
  if(mod(i,5000)==0) {grow,ff1,ff;ff=[];};
  };
  
  //make a nice array
  ind=where(ff1==6900.);
  dind=int((ind(dif))(avg));
//nab=int(numberof(ind));
  nw=numberof(wave);
  sp=array(0.,dind,nab);
  
  for(i=1;i<=nab-2;i++){sp(,i)=ff1(1+dind*(i-1):dind*i);};
    bsp(,:-2,j)=sp(,:-2);
};

bsp=bsp(,,::-1);
_a=indgen(nab);
ta=ages;
_x0=wave;
bloc=bsp
Res=2000; // that's 3 angstroms at 6000 angstroms

// now cut it to what is really useful by now
nab=219;
bloc=bloc(2:6901,:nab,);
ta=ta(:nab);
_a=indgen(nab)

f=createb("bc03_Pa94_Cha_raw.yor");
save,f,bloc,_x0,_a,_m,ta,Res;
close,f;

// make it smaller !!
bloc=bloc(,::3,);
ta=ta(::3);
_a=_a(::3);
f=createb("bc03_Pa94_Cha_raw_smaller.yor");
save,f,bloc,_x0,_a,_m,ta,Res;
close,f;





                                       


