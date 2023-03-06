// read synthetic spectra from gonzalez delgado 2005




#if 0
file=["SSPPadova.z004","SSPPadova.z008","SSPPadova.z019"](::-1);
nm=numberof(file);
_m=[0.004,0.008,0.019](::-1);
trax="padova";
#endif

#if 1
file=["SSPGeneva.z001","SSPGeneva.z004","SSPGeneva.z008","SSPGeneva.z020","SSPGeneva.z040"](::-1);
nm=numberof(file);
_m=[0.001,0.004,0.008,0.02,0.04](::-1);
trax="geneva";
#endif


for(i=1;i<=nm;i++){
  f=open("/raid/ocvirk/BASE/"+file(i));
  a=rdline(f,18000);
  nend=13365;
  nstart=45;
  iages=43;
  if(trax=="geneva") {nstart=46;nend=13366;iages=44;};
  
  ages=str2float(split2words(a(iages))(4:-1));
  nages=numberof(ages);
  nl=nend-nstart+1;
  if(is_void(bloc)) bloc=array(0.,nl,nages,nm);
  if(is_void(_x0)) _x0=array(0.,nl);
  
  for(k=1;k<=nl;k++){
    tra=str2float(split2words(a(k+nstart-1))(1:-1));
    _x0(k)=tra(1);
    //write,tra(1);
    bloc(k,,i)=tra(2::3);
  };
  write,file(i);
  
};


_a=indgen(nages);
ta=ages/1.e6;
Res=8300.; //?
Rdlambda=.3 // smapling in angstroms

f=createb("/raid/ocvirk/BASE/GD05_salpeter_"+trax+".yor");
save,f,bloc,_x0,_a,_m,ta,Res,Rdlambda;
close,f;






