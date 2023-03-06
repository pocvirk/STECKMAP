// reads jimenez2004 models
nage=40;
nlambda=1221;
nm=4;
_m=[0.001, 0.004, 0.02, 0.05](::-1);
bloc=array(0.,nlambda,nage,nm);
file=["speed0500_35_kroupa.spec","speed0200_27_kroupa.spec","speed0040_25_kroupa.spec","speed0010_25_kroupa.spec"];
for(i=1;i<=nm;i++){
f=open(file(i),"r");
a=rdline(f,3+nage*nlambda);
close,f; 
ages=str2float((split2words(a(2),sep=",: "))(2:-1))*1.e3;
tra=a(4:);
stra=split2words(tra,sep=" ");
fstra=str2float(stra);
rfstra=reform(fstra(,2),nlambda,nage);
wave=fstra(:nlambda,1);
bloc(,,i)=rfstra;
};

_x0=wave;
_a=indgen(nage);
ta=ages;
Res=250.; //?

f=createb("/raid/ocvirk/BASE/SPEED_kroupa.yor");
save,f,bloc,_x0,_a,_m,ta,Res;
close,f;



