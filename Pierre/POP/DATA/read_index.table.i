#include "Pierre/POP/sfit.i"

f=open("index.table");
a=rdline(f,100);
ta=split2words(a,sep=" ");
ta=ta(4:28,);
u=str2float(ta(,:8));
u=transpose(u);
// last column if u tells if index should be given adimensional or in magnitudes
cent=u(2:3,);
acent=cent(avg,);
ind=sort(acent);
u=u(,ind);

names=ta(,-1);
names=names(ind);
f=createb("Lickindex.yor");
save,f,u,names;


// create a smaller list without repeated indices

cent=u(2:3,);
acent=cent(avg,);
bc=u(4:5,);
rc=u(6:7,);
rc1=rc(1,);

gind=where(acent(dif)>=5.);
grow,gind,25;

u=u(,gind);
cent=u(2:3,);
acent=cent(avg,);
bc=u(4:5,);
rc=u(6:7,);
names=names(gind);
bands=bc*0.;
for(i=1;i<=numberof(acent);i++){bands(1,)=bc(avg,);bands(2,)=rc(avg,);};

func mergebands(bands){
  /*DOCUMENT
    bands is a mask, i.e. an array 2by the number of domains to mask
    if these bands are somehow overlapping, mergebands returns the domains
    merged so that none overlaps anymore.
  */
    
  nb=numberof(bands(1,));
  res=[];
  db=bands(*)(dif);
  ind=where(db<=0.);
  ind2=ind+1;
  grow,ind,ind2;
  ind2;

  //  res=bands(1);
  for(i=1;i<=2*nb;i++){ if ((is_void(dimsof(where(i==ind))))) grow,res,(bands(*))(i);
    write,bands(*)(i);
  };

  nbands=reform(res,2,int(numberof(res)/2));
  return nbands;
};

u1=mergebands(bands);
f=createb("merged_lick.yor");
save,f,u1;
close,f;
