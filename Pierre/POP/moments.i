

func moments(df,v,p=,N=){

  /* DOCUMENT moments(df,s)
   computes and prints (if p non void) mean, standard deviation, skewness and curtosis of distribution function df defined on the support v.
   recall that standard deviation (ecart type)=sqrt(Variance) and ecart type=RMS only if mean=0

  */

  local cur,mea,stand,skew,curtos;

  if(is_void(N)) N=1;

cur=df;

if(N==1) cur=cur/integ(cur,v,v(0));
 
mea=integ(cur*v,v,v(0));
stand=sqrt(integ(cur*(v-mea)^2,v,v(0)));
skew=integ(cur*(v-mea)^3,v,v(0))/stand^3;
curtos=integ(cur*(v-mea)^4,v,v(0))/stand^4-3.;

 if (!is_void(p)){ 
print,"mean=",mea;
print,"stand=",stand;
print,"skew=",skew;
print,"curtos=",curtos;
 };
 
 
 return [ mea,stand,skew,curtos];
};



