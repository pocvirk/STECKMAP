#include "STECKMAP/Pierre/POP/sfit.i"
//#include "Pierre/POP/pop2ls.i"
// lAST CHANGES: 4th aug 2006 added D4000

//upload,yodir+"STECKMAP/Pierre/POP/DATA/lickindex.yor";
upload,steckmapdir+"Pierre/POP/DATA/merged_lick.yor";


func D4000(sp,wave){
  /* DOCUMENT
     returns Dn4000 defined in Balogh et al. 1999
     "differential galaxy evolution in cluster"
  */

  bb=[3850.,3950.];
  rb=[4000.,4100.];
  ibb=where(wave>=bb(1)&wave<=bb(2));
  irb=where(wave>=rb(1)&wave<=rb(2));
  return (sp(irb)(sum)/(sp(ibb)(sum)));
};

func __index(i,sp,wave){
  /* DOCUMENT
     returns the spectral index number i measured in sp
     have to have u 
  */

  bb=u([4,5],i);bb=bb(sort(bb));ibb=where(wave>=bb(1)&wave<=bb(2));
  rb=u([6,7],i);rb=rb(sort(rb));irb=where(wave>=rb(1)&wave<=rb(2));
  cb=u([2,3],i);icb=where(wave>=cb(1)&wave<=cb(2));
  
  bpc=integ(sp(ibb),wave(ibb),wave(ibb)(0))/(wave(ibb(0))-wave(ibb(1)));
  rpc=integ(sp(irb),wave(irb),wave(irb)(0))/(wave(irb(0))-wave(irb(1)));

  // pseudo-continuum definition
  a=(rpc-bpc)/(rb(avg)-bb(avg));
  b=bpc-a*bb(avg);
  pc=a*wave(icb)+b;

  if(db==1) error;
  res=integ(1.-(sp(icb)/pc),wave(icb),wave(icb)(0));
  //if(u(8,i)==1) res=-2.5*log10(abs(res)/(wave(icb)(0)-wave(icb)(1)));
  return res;

};

func s_index(s,sp,wave){
  return __index(where(names==s),sp,wave);
};

func index(i,sp,wave){

  di=dimsof(sp);
  if(di(1)==1) return __index(i,sp,wave);
  if(di(1)==2) {
    res=array(0.,di(3));
    for(j=1;j<=di(3);j++){res(j)=__index(i,sp(,j),wave);};};
  if(di(1)==3){
    res=array(0.,di(3),di(4));
    for(j=1;j<=di(4);j++){res(,j)=index(i,sp(,,j),wave);};};
  return res;
};




  
  

  

