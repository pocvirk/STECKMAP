func cohm(a){
  /*DOCUMENT
    COHM: count how many: returns a vector containing the different values appear in a.
  */

  local na;
  na=numberof(a);
  res=a(1);
  for(i=2;i<=na;i++){
    if((a(i)!=a(i-1))&(is_void(dimsof(where(a(i)==res))))) grow,res,a(i);};
  return res;
};
    
