symbol_set,swrite(format="int%d",sizeof(  char)),char;
symbol_set,swrite(format="int%d",sizeof( short)),short;
symbol_set,swrite(format="int%d",sizeof(   int)),int;
symbol_set,swrite(format="int%d",sizeof(  long)),long;
symbol_set,swrite(format="float%d",sizeof( float)),float;
symbol_set,swrite(format="float%d",sizeof(double)),double;


func getSData(fstrm,&adr,nb,typ)
{
  return char2Str(getData(fstrm,adr,nb,typ));
}
  
func getData(fstrm,&adr,nb,typ)
{
  if(typ!=string)
    rslt=array(typ,nb);
  else
    rslt=array(char,nb);

  _read,fstrm,adr,rslt;
  adr+=sizeof(rslt);
  if(typ==string) rslt=string(&rslt);
  if((nb==1)||(typ==string)) rslt=rslt(1);
  return rslt;
}


func char2Str(flag)
{
  dgt=['0','1','2','3','4','5','6','7','8','9','A','B','C','D','E','F']
  rslt=array(char,numberof(flag)*2);
  rslt(1::2)=dgt((flag>>4)+1);
  rslt(2::2)=dgt((flag&15)+1);
  return string(&rslt);
}
