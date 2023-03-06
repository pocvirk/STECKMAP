#include "Bastien/utils_binary.i"

func readChars(fstrm,&adr,modu=)
{
  if(is_void(mod)) modu=4;
  size=getData(fstrm,adr,1,int4);
  rsize=((size+modu-1)/modu)*modu;
  rslt=getData(fstrm,adr,rsize,char);
  return rslt(1:size);
}

func readChunk(fstrm,&adr,&name,&data,only_name=)
{
  if(is_void(only_name)) only_name=0;
  
  _adr=adr;
  crap=getData(fstrm,adr,1,int4);
  if(crap!=2) return 1;
  off=getData(fstrm,adr,1,int4);
  //skip "0000 0000 0000 0000"
  adr+=8;
  name=string(&readChars(fstrm,adr,modu=4));
  type=getData(fstrm,adr,1,int4);
  crap=getData(fstrm,adr,2,int4);
  size=getData(fstrm,adr,1,int4);
  nbytes=getData(fstrm,adr,1,int4);
  ndata=getData(fstrm,adr,1,int4);
  start=off-nbytes;
  if(type==4) type=float;
  else if (type==5) type=double;
  if(!only_name) data=getData(fstrm,start,ndata,type);
  adr=off;
  return 0;
}

func skipHeader(fstrm)
{
  adr=0;
  crap=getData(fstrm,adr,8,char);
  adr =getData(fstrm,adr,1,int4);
  
  crap=getData(fstrm,adr,1,int4);
  offs=getData(fstrm,adr,1,int4);
  return offs;
}

func getvarsIDL(filename,vname,rev=)
{
  if(is_void(rev)) rev=1;
  if(is_void(  s)) s=0;
  fstrm=open(filename,"rb");
  if(rev) sun_primitives,fstrm;
  off=skipHeader(fstrm);
  names=[];data=[];
  while(1)
    {
      tmp=off;
      if(readChunk(fstrm,off,name,only_name=1)) break;
      if(is_void(vname)) grow,names,name;
      else
        {
          if(vname!=name) continue;
          readChunk,fstrm,tmp,name,data,only_name=0;
          break;
        }
    }
  close,fstrm;
  if(!is_void(vname)) return data;
  return names;
}

func restoreIDL(filename,rev=,s=)
{
  if(is_void(rev)) rev=1;
  if(is_void(  s)) s=0;
  fstrm=open(filename,"rb");
  if(rev) sun_primitives,fstrm;
  off=skipHeader(fstrm);
  data=[];
  while(1)
    {
      if(readChunk(fstrm,off,name,data,only_name=0)) break;
      if(!s) write,format="restoring %s\n",name;
      symbol_set,name,data;
    }
  close,fstrm;
}
