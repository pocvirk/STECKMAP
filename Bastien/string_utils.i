func getLastWord(strings,sep=)
{
  words=split2words(strings,sep=sep);
  dms=dimsof(words);
  if(dms(1)==1) dms=[1,1,dms(2)];
  plast=(strlen(words)!=0)(..,sum);
  idx=where(plast==0);if(is_array(idx)) plast(idx)=1;

  idx=(plast-1)*dms(2)+indgen(1:dms(2));
  return words(idx);
}

func pchar(value)
{
  return string(&char(value));
}

func plurial(nb,letter)
{
  if(is_void(letter)) letter="s";
  return swrite(format="%s",((nb>1)?letter:""));
}

func add_suffix(strings,one=)
{
  if(is_void(one)) one=0;
  rslt=strings;
  base=suniq(strings);
  for(i=1;i<=numberof(base);i++)
    {
      idx=where(strings==base(i));
      if(!one&&(numberof(idx)<2)) continue;
      rslt(idx)=swrite(format="%s_%03d",rslt(idx),indgen(numberof(idx)));
    }
  return rslt;
}

func strreplace1(&str,pattern,newstr,start=)
{
  if(is_void(start)) start=1;
  sptr=*pointer(str);
  pptr=(*pointer(pattern))(:-1);
  len=strlen(pattern);
  if(!strmatch(strpart(str,start:),pattern)) return 0;
  idx=where(sptr(start:)==pptr(1))+start-1;  
  i=1;
  while(len!=(sum(sptr(idx(i):idx(i)+len-1)==pptr))) ++i;
  idxinf=idx(i)-1;
  idxsup=idx(i)+len;
  if(idxinf==0)
    if (idxsup>strlen(str)) rslt=newstr;
    else rslt=newstr+strpart(str,idxsup:);
  else
    if(idxsup>strlen(str)) rslt=strpart(str,1:idxinf)+newstr;
    else rslt=strpart(str,1:idxinf)+newstr+strpart(str,idxsup:);
  str=rslt;
  return idxinf+strlen(newstr)+1;
}

func strreplace(str,pattern,newstr)
/* DOCUMENT strreplace(str,pattern,newstr)
      Remplace all occurence of pattern in str by newstr
   SEE ALSO:
 */
{
  rslt=str;
  for(i=1;i<=numberof(rslt);i++)
    {
      start=1;tmp=rslt(i);
      while((start=strreplace1(tmp,pattern,newstr,start=start)));
      rslt(i)=tmp;
    }
  return rslt;
}


func strglue(words,sep=)
/* DOCUMENT strglue(words,sep=)
   glue the array string along the last dimension. Put a separator sep
   between words.
   SEE ALSO:
 */
{
  if(is_void(sep)) sep="";

  rslt=words(..,1);
  nb_words=dimsof(words)(0);
  for(i=2;i<=nb_words;i++)
    rslt=rslt+sep+words(..,i);
  return rslt;
}

func split2words(str,sep=)
/* DOCUMENT split2words(str,sep=)

   return the array of words of string line str.
   Separators are listed in string keyword sep.
   By default sep=" \t" (space and tabulation).

   EXAMPLE  split2words("Hello world.")          return ["Hello","world."]
            split2words("Hello world.",sep="l")  return ["He","o wor","d."]
   SEE ALSO: splittok
 */
{
  local ptr,idx;
  local sptr,j,nbcol;

  if(is_void(sep)||strlen(sep)<1) sep=" \t";
  dms=dimsof(str);
  if(is_scalar(str))  dms=[1,1];
  nstr=numberof(str);
  ptrrslt=array(pointer,dms);
  nbwords=array(long,dms);
  sptr=*pointer(sep);
  for(i=1;i<=numberof(str);i++)
    {
      ptr=*pointer(str(i));
      j=1;
      while(sptr(++j))
        {
          idx=where(ptr==sptr(j));
          if(numberof(idx)) ptr(idx)=sptr(1);
        }
      tmp=splittok(string(&ptr),tok=strpart(sep,1:1));
      nbwords(i)=numberof(tmp);
      ptrrslt(i)=&tmp;
    }
  dms=grow(dms(1)+1,max(nbwords),dms(2:));
  rslt=array(string,dms);
  for(i=1;i<=numberof(str);i++)
    rslt(1:nbwords(i),i)=*ptrrslt(i);
  rslt=transpose(rslt);
  if(is_scalar(str)) return rslt(1,);
  return rslt;
}

func splittok(s,tok=)
/* DOCUMENT   splittok(s,tok=)
     split s into strings by using tok as separator (default " ")
 */
{
  local r;
  if(is_void(tok)) tok=" ";
  r=[];
  do {s=strtok(s,tok);if(s(1)) grow,r,s(1);else break;s=s(2);} while(1);  
  return r;
}

func getformat(n)
/* DOCUMENT getformat(n)
     return "%d" if n is a long or an int, and "%f" if n is a double
     of a float.
 */
{
  local tpf;
  tpf=typeof(n);
  if(tpf=="long"||tpf=="int") return "%d";
  if(tpf=="double"||tpf=="float") return "%f";
  return "";
}

func nb2str(l,fmt=)
/* DOCUMENT nb2str(l,fmt=)
   convert an array of number to an array of string.

   KEYWORDS: fmt is the format of the convertion.
             if is void, "%d" is used for long and int
             and "%f" for float and double.
     
   SEE ALSO: long2str,float2str;
 */
{
  if(is_void(fmt)) fmt=getformat(l);
  return swrite(format=fmt,l);
}

func long2str(l)
/* DOCUMENT long2str(l)
   convert an array of int or long to an array of string.
 */
{
  return swrite(format="%d",l);
}

func float2str(f)
/* DOCUMENT float2str(l)
   convert an array of float or double to an array of string.
*/
{
  return swrite(format="%f",f);
}

func str2nb(s,type,fmt=)
/* DOCUMENT str2nb(s,type,fmt=)
    Converts an array of string s to an array of number of type 'type'
    by using fmt as convertion format.
    If fmt is void, the function use the result of
    getformat(type(0)).
 */
{
  local rslt;
  if(is_void(fmt)) fmt=getformat(type(0));
  rslt=array(type,dimsof(s));
  sread,s,format=fmt,rslt;
  return rslt;
}

func str2int(s)
/* DOCUMENT str2long(s)
   convert an array of string to an array of int.
     
   SEE ALSO: str2nb,str2long,str2float,str2double;
 */
{
  return str2nb(s,int,fmt="%d");
}

func str2long(s)
/* DOCUMENT str2long(s)
   convert an array of string to an array of long.
     
   SEE ALSO: str2nb,str2int,str2float,str2double;
 */
{
  return str2nb(s,long,fmt="%d");
}

func str2float(s)
/* DOCUMENT str2long(s)
   convert an array of string to an array of float.
     
   SEE ALSO: str2nb,str2int,str2long,str2double;
 */
{
  return str2nb(s,float,fmt="%f");
}

func str2double(s)
/* DOCUMENT str2double(s)
   convert an array of string to an array of double.
     
   SEE ALSO: str2nb,str2int,str2long,str2float;
 */
{
  return str2nb(s,double,fmt="%f");
}


func suniq(sarray)
/* DOCUMENT suniq(strin_array)
   Return an array of string equal to string_array without doublon element.

   EXAMPLE: suniq(["A","B","C","A","AA"]);
            ["A","B","C","AA"]
 */
{
  local rslt;
  rslt=sarray(sort(sarray));
  grow,rslt,rslt(0)+"XX";
  idx=where(rslt(:-1)!=rslt(2:));
  return rslt(idx);
}


func scomplete(str,chr,nbc)
/* DOCUMENT scomplete(str,nbc,chr)
   complete the string str with character chr till its length
   is equal to nbc.
   SEE ALSO:
 */
{
  lgt=strlen(str);
  if(is_void(nbc)) nbc=lgt(max);

  if(typeof(chr)=="string") chr=(*pointer(chr))(1);
  chr=char(chr);
  nb2add=nbc-min(lgt);
  if(nb2add<=0) return str;
  
  suff=string(&array(chr,nb2add));
  rslt=str;
  idx=where(lgt<nbc);
  if(is_array(idx)) rslt(idx)=strpart(rslt(idx)+suff,1:nbc);
  return rslt
  
}
