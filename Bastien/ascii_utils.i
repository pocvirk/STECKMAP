func rreadAscii(filename)
{
  stride=1000;
  rslt=[];
  fstrm=open(filename,"r");
  while(1) {
    tmp=rdline(fstrm,stride);
    idx=where(tmp);
    if(is_array(idx)) grow,rslt,tmp(:idx(0));
    if(idx(0)!=stride) break;
  }
  return rslt;
}


func sreadAscii(filename,startl=,endl=,cmt=,verb=)
/* DOCUMENT sreadAscii(filename,startl=,endl=,cmt=,verb=)
   return an array of string that contains the lines of the file
   'filename' starting at line startl (default the first one), ending
   at line endl (default the last one).
   
   Comment lines  (lines containing the  string 'cmt'), if  any found,
   are displayed.
   
   Void  lines  (mean  with  no characters  except  carriage return)  are
   ignored, that's  mean you  must take into  account the number  of void
   lines when  you used the startl and  endl keywords. So if  you want to
   start to read your file from line 'sl' and if there is 'vl' void lines
   before this line, you should set keywords startl to 'sl-vl'.
     
 */
{
  local lines;

  if(is_void(  endl)) endl  = -1;
  if(is_void(startl)) startl=  1;
  if(is_void(   cmt)) cmt   ="#";
  if(is_void(  verb)) verb  =  1;
  
  wcr=wordcount(filename);

  if((wcr(1)==0)&&anyof(wcr(2:))) wcr(1)=1;
  if(wcr(1)==0)     error,"empty file!";
  if(wcr(1)<startl) error,"no line to read";

  nblines=wcr(1);
  if(endl<0) endl=nblines;
  endl  =min(nblines,  endl);
  startl=max(      1,startl);
  
  if(startl>endl) error,"no line to read";

  lines=array(string,nblines);
  f=open(filename,"r");
  read,f,format="%[^\n]",lines;
  close,f;

  lastg=0;
  while(strlen(lines(lastg))<=0) lastg--;
  if(lastg<0)
    {
      lines=lines(1:lastg);
      endl=min(endl,numberof(lines));
      if(startl>endl) error,"no line to read";
    }

  if(strlen(cmt)>0)
    mtc=strmatch(lines,cmt);
  else
    mtc=array(char(0),numberof(lines));
  idxc=where(mtc);
  if(numberof(idxc)&&verb) write,format="%s\n",lines(idxc);
  gidx=where(!mtc(startl:endl));
  if(numberof(gidx)) return lines(gidx+startl-1);
  return [];
  
}

func wordcount(name,startl=)
/* DOCUMENT return wc result apply on 'file' 
   [nb_lines,nb_words,nb_characters]
*/
{
  local sz;
  if(is_void(startl)) startl=1;
  if (!open(name, "", 1))
    error, "no such file or directory \""+name+"\"";
  sz= array(0L,3);
  f=popen(swrite(format="tail +%d %s | wc | awk '{print $1,$2,$3}'",startl,name), 0);
  if (sread(rdline(f), sz(1),sz(2),sz(3)) != 3)
    {close,f;error, "cannot get size of file \""+name+"\"";}
  close,f;
  return sz;
}


func readAsciiStream(fstrm,stride=)
/* DOCUMENT readAsciiStream(file_stream)
     read output of the ascii stream fstrm
   SEE ALSO: readAscii
 */
{
  if(is_void(stride)) stride=1000;
  rslt=[];
  i=0;
  do
    {
      grow,rslt,rdline(fstrm,stride);
    } while(rslt(0)!=string(0));
  idx=where(rslt!=string(0));
  if(!is_array(idx)) return [];
  return rslt(1:idx(0));
}



func setltf(&firstline,dims,type) 
/* DOCUMENT setltf(&firstline,dims,type)
    function used by writetf.

   EXAMPLE:
   data=array(int,[3,4,5,6])
   setltf,str,dimsof(data),typeof(data)
   str is equal to "int 3 4 5 6"

   SEE ALSO: writetf, readtf
 */
{
  sd=swrite(format=" %i",dims);
  firstline=type;
  for(i=1;i<=numberof(sd);i++)
    firstline+=sd(i);
  return;
}
  
func getltf(firstline,&dims,&type)
/* DOCUMENT getltf(&firstline,dims,type)
    function used by readtf. Do the inverse of setltf.

   EXAMPLE:
   getltf,"int 3 4 5 6",dims,type
   data=array(type,dims)
   info,data return array(int,[3,4,5,6])
     
   SEE ALSO: writetf, readtf
 */
{
  local _type;
  local _dims;

  sf=split2words(firstline,sep=" ");
  _type=sf(1);
  _dims=sf(2:);
  dims=str2long(_dims);
  type=symbol_def(_type);
}

func readtf(file)
/* DOCUMENT readtf(filename,&data)
            filename: name of the file

            read a array saved with writetf.
   SEE ALSO: writef
*/
{
  local dims,type,data;
  data=[]
 
  f=open(file,"r");
  firstline=rdline(f);
  getltf,firstline,dims,type;
  data=array(type,dims);
  st=swrite(format=",%i",dims);
  ss="";
  for(i=2;i<=numberof(st);i++)
    {
      ss+=st(i);
    }
  write,format="read : array(%s%s)\n",typeof(data),ss;
      
  read,f,data;
  return data;
}

func writetf(file,data)
/* DOCUMENT writetf(filename,data)
            filename: name of the file

            save the array data with its type and dimensions
            in ASCII file 'filename'. data can be recovered with
            readtf.
*/
{
  local firstline,f;

  f=open(file,"w");
  setltf,firstline,dimsof(data),typeof(data);
  write,format="%s\n",f,firstline;
  write,f,data;
  close,f;
  return;
}



func ReadAsciiFile(file,startl=,endl=)
/* DOCUMENT ReadAsciiFile(file,startl=,endl=)
   !!!! USE readAscii !!!!
     
 */
{
  write,format="[%s] **************************************\n","INFO";
  write,format="[%s] * You should use readAscii instead ! *\n","INFO";
  write,format="[%s] **************************************\n","INFO";
  return readAscii(file,startl=,endl=);
}

func readAscii(filename,startl=,endl=,cmt=,verb=)
/* DOCUMENT readAscii(file,startl=,endl=,verb=)

   read an  ASCII table from line  startl (1 by default)  to line endl
   (the last of the table by default).

   
   KEYWORDS  cmt : string that identify comment lines

   SEE ALSO: sreadAscii
   
     
 */
{
  lines=sreadAscii(filename,startl=startl,endl=endl,cmt=cmt,verb=verb);

  
  nblines=numberof(lines);
  if(nblines)
    {
      nbcolumns=numberof(split2words(lines(1),sep=" "));
      rslt=array(double,nbcolumns,nblines);
      n=sread(lines,rslt);
      return transpose(rslt);
    }
}

