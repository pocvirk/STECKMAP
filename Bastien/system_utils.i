require,"Bastien/ascii_utils.i"

func getVars(file)
{
  fstrm=openb(file);
  nvars=*get_vars(fstrm)(1);
  for(i=1;i<=numberof(nvars);i++)
    {
      write,format="%s ",nvars(i);
      info,get_member(fstrm,nvars(i));
    }
  close,fstrm;
}

func listFunc(file)
{
  local tmp;
  if(fileExist(Y_SITE+"i0/"+file)) tmp=Y_SITE+"i0/"+file;
  if(fileExist(Y_SITE+"i/"+file))  tmp=Y_SITE+"i/"+file;
  if(fileExist("~/Yorick/"+file))  tmp="~/Yorick/"+file;
  if(fileExist(file))              tmp=file;
  if(is_void(tmp)) error,"Cannot find file '"+file+"'";

  lines=sreadAscii(tmp,verb=0);
  idx=where(strmatch(lines,"func "));
  if(is_array(idx))
    {
      lines=strtok(lines(idx),"(")(1,);
      write,format="  %s\n",lines;
    }
}

func format_dir_name(dir)
{
  if(strpart(dir,0:0)=="/") return strpart(dir,1:-1);
  return dir;
}

func fileExist(filename)
{
  if(f=open(filename,"r",1)) {close,f;return 1;}
  else return 0;
}

func makeUniq(template,dir=,suffix=,crt=,nbr=)
/* DOCUMENT makeUniq(template,dir=,suffix=,crt=,nbr=)
   Return a list of file names created with template that don't already exist.

   the results name are " template+"XXX"+suffix
   where XXX is choose to make the file name uniq.
   SEE ALSO:
 */
{
  suffix=is_void(suffix)?"":suffix;
  crt   =is_void(   crt)? 1:crt;
  dir   =is_void(   dir)? 0:dir;
  nbr   =is_void(   nbr)? 0:1;
  
  if(dir)
    template+=swrite(format="%03d",long(random(dimsof(template))*1000));
  if(nbr&&!dir)
    rslt=swrite(format="%s%03d",template,array(1,numberof(template)))+suffix;
  else
    rslt=template+suffix;
    
  for(i=1;i<=numberof(template);i++)
    {
      if(fileExist(rslt(i)))
        {
          ref=(nbr&&!dir);
          do{
            ref++;
            rslt(i)=swrite(format="%s%03d",template(i),ref)+suffix;
          } while(fileExist(rslt(i)));
        }
      if(crt)
        if(dir)
          mkdir,rslt(i);
        else
          close,create(rslt(i));
    }
  return rslt;
}

func exec(order,&status)
{
  fstrm=popen(order+"; echo $?",0);
  rslt=readAsciiStream(fstrm);
  close,fstrm;
  status=0;
  sread,format="%d",rslt(0),status;
  if(numberof(rslt)==1) return [];
  return rslt(:-1);
}

func get_uniq_name(void)
/* DOCUMENT get_uniq_name()
   return a file name that not exist and then can be used safely
 */
{
  while(1)
    {
      filename=swrite(format=".crap_tmp_%5d",int(random([0])*100000));
      if(open(filename,"r",1)) continue;
      create,filename;
      break;
    } 
  return filename;
}

func eval_string(str)
/* DOCUMENT eval_string(str)
     launch the command in string str:

   EXAMPLES: eval_string("a=[1,2,3]")
   print,a;
   SEE ALSO:
 */
{
  filename=get_uniq_name();
  write,create(filename),format="%s\n",str;
  include,filename,1;
  remove,filename;
}

func screen_size(void)
{
  status=0;
  cmd="xrdb -symbols -screen | grep 'D%s' | awk -F \"=\" '{print $2}'";
  width =str2long(exec(swrite(format=cmd,"WIDTH"),status))(1);
  height=str2long(exec(swrite(format=cmd,"HEIGHT"),status))(1);
  return [width,height];
}
