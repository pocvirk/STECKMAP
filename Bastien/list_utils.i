func readList(filename)
/* DOCUMENT readList(filename)
     return the list saved by saveList in file "filename"
   SEE ALSO: saveList
 */
{
  local rslt;

  fstrm=openb(filename);
  vname=*get_vars(fstrm)(1);
  rslt=__readList(fstrm,vname,0,1);
  return rslt;
}

func saveList(filename,list)
/* DOCUMENT saveList(filename,list)
   save the list 'list' in file filename.
   SEE ALSO: readList
 */
{
  fstrm=createb(filename);
  __saveList,fstrm,list,0;
  close,fstrm;
}


func __readList(fstrm,vname,level,&idx)
{
  head_varname=swrite(format="___list_head_%d",level);
  tail_varname=swrite(format="___list_tail_%d",level);
  if(vname(idx)!=head_varname) error,"Bad file format !!";
  lst=[];nelmt=0;
  while(vname(++idx)!=tail_varname)
    {
      nelmt++;
      vnelmt=swrite(format="___list_%d_elmt_%d",level,nelmt);
      if(vname(idx)!=vnelmt)
        lst=_cat(lst,_lst(__readList(fstrm,vname,level+1,idx)));
      else
        lst=_cat(lst,get_member(fstrm,vnelmt));
    }
  return lst;
}
  
func __saveList(fstrm,sublist,level)
{
  head_varname=swrite(format="___list_head_%d",level);
  tail_varname=swrite(format="___list_tail_%d",level);
  headvalue="HEAD LIST FLAG";
  tailvalue="TAIL LIST FLAG";
  add_variable,fstrm,-1,head_varname,structof(headvalue),dimsof(headvalue);
  get_member(fstrm,head_varname)=headvalue;
  tmp=sublist;
  nelmt=0;
  while(!is_void(tmp))
    {
      elmt=_nxt(tmp);
      nelmt++;
      if(typeof(elmt)=="list")
        {
          __saveList,fstrm,elmt,level+1;
          continue;
        }
      vnelmt=swrite(format="___list_%d_elmt_%d",level,nelmt);
      add_variable,fstrm,-1,vnelmt,structof(elmt),dimsof(elmt);
      get_member(fstrm,vnelmt)=elmt;
    }
  add_variable,fstrm,-1,tail_varname,structof(tailvalue),dimsof(tailvalue);
  get_member(fstrm,tail_varname)=tailvalue;
}

