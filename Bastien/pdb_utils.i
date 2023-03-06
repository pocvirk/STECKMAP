require,"Bastien/idl_utils.i"
func upload(__fname,s=)
/* DOCUMENT  upload data file. Can be Yorick (.pdb) or IDL (.sav) file
*/
{
  if(is_void(s)) s=0;
  
  suffix=strpart(__fname,-2:);
  if(suffix=="sav")  __ff=idl_open(__fname);
  else __ff=openb(__fname);
    
  if(!s) write,"recovered", *get_vars(__ff)(1);
  restore,__ff;
  close,__ff;
}

func ssave(fstrm,var_name,var_value)
{
  if(is_void(var_value)) var_value=symbol_def(var_name);
  add_variable,fstrm,-1,var_name,structof(var_value),dimsof(var_value);
  get_member(fstrm,var_name)=var_value;
}


func duplicateb(fstrmi,fstrmo,vars,nvars)
{
  name_list=(is_void(vars)?(*get_vars(fstrmi)(1)):vars);
  if(!is_void(nvars))
    {
      idx=where((name_list(,-)!=nvars(-,))(,sum));
      if(is_array(idx))
        name_list=name_list(idx);
      else
        return;
    }
  name_list=suniq(name_list);
  for(i=1;i<=numberof(name_list);i++)
    ssave,fstrmo,name_list(i),get_member(fstrmi,name_list(i));
}
