require,"Bastien/string_utils.i"
require,"Bastien/system_utils.i"

struct __StructStruct {
  string name;
  long    nbm;
  pointer members; 
}

struct StructMember {
  string type;
  string name;
  long   nb;
}


struct_type_list=["char","short","int","long","float","double","complex","string","pointer"];
struct_type_rank=[     1,      2,    3,     4,      5,       6,        7,       0,        -1];

func get_struct_definition(x)
{
  tpx=typeof(x);
  if((tpx!="struct_definition")&&(tpx!="struct_instance")) error,"x is not a structure !";
  if(tpx=="struct_instance") x=structof(x);
  return x
}

func struct_get_member(x)
/* DOCUMENT struct_get_member(x)
     Return a StructStruct structure of the structure of x.
     x can be a structure instance or a structure definition.

   EXAMPLE: > struct strct {int a(10); string b;double d(4);}
            > s = struct_get_member(strct)
            > s.name
            "strct"
            > s.nbm
            3 // three members (a,b,d)
            > s.members(1)
            StructMember(type="int",name="a",nb=10)
            > s.members(2)
            StructMember(type="string",name="b",nb=1)
            > s.members(3)
            StructMember(type="double",name="d",nb=4)
            

 */
{
  x=get_struct_definition(x);
  decomp=print(x);
  nb_member=numberof(decomp)-2;

  rslt=__StructStruct(name="StructStruct");
  tmp=array(StructMember,3);
  tmp.type=["string","long","StructMember"];
  tmp.name=["name"  ,"nbm","members"];
  tmp.nb  =[       1,    1,nb_member];
  rslt.members=&tmp;
  
  struct_define,rslt;
  
  rslt=StructStruct(name=split2words(decomp(1))(2),nbm=nb_member);
  decomp=split2words(decomp(2:-1),sep=" ();");
  if(dimsof(decomp)(0)==3)
    {
      idx=where(strlen(decomp(,3))==0);
      if(is_array(idx)) decomp(idx,3)="1";
      snb=str2long(decomp(,3));
    }
  else
    snb=array(1,nb_member);
      
  tmp=array(StructMember,nb_member);
  tmp.type=decomp(,1);
  tmp.name=decomp(,2);
  tmp.nb  =snb;
  rslt.members=tmp;
  
  return rslt;
}

func struct_write_definition(s,file=)
/* DOCUMENT struct_write_definition(s)
   Write the structure definition describe by s (a StructStruct structure)

   KEYWORD: file : is set, the file name used, else a uniq random
                   file name is choose
   SEE ALSO:
 */
{
  if(is_void(file)) file=get_uniq_name();
  fstrm=create(file);
  write,fstrm,format="struct %s {\n",s.name;

  if(typeof(s.members)=="pointer") members=*s.members;
  else members=s.members;

  write,fstrm,format=" %s %s%s;\n",members.type,members.name,merge2(swrite(format="(%d)",members.nb),array("",numberof(members)),members.nb>1);
  write,fstrm,"}";
  close,fstrm;
  return file;
}


func struct_define(s,file=,keep=)
{
  if(is_void(keep)) keep=0;
  file=struct_write_definition(s,file=file);
  include,file,1;
  if(!keep) remove,file;
}

func struct_get_identical_members(strct_1,strct_2,compatible_type=)
/* DOCUMENT struct_get_identical_members(strct_1,strct_2)
    find the members of strct_1 and strct_2 with the same name.

   KEYWORD: compatible_type : perform a check on the type of the members
                              assuming members of strct_2 will be initia-
                              lised with the ones of strct_1
                               0 -> no check
                               1 -> simple  (numerical; string; pointer)
                       default 2 -> up cast (allow cast to upper numerical type
                                              ex float -> double; string; pointer)
                               3 -> strict  (same type)
     
   SEE ALSO:
 */
{
  if(is_void(compatible_type)) compatible_type=2;

  s_1=struct_get_member(strct_1);
  s_2=struct_get_member(strct_2);

  idx_same=where2(s_1.members.name(,-)==s_2.members.name(-,));

  if(!compatible_type) return s_1.members.name(idx_same(1,));

  type1=s_1.members.type(idx_same(1,));
  type2=s_2.members.type(idx_same(2,));

  r1=struct_type_rank((type1(,-)==struct_type_list(-,))(,mxx));
  r2=struct_type_rank((type2(,-)==struct_type_list(-,))(,mxx));
  
  if(compatible_type==1)
    {
      r1=(r1>0)-(r1<0);
      r2=(r2>0)-(r2<0);
      idx=where(r1==r2);
    }
  else if(compatible_type==2)
    idx=where( ((r1<=0)&(r1==r2)) | ((r1>0)&(r1<=r2)) );
  else
    idx=where(type1==type2);
    
  if(!is_array(idx)) return [];
  return s_1.members.name(idx_same(1,idx));
}

func struct_copy_members(dest_strct,src_strct,dest_fields,src_fields,copy_same=,compatible_type=)
/* DOCUMENT  struct_copy_members(dest_strct,src_strct,dest_fields,src_fields)
   copy members src_fields of src_strct to members dest_fields of des_strct.
   dest_strct can be an array of structure or a structure definition.
   src_strct must be an array of structure (the data source).
   dest_fields and src_fields are two string arrays with the same dimensions.
   If one of them is void, it is initialised to the values of the other.
   If they are both void, copy_same is set to 1 (except if the user set it to 0
   and then nothing append)

   KEYWORDS: copy_same       : if set (not the default) copy also
                               the fields with the same names
             compatible_type : see struct_get_identical_members
   SEE ALSO:
 */
{
  if(typeof(dest_strct)=="struct_definition")
    rslt=array(dest_strct,numberof(src_strct));
  else
    rslt=dest_strct;

  nbdf=numberof(dest_fields);
  nbsf=numberof(src_fields);
  if(nbdf&&nbsf&&(nbdf!=nbsf))
    error,"Dimension probleme with destination and source fields !";

  same_members=[];
  
  if(!nbdf&&!nbsf)
    {
      if(!is_void(copy_same)&&copy_same==0) write,format="[struct_utils.i] %s\n","nothing to do !";
      copy_same=1;
    }

  if(copy_same)
    same_members=struct_get_identical_members(src_strct,dest_strct);

  if(!nbsf&&nbdf) {
    src_fields=dest_fields;
    nbsf=nbdf;
  }
  if(!nbdf&&nbsf) {
    dest_fields=src_fields;
    nbdf=nbsf;
  }

  grow,src_fields,same_members;
  grow,dest_fields,same_members;

  for(i=1;i<=numberof(src_fields);i++)
    get_member(rslt,dest_fields(i))=get_member(src_strct,src_fields(i));

  return rslt;
}

func struct_cast(data,new_strct)
/* DOCUMENT struct_cast(data,new_strct)
     Cast data to the new_strct type. The structure of data must be
     identical to the one of new_strct (except for the size of array members).
   SEE ALSO:
 */
{
  s_type_new=get_struct_definition(new_strct);

  s_old=struct_get_member(data);
  s_new=struct_get_member(new_strct);

  rslt=array(s_type_new,dimsof(data));
  
  if((s_old.nbm!=s_new.nbm)||
     anyof(s_old.members.type!=s_new.members.type)||
     anyof(s_old.members.name!=s_new.members.name))
    error,"Structures not compatible for a cast";

  for(i=1;i<=s_new.nbm;i++)
    {
      member_old=s_old.members(i);
      member_new=s_new.members(i);
      nb_old=member_old.nb;
      nb_new=member_new.nb;
      
      tmp=get_member(data,member_old.name);

      nb_data=min(nb_old,nb_new);

      type=member_old.type;
      
      if(noneof(["char","int","long","float","double","string","pointer"]==type))
        tmp=struct_cast(tmp,symbol_def(member_new.type));

      if(nb_old==nb_new)
        get_member(rslt,member_new.name)=tmp;
      else if(nb_old!=1&&nb_new!=1)
        get_member(rslt,member_new.name)(1:nb_data,..)=tmp(1:nb_data,..);
      else if(nb_old==1)
        get_member(rslt,member_new.name)(1:nb_data,..)=tmp(-,..);
      else if(nb_new==1)
        get_member(rslt,member_new.name)=tmp(1,..);
    }
  return rslt;
}

func struct_change_member_size(s,members,new_sizes,copy=)
/* DOCUMENT struct_change_member_size(s,members,new_sizes,copy=)
   change the size of  members to new_sizes.

   KEYWORDS copy : if set, convert s to the new structure, else
                   just return the new structure type


   EXAMPLE : > struct x { int m1;long m2;double m3;}
             > v1=array(x,2)
             > v1.m1=1
             > v1.m3=3
             > v1.m1;v1.m3
             [1,1]
             [3,3]
             > v1=struct_change_member_size(v1,["m1","m3"],[2,3])
             > v1.m1;v1.m3
             [[1,0],[1,0]]
             [[3,0,0],[3,0,0]]
             > 
   SEE ALSO:
 */
{
  if(is_void(copy)) copy=1;
  s_strct=struct_get_member(s);

  for(i=1;i<=numberof(members);i++)
    {
      idx=where(s_strct.members.name==members(i));
      if(!is_array(idx)) write,format="cannot find member : %s !\n",members(i);
      idx=idx(1);
      s_strct.members(idx).nb=new_sizes(i);
    }

  struct_define,s_strct;
  if(copy) return struct_cast(s,symbol_def(s_strct.name));
  return symbol_def(s_strct.name);
}

func struct_set_array_member(s,member,idx,data)
/* DOCUMENT struct_set_array_member(s,member,idx,data)

   Simply do : s(idx).member=data
          idx must be a scalar, member is a string, and data
          can be an array. The data s will be resized if necessary

   EXMAPLE: > struct x {int m1;long m2;double m3;}
            > v1=array(x,2)
            > v1.m2
            [0,0]
            > v1=struct_set_array_member(v1,"m2",2,[3,1,4,1,5])
            > v1.m2
            [[0,0,0,0,0],[3,1,4,1,5]]
            
   SEE ALSO:
 */
{
  ndata=numberof(data);
  nmemb=numberof(get_member(s(1),member));
  
  if(ndata>nmemb) s=struct_change_member_size(s,member,ndata,copy=1);
  if(numberof(get_member(s(idx),member))==1)
    get_member(s(idx),member)=data(1);
  else
    get_member(s(idx),member)(1:ndata)=data;
  return s;
}



                    
