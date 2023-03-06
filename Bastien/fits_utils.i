func fits_duplicate_old(source_file,target_file,new_data,overwrite=)
/* DOCUMENT fits_duplicate(source_file,target_file,new_data,overwrite=)
   After calling fits_duplicate, target_file will have the same header than
   source_file but with data array new_data.

   KEYWORDS: overwrite : if set and target_file exists, overwrite it
 */
{
  fh=fits_open(source_file);
  header=_car(fh,1);
  fits_close,fh;

  fh=fits_open(target_file,"w",overwrite=overwrite);
  idx=where(fits_ids(header));
  header=header(idx);
  _car,fh,1,header;
  _car,fh,2,fits_ids(header);
  fits_write_header,fh;
  fits_write_array,fh,data;
  fits_close,fh;
}

func fits_duplicate(source_file,target_file,new_data,new_ext,overwrite=)
/* DOCUMENT fits_duplicate(source_file,target_file,new_data,overwrite=)
   After calling fits_duplicate, target_file will have the same header than
   source_file but with data array new_data.

   If the fits_file is made of multiplet extension, then new_data(i) should
   be an array of pointer
   KEYWORDS: overwrite : if set and target_file exists, overwrite it
 */
{
  if(typeof(new_data)!="pointer") new_data=[&new_data];
  
  next=numberof((ext=fits_list_card(source_file,"NAXIS")));
  ndext=numberof((dext=where(ext)));

  ndata=numberof(new_ext);
  
  //check the format of new_data
  if(ndext<ndata)
    error,swrite(format="To much data array (%d), there is only %d [%d]  extension) !",ndata,ndext,next);

  //check dimension of each extension
  fh_src=fits_open(source_file);
  j=0;
  for(i=1;i<=ndata;i++)
    {
      fits_goto_hdu,fh_src,new_ext(i);
      
      dms_fits=fits_get_dims(fh_src);
      dms_data=dimsof(*new_data(++j));

      if((dms_data(1)!=dms_fits(1))||anyof(dms_data!=dms_fits))
        error,swrite(format="Dimensions don't match for extension %d",i);
    }

  fits_rewind,fh_src;
  fh_trg=fits_open(target_file,"w",overwrite=overwrite);

  j=0;
  for(i=1;i<=next;i++)
    {
      fits_goto_hdu,fh_src,i;

      if(i>1) fits_new_hdu,fh_trg,fits_get_xtension(fh_src);
      fits_copy_header,fh_src,fh_trg;
      j=where(new_ext==i);
      fits_write_header,fh_trg;
      if(is_array(j))
        {
          data=*new_data(j(1));
          fits_set,fh_trg,"BSCALE",1.0;
          fits_set,fh_trg,"BZERO",0.0;
        }
      else
        data=fits_read_array(fh_src);

      if(is_array(data))
        fits_set,fh_trg,"BITPIX",fits_type_bitpix(structof(data));
      fits_write_array,fh_trg,data;
      
    }
  fits_close,fh_trg;
  fits_close,fh_src;
}




func fits_bitpix_type(bitpix)
{
  type=["char","short","int","long","float","double"];
  ntype=numberof(type);
  stype=array(long,ntype);
  for(i=1;i<=ntype;i++) stype(i)=8*sizeof(symbol_def(type(i)));
  stype(5:6)*=-1;
  
  idx=where(bitpix==stype);
  return symbol_def(type(idx(1)));
}

func fits_type_bitpix(type)
{
  types=["char","short","int","long","float","double"];
  ntypes=numberof(types);
  stypes=array(long,ntypes);
  for(i=1;i<=ntypes;i++) stypes(i)=8*sizeof(symbol_def(types(i)));
  stypes(5:6)*=-1;

  idx=where(typeof(type(0))==types);
  return stypes(idx(1));
}

func fits_copy_header(fh_in,fh_out)
{
  local header;
  header=_car(fh_in,1);
  for(i=1;i<=numberof(header);i++)
    {
      key=strpart(header(i),1:8);
      if(key=="        ") continue;
      rslt=fits_parse(header(i),safe=1);
      if(is_void(rslt)) continue;
      fits_set,fh_out,key,rslt,_fits_parse_comment;
    }
}                   
                        
func fits_list_card(fh,card)
/* DOCUMENT fits_list, fh, card;
       -or- fits_list(fh,card)

    Get the FITS card values in all  hdu in FH.  FH can be the name of
    a  FITS  file or  a  FITS  handle FH  (the  input  handle is  left
    unchanged).  When called  as a subroutine, the list  is printed to
    terminal;  when called  as a  function,  the returned  value is  a
    string array with the names of the FITS extensions in FH.

    if void, card is equal to "XTENSION", so fits_list_card do the same
    than fits_list.
    
   SEE ALSO: fits, fits_list, fits_read_header, fits_next_hdu. */
{
  if(is_void(card)) card="XTENSION";
  /* Get header of primary HDU. */
  if (structof(fh) == string) {
    /* open FITS file for reading */
    fh = fits_open(fh);
  } else {
    /* make private copy of FITS handle */
    if (typeof(fh) != "list" || _len(fh) != 4) error, "bad FITS handle";
    filemode = _car(fh,3)(5);
    stream = _car(fh,4);
    if (filemode != 'r') error, "FITS file not open for reading";
    fh = fits_read_header(_lst([], [], [1, 0, 0, 0, filemode], stream));
  }
  descr = [];
  for (;;) {
    if (is_void(_car(fh,1))) {
      if (! am_subroutine()) return descr;
      write, format="HDU=%-3d  %-8s=\"%s\"\n",
        indgen(numberof(descr)), array(card,numberof(descr)),descr;
      return;
    }
    tmp=fits_get(fh,card);
    grow, descr, tmp;
    fits_next_hdu, fh;
  }
}

func fits_get_data_hdu(file)
{
  next=numberof((ext=fits_list_card(file,"NAXIS")));
  return where(ext);
}

func fits_raw_get_adr(file,&dims,&dext)
{
  dext=fits_get_data_hdu(file);
  ndext=numberof(dext);

  addr=array(long   ,ndext);
  dims=array(pointer,ndext);

  fh=fits_open(file);
  for(i=1;i<=ndext;i++)
    {
      fits_goto_hdu,fh,dext(i);
      dims(i)=&fits_get_dims(fh);
      offset=_car(fh,3);
      addr(i)=offset(3);
    }
  return addr;
}


func fits_raw_get_array(file,&dext)
{
  local dims,dext;


  bitpix=fits_list_card(file,"BITPIX");
  
  addr=fits_raw_get_adr(file,dims,dext);
  ndata=numberof(addr);

  bitpix=bitpix(dext);

  rslt=array(pointer,ndata);

  fstrm=open(file,"rb");
  sun_primitives,fstrm;
  for(i=1;i<=ndata;i++)
    {
      data=array(fits_bitpix_type(bitpix(i)),*dims(i));
      _read,fstrm,addr(i),data;
      rslt(i)=&data;
    }
  close,fstrm;
  return rslt;
}

func fits_raw_set_array(file,data)
{
  local dims,dext;


  bitpix=fits_list_card(file,"BITPIX");
  
  addr=fits_raw_get_adr(file,dims,dext);
  ndata=numberof(addr);
  bitpix=bitpix(dext);
  
  //check consitency
  for(i=1;i<=ndata;i++)
    {
      dms=*dims(i);
      npt=1;for(j=1;j<=dms(1);j++) npt*=dms(j+1);
      size=abs(bitpix(i)>>3)*npt;
      if(sizeof(*data(i))!=size)
        error,swrite(format="Size doesn't match for extension %d",i);
    }

  fstrm=open(file,"r+b");
  sun_primitives,fstrm;
  for(i=1;i<=ndata;i++)
      _write,fstrm,addr(i),*data(i);
  close,fstrm;
}

