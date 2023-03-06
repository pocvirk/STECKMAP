// Cool tools for handling spectra

func rebin(a,x0,nr,&x1,&rmat,RM=){
  /* DOCUMENT
     tentative rebin with top-hat rebinning
     ok, gives correct chi^2, showing the cariance-covariance matrix remains diagonal
     NOTE: resampling is linear in wavelength
     if RM=1 then the matrix of the linear matrix representing the variable change is returned in rmat
  */

  nx1=int(numberof(a)/nr);
  res=array(0.,nx1);
  x1=res;
  if(RM==1) rmat=array(0.,nx1*nr,nx1);
  for(i=1;i<=nx1;i++){
    res(i)=(a(1+nr*(i-1):nr*i))(avg);
    x1(i)=(x0(1+nr*(i-1):nr*i))(avg);
    if(RM==1) rmat(1+nr*(i-1):nr*i,i)=1./nr;
  };
  if(RM==1) rmat=transpose(rmat);
  return res;
};

func fits_list_header(fh,pr=){
  /* DOCUMENT
     prints all the keywords and their values listed from fits file handle fh
  */
  li=fits_get_keywords(fh);
  res=array(string,2,numberof(li));
  if (pr==1) for(i=1;i<=numberof(li);i++) write,li(i)+" : "+pr1(fits_get(h,li(i)));
  for(i=1;i<=numberof(li);i++) {
    res(1,i)=li(i);
    res(2,i)=pr1(fits_get(h,li(i)));
  };
  return res;
};


func resample(y,x,xp){

  /*DOCUMENT
    this is similar to a resample by sum. We just find out which pixels of (y,x) fall into the xp bins and sum them
  */
  
  nx=numberof(x);
  nxp=numberof(xp);
  res=array(0.,nxp-1);
  for(i=1;i<=nxp-1;i++){
    ind=where((x>=xp(i))&(x<=xp(i+1)));
    if(!is_void(dimsof(ind))) res(i)=y(ind)(sum);
  };

  return res;
};
