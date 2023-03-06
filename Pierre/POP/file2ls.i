// SOME TOOLS TO HELP FILE CONVERSIONS OF OUTPUTS

//#include "sfit.i"

func geni_restore(s,file){
  /* DOCUMENT
     s is a string, file is a filename.
     geni_restore generates the crap_restore.i file which will contain something like:
     f=open(file);restore,f,s;
     and then includes it so that the var with name given in s AS A STRING will be restored.
  */

  f=create("crap_restore.i");
  write,f,"crapf=openb(\""+file+"\");restore,crapf,"+s+";";
  close,f;
  include,"crap_restore.i",1;
  return;
};

func geni_write(s,file){
  /* DOCUMENT
     s is a string, file is a filename.
     geni_write generates the crap_write.i file which will contain something like:
     geni_restore(s,file1);f=open(file);write,f,dimsof(s);write,f,s;close,f;
     and then includes it so that the var with name given in s AS A STRING will be appended to file.
     NB: file must exist beforehand
  */

  f=create("crap_write.i");
  //write,f,"crapf=open(\""+file+"\",\"a\");write,crapf,dimsof("+s+"),format=\"\\n\";write,crapf,"+s+";close,crapf;"
  write,f,"write,af,dimsof("+s+"),format=\"\\n\";write,af,"+s+";"
  close,f;
  include,"crap_write.i",1;
  return;
};


func pdb2asc(file){
  /* DOCUMENT
     converts a pdb file to ascii (file is the pdb file), formatted as follows:
     for each variable in file, the name is given followed on the same line by its dimensions and then follows the data

  */

  upload,file,s=1;
  f=openb(file);
  bla=get_vars(f);
  varlist=*bla(1);
  nvarlist=numberof(varlist);
  //nvarlist=2;
  ascfile=file+".asc";
  af=create(ascfile);
  write,af,numberof(varlist),format="\n";

  for(i=1;i<=nvarlist;i++){
    write,af,varlist(i),format="\n";
    geni_write(varlist(i),ascfile);
  };
     
  close,af;
};

func dumpfits(file1,file2){
  /* DOCUMENT
     to be used on .res files created by sfit
     dumps the wavelengths, the data, the best fit pmodel1, the flux correction, the residuals, and the weights/mask used
     also adds keywords LWAge, LWMet, the ki2, and a code for the SSP model used
     use with pfits described in plot2ls.i
  */

  upload,file1,s=1;
  nr=numberof(x0);
  fh = fits_open(file2, 'w',overwrite=1);      // create new file
  fits_set, fh, "SIMPLE", 'T';    //"true FITS file";
  fits_set, fh, "BITPIX", -32; //"bits per pixel";
  fits_set, fh, "NAXIS",  2;  //"number of dimensions";
  fits_set, fh, "NAXIS1", nr;   //"length of 1st dimension";
  fits_set, fh, "NAXIS2", 6;   //"length of 2nd dimension";
  fits_set, fh, "LWAge", LWAge; // LWAge of the best fit
  fits_set,fh,  "LWMet", LWMet; // LWMet of the best fit
  fits_set,fh,  "SSPmodel", basisname; // name of SSP basis used
  fits_set,fh,  "chi2", _ki(1);

  a=array(0.,nr,6);
  a(,1)=x0;
  a(,2)=d;
  a(,3)=pmodel1;
  a(,4)=npec;
  a(,5)=W;
  a(,6)=d-pmodel1;
  
  fits_write_header, fh;          // write header part of current HDU
  fits_write_array, fh, a;    // write data part of current HDU
  fits_close, fh;                   // close stream of FITS handle, the
  
  return [];
};
