// some tools specific to SDSS

func buildsdssspecname(plate,mjd,fiberid){
  // the basic assumption is that mjd is a 5-digit number, plate is 4 digits and fiberid is 3 digits
  if (mjd>9999.) smjd=pr1(mjd);
  if ((mjd<=9999)&(mjd>999.)) smjd="0"+pr1(mjd);
  if ((mjd>=99)&(mjd<=999.)) smjd="00"+pr1(mjd);

  if ((plate<=9999)&(plate>999.)) splate=pr1(plate);
  if ((plate>99)&(plate<=999.)) splate="0"+pr1(plate);
  
  if ((fiberid<=999)&(fiberid>99.)) sfiberid=pr1(fiberid);
  if ((fiberid>9)&(fiberid<=99.)) sfiberid="0"+pr1(fiberid);
  if ((fiberid<=9.)) sfiberid="00"+pr1(fiberid);

  return "spSpec-"+smjd+"-"+splate+"-"+sfiberid+".fit";
};

