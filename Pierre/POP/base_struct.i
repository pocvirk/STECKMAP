struct basStruct {
  string filename;
  double flux(di(1),di(2),di(3)); // spectral basis datacube
  double MsLratio(di(2),di(3)); // store some sort of M/L, right now its more a sort of renormalization rather than a true MsLratio
  double photo(nfilt,di(2),di(3)); // store photometry
  double filtmean(nfilt); // store mean wavelength of filters
  double filters_eff(nfilt); // store effective wavelength of filters
  double filters_range(nfilt,2); // store min and max wl of filters
  double wave(di(1));
  long nages(di(2));
  double met(di(3)); // Z/H .. i.e. total metals content, solar = 0.02
  double FeH(di(3)); // [Fe/H] here solar =0
  double alphaH(di(3)); // alpha elements ratio here solar =0
  double ages(di(2));
  string age_unit;
  double bab(di(2)+1);
  double R(1);
  string basistype;
  string filters(nfilt); // filter names
  
};
