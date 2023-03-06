//******************Constantes*******************
G=6.67e-11;
//c=3.e5;
c=299792.458;
func ns(l){
  /* DOCUMENT
     returns the air refraction index for wavelength lambda in Angstrom
  */
  sigma=1000./(l*10.);
  return (1 + 643.28e-7 + 294981.e-7 / (146. - sigma^2) + 2554.e-7 / (41. - sigma^2));
};

func salpIMF(m){
  /* DOCUMENT
     returns the number of star of mass in [m,m+dm] according to Salpeter's IMF, unnormalized, for m in solar masses
     
  */
  return 0.03*m^(-2.35);
  //return 1.;
};

Zsun=0.02; // SOLAR METALLICITY


 




