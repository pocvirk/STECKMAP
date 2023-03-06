func angular_separation(a1,d1,a2,d2,deg=)
/* DOCUMENT angular_separation(a1,d1,a2,d2,deg=)
     compute the angular separation between two object.
     a1,d1 : right ascension and declination of first object
     a2,d2 : right ascension and declination of second object

     if deg is not set, angles must be in radian else in degree.
 */
{
  local csep;
  
  if(deg)
    {
      deg2rad=pi/180;
      a1*=deg2rad;
      d1*=deg2rad;
      a2*=deg2rad;
      d2*=deg2rad;
    }
  csep=acos(sin(d2)*sin(d1)+cos(d2)*cos(d1)*cos(a2-a1));
  if(deg) csep=csep/pi*180;
  return csep;
}
