#include "string.i"

//print,"faire modif arnaud: w=(x1<max(binx))*(x1>=min(binx))";


func lbs(x,..)
/* DOCUMENT lbs 
   returns linear b spline
     in 1,2,3,4 D

     EXAMPLE
     > lbs([-1.,-0.5,0,0.5,1])
     [0,0.5,1,0.5,0]
   SEE ALSO:
 */
{
mm=more_args();
 if (mm==0)
   {
     w=where((abs(x)<=1));
     res=x*0;
     if (numberof(w)) res(w)= (1-abs(x(w)));
     return res;
   }
 else if (mm==1)
   {
     y=next_arg();
     w=where((abs(x)<=1)*(abs(y)<=1));
     res=x*0;
     if (numberof(w)) res(w)= (1-abs(x(w)))*(1-abs(y(w)));
     return res;
   }
 else   if(mm==2)
     {
     y=next_arg();
     z=next_arg();
     w=where((abs(x)<=1)*(abs(y)<=1)*(abs(z)<=1));
     res=x*0;
     if (numberof(w)) res(w)= (1-abs(x(w)))*(1-abs(y(w)))*(1-abs(z(w)));
     return res;
     }
 else  if(mm==3)
     {
     y=next_arg();
     z=next_arg();
     u=next_arg();
     w=where((abs(x)<=1)*(abs(y)<=1)*(abs(z)<=1)*(abs(u)<=1));
     res=x*0;
     if (numberof(w)) res(w)= (1-abs(x(w)))*(1-abs(y(w)))*
                                             (1-abs(z(w)))*(1-abs(z(w)));
     return res;
     }
}






 func histo1d(xx,binx,wght=)
/* DOCUMENT  histo1d(xx,binx)
   histogram binx is a vector defining the bins
 the keyword wght represents the weights to apply to each bincount
warning: check transpose
 */
{ local x1,y1,xi,yi,zi,u,nx,ny;
  x1=xx;
  w=(x1<max(binx))*(x1>=min(binx));
  if (numberof(w))  x1= x1(where(w));
  nx=dimsof(binx)(2);
  if (is_void(x1)) return binx(zcen)*0;
  xi= digitize(x1,binx);
  if (is_void(wght))  res=histogram(xi,top=numberof(binx));
  else res=histogram(xi,wght(where(w)),top=nx);
  res=res(2:);
  return res;
  
}


 
 func histo2d(xx,binx,biny,wght=)
/* DOCUMENT 2D histo2d(xx,binx,biny)
   histogram binx and biny are vectors defining the bins
 the keyword wght represents the weights to apply to each bincount
warning: check transpose
 */
{ local x1,y1,xi,yi,zi,u,nx,ny;
  x1=xx(,1);
  y1=xx(,2);
   w=(x1<max(binx))*(y1<max(biny))*(x1>=min(binx))*(y1>=min(biny));
   if (numberof(w))
     {
       x1= x1(where(w));
       y1= y1(where(w));
     }
  nx=dimsof(binx)(2);
  ny=dimsof(biny)(2);
  if (is_void(x1)) return array(0,[2,nx-1,ny-1]);
  xi= digitize(x1,binx);
  yi= digitize(y1,biny);
  zi=xi+ nx*(yi-1);
  if (is_void(wght))  res=reform(histogram(zi,top=nx*ny),[2,nx,ny]);
  else  res=reform(histogram(zi,wght(where(w)),top=nx*ny),[2,nx,ny]);
  res= res(2:,2:);
  return res;
  
}
/* pp=1000;nn=15;
x=random(pp);y=random(pp); z=random(pp)
u=histo2d([x,y],span(0,1,nn),span(0,1,nn));
pli,u
*/
 
 func histo3d(xx,binx,biny,binz,wght=)
/* DOCUMENT 3D histogram binx biny and binz are vectors defining the bins
 the keyword wght represents the weights to apply to each bincount
   warning: check transpose
 */
{ nx=dimsof(binx)(2);
  ny=dimsof(biny)(2);
  nz=dimsof(binz)(2);
  x1=xx(,1);
  y1=xx(,2);
  z1=xx(,3);
  w=(x1<max(binx))*(y1<max(biny))*(z1<max(binz))*
    (x1>=min(binx))*(y1>=min(biny))*(z1>=min(binz));
  if (numberof(w))
    {
      x1= x1(where(w));
      y1= y1(where(w));
      z1= z1(where(w));
    }
  if (is_void(x1)) return array(0,[3,nx-1,ny-1,nz-1]);
  xi= digitize(x1,binx);
  yi= digitize(y1,biny);
  zi= digitize(z1,binz);
  ui=xi+ nx*(yi-1)+ nx*ny*(zi-1);
  if (is_void(wght))  res=reform(histogram(ui,top=nx*ny*nz),[3,nx,ny,nz]);
  else  res=reform(histogram(ui,wght(where(w)),top=nx*ny*nz),[3,nx,ny,nz]);
  res= res(2:,2:,2:);
  return res;
  
}

/*
#include "random.i"
nn=25;pp=1000000;
x=random_n(pp);y=random_n(pp); z=random_n(pp)
u=histo3d([x,y,z],span(0,1,nn),span(0,1,nn),span(0,1,nn));
pli,u(,,1);
*/


 
 func histo4d(xx,binx,biny,binz,binw)
/* DOCUMENT 4D histogram binx biny  binz and binw are vectors defining the bins
 the keyword wght represents the weights to apply to each bincount
warning: check transpose
 */
{ nx=dimsof(binx)(2);
  ny=dimsof(biny)(2);
  nz=dimsof(binz)(2);
  nw=dimsof(binw)(2);
  x1=xx(,1);
  y1=xx(,2);
  z1=xx(,3);
  w1=xx(,4);
  w=(x1<max(binx))*(y1<max(biny))*(z1<max(binz))*(w1<max(binw))*
    (x1>=min(binx))*(y1>=min(biny))*(z1>=min(binz))*(w1>=min(binw));
    if (numberof(w))
     {
       x1= x1(where(w));
       y1= y1(where(w));
       z1= z1(where(w));
       w1= w1(where(w));
     }
 if (is_void(x1)) return array(0,[4,nx-1,ny-1,nz-1,nw-1]);
  xi= digitize(x1,binx);
  yi= digitize(y1,biny);
  zi= digitize(z1,binz);
  wi= digitize(w1,binw);
  ui=xi+ nx*(yi-1)+ nx*ny*(zi-1)+ nx*ny*nz*(wi-1);
   if (is_void(wght))  res=reform(histogram(ui,top=nx*ny*nz*nw),[4,nx,ny,nz,nw]);
else res=reform(histogram(ui,wght(where(w)),top=nx*ny*nz*nw),[4,nx,ny,nz,nw]);
  res= res(2:,2:,2:,2:);
  return res;
}

/*
nn=25;pp=1000000;
x=random_n(pp);y=random_n(pp); z=random_n(pp);w= random_n(pp);
u=histo4d([x,y,z,w],span(-3,3,nn),span(-3,3,nn),span(-3,3,nn),span(-3,3,nn));
pli,u(,,1,1);
*/


 
 func histo5d(xx,binx,biny,binz,binw,bint)
/* DOCUMENT 4D histogram binx biny  binz and binw are vectors defining the bins
 the keyword wght represents the weights to apply to each bincount
warning: check transpose
 */
{ nx=dimsof(binx)(2);
  ny=dimsof(biny)(2);
  nz=dimsof(binz)(2);
  nw=dimsof(binw)(2);
  nt=dimsof(bint)(2);
  x1=xx(,1);
  y1=xx(,2);
  z1=xx(,3);
  w1=xx(,4);
  t1=xx(,5);
  w=(x1<max(binx))*(y1<max(biny))*(z1<max(binz))*(w1<max(binw))*(t1<max(bint))*
    (x1>=min(binx))*(y1>=min(biny))*(z1>=min(binz))*(w1>=min(binw))*(t1>=min(bint));
  if (numberof(w))
     {
       x1= x1(where(w));
       y1= y1(where(w));
       z1= z1(where(w));
       w1= w1(where(w));
     }
  t1=t1(where(w));
 if (is_void(x1)) return array(0,[5,nx-1,ny-1,nz-1,nw-1,nt-1]);
  xi= digitize(x1,binx);
  yi= digitize(y1,biny);
  zi= digitize(z1,binz);
  wi= digitize(w1,binw);
  ti= digitize(t1,bint);
  ui=xi+ nx*(yi-1)+ nx*ny*(zi-1)+ nx*ny*nz*(wi-1)+ nx*ny*nz*nw*(ti-1);
  if (is_void(wght))res=reform(histogram(ui,top=nx*ny*nz*nw*nt),[5,nx,ny,nz,nw,nt]);
else res=reform(histogram(ui,wght(where(w)),top=nx*ny*nz*nw*nt),[5,nx,ny,nz,nw,nt]);
  res= res(2:,2:,2:,2:,2:);
  return res;
}

/*
nn=25;pp=10000;
x=random_n(pp);y=random_n(pp); z=random_n(pp);w= random_n(pp);t= random_n(pp);
u=histo5d([x,y,z,w,t],span(-3,3,nn),span(-3,3,nn),span(-3,3,nn),span(-3,3,nn),span(-3,3,nn));
pli,u(,,1,1,1);
*/



func dens1d(dat,X,s=) {
/* DOCUMENT 1D  non parametric  histogram X  is the  vector defining the bins
SEE ALSO histo
 */
  local x;
  x= dat;
  if(is_void(s)) s=X(2)-X(1);
  return (1./pi/s/numberof(x)*exp(-(X(,-)-x(-,))^2/s^2))(,sum);
}




func dens2d(dat,X,Y,sx=,sy=) {
/* DOCUMENT 2D non parametric histogram X,Y are the  vector defining the bins
SEE ALSO histo2d
 */
 local x,y,z,nx,px,py; 
 x= dat(,1); y= dat(,2);
  nx= numberof(x);    
  if(is_void(sx)) sx=(X(0)-X(1))/3./nx^(0.25); //sx=X(2)-X(1);
    if(is_void(sy)) sy=(Y(0)-Y(1))/3./nx^(0.25); //sy=Y(2)-Y(1);
   x= dat(,1); y= dat(,2);
  nx= numberof(x);   
  px = numberof(X);  py = numberof(Y);
 z=array(0.,px,py);
   for(i=1;i<=nx;i++){
     z+=exp(-((X(,-:1:py)-x(i))/sx)^2  -((Y(-:1:px,)-y(i))/sy)^2);
 }
 return 1./pi/sx/sy/numberof(dat)*z;
}


func dens3d(dat,X,Y,Z,sx=,sy=,sz=) {
/* DOCUMENT 3D non parametric histogram X,Y,Z are the  vector defining the bins
SEE ALSO histo3d
 */
  local x,y,z,u,nx,px,py,pz; 
  x= dat(,1); y= dat(,2); z= dat(,3);
  nx= numberof(x);   
  if(is_void(sx)) sx=(X(0)-X(1))/3./nx^(0.2); //sx=X(2)-X(1);
    if(is_void(sy)) sy=(Y(0)-Y(1))/3./nx^(0.2); //sy=Y(2)-Y(1);
      if(is_void(sz))sz=(Z(0)-Y(1))/3./nx^(0.2); //sz=Z(2)-Z(1);
  px = numberof(X);  py = numberof(Y);pz = numberof(Z);
 u=array(0.,px,py,pz);
   for(i=1;i<=nx;i++){
     u+=exp(-((X(,-:1:py,-:1:pz)-x(i))/sx)^2  
	    -((Y(-:1:px,,-:1:pz)-y(i))/sy)^2 
	    -((Z(-:1:px,-:1:pz,)-z(i))/sz)^2);
 }
 return 1./pi/sx/sy/sz/numberof(dat)*u;
}




 
func CIC2d(xx,binx,biny)
/* DOCUMENT 2D CIC2d(xx,binx,biny)
   histogram using cloud in cell procedure 
binx and biny are vectors defining the bins
EXEMPLE
#include "random.i"

  pp=1500;nn=20;
x1=random_n(pp);y1=random_n(pp); 
binx=span(-3,3,nn);
biny=span(-3,3,nn);

res= CIC2d([x1,y1],binx,biny);

 warning: check transpose */
{
     x1=xx(,1);
     y1=xx(,2);
 
w=(x1<max(binx))*(y1<=max(biny))*(x1>=min(binx))*(y1>=min(biny));
   x1= x1(where(w));
   y1= y1(where(w));
  nx=dimsof(binx)(2);
  ny=dimsof(biny)(2);
  xi= digitize(x1,binx);
  yi= digitize(y1,biny);
  zi=xi+ nx*(yi-1);

//u=histo2d([x1,y1],binx,biny);

  
aa= array(0.,nx,ny);
bb=aa(*);
ccx=binx(,-:1:ny);
ccy=biny(-:1:nx,);
ccx1=ccx(*);
ccy1=ccy(*);
 llx = binx(dif)(avg);
 lly = biny(dif)(avg); 
	    


 for(i=1;i<=(nx)*(ny);i++)
   {
       ww=where((zi==i) +(zi ==i-1)+(zi==i+1)+(zi==i-nx)+(zi==i+nx));
       if(is_array(ww)) {
         bb(i)=
             lbs((x1(ww)-ccx1(i))/llx,(y1(ww)-ccy1(i))/lly)(sum);
	      }
       else aa(i)=0;
   }
 aa(*)=bb(*);
   aa=aa(2:nx-1,2:ny-1);
   return aa;
   
}

func CIC3d(xx,binx,biny,binz)
/* DOCUMENT  CIC3d(xx,binx,biny,binz)
   histogram using cloud in cell procedure 
binx and biny are vectors defining the bins
 warning: check transpose 
x1=random_n(pp);y1=random_n(pp); z1=random_n(pp); 
binx=span(-3,3,nn);
biny=span(-3,3,nn);
binz=span(-3,3,nn);

res= CIC3d([x1,y1,z1],binx,biny,binz)

*
*/
{
  x1=xx(,1);
  y1=xx(,2);
  z1=xx(,3);
 
w=(x1<max(binx))*(y1<max(biny))*(x1>=min(binx))*
 (y1>=min(biny))*(z1>=min(binz))*(z1<max(binz));
   x1= x1(where(w));
   y1= y1(where(w));
   z1= z1(where(w));
  nx=dimsof(binx)(2);
  ny=dimsof(biny)(2);
  nz=dimsof(binz)(2);
  xi= digitize(x1,binx);
  yi= digitize(y1,biny);
  zi= digitize(z1,binz);
   ui=xi+ nx*(yi-1)+ ny^2*(zi-1);

  
aa= array(0.,nx,ny,nz);
bb=aa(*);
ccx=binx(,-:1:ny)(,,-:1:nz);
ccy=biny(,-:1:nx)(-:1:nz,);
ccz=binz(-:1:ny,)(-:1:nx,);

ccx1=ccx(*);
ccy1=ccy(*);
ccz1=ccz(*);
 llx = binx(dif)(avg);
	    lly = biny(dif)(avg); 
	     llz = binz(dif)(avg); 
	    

 for(i=1;i<=nx*ny*nz;i++)
 {
       ww=where((ui==i)+
  (ui ==i-1)+(ui==i+1)+
 (ui==i-nx)+(ui==i+nx)+
 (ui==i-ny^2)+(ui==i+ny^2));
       if(is_array(ww)) {
         bb(i)=
             lbs(
                 (x1(ww)-ccx1(i))/llx,
                 (y1(ww)-ccy1(i))/lly,
                 (z1(ww)-ccz1(i))/lly
                 )(sum);
  }
       else bb(i)=0;
  }
 aa(*)=bb(*);
   aa=aa(2:nx-1,2:ny-1,2:nz-1);
   return aa;
   
}

func CIC4d(xx,binx,biny,binz,binu)
/* DOCUMENT  CIC4d(xx,binx,biny,binz,binu)
   histogram using cloud in cell procedure 
binx and biny are vectors defining the bins
 warning: check transpose 
x1=random_n(pp);y1=random_n(pp); z1=random_n(pp);  u1=random_n(pp); 
binx=span(-3,3,nn);
biny=span(-3,3,nn);
binz=span(-3,3,nn);
binu=span(-3,3,nn);

res= CIC4d([x1,y1,z1,u1],binx,biny,binz,binu)

*
*/
{
  x1=xx(,1);
  y1=xx(,2);
  z1=xx(,3);
    u1=xx(,4);
 
w=(x1<max(binx))*(y1<max(biny))*(x1>=min(binx))*
 (y1>=min(biny))*(z1>=min(binz))*(z1<max(binz))*
 (u1>=min(binu))*(u1<max(binz));
   x1= x1(where(w));
   y1= y1(where(w));
   z1= z1(where(w));
   u1= u1(where(w));
  nx=dimsof(binx)(2);
  ny=dimsof(biny)(2);
  nz=dimsof(binz)(2);
 nu=dimsof(binu)(2);
  xi= digitize(x1,binx);
  yi= digitize(y1,biny);
  zi= digitize(z1,binz);
ui= digitize(u1,binu);				  
   ui=xi+ nx*(yi-1)+ ny^2*(zi-1)+ nz^3*(ui-1);
aa= array(0.,nx,ny,nz,nu);
bb=aa(*);
ccx=binx(,-:1:ny)(,,-:1:nz)(,,,-:1:nu);
ccy=biny(,-:1:nx)(-:1:nz,)(,,,-:1:nu);
ccz=binz(-:1:ny,)(-:1:nx,,)(,,,-:1:nu);
ccu=binu(-:1:nz,)(-:1:ny,)(-:1:nx,,,);				  

ccx1=ccx(*);
ccy1=ccy(*);
ccz1=ccz(*);
ccu1=ccu(*);				  
 llx = binx(dif)(avg);
	    lly = biny(dif)(avg); 
	     llz = binz(dif)(avg);
  llu = binu(dif)(avg);
	    

 for(i=1;i<=nx*ny*nz*nu;i++)
 {
       ww=where((ui==i)+
  (ui ==i-1)+(ui==i+1)+
 (ui==i-nx)+(ui==i+nx)+
 (ui==i-ny^2)+(ui==i+ny^2)+
 (ui==i-nz^3)+(ui==i+nz^3));
       if(is_array(ww)) {
         bb(i)= lbs((x1(ww)-ccx1(i))^2/llx,
                    (y1(ww)-ccy1(i))/lly,
                    (z1(ww)-ccz1(i))/llz,
 (u1(ww)-ccu1(i))/llu)(sum);
  }
       else bb(i)=0;
  }
 aa(*)=bb(*);
   aa=aa(2:nx-1,2:ny-1,2:nz-1,2:nu-1);
   return aa;
   
}









