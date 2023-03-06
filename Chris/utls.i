#include "histon.i"
#include "fits.i"
#include "ascii.i"
#include "utils.i"
#include "fft_utils.i"
#include "system.i"
#include "xplot.i"
__SetNameOfColors=1;
//#include "color_def.i"
#include "copy_plot.i"
#include "graphedit.i"
//#include "Chris/scatter3d.i"

cpp=copyWin;


if (!batch() & (__CHRISTOPHE__==1)) write,exec("date");
epsi=eps;// rename eps function since I use my own ...
// epsi has no bounding box pb...


   
func Colors(n) { return bytscl(span(1,n+1,n+1),cmin=0.5 ,cmax=n+1.5);}
func randCol(void) { return Colors(25)( 2+long(random()*15));}

func fft_indgenNyquist(nn) { tt= fft_indgen(nn); tt(nn/2+1)=0.;  return tt;}


func plv3(x0,y0,z0,bb,alt=,az=,width=,range=,color=,colorscheme=,type=,size=)
/* DOCUMENT
      plot color 3d vector fields;
   SEE ALSO:
 */
{
  if(is_void(width)) width=3;
  if(is_void(range)) range=1:dimsof(x0)(3)/2;
  if(is_void(color)) color=__green4;
  if(is_void(colorscheme)) colorscheme="length";
  if(is_void(type)) type=1;
  if(is_void(size)) size=0.125;

  _pl3xy, x0p1, y0p1,(x0)(,range,)(*),(y0)(,range,)(*),(z0)(,range,)(*);
  _pl3xy, x0p2, y0p2,(x0+bb(..,1))(,range,)(*),
    (y0+bb(..,2))(,range,)(*),(z0+bb(..,3))(,range,)(*);
  pldj,x0p1,y0p1,x0p2,y0p2,width=width,color=color,type=type;
  if(colorscheme=="length") cols=sqrt(abs2(bb(,range,,))(..,sum))(*);
  if(colorscheme=="depth") cols=y0(,range,)(*);
  plcolor,cols,y0p2(*),x0p2(*),   bar=0,size=size,npol=5;
}





func hhelp(kword,browser=)
/* DOCUMENT 
     launches html doc in yorick
     only works for yorick keywords for now
   SEE ALSO:
 */
{
  if(typeof(kword)!="string") kword=strpart(splittok(print(kword)(1))(0),:-2);
  k1=strcut(kword,1)(1);
  if (is_void(browser)) browser="";
  if(numberof(exec("ps -aux | grep mozilla"))>1)  system,"/usr/local/mozilla/mozilla -remote 'openURL(file:/home1/pichon/Yorick/doc/html_xref/global-"+k1+"s.html#"+kword+")'";  else  system,"/usr/local/mozilla/mozilla& ";
}




func pshft(x)
/* DOCUMENT 
     shifts plot by fraction x
   SEE ALSO:
 */
{
  tt=limits();
  x1=tt(1);
  x2=tt(2);
  limits,x1+(x2-x1)*x,x2+(x2-x1)*x,tt(3),tt(4);
}


func inform(i,n,msg,freq=)
/* DOCUMENT
   returns the progress of i out of n only if i=n*freq;
   FREQ= 10 % by default
   EXAMPLE
     inform inform(10,100)
   SEE ALSO:
 */
{
  if(is_void(freq)) freq=0.1; // 10 % 
  if(is_void(msg)) msg=""; 
  if((i%max(1,long(n*freq)))==0) write,"done "+msg+":"+pr1(long(i*100./n))+"\%";
  return;
}

func rootname(fname,name=,noext=)
/* DOCUMENT 
     return the directory name of a given file
     if NAME=1 it returns only the file name;
     EXAMPLE
rootname("../../tata",name=1)
"tata"

     SEE ALSO:
 */
{

  if (numberof(fname)>1){
    res=[]; for(i=1;i<numberof(fname);i++)  grow,res,rootname(fname(i),name=name,noext=noext); return res;
  }
  if(name && noext) return strjoin(splittok(splittok(fname,tok="/")(0),tok=".")(1:-1),".");
  if(name) return splittok(fname,tok="/")(0);
  if(strpart(fname,1:1)=="/")
    root="/"+strjoin(splittok(fname,tok="/")(1:-1),"/")+"/";
     else
         if(numberof(splittok(fname,tok="/"))>1)
           root=strjoin(splittok(fname,tok="/")(1:-1),"/")+"/";
         else root="./";
  return root;
     
}



func intersectIdx(lst1,lst2)
/* DOCUMENT 
     returns the index in lst2 of points in lst1;
   SEE ALSO:
 */
{
 return where2(abs(lst1(,-)-lst2(-,))==0)(1,);
}


func strIntersect(strlst1,strlst2)
/* DOCUMENT
      returns the intersect table of strlst1 and strlst2
      which are lists of strings
   EXAMPLE
   where(strIntersect(["a","b","c"],["b","c"])(sum,));
   SEE ALSO:
 */
{
LL=[];for(i=1;i<=numberof(strlst1);i++) grow,LL,[strmatch(strlst2,strlst1(i))];
 return LL; 
}



func __pidrawtick(x,y,up,width=,yaxis=)
{
  if(is_void(width)) width=2;
  if(is_void(yaxis))
    {
      pldj, x, y+0.015*up, x, y, width=width+1;
      pldj, x(:-1)+(x(2)-x(1))/2.,   y(:-1)+0.01*up  ,
        x(:-1)+(x(2)-x(1))/2.  , y(:-1), width=width;
      pldj, x(:-1)+(x(2)-x(1))*3/4., y(:-1)+0.0075*up,
        x(:-1)+(x(2)-x(1))*3/4., y(:-1), width=width;
      pldj, x(:-1)+(x(2)-x(1))*1/4., y(:-1)+0.0075*up,
        x(:-1)+(x(2)-x(1))*1/4., y(:-1), width=width;
    } else
      {//==== FIXME 
        pldj, x, y+0.015*up, x, y, width=width+1;
        pldj, x(:-1)+(x(2)-x(1))/2.,   y(:-1)+0.01*up  ,
        x(:-1)+(x(2)-x(1))/2.  , y(:-1), width=width;
        pldj, x(:-1)+(x(2)-x(1))*3/4., y(:-1)+0.0075*up,
          x(:-1)+(x(2)-x(1))*3/4., y(:-1), width=width;
        pldj, x(:-1)+(x(2)-x(1))*1/4., y(:-1)+0.0075*up,
          x(:-1)+(x(2)-x(1))*1/4., y(:-1), width=width;
      }
        
}




func encapsulatePlot(small,big,pos=)
/* DOCUMENT 
     adds interactively
     a plot in a box "small" into the frame "big";
     EXAMPLE
     encapsulatePlot,1,0,pos=[0.2,0.4,0.5,0.7]
   SEE ALSO:
 */
{
require, "Chris/draw.i";

   
  window,big;
  plsys,1;
  ll=limits();
  savePlot,"/tmp/aplot.gdb",big;
  draw_window;
  if(is_void(pos))  draw_manip; else draw_new,[pos(1),pos(2),pos(3),pos(4)];
  loadPlot,"/tmp/aplot.gdb",15;
  replot1Sys,15,1,big,1;
  window,big; plsys,1; limits,ll(1),ll(2),ll(3),ll(4);
  replot1Sys,small,1,big,2;
  window,small;
  ll=limits();
  window,big;
  plsys,2;limits,ll(1),ll(2),ll(3),ll(4);
  redraw;
  wk,15;
}



func pitick(off,asys=,height=,width=,notxt=,lessticks=,yaxis=)
/* DOCUMENT 
 replaces x label by pi fractions only works for
 [0,pi],[0,2pi] [-pi,pi] [-pi/2,pi/2]

      EXAMPLE
      x=span(-pi/2,pi/2,100);ws,0;plh,sin(x^2),x; limits;pitick;
      x=span(0,pi,100);ws,1;plh,sin(x^2),x; limits;pitick;
      x=span(0,2*pi,100);ws,2;plh,1+sin(x^2),x; limits;pitick;
      x=span(-pi,pi,100);ws,3;plh,1+sin(x^2),x; limits;pitick;

      SEE ALSO:
 */
{
   require, "style.i";
   local land, sys, legs, clegs;

   if (!is_void(asys)) plsys,asys;
   if (is_void(height)) height=24;
   get_style, land, sys, legs, clegs;
   if (off) {
     if (is_void(_retick_style)) {
       _retick_style = sys;
     } else {
       set_style, land, _retick_style, legs, clegs;
       _retick_style = [];
     }
   } else {
     if (is_void(_retick_style)) _retick_style = sys;
     nsys = plsys();
     if (!nsys) error, "no current coordinate system, use plsys";
     /* remove ticks, labels from bottom and left edges, see g/work.gs */
          sys(nsys).ticks.horiz.flags &= ~(0x001 | 0x020);
          sys(nsys).ticks.horiz.flags &= ~(0x002 | 0x020);

          if(yaxis)
            {
              sys(nsys).ticks.vert.flags &= ~(0x001 | 0x020);
              sys(nsys).ticks.vert.flags &= ~(0x002 | 0x020);
            }
          //     sys(nsys).ticks.vert.flags &= ~(0x001 | 0x020);
     set_style, land, sys, legs, clegs;
     lims = limits();
     plsys,asys;
     vp = viewport();
     plsys, 0;
     /* actually want to use lims, vp to compute correct tick and
      * label locations, don't bother for this example */

     if(lims(1)<-3)
       {
         limits,-pi,pi; 
         x = span(vp(1), vp(2), 5);
         y = 0.*x + vp(4);
         __drawtick,x,y,-1,width=width;
         y = 0.*x + vp(3);
         __pidrawtick,x,y,1,width=width;
        
 if(is_void(notxt)) plt1, ["-!p","-!p/2","0", "!p/2", "!p"], x, y-0.01, justify="CT",
           font="times", height=height;
       }
     else if(lims(1)<-1)
       { 
         x = span(vp(1), vp(2), 5);
         y = 0.*x + vp(4);
         __pidrawtick,x,y,-1,width=width;
         y = 0.*x + vp(3);
         __pidrawtick,x,y,1,width=width;
if(is_void(notxt))          plt1, ["-!p/2","-!p/4","0", "!p/4", "!p/2"], x, y-0.01, justify="CT",
           font="times", height=height;
 if(yaxis)
           {
        limits,-pi/2,pi/2,-pi/2,pi/2; 
         y = span(vp(1), vp(2), 5);
         x = 0.*x + vp(2);
         __pidrawtick,x,y,-1,width=width,yaxis=1;
         x = 0.*x + vp(1);
         __pidrawtick,x,y,1,width=width;
if(is_void(notxt))          plt1, ["-!p/2","-!p/4","0", "!p/4", "!p/2"], x, y-0.01, justify="CT",
           font="times", height=height;     
           } 

       }
     else if(lims(2)<4)
       {
         limits,0,pi; 
         x = span(vp(1), vp(2), 5);
         y = 0.*x + vp(4);
         __pidrawtick,x,y,-1,width=width;
         y = 0.*x + vp(3);
         __pidrawtick,x,y,1,width=width;
         if(is_void(notxt))  if(is_void(lessticks))
           plt1, ["0", "!p/4", "!p/2", "3!p/4", "!p"],
                               x, y-0.01, justify="CT",
           font="times", height=height;
         else plt1, ["0",  "!p/2",  "!p"],
                               x(::2), y(::2)-0.01, justify="CT",
           font="times", height=height;
       }
     else
       {
         limits,0,2*pi;
         x = span(vp(1), vp(2), 5);
         y = 0.*x + vp(4);
         __pidrawtick,x,y,-1,width=width;
         y = 0.*x + vp(3);
         __pidrawtick,x,y,1,width=width;
if(is_void(notxt))          plt1, ["0",  "!p/2",  "!p","2!p/3", "2!p"], x, y-0.01, justify="CT", font="times", height=height;
       }
     plsys, nsys;
   }
}








func inheritStruct(astruct,_new_struct_,savefields=,newfields=,debug=)
/* DOCUMENT
   exports a structure from an existing one;
   takes arrays of struct as arguments and does the inheritence of   the fields
   savefied= gives the index of the fields to retain
   new fields= gives a string describing the newfields
   EXAMPLE
   struct astruct{ int a; int b;}; qin=astruct(); qin.a=1; qin.b=2;
  inheritStruct([qin,qin],"newstruct",savefields=[1,2],newfields="double ss; int q;");
   SEE ALSO:
 */
{
  names=struct2name(astruct); // field name of structs
  if(is_void(savefields)) savefields=indgen(numberof(names));
  if(is_void(newfields)) newfields="";
  infields=names(savefields);
 types=struct2type(astruct)(savefields); // field type of struct
 name=transpose([types,infields,array(";",numberof(types))])(*);
 eval, "struct "+_new_struct_+ "{"+strjoin(name," ")+newfields+"}",debug=debug;
 _new_struct_ = symbol_def(_new_struct_);

 n = numberof(astruct);
 m = numberof(infields);
 q = array(_new_struct_, dimsof(astruct));
 r = _new_struct_();
 for (i=1;i<=n;++i) { // loop on array of struct
   t = astruct(i);
   for(j=1;j<=m;++j) { // loop on fields
     s = infields(j);
     get_member(r, s) = get_member(t, s); // carry identification
   }
   q(i) = r;
 }
 return q;
}




func replaceString(oldstring,newstring,filename,noconfirm=,noback=)
/* DOCUMENT
   
   EXAMPLE 
   replaceString("1.0\\)","1.5\\)","spectras.ps")
   SEE ALSO:
 */
{
   if(is_void(noconfirm))
    {
     write,"will do perl -pi.bak -e 's/"+oldstring+"/"+newstring+"/g' `ls -1 "+filename+"`";
  aaa="n";
      nn=read(prompt="confirm ?",aaa);
      if(aaa!="y") return;
    }
if(is_void(noback))  system, "perl -pi.bak -e 's/"+oldstring+"/"+newstring+"/g' `ls -1 "+filename+"`";
 else  system, "perl -pi -e 's/"+oldstring+"/"+newstring+"/g' `ls -1 "+filename+"`";
}


func rescaleboundingbox(fname,factor=)
/* DOCUMENT 
     rescaleboundingbox will try and
     make sure the labels are ok,factor=)
   SEE ALSO:
 */
{
  if (is_void(factor)) factor=0.1;
  input=exec("grep %%PageBoundingBox "+fname)(0);
  tt=str2nb(splittok(input,tok=" ")(2:),0.);
  output=splittok(input,tok=" ")(1)+" "+
    strjoin(pr01([tt(1:2)*(1-factor)*[0.9,1.05],tt([3,4])](*),ndigit=0)," ");
  replaceString(input,output,fname,noconfirm=1,noback=1);
  input=exec("grep %%BoundingBox "+fname)(0);
  output=splittok(input,tok=" ")(1)+" "+
    strjoin(pr01([tt(1:2)*(1-factor)*[0.9,1.05],tt([3,4])](*),ndigit=0)," ");
  replaceString(input,output,fname,noconfirm=1,noback=1);
 replaceString(input,output,fname,noconfirm=1,noback=1);
}


func circul(x)
/* DOCUMENT 
     returns circulent matrix corresponding to
     vector x
   SEE ALSO:
 */
{
  res=[];
  n=dimsof(x)(0);
  for(i=1;i<=n;i++)
    {
      grow,res,[roll(x,i-1)];
    }
  return res;
}


func incl(fname)
/* DOCUMENT 
     includes the file and displays functions
   SEE ALSO:
 */
{
  include,fname;
  showfun,fname;
  return;
}



func showfit(fname,aladin=)
/* DOCUMENT 
     display fits using panoramix  or aladin
   SEE ALSO:
 */
{
  require,"Eric/fits2.i";
 if(dimsof(fname)(0)>1)
   {
     fits_write,"/tmp/tmp.fits",fname,overwrite=1;
     fname="/tmp/tmp.fits";
   }
if(is_void(aladin)) system,"/usr/local/bin/panoramix "+ fname+" &";
else  system,"/usr/local/bin/AladinJava "+ fname+" &";
}


func cplim(src, dst)
/* DOCUMENT cplim, src..., dst...;
     Copy plot limits from source graphic window SRC to destination graphic
     window DST.  SRC and DST are  both specified by keyword and default to
     the current graphic window if not set.

   SEE ALSO: limits, window. */
{
  if ((cur = current_window()) < 0) return; /* no current window */
  if (is_void(dst)) dst = cur;
  if (is_void(src)) src = cur;
  if (dst == src) return;
  window, src;
  lim = limits();
  window, dst;
  limits, lim;
  window, cur;
}



func circle2seg(x0, y0, r, dim, rnd=,periodic=)
 /* DOCUMENT returns a mask describing the segment of length
    r at  centre (x0,y0) in units of dim

   EXAMPLE mask= circle2seg(256*random(64),256*random(64),256*random(64)/15.,256)
   SEE ALSO:
 */
{
  x = indgen(dim);
  y = indgen(dim)(-,);
  n = numberof(x0);
  mask = array(0n, dim, dim);
  r2 = (rnd ? (r + 0.5)^2 : r*r);
  for (i=1 ; i<=n ; ++i) {
    if(is_void(periodic))
      {
        u = x - x0(i);
        v = y - y0(i);
      }
    else
      {
        u = min(abs(x - x0(i)),dim-abs(x-x0(i)));
        v = min(abs(y - y0(i)),dim-abs(y-y0(i)));
      }
   
  }
  return mask;
}


func circle2mask(x0, y0, r, dim, rnd=,periodic=,seg=)
/* DOCUMENT returns a mask describing the inner region of cicles
   with radii r at  centre (x0,y0) in units of dim
   keyword seg will return line segment corresponding to the width of
   the intersection of the circle centered at the center.
   EXAMPLE mask= circle2mask(256*random(64),256*random(64),256*random(64)/15.,256)
   SEE ALSO:
 */
{
  x = indgen(dim);
  y = indgen(dim)(-,);
  n = numberof(x0);
  mask = array(0n, dim, dim);
  r2 = (rnd ? (r + 0.5)^2 : r*r);
  for (i=1 ; i<=n ; ++i) {
if(is_void(periodic))
  {
    u = x - x0(i);
    v = y - y0(i);
  }
 else
  {
    u = min(abs(x - x0(i)),dim-abs(x-x0(i)));
    v = min(abs(y - y0(i)),dim-abs(y-y0(i)));
  }
    if(is_void(seg))
      mask |= u*u + v*v < r2(i); // cicles
    else 
      mask |= (u*u  < r2(i))*(v*v<0.5); // segments
  }
  return mask;
}



func sphere2mask(x0, y0,z0, r, dim, rnd=,periodic=)
/* DOCUMENT returns a mask describing the inner region of cicles
   with radii r at  centre (x0,y0,z0) in units of dim
   rd= sets the rounding at the edge
   periodic= assumes periodic conditions
   EXAMPLE mask= sphere2mask(64*random(64),64*random(64),64*random(64),64*random(64)/15.,64)
   SEE ALSO:
 */
{
  x = indgen(dim);
  y = indgen(dim)(-,);
  z = indgen(dim)(-,-,);
  n = numberof(x0);
  mask = array(0n, dim, dim,dim);
  r2 = (rnd ? (r + 0.5)^2 : r*r);
  for (i=1 ; i<=n ; ++i) {
if(is_void(periodic))
  {
    u = x - x0(i);
    v = y - y0(i);
    w = z - z0(i);
  }
 else
  {
    u = min(abs(x - x0(i)),dim-abs(x-x0(i)));
    v = min(abs(y - y0(i)),dim-abs(y-y0(i)));
    w = min(abs(z - z0(i)),dim-abs(z-z0(i)));
  }
    mask |= u*u + v*v + w*w < r2(i);
  }
  return mask;
}



func pr01(x,ndigit=)
/* DOCUMENT like pr1 but with heading zeros
   if ndigit < 0 it will use scientific notation beyond n pts
   SEE ALSO: pr1
 */
{
  if(is_void(ndigit)) ndigit=2;
  if(ndigit<0) return swrite(format="%0."+pr1(abs(ndigit))+"g",x);
  if(is_integer(x))  return swrite(format="%0"+pr1(ndigit)+"d",x);
  if(is_real(x))  return swrite(format="%0."+pr1(ndigit)+"f",x);
}


func splitplot(y,x,cut=,color=,width=,type=,noerase=,rg=,nmax=)
/* DOCUMENT 
     splitplot(y,x,cut=,color=,width=,type=,noerase=)
   will split over numberof(y)/cut subplot the curve y(x)
   (at most 7x4 plots)
   SEE ALSO:
 */
{
  if(is_void(cut)) cut=500;
  if(is_void(nmax)) nmax=20;
  nn=numberof(y);
 if(is_void(x)) x=indgen(nn);
   qq=min(long((nn/cut/4))+1,nmax);
  if(nn>20*4*cut) write,"WARNING: the plot is truncated\n";
  for(j=0;j<=qq-1;j++)
    {
	if(j*4*cut==nn) continue;
     if(is_void(noerase)) window,j,style="win41.gs"; else   window,j;
         for(i=0;i<=min(4,nn/cut-j*4)-1;i++)
        {
          plsys,4-i;
      if((i+j*4)*cut==min(((i+j*4)+1)*cut,nn)) continue;
    	  plh,y((i+j*4)*cut+1:min(((i+j*4)+1)*cut,nn)),
      x((i+j*4)*cut+1:min(((i+j*4)+1)*cut,nn)),color=color,width=width,type=type;
  if(!is_void(rg)) range,rg(1),rg(2)   ;
  
  }
	 redraw; 
    }
    
}



func cube2DVR(mtx)
/* DOCUMENT  cube2DVR(mtx) converts the cube into an
   avf file compatible with dvr and displays it
   SEE ALSO: cube2avf
...      
*/
{
cube2avf("~/.tmp",mtx);
system,"rsh pleiades12 dvrinit &";
system,"xon pleiades12 'dvrf -file ~/.tmp.xvf -script $DVRHOME/runtime/startAndOpenFile.tcl &' & ";
}


func cube2avf(mfile,mtx)
/* DOCUMENT  cube2avf(mfile,mtx) converts the cube into an
   avf file compatible with dvr
   SEE ALSO: cube2DVR
...      
*/
{
file=open(mfile+".avf","w");
if(typeof(mtx)!="double") mtx=double(mtx);
nd = dimsof(mtx);
if(nd(1) !=3) error,"wrong dimension: should be a cube";
write,file,format="WIDTH %d\n",nd(2);
write,file,format="HEIGHT %d\n",nd(3);
write,file,format="SLICES %d\n",nd(4);
write, file,format="%s\n","FRAMES  1";
write,file,format="MIN %f\n",min(mtx);
write,file,format="MAX %f\n",max(mtx);
write, file,format="%s\n","FORMAT  SCALAR8    # 8 bit data";
write, file,format="%s\n","XDIST 1.0";
write, file,format="%s\n","YDIST 1.0";
write, file,format="%s\n","ZDIST 1.0";
write, file,format="%s\n","TIME  1.0";
write,file,format="%f\n",(mtx)(*),linesize=15;
 close ,file;
 ll=exec("ls "+mfile+".xvf");
 if(numberof(ll)==1) system,"rm "+mfile+".xvf";
 system,"rsh pleiades12  vconv "+mfile+".avf "+mfile+".xvf";
 system,"rm "+mfile+".avf";
 // system," xon pleiades12 dvrf -file "+mfile+".xvf"
return mfile+".xvf";
}


func lst2arr(x)
/* DOCUMENT  lst2arr(x) converts the list into an array when applicable
   SEE ALSO: arr2lst
*/
{ local res;
  res=[];
  do {
    res=grow(res, _car(x));   /* recurse */
    x= _cdr(x);
  } while (!is_void(x));
  return res;
}


func arr2lst(x)
/* DOCUMENT      arr2lst(x) converts the  array into a list
   SEE ALSO: lst2arr
*/
{
  nn=dimsof(x)(0);
  lst=[];
  for(i=1;i<=nn;i++)
    lst=_cat(lst,_lst(x(,i)));
  return lst;
}




func splitstat(xx,binx,&wghts,xidx=,yidx=,dummy=,ywghts=)
/* DOCUMENT  splitstat([x,y,...],binx)
   return list of y's falling in binx i
   dummy contains the number to return for empty bins
   e.g.
   pp=1000;nn=15;
   x=random(pp);y=random(pp);
   u=splitstat([x,y],span(0,1,nn));
   u1=lst2arr(_map(kurtosis,u));

   will give  the kurtosis of y for each bin of binx

   tt=splitstat([x,y],span(0,1,15),wghts,ywghts=dy);
   will also put in wghts the list of dy of points which fall in each bin.
   
   SEE ALSO: pler
*/
{ local x1,y1,xi,w;
  wws=wghts=[]; 
 if (is_void(xidx)) xidx=1; 
 if (is_void(yidx)) yidx=2; 
 if (is_void(dummy)) dummy=0.; else dummy=double(dummy);
 
  x1=xx(,xidx);
  y1=xx(,yidx);
  xi= digitize(x1,binx);
  res=[];
  for(i=1;i<=numberof(binx)+1;i++){
    w=where(xi==i);
    if(numberof(w)) {
      res=_cat(res,y1(w));
      //wws=_cat(wws,w);
      if(numberof(ywghts))      wghts=_cat(wghts,ywghts(w));
    } else {
      res=_cat(res,[dummy]);
      //wws=_cat(wws,void);
      if(numberof(ywghts))      wghts=_cat(wghts,[0.]);
    }
  }
  
  return res;
  
}



func pow(a,b) {
/* DOCUMENT  pow(a,b) returns a^b
     SEE ALSO: sqr
*/
  return a^b;}

func showfun(fname,v=) {
/* DOCUMENT  showfun(fname,v=)  returns list of function
     defined in fname
     SEE ALSO: HELP,about
*/
  require, "Eric/system.i";
  if (numberof(strtok(fname,"/"))<2) fname=get_cmd()+fname;
  if (strtok(fname,".")(0)!="i") fname=fname+".i";
 if(is_void(v))  ll =exec("\grep -hi 'func ' "+fname); else  ll =exec("\grep  -in 'func ' "+fname+ " ~/Yorick/"+fname);
  if (numberof(ll) == 1000) ll="NO FUNC FOUND";
  ll2=[];
  for(i=1;i<=numberof(ll);i++) grow,ll2,splittok(splittok(ll(i))(2),tok="(")(1);
  if(is_func(about))
    {
      //     include,fname,1;
      ll=select_item_in_string_list(ll2);
      help,symbol_def(ll); info,ll;
     } else write,format="%s %s %s %s \n ",linesize=40,ll2;
}



func OpenSpeedWindow(nb,win=,dpi=)
{
 require,"Eric/jgraph.i";

  local height,width,vwpl,yp,vwp;
  local xopt,yopt,s1,s2,s3,s4;

  if(is_void(win))
    win=0;
  if(is_void(dpi))
    dpi=75;

  height=long(11.04*dpi)+1;
  width=long(8.54*dpi)+1;
  vwpl=[0.08,0.72,0.03,1.01];
  yp=span(vwpl(3),vwpl(4),nb+1);
  vwp=array(double,[2,4,nb]);
  vwp(1,)=vwpl(1);
  vwp(2,)=vwpl(2);
  vwp(3,)=yp(1:-1);
  vwp(4,)=yp(2:0);

  xopt=array("tTiz",nb);
  xopt(1)="tTilz";
  yopt=array("tTil",nb);

  JWindow,win,width=width,height=height,color=254,viewport=vwp,xopt=xopt(*),yopt=yopt(*),xmargin=0,ymargin=0;
  get_style,s1,s2,s3,s4;
  s2(1).ticks.horiz.labelOff=0.006;
  s2(*).ticks.vert.labelOff=0.008;
  set_style,s1,s2,s3,s4;
  return;
}                    



func mod(x,n) {
/* DOCUMENT mod(x,n)  returns x modulus n
     SEE ALSO: %
*/
  return long(x-(x/n)*n);}

func ydoc { system, "netscape /home1/pichon/Yorick/doc/index.html &"; return 0;}


func yref { system, "gv -scale 2 -seascape  /usr/local/share/yorick/1.5/doc/refs.ps &"; return 0;}

func yman { system, "gv   /usr/local/share/yorick/1.5/doc/yorick.ps &"; return 0;}

func c_v { palette, "idl-15.gp"; return 0;}
func c_q { palette, "idl-17.gp"; return 0;}
func c_b { palette, "idl-01.gp"; return 0;}
func c_r { palette, "idl-03.gp"; return 0;}
func c_g { palette, "idl-08.gp"; return 0;}
func c_p { palette, "idl-07.gp"; return 0;}
func c_k { palette, "white-purple.gp"; return 0;}
func c_w { palette, "white-green.gp"; return 0;}
func c_e { palette, "earth.gp"; return 0;}
func c_a { palette, "idl-13.gp"; return 0;}
func c_n { palette, "yarg.gp"; return 0;}
func c_y { palette, "gray.gp"; return 0;}
func c_s { palette, "stern.gp"; return 0;}

func c_qq { palette, "ridl-17.gp"; rvp; return 0;}
func c_vv { palette, "ridl-15.gp"; return 0;}
func c_bb { palette, "ridl-01.gp"; return 0;}
func c_rr { palette, "ridl-03.gp"; return 0;}
func c_gg { palette, "ridl-08.gp"; return 0;}
func c_pp { palette, "ridl-07.gp"; return 0;}
func c_aa { palette, "ridl-13.gp"; return 0;}
func c_ss { palette, "stern.gp"; rvp; return 0;}



func splittok(s,tok=)
/* DOCUMENT 
      splittok(s,tok=)
      splits string s into list using tok as separator default ""
   SEE ALSO:
 */
{
  local r;
  if(is_void(tok)) tok=" ";
  r=[];
  do {s=strtok(s,tok);grow,r,s(1);s=s(2);} while(s);  
  return r;
}
  

func _filterdot(s,tok=)
/* DOCUMENT _filterdot(s) replaces "." by  tok (default ",") in string
     
   SEE ALSO:
 */
{
  tt=splittok(s,tok=".");
  if (is_void(tok)) tok =",";
  if(numberof(tt)>2) tt=strjoin(tt(1:-1),tok)+"."+tt(0); else tt=s;
  return tt;
}


func exportfig(ll,fname=,pp=,ps=)
/* DOCUMENT exportfig(ll,fname=,pp=,ps=) creates a latex file (default crap.tex)
   with a table containing the figures in the list ll
   EXAMPLE
   #include "Eric/system.i"
   ll =exec("ls clusters/*.jpg");
   exportfig(ll(1:21));
   or
   exportfig(ll(1:21),pp=3);      
   KEYWORDS
   fname = tex file name
   pp = number of fig per line (default 4 )
   SEE ALSO: dumppdf,eps,pdf,jpg
*/
{
if (is_void(fname))  fname="crap.tex";
  ff =open(fname,"w")

    nn= numberof(ll);
  //  for(i=1;i<=nn;i++) ll(i)=_filterdot(ll(i));
  if (is_void(pp))  pp=4;
 if(pp==3){ qq=9; q1=0.40; q2=0.25;}
 if(pp==4){ qq=24; q1=0.28; q2=0.15;}
 if(pp==2){ qq=6; q1=0.6; q2=0.32;}
 
write,ff,"\\documentclass{article}";
write,ff,"\\usepackage[T1]{fontenc}";
write,ff,"\\usepackage[latin1]{inputenc}";
write,ff,"\\usepackage{graphics}";
write,ff,"\\textwidth=16.5cm";
write,ff,"\\textheight=27.5cm"
write,ff,"\\topmargin=-3.cm";
write,ff,"\\evensidemargin=-2.5cm";
write,ff,"\\oddsidemargin=-2.5cm";    
write,ff,"\\begin{document}";

 for(k=0;k<=nn/qq;k++)
   {write,ff,"{\\centering \\begin{tabular}{c c c c}";
     tt=((k<nn/qq)?qq:(nn%qq));
 for(i=1;i<=tt ;i++)
   { 
     if(is_void(ps)) tt1=strpart(ll(i+k*qq),1:-4); else tt1=ll(i+k*qq);
  write,ff,
        format="\\resizebox*{%2f\\columnwidth}{%2f\\textheight}{\\includegraphics{%s}} ",
                       q1,q2,tt1;
      if (i % pp ) {write,ff,"& ";} else if (i<nn) {  write,ff,"\\\\  ";} else {write,ff,"";} 
    }

write,ff,"\\end{tabular}\\par}";
   }
 write,ff,"\\vspace{0.3cm}";
write,ff,"\\end{document}";
close,ff;
    
 return fname;
}



func xmgrshow(x,y)
/* DOCUMENT xmgrshow(x,y)uses xmgr to display cluster of points
   SEE ALSO:
*/
{
  if (structof(x)!=double) x= (*x.pos);
  ff=open("~/.crapxmgr.dat","w");
  write,ff,format="%f %f \n",y,x;
  close,ff;
  system,"xmgr ~/.crapxmgr.dat &";
}



func __xgobcol(x)
{
  y=array("white",numberof(x));
  w1= where(x==1);
  w2= where(x==2);
  w3= where(x==3);
  w4= where(x==4);
  w5= where(x==5);
  w6= where(x==6);
if (numberof(w1)>0)  y(w1)=array("red",numberof(w1));
if (numberof(w2)>0)  y(w2)=array("blue",numberof(w2));
if (numberof(w3)>0)  y(w3)=array("orange",numberof(w3));
if (numberof(w4)>0)  y(w4)=array("lightblue",numberof(w4));
if (numberof(w5)>0)  y(w5)=array("yellow",numberof(w5));
if (numberof(w6)>0)  y(w6)=array("purple",numberof(w6));
 return y; 
}


func xgobishow(x,range=,center=,cflag=,glyphs=)
/* DOCUMENT uses xgobi to display cluster of points
   xgobishow(x,range=,center=,cflag=,glyphs=)

   glyphs correspond to the symbol style coded as 
   1 through 5: Five sizes of '+'
   6 through 10: Five sizes of 'X'
   11 through 15: Five sizes of open rectangle
   16 through 20: Five sizes of filled rectangle
   21 through 25: Five sizes of open circle
   26 through 30: Five sizes of filled circle
   31: A single-pixel point              
   EXAMPLE
   x=random(150);y=random(150);z=random(150);
   xgobishow,transpose(100*[x,y,z]),glyphs=long(x*10)+10

   SEE ALSO:
 */
{
  if (typeof(x)=="struct_instance" ) x= (*x.pos); else x =double(x);
  
  ff=open("~/.crapxgobi.dat","w");
 ffc=open("~/.crapxgobi.dat.colors","w");
 ffg=open("~/.crapxgobi.dat.glyphs","w");   
 // if (dimsof(x)(3)> dimsof(x)(2))  x=transpose(x);
 x=double(x);
  if( dimsof(x)(2)==2)  write,ff,x,format= "%8.4e %8.4e  \n",linesize=30;
  if( dimsof(x)(2)==3)  write,ff,x,format= "%8.4e %8.4e %8.4e \n",linesize=40;
  if( dimsof(x)(2)==4)  write,ff,x,format= "%8.4e %8.4e %8.4e %8.4e \n",linesize=50;
  if( dimsof(x)(2)==5)  write,ff,x,format= "%8.4e %8.4e %8.4e %8.4e %8.4e \n",linesize=60;
  if( dimsof(x)(2)==6)  write,ff,x,format= "%8.4e %8.4e %8.4e %8.4e %8.4e %8.4e \n",linesize=70;
  if( dimsof(x)(2)==7)  write,ff,x,format= "%8.4e %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e \n",linesize=80;
    if (!is_void(cflag)) {
      if (numberof(cflag)==1) cflag=array(cflag,dimsof(x)(0));
      write,ffc,format="%s\n",__xgobcol(cflag);
    }
    if (!is_void(glyphs)) {
      if (numberof(glyphs)==1) glyphs=array(glyphs,dimsof(x)(0));
      write,ffg,format="%d\n",glyphs;
    }
    else write,ffc,"\n";
    close,ff;
    close,ffc;
    if (is_void(cflag)) system,"rm ~/.crapxgobi.dat.colors"; 
    if (is_void(glyphs)) system,"rm ~/.crapxgobi.dat.glyphs"; 
 system,"xgobi ~/.crapxgobi.dat &";
}


func split_palette(name)
/* DOCUMENT split_palette
         or split_palette, "palette_name.gp"
     split the current palette or the specified palette into two
     parts; colors 0 to 99 will be a compressed version of the
     original, while colors 100 to 199 will be the reversed version.
   SEE ALSO: pl3tree, split_bytscl
 */
{
  if (!is_void(name)) palette, name;
  local r,g,b;
  palette,query=1, r,g,b;
  n= numberof(r);
  if (n<100) {
    palette, "idl-03.gp";
    palette,query=1, r,g,b;
    n= numberof(r);
  }
  r1= char(interp(r,indgen(n),span(1,n,100)));
  g1= char(interp(g,indgen(n),span(1,n,100)));
  b1= char(interp(b,indgen(n),span(1,n,100)));
  r=grow( r1,r1(::-1));
  g=grow( g1,g1(::-1));
  b=grow( b1,b1(::-1));

  palette, r,g,b;
}






func merge_palette(name1,name2)
/* DOCUMENT split_palette
         or split_palette, "palette_name.gp"
     split the current palette or the specified palette into two
     parts; colors 0 to 99 will be a compressed version of the
     original, while colors 101 to 200 will be frm the second palette reversed.
   SEE ALSO: pl3tree, split_bytscl
 */
{
  if (is_void(name2)) name2="idl-08.gp";
  if (is_void(name1)) name1="idl-01.gp";
  
  palette, name1;
  local r1,g1,b1;
  palette,query=1, r1,g1,b1;
  n= numberof(r1);
   palette, name2;
  local r2,g2,b2;
  palette,query=1, r2,g2,b2;
  n= numberof(r2);
  r1= char(interp(r1,indgen(n),span(1,n,100)));
  g1= char(interp(g1,indgen(n),span(1,n,100)));
  b1= char(interp(b1,indgen(n),span(1,n,100)));
  r2= char(interp(r2,indgen(n),span(1,n,100)));
  g2= char(interp(g2,indgen(n),span(1,n,100)));
  b2= char(interp(b2,indgen(n),span(1,n,100)));
  r=grow( r1,r2(::-1));
  g=grow( g1,g2(::-1));
  b=grow( b1,b2(::-1));
  palette, r,g,b;
}

     func diag(x)
/* DOCUMENT diag(x) return diagonal matrix with diag elements given by x
   SEE ALSO: unit
*/
{ n=dimsof(x)(0);
 res= array(0.,[2,n,n]);
 res(indgen(1:n*n:n+1))=x(*);
 return res;
    }


func normalise(x,var=) 
/* DOCUMENT normalise(x,var=) returns (x-min(x))/(max(x)-min(x))
        or  (x-avg(x))/(rms(x)) if var=1
*/
{  local tt;
  if(is_void(var))  tt= (x-min(x))/(max(x)-min(x));
  else tt = (x-avg(x))/(x(*)(rms));
  return tt;
}

func lighter(n) 
/* DOCUMENT lighter(n) ligthens current palette;
     
   SEE ALSO: darker
*/
{
  require,"color.i"; if (is_void(n)) n=1;brighten,n*1.5; return 0;}

func darker(n)  
/* DOCUMENT darker(n) darkens current palette;
     
   SEE ALSO: ligther
*/
 {
  require,"color.i"; if (is_void(n)) n=1;brighten,.75/n; return 0;}

func skewness(x)  
/* DOCUMENT skewness(x)  returns   statistical skewness (3rd cumulant)
   SEE ALSO: kurtosis, AVG,median
*/
{ if (numberof(x)>1) {return ((x(*)-avg(x))^3)(avg)/x(*)(rms)^3;} else
    return 0.;}

func kurtosis(x)  
/* DOCUMENT    kurtosis(x) returns the statistical  kurtosis excess (4rth cumulant)
     SEE ALSO: kurtosis, AVG,median
*/
{ if (numberof(x)>1) {return (((x(*)-avg(x))^4)(avg)-3*(x(*)(rms))^4)/x(*)(rms)^4;} else
    return 0.;}

func RMS(x)  
/* DOCUMENT   RMS(x)  returns the  statistical  rms (2nd cumulant)
   SEE ALSO: kurtosis, AVG,median
*/
{
  if (numberof(x)>1){ return x(*)(rms);} else
    return 0.;
}

func AVG(x)  
/* DOCUMENT  AVG(x)  returns   statistical  avg (1st cumulant)
   SEE ALSO: kurtosis, RMS,median
*/
{ return x(*)(avg);}


func rvp(void)
/* DOCUMENT rvp  reverses current color palette
     SEE ALSO: lighter,darker
*/
     {
       palette,red,green,blue,query=1;
       palette,red(::-1),green(::-1),blue(::-1);
       return 0
       }



func im(x)
/* DOCUMENT  im(x) returns the  imaginary part of x
   SEE ALSO: x.im;
   */
{
  return float(-1.i*(x-conj(x))/2.);}      



func uniqstr(x)
/* DOCUMENT removes duplicates in vector x of strings
   
   SEE ALSO: unique
 */
{
  local flag,backup;

  nline=numberof(x);
  flag=array(0,dimsof(x));
  backup="";
  for(i=1;i<=numberof(x);i++)
    {
      name=x(i);
      if(anyof(strmatch(backup,name))) continue;
      flag(i)=1;
      backup+="@"+name;
    }
  return x(where(flag));
}


func unique(x) 
/* DOCUMENT unique(x)   removes duplicates in vector x
     
 */
  {
local x1,t; x1=x(sort(x));
if (numberof(x)==1) return x;  t= x1(dif);  t=grow(t,(x1(0) != x1(1))); 
 if (!t(sum)) return x1(1);
return  x1(where(t !=0));}



func load(fname)
/* DOCUMENT load(fname) loads an ascii file 
   and returns the table corresponding to its content. 
   WARNING it should contain only numerical characters
   SEE ALSO: asciiRead,smread
*/
{
 local ll;
 if (!open(fname, "", 1))
        error, "no such file or directory \""+fname+"\"";
 sz= 0L; sz2=0L;
    if (sread(rdline(popen("wc -l " + fname, 0)), sz) != 1)
        error, "cannot get nol of file \""+fname+"\"";
    if (sread(rdline(popen("wc " + fname+" | awk '{print $2}'", 0)), sz2) !=1)
        error, "cannot get noc of file \""+fname+"\"";
if(sz2/sz>1)  {  mat =array(0.,sz2/sz,sz);
    write,format="file %s contains %d x %d  data \n",fname,sz,sz2/sz;
   }
else {mat= array(0.,sz);
    write,format="file %s contains %d   data \n",fname,sz;
   }
    read,open(fname,"r"),mat;
    return transpose(mat);
}


func _expstr(arg){
  require, "string.i";
      r=arg;b=[];
      while ((r)!= string(0))
        {
          tt=strtok(r);
          if(tt(1) !=string(0)) b=grow(b,tt(1));
          r=tt(2);
        }
     res="";
 for(j=1;j<=dimsof(b)(2);j++){
      u=b(j);
      if (strchr(u,':')){
	if (strlen(u)==3){
	x1=0; er=sread(strpart(u,strchr(u,':')-1:strchr(u,':')-1),x1); 
	x2=0; er=sread(strpart(u,strchr(u,':')+1:strchr(u,':')+1),x2);
	} else
	if (strlen(u)==4){
	x1=0; er=sread(strpart(u,strchr(u,':')-1:strchr(u,':')-1),x1); 
	x2=0; er=sread(strpart(u,strchr(u,':')+1:strchr(u,':')+2),x2);
	} else
	if (strlen(u)==5){
	x1=0; er=sread(strpart(u,strchr(u,':')-2:strchr(u,':')-1),x1); 
	x2=0; er=sread(strpart(u,strchr(u,':')+1:strchr(u,':')+2),x2);
	}
 x =indgen(x1:x2);
xx=""; for(i=1;i<=dimsof(x)(2);i++){xx+=" "; xx+=pr1(x(i));}
 res+=xx;  
      } else {res+=" "; res+=u;}
}
  return res;
}


func smwrite(fname,x,head=)
/* DOCUMENT   smwrite(fname,head=)
     write a la sm;
     SEE ALSO: asciiRead,load
*/
{
  ff=open(fname,"w");
  x=transpose(x);
  nn=dimsof(x)(2); nn;
  if(nn>50) error,"can't deal with more than 15 column";
  form=strjoin(array("%1.6e\t",nn));
  if(!is_void(head))  write,ff,format="# %8s\n",head;
  write,ff,format=form+"\n",double(x),linesize=14*nn;
  close,ff;
  return fname;
}

func smread(fname,arg=,delimiters=,verb=)
/* DOCUMENT   smread(fname,arg)
     read a la sm; arg is a string of positions in the table.
     eg 1 2 or 1:4 8 etc ...
     SEE ALSO: asciiRead,load
*/
{
  require, "string.i";
ff=open(fname,"r");
 if(is_void(verb)) verb=0; 
 if (is_void(arg) || arg=="-") arg="1:19";
 arg =_expstr(arg); 
if (is_void(delimiters)) delimiters=": \t\n";
a=1;
b1=[];idx=0;
 iii=0; 
while (a!= string(0))
{
  iii++;
  if( !(iii %500)) write,"read 500*",(iii/500);
  a= rdline(ff);
  if( strpart(strtok(a)(1),1:1) =="#") {if(verb!=0) write,a;}
  else 
    if( strtok(a)(1) !=string(0))  {idx++;
      r=a;b=[];
      while ((r)!= string(0))
        {
          tt=strtok(r,delimiters);
          if(tt(1) !=string(0)) b=grow(b,tt(1));
          r=tt(2);
        }
      b1=grow(b1,b);
    }
}
close,ff;
b2=array(string,dimsof(b1)(2)/idx,idx);
b2(*)=b1(*);
pos =array(long,128);
er=sread(arg,pos);
pos=unique(pos)(2:);
 if (min(pos)>dimsof(b2)(2)) {error,"error in selection";}
pos=min(pos,dimsof(b2)(2));
 pos=unique(pos);
if (dimsof(b2)(2)>1){
  res=array(0.,dimsof(b2)(3),dimsof(pos)(2));
  x=array(0.,dimsof(b2)(3));
  for(i=1;i<=dimsof(pos)(2);i++) {
    err= sread(b2(pos(i),),x);
    if (err==0) { x=b2(pos(i),); write,
                    format=" ascii col %i\n",pos(i); res=x;  return res; }
    else {res(,i)=x; if ((verb>0)) write, format=" col %i\n",pos(i);}
  }
  res= transpose(res);
}
else
{
  res=array(0.,dimsof(b2)(3));
  x=array(0.,dimsof(b2)(3));
  err= sread(b2(*),x); 
    if (err<dimsof(b2)(3)) {if (verb>0)  write,
                    format=" ascii col %i\n",1; res=b2(pos,); info,res; return res; }
    else res=x;
}
res=transpose(res);
if(verb!=0)info,res;
return res;
}


func get_filelength(name)
/* DOCUMENT nline = get_filelength(name)
     Returns number of line of NAME.

   RESTRICTIONS:
     Makes use of UNIX command "wc" which may not exists on other systems.

   SEE ALSO: exec.
*/
{
    if (!open(name, "", 1))
        error, "no such file or directory \""+name+"\"";
    sz= 0L;
    if (sread(rdline(popen("wc -l " + name, 0)), sz) != 1)
        error, "cannot get nol of file \""+name+"\"";
    return sz;
}

func sqr(x)
/* DOCUMENT sqr(x) returns x^2
     SEE ALSO: pow;
*/
{
return x*x;}

func  spanf(a,b,n,f,g)
/* DOCUMENT  spanf(a,b,n,f,g) generalises spanl to f and g= f-1;
 EXAMPLE spanf(1,10,10,exp,log)-spanl(1,10,10) yields zero
 SEE ALSO: span,spanl;
*/
  {
  if (is_void(f)) { f=log; g=exp; }
 return f(span(g(a),g(b),n));} 


func thread(f,x,inter=)  
/* DOCUMENT    thread f over x  a la MMA 
   SEE ALSO: map,_map
*/
{
 local i,n;
 n= dimsof(x)(0);
   for(i=1;i<=n;i++) { if(!is_void(inter)) {
     i;  rdline,prompt="hit RET or Enter to continue";}
                       res=f(x(..,i));}
 
} 




func map(f,x,repeat=,rk=,pos=) 
/* DOCUMENT    maps f over x  a la MMA 
  if rk not void map assumes f is contracting dimsof(f(x(1,..))(0)  is dimsof(x(1,..))(0) +rk;
   if repeat is non zero integer >=2 f is map repeat times.
   if pos is non zero it applies f at position pos in x;
   EXAMPLE  func g(x) { return [x(1),x(2),x(3),x(3)];}
    map(g,[[1,2],[1,1],[2,2]],rk=1)=[[1,2],[1,1],[2,2],[2,2]]
  SEE ALSO: lst2arr(_map(f,arr2lst(x))),thread
*/
{
  local i,n;  
  // if (!is_void(pos)) x=transpose(x,[q,pos]);
 if (!is_void(pos)) error,"not implemented properly";
   n= dimsof(x)(0);
   if (numberof(x)==1) return f(x);
 if (is_void(rk)) { res= array(structof(f(x(1))),dimsof(x)); }
 else {res=array(structof(f(x(1,))),dimsof(x))(1,..);}
 if(is_void(repeat) || repeat <=1)
    {
    if(is_void(rk))  for(i=1;i<=n;i++) { res(..,i)= f(x(..,i));}
 else {
   
      d=long(exp(log(dimsof(x)(2:-1))(sum))+1e-5);
      x1=reform(x,[2,d,dimsof(x)(0)]);
      res=lst2arr(_map(f,arr2lst(transpose(x1))));
      res= transpose(reform(res,[2,numberof(res)/dimsof(x1)(2),dimsof(x1)(2)]));
      d=dimsof(x);
      dd=grow(d(1)-1,d(2:-1));
      res=reform(res,dd);
 }
      return res;
    }
 // if (!is_void(pos)&& p==1) x=transpose(x,[pos,q]);
 return map(f,map(f,x),repeat=repeat-1);
 } 



func evalfunc(f,x)
/* DOCUMENT   evalfunc(f,x) evaluate f(x) even when f is a string  
   SEE ALSO: map;
*/
{
  if (structof(f)==string) f=symbol_def(f);
  return f(x);
}



func INFO(x) 
/* DOCUMENT    info data with min  
   SEE ALSO: stat,structof
*/
{
  stat,x;
}

func upload(__fname,s=)
/* DOCUMENT    upload data file 
s=1 silent mode //=== ATTENTION PB NOM s pris a reglet
// save,createb("toto.yor"),a,b;
// restore,openb("toto.yor"),a; //loads only a
// save,updateb("toto.yor"),b; //adds only a 
 SEE ALSO: createb,restore,save
*/
{

  __ff=openb(__fname); 
 if(is_void(s)  ) write,"recovered", *get_vars(__ff)(1);
 s=[];
 restore,__ff; close,__ff;
}

func cinterp3(z,n=,s=)
/* DOCUMENT cinterp3(z,n=,s=) reinterpolate cube
      keywords n and s correspond to n-fold and s offset
*/
{
  local x,xp;
 if (is_void(n)) n=2;
 if (is_void(s)) s=0;
 x=span(-1,1,dimsof(z)(2));
 xp=span(-1,1,n*dimsof(z)(2)+s);
  return interp(interp(interp(z,x,xp,1),x,xp,2),x,xp,3);
}

func cinterp2(z,n=,s=)
/* DOCUMENT  cinterp2 reinterpolate square
      keywords n and s correspond to n-fold and s offset
*/
{
  local x,xp;
 if (is_void(n)) n=2;
 if (is_void(s)) s=0;
 x=span(-1,1,dimsof(z)(2));
 xp=span(-1,1,n*dimsof(z)(2)+s);
  return interp(interp(z,x,xp,1),x,xp,2);
}



func HELP2(x,v=,l=,a=,h=)  
/* DOCUMENT    HELP(x,v=,l=,a=,h=)
   scans the yorick and user defined routines in   the $YHOME directory
   flags
   v = verbose on file name
   l = looks also in ~/Project/ for *.i files
   a = standard grep in ~/Project/ for *.i files 
   h = standard grep in HTML doc
   SEE ALSO: stat,INFO
*/
{
local flag;
if (is_void(v))  flag=" -h ";  else flag ="";
 tt=[".",  "/usr/local/share/yorick/1.5/i/",
   " ~/Yorick/include", "~/Yorick/include0",
   " ~/Yorick/Eric/"," ~/Yorick/Chris/",
   " ~/Yorick/yeti/"," ~/Yorick/Bastien"," ~/Yorick/Chris/devel/"];
 for(i=1;i<=numberof(tt);i++)
   { 
    if(i==1) ll=exec( " grep -i "+flag+pr1(x)+" "+tt(i)+"/*.i | grep extern");
    ll=exec( " grep -i "+flag+pr1(x)+" "+tt(i)+"/*.i | grep func");
   if (numberof(ll)<1000)   error,"ah";
   }
if (!is_void(l)) {ll=exec(  " grep -i "+flag+pr1(x)+"  `find /home1/pichon/Project/ -name '*.i' -print ` | grep func" ); if (numberof(ll)<1000) print,ll; }
if (!is_void(a)) {ll=exec(  " grep -i "+flag+pr1(x)+"  `find /home1/pichon/Project/ -name '*.i' -print ` " ); if (numberof(ll)<1000) print,ll; }
if (!is_void(h)) {ll=exec(  " grep -i "+flag+pr1(x)+"  `find /home1/pichon/Yorick/doc/ -name '*.html' -print ` " );if (numberof(ll)<1000) print,ll; }
}



func HELP(x,v=,l=,a=,h=)  
/* DOCUMENT    HELP(x,v=,l=,a=,h=)
   scans the yorick and user defined routines in   the $YHOME directory
   flags
   v = verbose on file name
   l = looks also in ~/Project/ for *.i files
   a = standard grep in ~/Project/ for *.i files 
   h = standard grep in HTML doc
   SEE ALSO: stat,INFO
*/
{
local flag;
if (is_void(v))  flag=" -h ";  else flag ="";
system, " grep -i "+flag+pr1(x)+" /usr/local/share/yorick/1.4/startup/*.i | grep extern";
system, " grep -i "+flag+pr1(x)+" ~/Yorick/include/*.i | grep func";
system, " grep -i "+flag+pr1(x)+" ~/Yorick/include0/*.i | grep func";
system, " grep -i "+flag+pr1(x)+" ~/Yorick/Eric/*.i | grep func";
system, " grep -i "+flag+pr1(x)+" ~/Yorick/Chris/*.i | grep func";
system, " grep -i "+flag+pr1(x)+" ~/Yorick/yeti/*.i | grep func";
system, " grep -i "+flag+pr1(x)+" ~/Yorick/Bastien/*.i | grep func";
 system, " grep -i "+flag+pr1(x)+" ~/Yorick/Chris/devel/*.i | grep func";
system, " grep -i "+flag+pr1(x)+" ./*.i | grep func";
if (!is_void(l)) {system,  " grep -i "+flag+pr1(x)+"  `find /home1/pichon/Project/ -name '*.i' -print ` | grep func" ;}
if (!is_void(a)) {system,  " grep -i "+flag+pr1(x)+"  `find /home1/pichon/Project/ -name '*.i' -print ` " ;}
if (!is_void(h)) {system,  " grep -i "+flag+pr1(x)+"  `find /home1/pichon/Yorick/doc/ -name '*.html' -print ` " ;}
}

func ymread(mfile)
/* DOCUMENT   ymread(mfile) tensor read in an ascii format compatible
   with mathematica.
   note that mma and yorick have oposite convention for
   dimensions
ie [[1],[2]] has dimensions
 1,2 in yorick and
 2,1 in mma
 SEE ALSO: ymwrite
*/
{
nd=1; 
file=open(mfile,"r"); read,file,nd;
dd =array(0,nd); 

read,file,dd; print,dd;
a =array(0.,grow(nd,dd(::1)));
read,file,a;
return a; // this transpose accounts for the convention mma <> yorick
}





func ymwrite(mfile,mtx,s=)
/* DOCUMENT 
   ymwrite(mfile,mtx,s=) tensor write
   note that mma and yorick have oposite convention for
   dimensions
   ie [[1],[2]] has dimensions 1,2 in yorick and
   2,1 in mma
 SEE ALSO: ymread
*/
{
file=open(mfile,"w");
nd = dimsof(mtx)(1);
write,file,nd;
for(i=1;i<=nd; i++) write,file,dimsof(mtx)(i+1);
write,file,(mtx)(*),format="%4.16g ";
return mfile;
}

func prm(x)
/* DOCUMENT 
prm(x) displays a matrix
*/
{
  n =dimsof(x)(2);
  for(i=1;i<=n;i++){ 
    write, x(i,),linesize=500;
  }
}

func pprt(x)
/* DOCUMENT 
  pprt(x) displays a list
*/
{
    pretty_print_string_list,map(pr1,x(*));
}

func plb(z,y,x,bar=,every=, legend=, hide=, type=,
         width=, color=, closed=, smooth=,
	 marks=, marker=, mspace=, mphase=, rays=,
         arrowl=, arroww=, rspace=,
	 rphase=,fun=,ofset=,inter=,delay=,labs=,vert=,noerase=)
/* DOCUMENT plb(z,y,x,bar=,every=, legend=, hide=, type=,
         width=, color=, closed=, smooth=,
	 marks=, marker=, mspace=, mphase=, rays=,
         arrowl=, arroww=, rspace=,
	 rphase=,fun=,ofset=,inter=,delay=)
 will plot the bundle z versus y and x
 KEYWORDS fun alows the user to specify another function instead of plg
         like pla but in color
 SEE ALSO: PLB,pla,thread
*/
{ local n,colors,levs;
 if (is_array(z)) n =dimsof(z)(3); else n=_len(z);
 if(is_void(x)){cmax=n; cmin=1;} else {cmax=x(*)(max); cmin=x(*)(min);}
 if(is_void(ofset)) ofset=0.;
 colors= bytscl(span(1,n+1,n+3),cmin=0.5 ,cmax=n+1.5); colors=colors(1:n+1);
if (is_void(width)) width=3;
 for(i=1;i<=n;i++){
   if (is_array(z)) zz= z(,i); else zz=_car(z,i);
   if(is_void(fun)){ if(type){ if(type<0) type1=i; else type1=type; } 
  plg,zz+ofset*i,y,color=colors(i),width=width,legend=legend, hide=hide, type=type1,closed=closed, smooth=smooth, marks=marks,
      mspace=mspace, mphase=mphase, rays=rays, arrowl=arrowl, arroww=arroww,
                     rspace=rspace, rphase=rphase;
if(marker) PL,zz+ofset*i,y,color=colors(i),marker=i % 5;
 if (!is_void(inter)) { plt,pr1(i)+"/"+pr1(n),0.5,0.5,tosys=0,color=Colors(n)(i);   tt=mouse(,,"");if (tt(-1)==3) i-=2;  if (tt(-1)==2) exit; if((i<n) & !noerase) fma;}
 if (!is_void(delay)) {plt,pr1(i)+"/"+pr1(n),0.5,0.5,tosys=0,color=Colors(n)(i); pause,delay;  if((i<n) & !noerase) fma;}
} else {
   if(type){ if(type<0) type1=i; else type1=type; } 
         fun,z(,i)+ofset*i,y,color=colors(i),width=width,type=type1;
         if (!is_void(inter)){plt,pr1(i)+"/"+pr1(n),0.5,0.5,tosys=0;
         tt=mouse(,,"");
         if (tt(-1)==3) i-=2;if (tt(-1)==2) exit;if(i<n) fma;}
         if (!is_void(delay)) {plt,pr1(i)+"/"+pr1(n),0.5,0.5,tosys=0;  pause,delay;  if((i<n) & !noerase) fma;}

}
}
 if(is_void(x)) levs= span(cmin,cmax,n); else levs=x;
if (is_void(bar))    color_bar,levs,colors,vert=1,font="timesI", 
                       height=12,format="%3.2f",labs=labs;
 
}

func gpr(file,&x,&y)
/* DOCUMENT 
   gnuplot style read
*/
{
z= gpRead(file,x,y,msbf=1);
fma;
pl_fc,z,y,x;
color_bar,vert=1
return z;
}

func plct(z,y,x,type=,levs=,color=,width=)
/* DOCUMENT plct(z,y,x,type=,levs=,color=,width=)
     contour plot without fill
     SEE ALSO: plcont;
*/
{
  return plcont(z,y,x,mono=1,bar=0,type=type,width=width,levs=levs,color=color);
}

func plcont(z,y,x,levs=, mono=,width=,color=,triangle=, region=,bar=,cmin=,cmax=,type=) 
/* DOCUMENT plcont(z,y,x,   levs=, mono=,width=,color=,
              triangle=, region=,bar=,cmin=,cmax=,type=) 
              contour plot with legend in colour with filling
   SEE ALSO: plct,plcl;
*/
{
  if(is_void(bar)) { fma; dobar=1;} 
  if(mono) {

if(is_void(x)) x=indgen(dimsof(z)(2))(,-:1:dimsof(z)(3));else x=x(,-:1:dimsof(z)(3)); 
if(is_void(y)) y=indgen(dimsof(z)(3))(-:1:dimsof(z)(2),); else y=y(-:1:dimsof(z)(2),); 
    plc,z,y,x,levs=levs,type=type,color=color,width=width;
  } else
{
  pl_fc,z,y,x,levs=levs,colors=colors,triangle=triangle,
      region=region,type=type,width=width; //,font="timesI", height=12;
 if(dobar){    color_bar,vert=1,font="timesI", 
		 height=12,format="%3.2f";}
}
}

func plcl(z,y,x,levs=, opaque=,height=,mono=) 
/* DOCUMENT plcl(z,y,x,levs=, opaque=,height=,mono=) 
     contour plot with label
     SEE ALSO: plct,plcl;
*/
{  require, "plclab.i"; 
if(is_void(levs))   levs =span(min(z),max(z),15);
if(is_void(x)) {
  x =span(-1,1,dimsof(z)(2))(,-:1:(dimsof(z)(3)));
  y =span(-1,1,dimsof(z)(3))(-:1:(dimsof(z)(2)),);
} 
if(is_void(mono))  	   plfc,z,y,x,levs=levs;
     	   plc,z,y,x,levs=levs;
      	   plc_label,z,y,x,levs,opaque=opaque,height=height;
}




func plcol(c,y,x,pal=,marker=,width=,type=)
/* DOCUMENT plcol(c,y,x,pal=,marker=,width=,type=)
    displays a cluster of points (x,y) colored by vector c with marker
    SEE ALSO: plcolor,plmk,plp;
*/
{  // see bytscl for palette;
local c1,ll,i;
  if(is_void(pal)) pal="idl-03.gp";
  if(is_void(marker)) marker=2;  
if(is_void(width)) width=5;
if(is_void(type)) type=0;
  palette,pal;
c1 =long(100+100*(c-c(min))/(c(max)-c(min)));
 ll=dimsof(c)(2);
 for(i=1;i<ll;i++){ plg,y(i),x(i),type=type,marker=marker,width=width,color=c1(i);}
}




func plcir(r,y,x,color=,width=,number=,type=,legend=)
{
  require, "Bastien/jgraph.i";
  for(i=1;i<=numberof(x);i++)
    JDrawCircle,x(i),y(i),r(i),color=color,width=width,type=type,legend=legend,number=number;
}

func plcolor(c,y,x,npol=,cmin=,cmax=,size=,bar=,scale=)
/* DOCUMENT plcolor(c,y,x,npol=,cmin=,cmax=,size=,bar=,scale=)
   cluster of points (x,y) colored by vector c
    SEE ALSO: plcol,plmk,plp;
*/
{ local n,list,phase,theta,x1,y1,sx,sy,s,n0,n1,levs,dobar;
 if(is_void(npol)) npol=3;
 if(is_void(size)) size=1;
 if(is_void(cmin)) cmin=c(min);
 if(is_void(bar)) dobar=1;
  if(is_void(cmax)) cmax=c(max);
 if (npol < 3) {write,"npol should be larger than 3"; return;}
  sx=(x(max)-x(min))/15.*size;
  sy=(y(max)-y(min))/15.*size;
  n0 =numberof(x);
 n= (indgen(n0)*0+1)*npol;
 list= histogram(n(cum)+1)(psum:1:-1);
 phase= indgen(numberof(list))-n(cum:1:-1)(list);
 theta= 2.*pi*phase/n(list);
  if(!is_void(scale)) 
    { scl= (c-cmin)/(cmax-cmin);
      scl =exp(1.5*scl);
  } else scl=c*0+1;
x1= scl(list)*sx*cos(theta) + x(list);
    y1= scl(list)*sy*sin(theta) + y(list);
    plfp, c, y1, x1, n, cmin=cmin,cmax=cmax,edges=0;
      n1= 10;
    colors= bytscl(span(1,n1+1,n1+1),cmin=0.5,cmax=n1+1.5);
 levs= span(cmin,cmax,n1);
 if(dobar){ color_bar,levs,colors,vert=1,font="timesI", 
		 height=12,format="%3.2f";}
}





func plii(z,y,x)
/* DOCUMENT plii(z,y,x)
   like pli but with proper scaling
*/
{ pli, z, min(x),min(y),max(x),max(y); 
}






func plia(z, ..,legend=,hide=,top=,factor=)
/* DOCUMENT plia
   like pli but with cuts in z
*/
{
  zm=z(avg);
  zs=z(*)(rms);
  
  if(is_void(factor))
    {
      factor=3.;
    }

  ma=more_args();
  if(ma!=0&&ma!=2&&ma!=4) error,"plia needs either 0, 1, or 2 corner (x,y) points";

  cmin=zm-zs/factor;
  cmax=zm+zs/factor;
  if(ma==0)
    pli,z,legend=legend,hide=hide,top=top,cmin=cmin,cmax=cmax;
  
  if(ma==2)
    {
      x1=next_arg();
      y1=next_arg();
      pli,z,x1,y1,legend=legend,hide=hide,top=top,cmin=cmin,cmax=cmax;
    }
  if(ma==4)
    {
      x0=next_arg();
      y0=next_arg();
      x1=next_arg();
      y1=next_arg();
      pli,z,x0,y0,x1,y1,legend=legend,hide=hide,top=top,cmin=cmin,cmax=cmax;
    }

  return [cmin,cmax];
}


func scale(data,..)
{
  if (more_args()!=0&&more_args()!=2)
    error,"need 1 or 3 arguments !";
  mini=0.;
  maxi=1.;
  if (more_args())
    {
      mini=double(next_arg());
      maxi=double(next_arg());
    }
  mn=data(*)(min);
  mx=data(*)(max);
  mx-=mn;
  maxi-=mini;
  return (data-mn)/mx*maxi+mini;
}         



func plk(z,y,x,cmin=,cmax=)
/* DOCUMENT plk(z,y,x,cmin=,cmax=) displays density map
   z versus y and x together with a colour bar
   like pli but with colorbar
   SEE ALSO: pli,plii,plcont,plb
*/
{ local n,colors,levs;
 n =min(15,dimsof(z)(3));
colors= bytscl(span(1,n+1,n+1),cmin=0.5 ,cmax=n+1.5);
 if(is_void(x))   {pli,z,cmin=cmin,cmax=cmax;} else
   {pli,z, min(x),min(y),max(x),max(y),cmin=cmin,cmax=cmax; }
 if(is_void(cmin) & is_void(cmax)) levs= span(z(*)(min),z(*)(max),n);
 else if (is_void(cmax)) levs= span(cmin,z(*)(max),n);
          else if (is_void(cmin)) levs= span(z(*)(min),cmax,n);
 else levs= span(cmin,cmax,n);
    color_bar,levs,colors,vert=1,font="timesI", 
      height=12,format="%3.2f";
}




func sh(z,y,x)
/* DOCUMENT sh(x,y,x) returns the 
   3D curve z=z(x,y)
   in colour
   SEE ALSO: sh3
*/
{
require, "plwf.i"
 
  if (!is_void(y)) plwf,z,y,x, shade=1, ecolor="red";
  else  plwf,z, shade=1, ecolor="red";
  orient3; cage3,1;
}





func plp3(z,y,x,alt=,az=, legend=, hide=, type=, width=, color=,nocage=)
/* DOCUMENT plp3
   plots clouds of points in 3d
*/
{ local eps;
 if(is_void(alt)) alt=45;
 if(is_void(az)) az=30;
 if(is_void(color)) color="black";
 require, "Eric/plot.i"
   // window,style="nobox.gs";
  eps =1.e-5*max(x,y,z);
x-=min(x);
 y-=min(y);
 z-=min(z);
 pl3dj,x+eps,y+eps,z+eps,x,y,z,alt=alt,az=az,
    legend=legend, hide=hide, type=type, width=width, color=color;
 if(is_void(nocage)){
 LL=max([x,y,z](*));
 pl3dj,LL,0,0,0,0,0,width=3,type=2,color="blue";
    pl3dj,0,LL,0,0,0,0,width=3,type=2,color="blue";
    pl3dj,0,0,LL,0,0,0,width=3,type=2,color="blue";

    pl3dj,LL,0,0,LL,LL,0,width=3,type=2,color="blue";
    pl3dj,LL,0,0,LL,0,LL,width=3,type=2,color="blue";
    pl3dj,LL,LL,0,LL,LL,LL,width=3,type=2,color="blue";

    pl3dj,LL,LL,LL,LL,0,LL,width=3,type=2,color="blue";
    pl3dj,LL,LL,LL,0,LL,LL,width=3,type=2,color="blue";

    pl3dj,LL,0,LL,0,0,LL,width=3,type=2,color="blue";
    pl3dj,0,LL,LL,0,0,LL,width=3,type=2,color="blue";
    pl3dj,0,LL,0,0,LL,LL,width=3,type=2,color="blue";
 }
}




func sh3(z)
/* DOCUMENT sh(z)
   3D curve in colour
*/
{
window, max(0, current_window()), wait=1, style="nobox.gs";
nx = dimsof(z)(2);
ny = dimsof(z)(3);
 x= indgen(nx)*1.;
 y= indgen(ny)*1.;
        pl3s, z, y, x, axis=1, fill=2, edges=1, font="times", height=10,box=1;
}

func _spin3(i)
{
  if (i>=nframes) return 0;
  rot3,,,-phi
  rot3,,-theta,dtheta;
  rot3,,theta,phi;
  draw3; limits,-1,1,-1,1;
  
  return 1;
}


func slicer(w,z,y,x,val=)
/* DOCUMENT slicer
   3D isocontour slicer(w,z,y,x,val=)
EXAMPLE
x0 =span(-1,1,15);
x =x0(,-:1:15)(,,-:1:15);
y =x0(-:1:15,)(,,-:1:15);
z =x0(-:1:15,)(-:1:15,,);
w= 2.*exp(-x^2-y^2-2*z^2);
slicer(w,z,y,x,val=0.5);
SEE ALSO: slice
*/
{
  require, "plwf.i";
  require, "slice3.i";
  local nx,ny,nz,xyz,m3,nv,xyzv,lims,dx,dy;

 window,style="nobox.gs"; 

if(is_void(val)) val=abs(w(ptp))/2.;
  print,"cut at ",val;
if(is_void(z)) {
  x =span(-1,1,dimsof(w)(2))(,-:1:(dimsof(w)(3)))(,,-:1:(dimsof(w)(4)));
  y =span(-1,1,dimsof(w)(3))(-:1:(dimsof(w)(2)),)(,,-:1:(dimsof(w)(4)));
  z =span(-1,1,dimsof(w)(4))(-:1:(dimsof(w)(2)),)(-:1:(dimsof(w)(3)),,);
}

nx=dimsof(x)(2);ny=dimsof(x)(3);nz=dimsof(x)(4);
xyz= array(0.0, 3, nx,ny,nz);
xyz(1,..)= x(,1,1);
xyz(2,..)= (y(1,,1))(-,);
xyz(3,..)= (z(1,1,))(-,-,);
m3= mesh3(xyz, w);
slice3, m3, 1,value=val, nv,xyzv;  /* inner isosurface */
fma;
pl3tree, nv,xyzv;
orient3;
light3,diffuse=.2,specular=1;
limits;
lims= limits();
dx= 1.1*max(lims(2),-lims(1));
dy= 1.1*max(lims(4),-lims(3));
limits, -dx,dx,-dy,dy;
 palette,"idl-03.gp";limit3,-1,1,-1,1,-1,1; cage3,1;
spin3,50,[1,0,0];
//limits; 
}






func ple(x,pal1=,pal2=,fun=,eps=)
/* DOCUMENT ple(x,pal1=,pal2=,fun=,eps=)
   uses graphics function fun (default plk) to display x with
   a split palette ( default plk "idl-03.gp" and "idl-01.gp";)
   with a log scale shifted by eps around zero
   KEYWORDS ,pal1=,pal2=,fun=, eps=
   require,"legndr.i"
   dims=[2,30,60];
   lat=span(0,pi,dims(2))(,-:1:dims(3));lon=span(-pi,pi,dims(3))(-:1:dims(2),); 
   y53=legndr(5,3,cos(lat))*sin(3*lon);
   ple,float(y53),fun=plsphere,eps=10.
   SEE ALSO: pli,plk
*/
{
  if(is_void(pal1)) pal1="idl-03.gp";
  if(is_void(pal2)) pal2="idl-01.gp";
  if(is_void(fun)) fun=plk;
  if(is_void(eps)) eps=0.1;
  
  merge_palette,pal1,pal2;
  wp= (x>=0);wn= (x<0);
  M= max(x);
  m= min(x);
  x1= log(x/M*wp+eps)-log(x/m*wn+eps);
  fun ,x1;
  
  return x1;
  
  }
     


func slice(w,val=,save=,opengl=,noerase=)
/* DOCUMENT slice   3D isocontour slice(w,z,y,x,val=)
   EXAMPLE
   x0 =span(-1,1,15);
   x =x0(,-:1:15)(,,-:1:15);
   y =x0(-:1:15,)(,,-:1:15);
   z =x0(-:1:15,)(-:1:15,,);
   w= 2.*exp(-x^2-y^2-2*z^2);
   slice(w,z,y,x,val=0.5);
   SEE ALSO: slicer
*/
{
require, "plwf.i"
require, "slice3.i"
  local nx,ny,nz,xyz,m3,nv,xyzv,lims,dx,dy;

 if(is_void(val)) val =AVG(w)+RMS(w)/2; 
  print,"cut at ",val;
  x =span(-1,1,dimsof(w)(2))(,-:1:(dimsof(w)(3)))(,,-:1:(dimsof(w)(4)));
  y =span(-1,1,dimsof(w)(3))(-:1:(dimsof(w)(2)),)(,,-:1:(dimsof(w)(4)));
  z =span(-1,1,dimsof(w)(4))(-:1:(dimsof(w)(2)),)(-:1:(dimsof(w)(3)),,);

 window,style="nobox.gs"; palette,"idl-03.gp";
nx=dimsof(x)(2);ny=dimsof(x)(3);nz=dimsof(x)(4);
xyz= array(0.0, 3, nx,ny,nz);
xyz(1,..)= x(,1,1);
xyz(2,..)= (y(1,,1))(-,);
xyz(3,..)= (z(1,1,))(-,-,);
m3= mesh3(xyz, w);


  if (is_void(noerase)){
   light3, ambient=.1,diffuse=.1,specular=2.,
     sdir=[[0,0,-1],[1,.5,1]],spower=[4,4];
   clear3;fma;
 }

  if(opengl) save=1;
 if(save)
   {
     ff=open(".crap.ts","w");
     write,ff,"% Graphics3D objects";
     write,ff,"boundingbox";
     write,ff,"-0.03 -0.03 -0.03";
     write,ff,"1.03 1.03 1.03";
     write,ff,"viewpoint";
     write,ff,"1.3 -2.4 2.";
     write,ff,"ambientlight";
     write,ff,"0 0.5 0";
     write,ff,"lightsources";
     write,ff,"1.0 0. 1.";
     write,ff,"1 0 0";
     write,ff,"1. 1. 1.";
     write,ff,"0 1 0";
     write,ff,"0. 1. 1.";
     write,ff,"0 0 1";
   }
if (dimsof(val)(0)==0) 
{ 
slice3, m3, 1,value=val, nv,xyzv;  /* inner isosurface */
pl3tree, nv,xyzv;
if(save){     for(i=1;i<=numberof(nv);i++)
       {
         write,ff,format="\npolygon%s\n"," ";
         form=strjoin(array("%1.6e\t",3));
         //form=strjoin(array("%1.8f\t",3));
         idx=(nv(1:(i-1))(sum)+1);
         idx=((i==1)?1:idx);
         xxx=xyzv(,idx:(idx-1+nv(i)));
         write,ff,format=form+"\n",double(xxx),linesize=14*3;
       }
} 
} else
  for(i=1;i<=dimsof(val)(0);i++)
    { 
      slice3, m3, 1,value=val(i), nv,xyzv;  /* inner isosurface */
      pl3tree, nv,xyzv;
      if(save){
        for(ii=1;ii<=numberof(nv);ii++)
        {
          write,ff,format="\npolygon%s\n"," ";
          form=strjoin(array("%1.6e\t",3));
          //form=strjoin(array("%1.8f\t",3));
          idx=(nv(1:(ii-1))(sum)+1);
          idx=((ii==1)?1:idx);
          xxx=xyzv(,idx:(idx-1+nv(ii)));
          write,ff,format=form+"\n",double(xxx),linesize=14*3;
        }
      } 
    }
 if(save)  close,ff;
orient3;
light3,diffuse=.2,specular=1;
limits;
lims= limits();
dx= 1.1*max(lims(2),-lims(1));
dy= 1.1*max(lims(4),-lims(3));
limits, -dx,dx,-dy,dy;
//spin3,100;
limit3,-1,1,-1,1,-1,1;
cage3,1;


  if(opengl) system,"mathview3d -mvps 20 -mvfile .crap.ts &";
}




func plaitof(z, cmin=, cmax=,smooth=,minr=,maxr=,bar=,color=,width=,type=,nlevs=,nofill=,cont=,nogrid=)
/* DOCUMENT plaitof(w,az=,anim=, cmin=, cmax=)  
   projects z onto a sphere
   KEYWORDS cmin=, cmax=,smooth=,bar=,minr=,maxr=
   EXAMPLE  require,"legndr.i"
   dims=[2,30,60];
   lat=span(0,pi,dims(2))(,-:1:dims(3));lon=span(-pi,pi,dims(3))(-:1:dims(2),); 
   y21=legndr(2,1,cos(lat))*sin(1*lon);
   plaitof(y21)
   SEE ALSO: ploctan,plsphere
*/
{     require, "Chris/aitof.i";

 if (is_void(rule)) rule=1;
 er= aitof_zdens(transpose(z)(::-1,::-1), cmin=cmin,cmax=cmax,smooth=smooth,rule=bar,
                 minr=minr,maxr=maxr,color=color,width=width,type=type,nlevs=nlevs,nofill=nofill,cont=cont);
 limits;
 
    }

func plsphere(w1,az=,anim=, cmin=, cmax=,border=)
/* DOCUMENT plsphere(w,az=,anim=, cmin=, cmax=)  
projects the colors w1 onto a sphere
   KEYWORDS ,az=,alt=,anim=, cmin=, cmax=, anim=
 require,"legndr.i"
 dims=[2,30,60];
 lat=span(0,pi,dims(2))(,-:1:dims(3));lon=span(-pi,pi,dims(3))(-:1:dims(2),); 
 y53=legndr(5,3,cos(lat))*sin(3*lon);
 SEE ALSO: plsphere;
*/
{
require, "plwf.i"

  if (is_void(az)) az=0;
 dims=dimsof(w1);
 lat=span(0,pi,dims(2))(,-:1:dims(3));
 lon=span(-pi,pi,dims(3))(-:1:dims(2),);
 
 xyz=[sin(lat)*cos(lon+az), sin(lat)*sin(lon+az), cos(lat)];
 
 xx=[1,1,0](+)*xyz(,,+)/abs(1,1,0);
 yy=[1,-1,1](+)*xyz(,,+)/abs(1,-1,1);
 zz=[-1,1,2](+)*xyz(,,+)/abs(-1,1,2);



plf,w1,yy,xx,zz>0, cmin=cmin, cmax=cmax;
 if (!is_void(border)) { tt= span(0,2*pi,400); plg,sin(tt),cos(tt);}


 if(!is_void(anim) && anim >0)
  {
    animate,1;
    for(i=1;i<=150;i++)
   {   az=(i-1)*pi/149.; 
lat=span(0,pi,dims(2))(,-:1:dims(3));lon=span(-pi,pi,dims(3))(-:1:dims(2),);
xyz=[sin(lat)*cos(lon+az), sin(lat)*sin(lon+az), cos(lat)];
xx=[1,1,0](+)*xyz(,,+)/abs(1,1,0);
yy=[1,-1,1](+)*xyz(,,+)/abs(1,-1,1);
zz=[-1,1,2](+)*xyz(,,+)/abs(-1,1,2);
    plf,w1,yy,xx,zz>0, cmin=cmin, cmax=cmax;    
   } }
animate,0;
}





func showsec(z,y,x,sty=,vmin=,vmax=,fun=,color=,noerase=,rrange=,all=) 
/* DOCUMENT  showsec(z,y,x,sty=,vmin=,vmax=,fun=,color=,noerase=,rrange=)
     displays n sections of nxnxn cube
     
     EXAMPLE showsec(z,sty="win33.gs");
     KEYWORDS sty= the style sheet to use
              fun= the function to use for plot (default pli)
     SEE ALSO: slice,sec3
*/
{
  require,"style.i";

  if (all)
    {
      pp= dimsof(z)(0);
      nn=pp/16;
      for(i=0;i<=nn-1;i++) {window,10+i; showsec,z(,,i*16+1:(i+1)*16),y,x,sty="win44.gs",
                                        vmin=vmin,vmax=vmax,fun=fun,color=color,noerase=noerase,rrange=rrange;
     
      } 
      qq= pp -16*nn; if(qq>0){ window,10+nn+1; showsec,z(,,nn*16+1:),y,x,sty="win44.gs",
                                 vmin=vmin,vmax=vmax,fun=fun,color=color,noerase=noerase,rrange=rrange;
      }
      return;
    }
  local i,vmax,tt,ll1,sys1;
  tt=get_style(ll1,sys1);
  if(is_void(sty)) { window,style="win22.gs"; imax=4;}
 else {
   if(is_void(noerase)) window,style=sty;
   imax=min(numberof(sys1),dimsof(z)(0));
 }
 animate,0;
 if (is_void(vmax)) vmax=max(z);
 if (is_void(vmin)) vmin=min(z);
 if (is_void(fun)) for(i=1;i<=imax;i++){
   plsys,i;
   pli,z(,,i),cmax=vmin,cmax=vmax; 
 }
 else for(i=1;i<=imax;i++)
   {
     plsys,i;
     fun,z(..,i),y,x;
     if (!is_void(rrange)) range,rrange;
   }
 return void;
}



func showanim(x,y,fun=,lim=,inter=,delay=)
/* DOCUMENT 
     displays 2*p cloud of point in a sequence
    KEYWORDS inter= keyboard interrupt ?
             sec= dimension to cut (default 3)
             fun= the function to use for plot (default PL)
             delay= 
    SEE ALSO: showsec,sec
 */
{
  if(is_void(fun)) fun=PL;
  if (is_void(delay)) delay=25;
  if(is_void(lim)) lim=[min([x,y](*)),max([x,y](*))];
  animate,1;
  for(i=1;i<=dimsof(x)(0);i++)
    {
      fma; limits,lim(1),lim(2),lim(1),lim(2),square=1; 
      fun,y(..,i),x(..,i);
      if(!is_void(inter)) { print,i;  rdline,prompt="hit RET or Enter to continue";}
  pause,delay;
    }
  animate,0; 
}







func sec3(z,y,x,sec=,inter=,vmin=,vmax=,delay=,fun=,noerase=,range=) 
/* DOCUMENT  sec3(z,y,x,sec=,inter=,vmin=,vmax=,delay=,fun=,noerase=,range=)
    displays n sections of nxnxn cube  
    KEYWORDS inter= keyboard interrupt ?
             sec= dimension to cut (default 3)
             fun= the function to use for plot (default pli)
             delay= 
    SEE ALSO: showsec
*/
  {
    local i,vmax;
    //window,style="style1.gs";
 animate,0;
if (is_void(vmax)) vmax=max(z);
 if (is_void(vmin)) vmin=min(z);
if (is_void(delay)) delay=25;
 iii=1;
if(is_void(sec)) sec=3; 
if (is_void(fun)) for(i=1;i<=dimsof(z)(sec+1);i++){ fma;
 if(sec==3) pli,z(,,i),cmax=vmin,cmax=vmax; 
                                  if(sec==2) pli,z(,i,),cmax=vmin,cmax=vmax;
                                  if(sec==1) pli,z(i,,),cmax=vmin,cmax=vmax; 
if(!is_void(inter)) { print,i;  rdline,prompt="hit RET or Enter to continue";}
 pause,delay; fma;}
 else for(i=1;i<=dimsof(z)(sec+1);i++){ 
 if(sec==3) b=z(,,i);
 if(sec==2) b=z(,i,);
 if(sec==1) b=z(i,,); 
 if (is_void(noerase)) { fma;  fun,b,y,x;
 } else { fun,b,y,x,type=iii++;} 
 if (!is_void(range)) range,range;
 if(!is_void(inter)) { print,i;  rdline,prompt="hit RET or Enter to continue";}
 pause,delay; }
 return void;
}


func wk(win)
/* DOCUMENT wk(win) wills window win (default current)
     
   SEE ALSO:
 */
{
  if (is_void(win) && (win= current_window())<0) win= 0;
window, win;
winkill;
}


func wkl(nmax=)
{
for(i=0;i<=(is_void(nmax)?15:nmax);i++) wk,i;
}


func ws(win)
/* DOCUMENT ws(win) cleans function win
   include logxy+animate
     
   SEE ALSO:ws1,WS,wsl 
 */
{
  if (is_void(win) && (win= current_window())<0) win= 0;
  window, win,wait=0;   fma;  animate, 0; 
  window, win,wait=0;  limits;  fma; logxy,0,0;
}



func pl3(x,width=,save=,opengl=,color=,erase=,delay=) {
  require,"pl3d.i"
  if(is_void(delay)) delay=10;
  if(is_void(width)) width=1;
  dd=dimsof(x);
  if(dd(0)>dd(-1)) x=transpose(x);
  get3_xy,transpose(x),x1,x2;
  dd=dimsof(x);
  if(erase)
    {
      limits,min(x),max(x),min(x),max(x);
      for(i=1;i<=dd(2);i++)
        {
          PL,x2(i),x1(i),msize=0.3*width,incolor=color,color=color;
          pause,delay;
          fma;
        }
    }
  else if(numberof(color)>1)
    plcolor,color,x2(*),x1(*),size=0.3*width,npol=15;
    else
      PL,x2(*),x1(*),msize=0.3*width,incolor=color,color=color;
    if(opengl) save=1;
  if(save) 
   {
     ff=open(".crap.ts","w");
     write,ff,"% Graphics3D objects";
     write,ff,"boundingbox";
     write,ff,"-0.03 -0.03 -0.03";
     write,ff,"1.03 1.03 1.03";
     write,ff,"viewpoint";
     write,ff,"1.3 -2.4 2.";
     write,ff,"ambientlight";
     write,ff,"0 0.5 0";
     write,ff,"lightsources";
     write,ff,"1.0 0. 1.";
     write,ff,"1 0 0";
     write,ff,"1. 1. 1.";
     write,ff,"0 1 0";
     write,ff,"0. 1. 1.";
     write,ff,"0 0 1";
     for(i=1;i<=numberof(x(,3));i++)
       {
         form="\npoint\n";
         write,ff,form;
         form=strjoin(array("%1.6e\t",3));
         write,ff,format=form+"\n",double(x(i,)),linesize=14*3+5;
       }
     close,ff;
   }
  if(opengl) system,"mathview3d -mvps 20 -mvfile .crap.ts &";
}


func pl2(x,width=,line=) {
 dd=dimsof(x);
  if(dd(0)>dd(-1)) x=transpose(x);
  if(line) plg,x(,2),x(,1),width=width;
  else pl,x(,2),x(,1),width=width;
}


func xytitle(x,y)
{
  if (is_void(y)) y="";
  return xytitles(x,y,[0.03,0.03]);
}



func WS(win,dpi=)
/* DOCUMENT WS(win) cleans function win and sets its dpi to dpi
   (default 100)
   include logxy+animate
     
   SEE ALSO:ws1,WS,wsl 
 */
{
  if (is_void(win) && (win= current_window())<0) win= 0;
  window,win,wait=0; winkill;
if (is_void(dpi))  dpi=100;
  width=long(6.54*dpi)+1;
  height=long(6.54*dpi)+1;
  window, win,wait=0,dpi=dpi,height=height,width=width;    fma;  animate, 0; 
  window, win,wait=0,height=height,width=width;    limits;  fma; logxy,0,0;
}

func WSL(win,dpi=)
/* DOCUMENT WS(win) cleans function win and sets its dpi to dpi
   (default 100) in portrait mode
   include logxy+animate
     
   SEE ALSO:ws1,WS,wsl 
 */
{
  if (is_void(win) && (win= current_window())<0) win= 0;
 if (is_void(dpi))  dpi=100;
  window,win,wait=0; winkill;


  width=long(10.54*dpi)+1;
  height=long(7.04*dpi)+1;
  
window, win,dpi=dpi,style=sdir+"Gist/large.gs",height=height,width=width;  limits; 
  fma;  animate, 0; 
  window, win,wait=0;  limits;  fma; logxy,0,0;
}


func ws1(win,dpi=)
/* DOCUMENT WS(win) cleans function win and sets its dpi to dpi
   (default 75) in portrait mode
   include logxy+animate
     
   SEE ALSO:ws1,WS,wsl 
 */
{
  while(catch(-1))
    {
      winkill,win;
    }
  if(is_void(dpi)) dpi=75;
  width=long(6.54*dpi)+1;
  height=long(6.54*dpi)+1;
  
  if (is_void(win) && (win= current_window())<0) win= 0;
window, win,dpi=dpi,width=width,height=height,style="bboxed.gs";  limits;  fma;
}

func ws0(win,dpi=)
/* DOCUMENT WS(win) cleans function win and sets its dpi to dpi
   (default 50) in portrait mode
   include logxy+animate
     
   SEE ALSO:ws1,WS,wsl 
 */
{
  while(catch(-1))
    {
      winkill,win;
    }
  if(is_void(dpi)) dpi=50;
  width=long(8.54*dpi)+1;
  height=long(8.54*dpi)+1;
  
  if (is_void(win) && (win= current_window())<0) win= 0;
window, win,dpi=dpi,width=width,height=height,style="bboxed.gs";  limits;  fma;
}


func ws12(win,dpi=)
/* DOCUMENT WS(win) cleans function win and sets its dpi to dpi
   (default 100) in 2 plsys per page mode
   include logxy+animate
     
   SEE ALSO:ws1,WS,wsl 
 */
{
  while(catch(-1))
    {
      winkill,win;
    }
  if(is_void(dpi)) dpi=75;
  width=long(6.54*dpi)+1;
  height=long(6.54*dpi)+1;
 if (is_void(win) && (win= current_window())<0) win= 0;
 // window, win;   fma;  animate, 0; 
window, win,dpi=dpi,width=width,height=height,style="win12.gs";  limits;  fma;
}


func ws21(win,dpi=)
/* DOCUMENT WS(win) cleans function win and sets its dpi to dpi
   (default 100) in 2 plsys per page mode
   include logxy+animate
     
   SEE ALSO:ws1,WS,wsl 
 */
{
  while(catch(-1))
    {
      winkill,win;
    }
  if(is_void(dpi)) dpi=75;
  width=long(6.54*dpi)+1;
  height=long(6.54*dpi)+1;
 if (is_void(win) && (win= current_window())<0) win= 0;
 // window, win;   fma;  animate, 0; 
window, win,dpi=dpi,width=width,height=height,style="win21.gs";  limits;  fma;
}


func ws22(win,dpi=)
/* DOCUMENT WS(win) cleans function win and sets its dpi to dpi
   (default 75) in 4 plsys per page mode
   include logxy+animate
     
   SEE ALSO:ws1,WS,wsl 
 */
{
  while(catch(-1))
    {
      winkill,win;
    }
  if(is_void(dpi)) dpi=75;
  width=long(6.54*dpi)+1;
  height=long(6.54*dpi)+1;
 if (is_void(win) && (win= current_window())<0) win= 0;
 // window, win;   fma;  animate, 0; 
window, win,dpi=dpi,width=width,height=height,style="win22.gs";  limits;  fma;
}


func ws32(win,dpi=)
/* DOCUMENT WS(win) cleans function win and sets its dpi to dpi
   (default 100) in 6 plsys per page mode
   include logxy+animate
     
   SEE ALSO:ws1,WS,wsl 
 */
{
  while(catch(-1))
    {
      winkill,win;
    }
  if(is_void(dpi)) dpi=100;
  width=long(6.54*dpi)+1;
  height=long(6.54*dpi)+1;
 if (is_void(win) && (win= current_window())<0) win= 0;
window, win,dpi=dpi,width=width,height=height,style="win32.gs";  limits;  fma;
}



func ws33(win,dpi=)
/* DOCUMENT WS(win) cleans function win and sets its dpi to dpi
   (default 100) in 6 plsys per page mode
   include logxy+animate
     
   SEE ALSO:ws1,WS,wsl 
 */
{
  while(catch(-1))
    {
      winkill,win;
    }
  if(is_void(dpi)) dpi=100;
  width=long(6.54*dpi)+1;
  height=long(6.54*dpi)+1;
 if (is_void(win) && (win= current_window())<0) win= 0;
window, win,dpi=dpi,width=width,height=height,style="win33.gs";  limits;  fma;
}



func wsn(win,dpi=)
/* DOCUMENT WS(win) cleans function win and sets its dpi to dpi
   (default 100) in 9 plsys per page mode
   include logxy+animate
     
   SEE ALSO:ws1,WS,wsl 
 */
{
  while(catch(-1))
    {
      winkill,win;
    }
  if(is_void(dpi)) dpi=75;
  width=long(6.54*dpi)+1;
  height=long(6.54*dpi)+1;
 if (is_void(win) && (win= current_window())<0) win= 0;
window, win,dpi=dpi,width=width,height=height,style="nobox.gs";  limits;  fma;
}



func wsl(n,display=,dpi=,wait=,private=,hcp=,dump=,legends=)
/* DOCUMENT WS(win) cleans function win and sets its dpi to dpi
   (default 75) in portrait mode
   include logxy+animate
     
   SEE ALSO:ws1,WS,wsl 
 */
{
  local height,width;
  if(is_void(dpi)) dpi=75;
  if (dpi==75)
    {
      width=828;
      height=640;
    }
  else
    {
      width=1103;
      height=874;
    }
  width=long(11.04*dpi)+1;
  height=long(8.54*dpi)+1;

  while(catch(-1))
    {
      winkill,n;
    }

window,n,style="large.gs",width=width,height=height,display=display,dpi=dpi,wait=wait,private=private,hcp=hcp,dump=dump,legends=legends;
  return;
}


func fman(nmax=)
{
  if(is_void(nmax)) nmax=30;
  for(i=0;i<=nmax;i++)
    {
      window,i; fma;
    }

}



func pler(y,x, dy=,dx=,marker=, width=, color=,incolor=,msize=,
          mean=,skew=,kurto=,med=,yweight=,xweight=,save=,noerror=,line=,logx=,logy=,type=,trimmed=,errorweight=
          ,shading=)
  /* DOCUMENT
   pler(y,x, marker=, width=, color=,incolor=,
   msize=,mean=,skew=,kurto=,med=)
   plot without fuss with error bars x and y
   can be either arrays or sorted lists
   mean= divides by sqrt(number of data points in bin)
   skew= computes skewness as well
   kurto= computes kurtosis as well
   med= computes median and quartile instead of mean and rms
   trimmed= trims 100*trimmed points on both side of the distribution
   yweight= can be a list a vector or an array of dimsof(y)
   xweight= itou
   noerror= doesn't plot error bar
   line= uses lines to join points
   logx= plots in log  without logxy
   logy=  itou
   note color=-1   and/or marker will pick random colors resp. marker
   SEE ALSO: splitstat,median,quartile
 */
{
  
  if(trimmed)
    {
      yy=xx=[];
      wwy=wwx=[];
      if(typeof(y)=="list")
          for(i=1;i<=_len(y);i++)
       {       yy=_cat(yy,trim(_car(y,i),wy,lcut=trimmed,cutedges=1));
       if(typeof(yweight)=="list") wwy=_cat(wwy,trim(_car(yweight,i)(wy),lcut=trimmed,cutedges=1));}
      else
      if(typeof(x)=="list")
          for(i=1;i<=_len(x);i++)
       {     xx=_cat(xx,trim(_car(x,i),wx,lcut=trimmed,cutedges=1)); 
       if(typeof(xweight)=="list") wwx=_cat(wwx,trim(_car(xweight,i)(wx),lcut=trimmed,cutedges=1));}
      else
      for(i=1;i<=dimsof(y)(2);i++)
        {
          yy=_cat(yy,trim(y(i,),wy,lcut=trimmed,cutedges=1)); 
        if(is_array(yweight)) wwy=_cat(wwy,yweight(wy));
        if (numberof(x))    if (dimsof(x)(1)==2)
          {
            xx=_cat(xx,trim(x(i,),wx,lcut=trimmed,cutedges=1));
            if(is_array(xweight)) wwx=_cat(wwx,xweight(wx));
          }
        }
      
      if (numberof(x)) if (dimsof(x)(1)==2)  x=xx;
      ss=pler(yy,x,marker=marker, width=width, color=color,incolor=incolor,msize=msize,
        mean=mean,skew=skew,kurto=kurto,med=med,yweight=wwy,xweight=wwx,save=save,noerror=noerror,
        line=line,logx=logx,logy=logy,type=type,errorweight=errorweight);
      
     if(save) return ss;
     return;
    }

  if( is_void(dy) && is_void(dx))
    {

      if (typeof(y)=="list")
        {
          if(typeof(yweight)=="list")
            {
              
              yw=yweight;  sw=[]; nn=_len(yw); for(i=1;i<=nn;i++){grow,sw,(1.*_nxt(yw))(sum); }
              w1=where(!sw); if(is_array(w1)) sw(w1)=1.;
              yy=y; yw=yweight; 
            y1=[]; for(i=1;i<=nn;i++){ grow,y1,(1./sw(i)*_nxt(yw)*_nxt(yy))(sum);}
             yy=y;yw=yweight;
            ye=[]; for(i=1;i<=nn;i++){ grow,ye,sqrt((1./sw(i)*_nxt(yw)*(_nxt(yy)-y1(i,-))^2)(sum));}
            if(!is_void(mean)) {
              yw=yweight;  yn=[]; for(i=1;i<=nn;i++){ grow,yn,numberof(_nxt(yw));}
              ye = ye/sqrt(yn); }
            }
          else
            {
          if(is_void(med)) y1=lst2arr(_map(AVG,y)); else {  y1=lst2arr(_map( median,y));};
          if(is_void(med)) ye=lst2arr(_map(RMS,y)); else ye=lst2arr(_map(quartile,y));
          yn=lst2arr(_map(numberof,y));
          if(!is_void(skew))     ys=lst2arr(_map(skewness,y));
          if(!is_void(kurto))     yk=lst2arr(_map(kurtosis,y));
          if(!is_void(mean)) ye = ye/sqrt(yn);
            }
        }
      else
        {
          if(!is_void(yweight))
            {
              if (dimsof(yweight)(1) == 1)
                {
                  yweight *= 1/sum(double(yweight));
                  if(is_void(med))
                    {
                      y1=(y*yweight(-,))(,sum);
                      ye=sqrt((yweight(-,)*((y-y1(,-))^2))(,sum));
                      if(!is_void(mean)) ye = sqrt(((yweight(-,))^2*((y-y1(,-))^2))(,sum));
                    }
                }
              else
                {
                  yweight *=1. /double(yweight(,sum));
                  y1=(y*yweight)(,sum);
                  ye=sqrt((yweight*((y-y1(,-))^2))(,sum));
                  if(!is_void(mean)) ye = sqrt((yweight^2*((y-y1(,-))^2))(,sum));
                }
            }
          else
            {
              if(is_void(med))
                {
                  ye=y(,rms);    y1=y(,avg);
                  if(!is_void(skew)) ys=lst2arr(_map(skewness,arr2lst(transpose(y))));
                  if(!is_void(kurto)) yk=lst2arr(_map(kurtosis,arr2lst(transpose(y))));

                }
              else
                {
                  y1=median(y,2);ye=quartile(transpose(y));
                }
              if(!is_void(mean)) ye = ye/sqrt(dimsof(y)(0));
            }
        }
      if (is_void(x)) x=1.*indgen(dimsof(y1)(0));

      if (typeof(x)=="list")
        {
          if(typeof(xweight)=="list")
            {
              xw=xweight;  sw=[]; nn=_len(xw); for(i=1;i<=nn;i++){ grow,sw,(1.*_nxt(xw))(sum);}
             w1=where(!sw); if(is_array(w1)) sw(w1)=1.;
             xx=x; xw=xweight;
            x1=[]; for(i=1;i<=nn;i++){ grow,x1,(1./sw(i)*_nxt(xw)*_nxt(xx))(sum);}
            xx=x;xw=xweight;
            xe=[]; for(i=1;i<=nn;i++){ grow,xe,sqrt((1./sw(i)*_nxt(xw)*(_nxt(xx)-x1(i,-))^2)(sum));}
            if(!is_void(mean)) {
              xw=xweight;  xn=[]; for(i=1;i<=nn;i++){ grow,xn,numberof(_nxt(xw));}
              xe = xe/sqrt(xn); }
            }
          else
            {
              if(is_void(med))       x1=lst2arr(_map(AVG,x)); else x1=lst2arr(_map(median,x));
              if(is_void(med))       xe=lst2arr(_map(RMS,x)); else xe=lst2arr(_map(quartile,x));
              xn=lst2arr(_map(numberof,x));
              if(!is_void(skew))     xs=lst2arr(_map(skewness,x));
              if(!is_void(kurto))     xk=lst2arr(_map(kurtosis,x));
            if(!is_void(mean)) xe = xe/sqrt(xn);
            }
        }
      else  if (dimsof(x)(1) > 1) 
        {
          if(!is_void(xweight)){
            if (dimsof(xweight)(1) == 1)
              {
                xweight *= 1/sum(double(xweight));
                if(is_void(med))
                  {
                    x1=(x*xweight(-,))(,sum);
                    xe=sqrt((xweight(-,)*((x-x1(,-))^2))(,sum));
                    if(!is_void(mean)) xe = sqrt(((xweight(-,))^2*((x-x1(,-))^2))(,sum));
                  }
              }
            else
              {
                xweight *=1. /double(xweight(,sum));
                x1=(x*xweight)(,sum);
                xe=sqrt((xweight*((x-x1(,-))^2))(,sum));
                if(!is_void(mean)) xe = sqrt((xweight^2*((x-x1(,-))^2))(,sum));
              }
          }
          else
            {
              if(is_void(med))
                {
                  x1=x(,avg);  xe=x(,rms);
                   if(!is_void(skew)) xs=lst2arr(_map(skewness,arr2lst(transpose(x))));
                  if(!is_void(kurto)) xk=lst2arr(_map(kurtosis,arr2lst(transpose(x))));
                }
              else
                {
                  x1=median(x,2);
                  xe=quartile(transpose(x));
                }
              if(!is_void(mean)) xe =xe/sqrt(dimsof(x)(0)); 
            }
        }
      else {
        x1=x; xe=0.;
      }
    } else{
      if(is_void(x)) x=indgen(numberof(y));
      if(!is_void(dx))  xe=dx; else xe=x*0;
      if(!is_void(dy))  ye=dy;  else ye=y*0;
      y1=y; x1=x;
    }
  if(errorweight) {ye *= errorweight; xe *= errorweight;}
  
  if (is_void(color)) color=incolor;
if(numberof(color)==1) if(color==-1) { color=Colors(8)(1+long(random()*6)); incolor=color;}

  if(!is_void(logx)) {
    ww=where((x1>0)*(x1>xe)*(x1+xe>0));
    if(is_array(xe>0))
      {xlo=log10((x1-xe)(ww));
      xhi=log10((x1+xe)(ww));} else
        {xlo=xhi=log10(x1);}
    x1=log10(x1(ww));
    y1=y1(ww);
  }
  if(!is_void(logy)) {
    ww=where((y1>0)*(y1>ye)*(y1+ye>0));
    if(is_array(xe>0) ){
      ylo=log10((y1-ye)(ww));
      yhi=log10((y1+ye)(ww));
    } else {ylo=yhi=log10(y1);}
    y1=log10(y1(ww));
    x1=x1(ww);
  }

  if(is_void(xlo)){ xlo=x1-xe; xhi=x1+xe;}
  if(is_void(ylo)){ ylo=y1-ye; yhi=y1+ye;}
  if(is_void(noerror)&is_void(shading)) 
 {   plp,y1,x1,xlo=xlo,xhi=xhi,ylo=ylo,yhi=yhi,ticks=1,color=color,width=width,type=1,symbol=0;
  PL,y1,x1,marker=marker,incolor=incolor,width=width,color=color,msize=msize;}
else
  if(is_void(noerror)&!is_void(shading))  plshade,[ylo,yhi],x1,color=color,edge=1;
 else
   PL,y1,x1,marker=marker,incolor=incolor,width=width,color=color,msize=msize;
  if(!is_void(line)) plg,y1,x1,width=width,color=color,type=type;
  if(!is_void(skew))
    {
      PL,y1+ye*ys,x1,marker=5,incolor=incolor,width=width,color=color,msize=0.4;
      if(numberof(xs))  PL,y1,x1+xe*xs,marker=5,incolor=incolor,width=width,color=color,msize=0.4;
      if(!is_void(save)& is_void(kurto)) return [x1,xe,y1,ye,ys];
    }
  if(!is_void(kurto))
    {
      PL,y1+yk*ye,x1,marker=6,incolor=incolor,width=width,color=color,msize=0.3;
      PL,y1-yk*ye,x1,marker=6,incolor=incolor,width=width,color=color,msize=0.3;
      if(numberof(xk))
        {
          PL,y1,x1+xe*xk,marker=6,incolor=incolor,width=width,color=color,msize=0.4;
          PL,y1,x1-xe*xk,marker=6,incolor=incolor,width=width,color=color,msize=0.4;
        }
      if(!is_void(save) & is_void(skew)) return [x1,xe,y1,ye,yk];
      if(!is_void(save)) return [x1,xe,y1,ye,ys,yk];
    }
  if(!is_void(save) & !is_void(noerror)) return [x1,y1];
  if(!is_void(save)) return [x1,xe,y1,ye];
}



func plshade(y,x,color=,edge=)
/* DOCUMENT 
     plots a shaded region between y(,1) and y(,2) as a function of x

     EXAMPLE
      x=span(-5,5,25);y1=sin(x); y2=sin(x)+1*cos(x); y=[y1,y2];
     plshade,y,x,color=__rgb(,20)
   SEE ALSO:
 */
{
  n=dimsof(y)(2);
  if(is_void(x)) x=indgen(n);
  if(is_void(color)) color=__gold2;
  color=char(color);
  xx=x(,-:1:n);
  y1=y(,1);
  y2=y(,2);
  yy=transpose(span(0,1,n)(,-:1:n))*((y2-y1+0*(y2-y1)(ptp)/2.)(,-))+y1(,-);
  plf,array(color,dimsof(xx)),yy,xx,edges=0;
  // if(edge){ plg,y1,x,color=__black,type=3; plg,y2,x,color=__black,type=3;}
}



func pl(y,x,legend=, hide=, marker=, width=, color=)
/* DOCUMENT  pl(y,x,legend=, hide=, marker=, width=, color=)
     point plot without fuss
     SEE ALSO: plp,plg,plh
*/
{ if(is_void(marker)) marker=1;
if(is_void(width))
  { plg,y,x, marker=marker,type=0,color=color,legend=legend,hide=hide;}
else
{  plmk,y,x, color=color,msize=width/10.,width=width;}
}


func PL(y,x, marker=, width=, color=,incolor=,msize=,line=,type=)
/* DOCUMENT
   PL(y,x, marker=, width=, color=,incolor=,msize=) point plot without
   in colour
   SEE ALSO: pl
*/
{ if(is_void(incolor)) incolor="white";
   if(is_void(width)) width=3;
   if(is_void(msize)) msize=0.75;
   if(is_void(marker)) marker=4;
   if(marker<0) marker=[];

   if(is_void(color)) color="red";
   if(numberof(color)==1)if(color==-1) { color=Colors(8)(2+long(random()*6)); incolor=color;}
   plmk2,y,x, marker=marker,incolor=incolor,width=10,color=color,msize=msize;
   if(line)
     {
       if(type){ if(type<0) type1=1+long(random()*6); else type1=type; } 
     plg,y,x, width=width,color=color,type=type1;
     }
}



func PLB(y,x, marker=, width=, color=,incolor=,msize=,ofset=)
/* DOCUMENT
   PLB(z,y, marker=, width=, color=,incolor=,msize=,ofset=)
   will plot the bundle z versus y and x using different symbols for
   each y
   KEYWORDS fun alows the user to specify another function instead of plg
   like pla but in color and symbols
   SEE ALSO: plb,pla,thread
*/
{ if(is_void(incolor)) incolor=-5;
   if(is_void(width)) width=10;
   if(is_void(ofset)) ofset=0.;
   if(is_void(msize)) msize=0.75;
   if(is_void(color)) color="blue";
   n=dimsof(y)(0); 
   ii=0;
   if (!is_void(x))
   for(i=1;i<=n;i++){
  plmk2,y(,i)+i*ofset,x(,i),incolor=incolor-i,width=width,color=color,msize=msize;
   } else
   for(i=1;i<=n;i++){
  plmk2,y(,i)+i*ofset,incolor=incolor-i,width=width,color=color,msize=msize;
   }}


func __rgbPlot(void,odd=)
/* DOCUMENT __rgbPlot(void) displays the 40 colours
     
   SEE ALSO:
 */
{
  ws;
  for(i=1;i<=40;i++)
    {
      odd=(is_void(odd)?0:odd);
      if(odd==1) if (i%2) continue;
      if(odd==2) if ((i+1)%2) continue;
      plh,[0,1],[i,i],color=__rgb(,i),width=20;
    }
}


func cutRegion(X)
/* DOCUMENT cutRegion(X) cuts selected subregion  from array X
   via the mouse 
*/
{ local fl,pts,x,y,Y;
 fl=1; x= []; y=[];
 pli,X;
  write,"press 1st button to define rec"+
    " other button to accept" ;
  while (fl) {
     res= mouse(1,1);
    grow,x,res(1);
    grow,y,res(2);
    grow,x,res(3);
    grow,y,res(4);
     //     if (res(10)>1) fl =0;   
  Y= X(long(min(x)):long(max(x)),long(min(y)):long(max(y)));
   fl=0;
   }
  fma;  limits;
  pli,Y;
  print,long(min(x)):long(max(x)),long(min(y)):long(max(y));
  return Y;
}



func centile(x,lcut=,ucut=)
/* DOCUMENT  centile(x,lcut=,ucut=)
   returns upper and lower centile (default 5 %)
   modified by Pierre 16/02/2011
   SEE ALSO: quartile,trim
*/
{
  if (is_void(lcut)) lcut=0.05;
  if (is_void(ucut)) ucut=lcut;
  if (lcut > 1-ucut) {write,"error in cut"; return 0;
}
  s =sort(x,1);
  n=dimsof(x)(2);
  //  s1=s(1+long(lcut*n):n-long(ucut*n),);
  s1=s(1+int(lcut*n):n-int(ucut*n),);

  return transpose([x(s1)(min,),x(s1)(max,)]);
}


func quartile(x)
/* DOCUMENT quartile(x) quartile of x
   SEE ALSO: centile
   //# note that the interquartile equals twice the quartile here #//
*/
{
  tt= centile(x,lcut=0.25,ucut=0.25);
  return (tt(2,)-tt(1,))/2.;
}

func trim(x,&ww,lcut=,ucut=,cutedges=)
/* DOCUMENT trim(x,lcut=,ucut=,cutedges=)
   returns  x trimmed at upper and lower band (default 5 %)
   //#  NOTE that numberof(ftrim(x))=numberof(x) by default #//
   use cutedges=1 to avoid this effect.
   SEE ALSO: centile,quartile
*/
{
  local x0,x1,s,n,m1,m2,s1;
 if (is_void(lcut)) lcut=0.05;
  if (is_void(ucut)) ucut=lcut;
  if (lcut >= 1-ucut){ error,"bound in cut incorrect";  }
  x0= x*0;
  x1=x(*);
  s =sort(x1);
  n=dimsof(x1)(2); 
s1=s(1+long(lcut*n):n-long(ucut*n));
 m2=max(x1(s1)); 
 m1=min(x1(s1));
  // m=s(long(n/2));
  x1= max(min(x,m2),m1);
  x0(*)=x1(*);
  if(is_void(cutedges))  return x0; else
    { ww=where((x>=m1)*(x<=m2));
      return x(ww);
    }
}


func ftrim(x,lcut=,ucut=,cutedges=)
/* DOCUMENT
   ftrim(x,lcut=,ucut=,cutedges=)
   returns  x statistically trimmed at upper and lower band (default 0.5 %)
   //#  NOTE that numberof(ftrim(x))=numberof(x) by default #//
   use cutedges=1 to avoid this effect.
   SEE ALSO: trim,quartile,centile
*/
{  require,"Eric/histo.i";
  
  local x0,x1,s,n,m1,m2,s1;
 if (is_void(lcut)) lcut=0.005;
  if (is_void(ucut)) ucut=lcut;
  if (lcut >= 1-ucut){ error,"bound in cut incorrect";  }
  x0= x*0;
  x1=double(x(*));
  // x1= x1(where(x1>0)); // positivity ?
  h= histo2(x1,h1,binsize=(max(x1)-min(x1))/200.);
  ch=integ(h,h1,h1);
  ch/=float(max(ch));
  c1=where(long((ch-lcut)*100)>=0)(1);
  c2=where(long((ch-(1-ucut))*100)<=0)(0);
  //   plh,h(c1:c2),h1(c1:c2);
  m2=h1(c2);
  m1=h1(c1);   
  // m=s(long(n/2));
  x1= max(min(x,m2),m1);
  x0(*)=x1(*);
  if(is_void(cutedges))  return x0; else
    { ww=where((x>=m1)*(x<=m2));
      return x(ww);
    }
}






func spinterp(zz,r=,th=,ph=,n=)
/* DOCUMENT  spinterp(zz,r=,th=,ph=,n=) sphere cut in unit cube
   KEYWORDS r= radius at which to cut
*/
{  require, "linterp.i";
  
  local x1,y1,z1,v;
 if (is_void(r)) r=1/2.;
if (is_void(n))  n=dimsof(zz)(2);
 p=dimsof(zz)(2);
  if (is_void(th)) th=span(0,pi,n)(-:1:n,);
if (is_void(ph))  ph= span(0,2*pi,n)(,-:1:n);

  x1=r*cos(ph)*sin(th);
 y1=r*sin(ph)*sin(th);
 z1=r*cos(th); 
v=LInterp(zz,x1,y1,z1,x=span(-1,1,p),y=span(-1,1,p),z=span(-1,1,p));
return v;
}






func select(Y,X,&pth,path=,cont=,rest=,xonly=,yonly=,logs=,color=)
/* DOCUMENT    w=select (Y,X,path=) 
   returns w, the index of (X,Y) contained in the path 
   which is selected via the mouse 
   if path=1 the path need not be rectangular
   EXAMPLE
   Y= random(250); X= random(250);
   fma; plp,Y,X;
   w= select(Y,X,path=1); plp,Y(w),X(w);
   ww=select(yy,xx,logs="1") to select in logxy,1,1
   ww=select(yy,xx,logs="x") to select in logxy,1,0
   ww=select(yy,xx,logs="y") to select in logxy,0,1
   SEE ALSO: cut
*/
{ require, "Chris/scan-button.i";

  local fl,pts,x,y,x1,y1,mx,my,Mx,My,cx,cy,p,z,Z,c2x,c2y,w;
 fl=1; pts= []; w=[];
 ll=limits();
 if (is_void(color)) color=-6;
 
 /*
 if(ll(0)==15) logs=[];
 if(ll(0)==399) logs=1;
 if(ll(0)==271) logs="y";
 if(ll(0)==143) logs="x";
 to get automatically the logxy flag
 */
 if(logs==1)
   {
     curwin=window();
     cpp,curwin,35;
     window,35; w1= where(( Y>0)&(X >0));
     ws; logxy,0,0; pl,log10(Y(w1)),log10(X(w1)),color=color,width=3; pause,0;
      ww=select(log10(Y(w1)),log10(X(w1)),&pth,path=path,cont=cont,rest=rest,xonly=xonly,yonly=yonly);
      wk,35;  window,curwin; pl,Y(w1(ww)),X(w1(ww)), width=3,color=-6; return w1(ww);
   }
 if(logs=="x")
   {
     cpp,window(),35;
     window,35;w1= where((X >0));
     ws; logxy,0,0; pl,Y,log10(X),color=color,width=3; pause,0;
      ww=select(Y(w1),log10(X(w1)),&pth,path=path,cont=cont,rest=rest,xonly=xonly,yonly=yonly);
      wk,35; window,curwin; pl,Y(w1(ww)),X(w1(ww)), width=3,color=color; return w1(ww);
   }
 if(logs=="y")
   {
     cpp,window(),35;
     window,35;w1= where(( Y>0));
     ws; logxy,0,0; pl,log10(Y),X,color=color,width=3; pause,0;
      ww=select(log10(Y(w1)),X(w1),&pth,path=path,cont=cont,rest=rest,xonly=xonly,yonly=yonly);
      wk,35; window,curwin; pl,Y(w1(ww)),X(w1(ww)), width=3,color=color; return w1(ww);
   }

 if(!numberof(cont)){
   write,"press shift 1st button to select point shift 3nd button to close path";
   XX=YY=[];
   for (;;) {
   ms = mouse(-1,2,"");
         x1=ms(1);
         y1=ms(2);
         if (is_void(ms)) break;
         if ((ms(-1)==3)&(ms(0)==1)) break;
         if ((ms(0)==1))
           {
             grow,XX,x1;
             grow,YY,y1;
             pl,y1,x1,width=3,marker=2,color=-8; pause,0;
           }
        else   gg_zoom, ms;
 }
 pts=transpose([XX,YY]);
 n=numberof(pts)-2;
  }
 else{ pts=cont; n=numberof(pts);}
 pts =pts(1:n);
 grow,pts,pts(1);
 grow,pts,pts(2);
 y= pts(2:n+2:2); x= pts(1:n+2:2);
 mx= min(x); Mx= max(x);
 my= min(y); My= max(y);
 pth=[[mx,Mx],[my,My]]; 
   write,mx,Mx,my,My,format= "range= [[%f,%f],[%f,%f]] \n";
x=interp(x,indgen(numberof(x)),span(1,numberof(x),max(500,n*3)));
y=interp(y,indgen(numberof(y)),span(1,numberof(y),max(500,n*3)));
 if(xonly) w =where( (X>mx)*(X<Mx)); else
 if(yonly) w =where( (Y>my)*(Y<My)); else
 if(rest) w=where( !((X>mx)*(X<Mx)*(Y>my)*(Y<My))); else
 w =where( (X>mx)*(X<Mx)*(Y>my)*(Y<My));
 if (!numberof(w)) {write,"found no match"; return w;}
 if (!path) {
   plg,[my,my,My,My,my],[mx,Mx,Mx,mx,mx],color="red";
   c2y=Y(w); c2x=X(w);
plg,c2y,c2x,type=0,marker=2,color=color;
return w;
}
plg, pts(2:n+2:2),pts(1:n+2:2),color=color;
cx =X(w); cy =Y(w);
p=numberof(w);
c =cy*1.i+cx;
z =y*1.i+x;
c2x=[]; c2y=[]; 
for(i=1;i<=p;i++)
{
Z=(1/(c(i)-z));
Z=(Z(1:-1)*(z(dif)))(sum); 
if(abs(abs(im(Z))-2*pi) <0.15) { grow,c2y,cy(i); grow,c2x,cx(i);  }
}
write,numberof(c2x);
 if (numberof(c2x)) plg,c2y,c2x,type=0,marker=2,color=color; 
 w= []; 
 for(i=1;i<=numberof(c2x);i++){
 tt =where(abs(c2x(i)-X,c2y(i)-Y)==0)(1);
if (!is_void(tt)) grow,w,tt;
 } 
 return w;

 

}



func pdf(str)
/* DOCUMENT    pdf(str) 
   produces a pdf output and displays it
*/

{
hcp_file,str;
hcp;
hcp_finish;
system,"ps2pdf "+str;
system,"acroread "+str+":r.pdf &";
}




func eps(str,a=,dir=,nodisp=,norescale=)
/* DOCUMENT    eps(str) 
   produces a ps output and displays it

   to produce multiple output per page use a flag
       psnup -pa4 -8 -l -Pb5 -l test.ps > !test2.ps ; gv -a4 test2.ps 
*/
{
   if (is_void(dir)) dir ="./";
   dir=dir+"/";
   if (is_void(a)){
     hcp_file,dir+str;
     hcp;
hcp_finish;
if (is_void(nodisp)) system,"gv "+dir+str+"&";
  } 
  else {hcp_file,dir+str; for(i=0;i<=a;i++){ window,i; hcp; } hcp_finish;
if (is_void(nodisp))  system,"gv "+dir+str+"&";
  }
if (is_void(norescale)) rescaleboundingbox(dir+str);
}





func jpegread(fname,rgb=)
/* DOCUMENT    jpegread(str) 
   reads a jpg file and displays it
   returns a 2D image or 3 2D images
   corresponding to the RGB colors
*/
{       require,"pnm.i";
       system,"convert "+fname+" "+fname+".pnm";
       image=pnm_read(fname+".pnm");
       v=pnm_display(image);
       system, "touch "+fname+".pnm; rm  "+fname+".pnm";
       if (is_void(rgb))       return v;
       else return image;
     }


func jpegwrite(image,fname)
/* DOCUMENT    jpegwrite(image,fname) 
   write image into jpg file fname image is
   either a 2D image or 3 2D images
   corresponding to the RGB colors
*/
{
       require,"pnm.i";
         er=pnm_write(image,".tmp.pnm");
       system,"convert .tmp.pnm "+fname;
       system, "touch .tmp.pnm; rm .tmp.pnm";
       return er;
}



func jpeg(str,nodisp=,hires=,dir=)
/* DOCUMENT    jpeg(str,nodisp=,hires=,dir=)
   produces a jpg output of current windows
   saves it into str and and displays it
*/
{
local res;
  if (is_void(dir)) dir ="./";
  dir=dir+"/";
  str=_filterdot(str);
  hcp_file,dir+str+".ps";
hcp;
hcp_finish;
if (!is_void(hires)) 
{
if(hires==1) system,"convert -density 144x144 "+dir+str+".ps "+dir+str+"";
if(hires==2) system,"convert -density 288x288 "+dir+str+".ps "+dir+str+"";
if(hires==3) system,"convert -density 576x576 "+dir+str+".ps "+dir+str+"";
}
else system,"convert  "+dir+str+".ps "+dir+str+""; 
system," touch "+dir+str+".ps; rm "+dir+str+".ps "
if (is_void(nodisp)) system,"ee "+dir+str+" &";
return ;
}







func add_txt(str,pos=, color=,sys=, font=, height=, opaque=, orient=, justify=)
/* DOCUMENT    add_txt(str) 
   puts text where clicked
   xytitles,"hello","hello",[0.03,0.03]
   pltitle,"hello"
   SEE ALSO:putt_txt,xytitles,pltitle
*/
{
local pos;
if(is_void(sys)) sys= 0;
if(is_void(height)) height= 18;
 plsys,sys; 
 if(is_void(pos) &(sys==0)) pos= mouse()(5:6);
 if(is_void(pos) &(sys>0)) pos= mouse()(1:2);
 plt,str,pos(1),pos(2),tosys=sys,color=color,
  font=font, height=height, opaque=opaque, orient=orient, justify=justify;
return [pos(1),pos(2)];
}


func add_label(str,color,pos=,sys=, font=, height=, opaque=, orient=, justify=,marker=,type=,line=,msize=,width=,length=)
/* DOCUMENT    add_label(str,color) 
   puts text where clicked
   xytitles,"hello","hello",[0.03,0.03]
   pltitle,"hello"

   if MARKER <0 no marker is display
   
   SEE ALSO:putt_txt,xytitles,pltitle
*/
{
local pos;

 
if(is_void(sys)) sys= 0;
if(is_void(height)) height= 18;
if(is_void(pos)) pos= mouse()(5:6);
if(is_void(msize)) msize= 0.75;
 if(is_void(marker)) marker= 4;
 oldsys=plsys(); 
 plsys,sys;
 pos=add_txt(str,sys=sys,color=color,pos=pos,
  font=font, height=height, opaque=opaque, orient=orient, justify=justify);
 
 if(length)   dlx=length;
 if(sys==0) {dlx=0.05; dly=0.0075;} else { dlx=limits()(1:2)(dif)*0.1; dly=limits()(3:4)(dif)*0.005;}
 plsys,sys;
 if(!is_void(line)) plg,[pos(2)+dly,pos(2)+dly],[pos(1)-dlx,pos(1)-dlx/5.],width=5,color=color,type=type;
 if(marker>0) PL,pos(2)+dly,pos(1)-dlx/5.,
              color=color,marker=marker,incolor=color,msize=msize,width=width;

 plsys,oldsys;
}



func get_pts(void)
/* DOCUMENT    xy=get_pts() 
   returns xy, the points defining the path
   SEE ALSO: select,cut
*/
{
local fl,pts,x,y,res,n;
   fl=1; pts= []; w=[];
 write,"press 1st button to select point 2nd button to close path";
while (fl) {
res= mouse();
grow,pts,transpose([res(1),res(2)]);
plg,[res(2),res(2)*1.0001],[res(1),res(1)*1.0001],width=5,color="blue";
if (res(10)==2) fl =0; 
}
n=numberof(pts)-2;
pts =pts(1:n);
//grow,pts,pts(1);
//grow,pts,pts(2);
y= pts(2:n:2); x= pts(1:n:2);

 return transpose([x,y]);
}


func blockmat(a,b,c,d){
/* DOCUMENT blockmat(a,b,c,d)
   defines matrix by block
*/
  n =dimsof(a)(2);
  p =dimsof(a)(3);
  res =array(a(1,1)*0,2*n,2*p);
  for(i=1;i<=n;i++)
    {
      for(j=1;j<=n;j++)
        {
          res(i,j)=a(i,j);
          res(i+n,j+p)=d(i,j);
          res(i+n,j)=b(i,j);
          res(i,j+n)=c(i,j);
        }
    }
return transpose(res);
}                      


func plmk2(y,x,marker=,width=,color=,msize=,incolor=)
/* DOCUMENT plmk, y,x

     Make a scatter plot of the points Y versus X.  If X is nil,
     it defaults to indgen(numberof(Y)).  By default, the marker
     cycles through 7 predefined marker shapes.  You may specify a shape
     using the marker= keyword, line width using the width= keyword (you
     get solid fills for width>=10), color using the color= keyword.
     You can also use the msize= keyword to scale the marker (default
     msize=1.0).  You can change the default width, color, or msize
     using the plmk_default function.

     The predefined marker= values are:

     marker=
       1	square
       2	cross
       3	triangle
       4	circle
       5	diamond
       6	cross (rotated 45 degrees)
       7	triangle (upside down)

     You may also put marker=[xm,ym] where xm and ym are vectors
     of NDC coordinates to design your own custom marker shapes.

   SEE ALSO: plmk_default, plg (type=0 keyword)
*/
{
  if (is_void(marker)) {
    marker= (_plmk_count-1)%7 + 1;
    _plmk_count++;
  }
  if (numberof(marker)==1) {
    marker= *_plmk_markers(marker);
  } else if (dimsof(marker)(1)!=2 || dimsof(marker)(3)!=2 ||
	     dimsof(marker)(2)<=2) {
    error, "illegal marker= keyword value";
  }
  xm= marker(,1);
  ym= marker(,2);
  if (is_void(msize)) msize= _plmk_msize;
  if (!is_void(msize)) {
    xm*= msize;
    ym*= msize;
  }
  if (is_void(color)) color= _plmk_color;
  ecolor= color;
  if (structof(color)==string) {
    n= where(color==["bg","fg","black","white",
		     "red","green","blue","cyan","magenta","yellow"]);
    if (numberof(n)!=1) error, "unrecognized color name: "+color;
    color= char(-n(1));
  }

  if (is_void(incolor)) incolor= ecolor;
  if (structof(incolor)==string) {
    n= where(incolor==["bg","fg","black","white",
		     "red","green","blue","cyan","magenta","yellow"]);
    if (numberof(n)!=1) error, "unrecognized color name: "+incolor;
    incolor= char(-n(1));
  }
  incolor=char(incolor);


  if (is_void(width)) width= _plmk_width;
  if (!is_void(width)) {
    if (width>=10) {
      solid= 1;
      z= array(incolor, 1+numberof(y));
      width= [];
    }
  }
  n= array(1,1+numberof(y));
  n(1)= numberof(ym);
  if (is_void(x)) x= indgen(numberof(y));
  plfp, z,grow(ym,y),grow(xm,x),n,edges=1,ewidth=width,ecolor=ecolor;
}



func zload(fname)
/* DOCUMENT  zload(fname)
   loads ascii compressed files
*/
{
  if (sread(rdline(popen("zcat "+ fname+"; ", 0)), sz) != 1)
    error, "cannot decompress tmp file ";
  return res;
}





func strjoin(str, glue)
/* DOCUMENT strjoin(str)
       -or- strjoin(str, glue)
     Join strings from array STR into a single long string.  The string GLUE
     (default "") is used between each pair of element from STR.

   SEE ALSO: strcut
*/
{
  if ((n= numberof(str)) >= 1) {
    s= str(1);
    if (glue) for (i=2 ; i<=n ; ++i) s+= glue+str(i);
    else      for (i=2 ; i<=n ; ++i) s+= str(i);
    return s;
  }
}



func struct2string(struc,separ=)
/* DOCUMENT struct2string(struc)
      return a string array corresponding to the content of
      the structure struc
   SEE ALSO: struct2name,struct2type
*/
{
  print_format,10000;
  if (typeof(struc) != "struct_instance") { u=print(struc);
  if ( typeof(struc) == "list") u="list";
  if ( typeof(struc) == "pointer") u="*";
    return u; }
  ww=(lst2arr(_map(typeof,struct2list(struc)))=="list")(sum);
 if (!ww)
    {
      nn=struct2name(struc);
      v=[];
      u=print((struc)); u=strjoin(u);
      u=strreplace(u,"\"","");
      u=strreplace(u,"(nil)","nil");
      u=strreplace(u,"_","-");
      uu=splittok(u,tok="=");
      for(i=2;i<=numberof(uu)-1;i++)
        grow,v,splittok(uu(i),tok=",")(1);
          //strpart(strtok(uu(i),strpart(nn(i),1:1))(1),:-1);
      grow,v,strpart(strtok(uu(i),")")(1),:0);
    }
  else
    {
      v=_map(struct2string,struct2list(struc));
      v=lst2arr(v);
     if(is_void(separ)) v=strjoin(v," "); else  v=strjoin(v,separ);
    }
  print_format,100;
  return v;
}


func splittok(s,tok=)
/* DOCUMENT      splittok(s,tok=)
   makes a list out of splitting in bits the string s
   SEE ALSO:
*/
{
  local r;
  if(is_void(tok)) tok=" ";
  r=[];
  do {s=strtok(s,tok);grow,r,s(1);s=s(2);} while(s);  
  return r;
}

func struct2list(str)
/* DOCUMENT struct2list(str)
     return a list  corresponding to the content of
     the structure str
      SEE ALSO: struct2name,struct2type,struct2string
*/
{
  nn=struct2name(str);
  v=[];
  for(i=1;i<=numberof(nn);i++)
    {
      tt=get_member(str,nn(i));
      if(typeof(tt)!="struct_instance") v=_cat(v,_lst(tt)); else
        {
           v=_cat(v,_lst(struct2list(tt)));
        }
    }
  return v;
}



func struct2name(str)
/* DOCUMENT struct2name(str)
      return a string array corresponding to the name of
      the fields of      the structure str

      to be used in conjunction with e.g.
      uu=struct2name(str);
      write,"%"+strjoin(uu,"\n %")
      
   SEE ALSO: struct2name,struct2type,struct2list
*/
{
  v=[];
  u=print(structof(str));
  for(i=2;i<=numberof(u)-1;i++)
    grow,v,strtok(strpart(strtok(u(i))(2),:-1),"(")(1);
  return v;
}


func struct2type(str)
/* DOCUMENT struct2type(str)
      return a string array corresponding to the type of
      the structure str
   SEE ALSO: struct2name,struct2list
*/
{
  v=[];
  u=print(structof(str));
  for(i=2;i<=numberof(u)-1;i++)
    grow,v,strpart(strtok(u(i))(1),:0);
  return v;
}









func struct2html(H,fname=,dir=,add=)
/* DOCUMENT struct2html(fname,H,dir=)
	saves in fname a html table describing the structure H.
  SEE ALSO: struct2name,struct2list
*/

{
  if (is_void(dir)) dir ="./"; else dir+="/";
  if (is_void(fname)) fname ="crap.html";
  if(is_void(add))  ff=open(dir+fname,"w"); else ff=open(dir+fname,"a");
  write,ff,"<BODY>\n";
  write,ff,"<TABLE FRAME=BELOW CELLSPACING=0 RULES=GROUPS BORDER=1>\n ";
  write,ff,"<TBODY><TR>\n";
  write,ff,"<TD>\n "+strjoin(struct2name(H(1)),"</TD>\n <TD>")+"</TD>";
  write,ff,"</TBODY></TR>\n";
    for(i=1;i<=numberof(H);i++){
      write,ff,"<TBODY><TR>\n";
            write,ff,"<TD>\n "+strjoin(struct2string(H(*)(i),separ="</TD>\n <TD>"),"</TD>\n <TD>")+"</TD>";
      write,ff,"</TBODY></TR>\n";
    }

    
  write,ff,"\n</TABLE>";
  write,ff,"\n</BODY>";
  close,ff; return fname;
}





func struct2txt(H,fname=,dir=,add=)
/* DOCUMENT struct2txt(fname,H,dir=)
	saves in fname a txt  describing the structure H.
  SEE ALSO: struct2name,struct2list,smwrite
*/

{
  if (is_void(dir)) dir ="./"; else dir+="/";
  if (is_void(fname)) fname ="crap.txt";
  if(is_void(add))  ff=open(dir+fname,"w"); else ff=open(dir+fname,"a");
  write,ff,"#"+strjoin(struct2name(H(1)),"\t"),linesize=10000;
    for(i=1;i<=numberof(H);i++){
            write,ff,strjoin(struct2string(H(*)(i),separ="\t"),"\t")+"",linesize=10000;
    }

    
  write,ff,"";
  write,ff,"";
  close,ff; return fname;
}



func struct2latex(H,fname=,dir=,add=,caption=)
/* DOCUMENT struct2latex(H,fname=,dir=,add=)
	saves in fname a latex table describing the structure H.
        dir=
        add=
        caption= specifies a string caption 
  SEE ALSO: struct2name,struct2list,struct2html
*/

{
  if (is_void(dir)) dir ="./"; else dir+="/";
  if (is_void(fname)) fname ="crap.tex";
  if(is_void(add))   {
    ff=  open(dir+fname,"w"); 
write,ff,"\\documentclass{article}";
write,ff,"\\usepackage[T1]{fontenc}";
write,ff,"\\usepackage[latin1]{inputenc}";
write,ff,"\\usepackage{graphics}";
write,ff,"\\usepackage{times}";
write,ff,"\\textwidth=16.5cm";
write,ff,"\\textheight=27.5cm"
write,ff,"\\topmargin=-3.cm";
write,ff,"\\evensidemargin=-2.5cm";
write,ff,"\\oddsidemargin=-2.5cm";
write,ff,"\\begin{document}";
  }
    else ff=  open(dir+fname,"a"); 
		  write,ff,"\\begin{table}[!hbp] \\begin{tabular}{"+strjoin(array("c",numberof(struct2name(H(1))))," ")+"}";
write,ff,""+strjoin(struct2name(H(1))," &")+"\\cr";
      write,ff,"\\hline"
    for(i=1;i<=numberof(H);i++){
             write,ff,""+strjoin(struct2string(H(*)(i),separ="&"),"& ")+"";
      write,ff,"\\cr \n";
      if (((i %45) ==0) && (i != numberof(H)))
        {
        write,ff,"\\end{tabular} ";
        if(!is_void(caption)) {
          write,ff,"\\caption{"+caption+"$\\ldots$ continued }";
        }
        write,ff,"\\end{table} \\eject";
        write,ff,"\\begin{table}[!hbp] \\begin{tabular}{"+strjoin(array("c",numberof(struct2name(H(1))))," ")+"} ";
write,ff,""+strjoin(struct2name(H(1))," &")+"\\cr";
      write,ff,"\\hline"
        }
    }
write,ff,"\\end{tabular}"
         if(!is_void(caption)) {
          write,ff,"\\caption{"+caption+"}";
        }
write,ff, "\\end{table}";

  write,ff,"\\vspace{0.3cm}";
write,ff,"\\end{document}";
   close,ff; return fname;
}


		  
func arr2latex(arr,fname=,dir=,add=)
/* DOCUMENT arr2latex(arr,fname=,dir=,add=)
	saves in fname a latex table arr.
  SEE ALSO: struct2name,struct2list,struct2html
*/

{
  if (is_void(dir)) dir ="./"; else dir+="/";
   if (is_void(fname)) fname ="crap.tex";
 if(is_void(add))   {
    ff=  open(dir+fname,"w"); 
write,ff,"\\documentclass{article}";
write,ff,"\\usepackage[T1]{fontenc}";
write,ff,"\\usepackage[latin1]{inputenc}";
write,ff,"\\usepackage{graphics}";
write,ff,"\\usepackage{times}";
write,ff,"\\textwidth=16.5cm";
write,ff,"\\textheight=27.5cm"
write,ff,"\\topmargin=-3.cm";
write,ff,"\\evensidemargin=-2.5cm";
write,ff,"\\oddsidemargin=-2.5cm";
write,ff,"\\begin{document}";
  }
    else ff=open(dir+fname,"a");
write,ff,"\\begin{tabular}{"+strjoin(array("r",dimsof(arr)(2))," ")+"}";
for(i=1;i<=dimsof(arr)(3);i++){
             write,ff,strjoin(map(pr1,arr(,i)),"& ");
      write,ff,"\\cr \n";
    }
write,ff,"\\end{tabular}";

  write,ff,"\\vspace{0.3cm}";
 write,ff,"\\end{document}";
 close,ff; return fname;
}

		  
 func dumppdf(lst,fname=,pp=,dir=)
/* DOCUMENT dumppdf((lst,fname=,pp=)
     makes a pdf file of the graphics windows in the list
   Keywords 
if lst is a 2D array or a structure it will produce a table in latex
		    default fname is /tmp/crap.pdf
  pp represents the number of widow per line
   SEE ALSO: eps,pdf,jpg
*/
{
  dd=get_cwd();
  if(is_void(dir)) dir="/tmp/"; 
 cd,dir;
  if(is_void(fname)) fname="crap";
 tt=splittok(fname,tok=".");
 if(numberof(tt)>1) fname=tt(-1);
if(typeof(lst)=="struct_instance") {
  
  struct2latex(lst,fname=fname+".tex",dir=dir);
} else if (dimsof(lst)(1)<=1)
  {
    if(is_void(pp)) pp=2;
    nnl=[];	
    if (numberof(lst)==1)  lst=[lst];
   for(i=1;i<=numberof(lst);i++)
      {
        nn=fname+pr1(lst(i))+".jpg"; 
        grow,nnl,nn;
        window,lst(i); jpeg,nn,nodisp=1,dir="/tmp/",hires=1;
        err=exportfig(nnl,pp=pp,fname=dir+fname+".tex");
      }
  } else if(dimsof(lst)(1)==2)
    {
      arr2latex(lst,fname=fname+".tex",dir=dir);
    } 
	
 system,"pdflatex "+dir+"/"+fname+".tex" ;
 system,"acroread -geometry 1200x950 "+dir+"/"+fname+".pdf &" ;
 cd,dd;
}


func eigs3D(a,&vec)
/* DOCUMENT
 returns 
 eigenvalues  and eigenvectors of 
 symmetric A=[[a11,...]

 EXAMPLE 
 tt=random(25,25,3,3); tt+=transpose(tt,[3,4]);ee=eigs3D(tt,vec);
 qq=(((tt*vec(,,,1))(,,sum,))/vec(,,,1))/ee(,,1,-); // first eigen
 stat,qq
SEE ALSO: 
 */
{
res= _eigs3D(a(..,1,1),a(..,1,2),a(..,1,3),a(..,2,2),a(..,2,3),a(..,3,3));
 
 vec= res(..,4:);
 dv= dimsof(vec);
 dv=dv(1:-1);
 dv(1)+=1;  grow,dv,[3,3];
 vec=reform(vec,dv);
 return res(..,:3);
}




 func _eigs3D(a11,a12,a13,a22,a23,a33)
/* DOCUMENT
 returns 
eigenvalues  and eigenvectors of 
symmetric A=[[a11,...]
 	SEE ALSO: 
 */
{
  a11+=0.;
  a12+=0.;
  a13+=0.;
  a22+=0.;
  a23+=0.;
  a33+=0.;
 out =(a33*0.)(..,-:1:12);
z1 = a11 + a22 + a33;
z2 = 0.333333333333333*z1;
z3 = pow(a11,2.);
z4 = -z3;
z5 = pow(a12,2.);
z6 = -3.*z5;
z7 = pow(a13,2.);
z8 = -3.*z7;
z9 = a11*a22;
z10 = pow(a22,2.);
z11 = -z10;
z12 = pow(a23,2.);
z13 = -3.*z12;
z14 = a11*a33;
z15 = a22*a33;
z16 = pow(a33,2.);
z17 = -z16;
z18 = z11 + z13 + z14 + z15 + z17 + z4 + z6 + z8 + z9;
z19 = pow(a11,3.);
z20 = 2.*z19;
z21 = 9.*a11*z5;
z22 = 9.*a11*z7;
z23 = -3.*a22*z3;
z24 = 9.*a22*z5;
z25 = -18.*a22*z7;
z26 = -3.*a11*z10;
z27 = pow(a22,3.);
z28 = 2.*z27;
z29 = 54.*a12*a13*a23;
z30 = -18.*a11*z12;
z31 = 9.*a22*z12;
z32 = -3.*a33*z3;
z33 = -18.*a33*z5;
z34 = 9.*a33*z7;
z35 = 12.*a11*a22*a33;
z36 = -3.*a33*z10;
z37 = 9.*a33*z12;
z38 = -3.*a11*z16;
z39 = -3.*a22*z16;
z40 = pow(a33,3.);
z41 = 2.*z40;
z42 = pow(z18,3.);
z43 = 4.*z42;
z44 = z20 + z21 + z22 + z23 + z24 + z25 + z26 + z28 + z29 + z30 + z31 + z32 + z33 + z34 + z35 + z36 + z37 + z38 + z39 + z41;
z45 = pow(z44,2.);
z46 = z43 + z45;
z47 = sqrt(z46+0.i);
z48 = z20 + z21 + z22 + z23 + z24 + z25 + z26 + z28 + z29 + z30 + z31 + z32 + z33 + z34 + z35 + z36 + z37 + z38 + z39 + z41 + z47;
 z49 = pow(z48+(z48==0),-1./3.);//====
z50 = pow(z48,1./3.);
z51 = 1/a13;
z52 = -a11;
z53 = -a22;
z54 = -a33;
z55 = z52 + z53 + z54;
z56 = 0.333333333333333*z55;
z57 = 0.419973683298291*z18*z49;
z58 = -0.2645668419947*z50;
z59 = a33 + z56 + z57 + z58;
z60 = -(a12*a23);
z61 = a22 + z56 + z57 + z58;
z62 = a13*z61;
z63 = z60 + z62;
z64 = 1/z63;
z65 = a13*a23;
z66 = -(a12*z59);
z67 = z65 + z66;
z68 = (-0.20998684164914555 -0.3637078786572405i)*z18*z49;
z69 = (0.13228342099734997-0.22912160616643376i)*z50;
z70 = a33 + z56 + z68 + z69;
z71 = a22 + z56 + z68 + z69;
z72 = a13*z71;
z73 = z60 + z72;
z74 = 1/z73;
z75 = -(a12*z70);
z76 = z65 + z75;
z77 = (-0.20998684164914555+0.3637078786572405i)*z18*z49;
z78 = (0.13228342099734997+0.22912160616643376i)*z50;
z79 = a33 + z56 + z77 + z78;
z80 = a22 + z56 + z77 + z78;
z81 = a13*z80;
z82 = z60 + z81;
z83 = 1/z82;
z84 = -(a12*z79);
z85 = z65 + z84;
out(..,1)= z2 - 0.419973683298291*z18*z49 + 0.2645668419947*z50;
out(..,2)= z2 + (0.20998684164914555+0.3637078786572405i)*z18*z49 - (0.13228342099734997-0.22912160616643376i)*z50;
out(..,3)= z2 + (0.20998684164914555-0.3637078786572405i)*z18*z49 - (0.13228342099734997+0.22912160616643376i)*z50;
out(..,4)= -(z51*z59) + a23*z51*z64*z67;
out(..,5)= -(z64*z67);
out(..,6)= 1.;
out(..,7)= -(z51*z70) + a23*z51*z74*z76;
out(..,8)= -(z74*z76);
out(..,9)= 1.;
out(..,10)= -(z51*z79) + a23*z51*z83*z85;
out(..,11)= -(z83*z85);
out(..,12)= 1.;
return out; }


func eigs3DG(a,&vec)
/* DOCUMENT
 returns 
 eigenvalues  and eigenvectors of  A=[[a11,...]
  need not be symmetric;
 EXAMPLE 
 tt=[[11,12,13],[21,22,23],[31,32,33.-1.i]];
  ee=eigs3D(tt,vec);
 qq=(((tt*vec(,,,1))(,,sum,))/vec(,,,1))/ee(,,1,-); // first eigen
 stat,qq
SEE ALSO: 
 */
{
  a+=0.i;
  res= _eigs3D2(a(..,1,1),a(..,1,2),a(..,1,3),
                a(..,2,1),a(..,2,2),a(..,2,3),
               a(..,3,1),a(..,3,2),a(..,3,3));
 
 vec= res(..,4:);
 dv= dimsof(vec);
 dv=dv(1:-1);
 dv(1)+=1;  grow,dv,[3,3];
 vec=reform(vec,dv);
 return res(..,:3);
}


 func _eigs3DG(a11,a12,a13,a21,a22,a23,a31,a32,a33)
/* DOCUMENT

   returns eigenvalues and vector of arbitrary matrix

 	SEE ALSO: 
 */
{ 
  a11+=0.;
  a12+=0.;
  a13+=0.;
  a21+=0.;
  a22+=0.;
  a23+=0.;
  a31+=0.;
  a32+=0.;
  a33+=0.;
 out =(a33*0.+0.i)(..,-:1:12);
z1 = a11 + a22 + a33;
z2 = 0.333333333333333*z1;
z3 = pow(a11,2.);
z4 = -z3;
z5 = -3.*a12*a21;
z6 = a11*a22;
z7 = pow(a22,2.);
z8 = -z7;
z9 = -3.*a13*a31;
z10 = -3.*a23*a32;
z11 = a11*a33;
z12 = a22*a33;
z13 = pow(a33,2.);
z14 = -z13;
z15 = z10 + z11 + z12 + z14 + z4 + z5 + z6 + z8 + z9;
z16 = pow(a11,3.);
z17 = 2.*z16;
z18 = 9.*a11*a12*a21;
z19 = -3.*a22*z3;
z20 = 9.*a12*a21*a22;
z21 = -3.*a11*z7;
z22 = pow(a22,3.);
z23 = 2.*z22;
z24 = 9.*a11*a13*a31;
z25 = -18.*a13*a22*a31;
z26 = 27.*a12*a23*a31;
z27 = 27.*a13*a21*a32;
z28 = -18.*a11*a23*a32;
z29 = 9.*a22*a23*a32;
z30 = -3.*a33*z3;
z31 = -18.*a12*a21*a33;
z32 = 12.*a11*a22*a33;
z33 = -3.*a33*z7;
z34 = 9.*a13*a31*a33;
z35 = 9.*a23*a32*a33;
z36 = -3.*a11*z13;
z37 = -3.*a22*z13;
z38 = pow(a33,3.);
z39 = 2.*z38;
z40 = pow(z15,3.);
z41 = 4.*z40;
z42 = z17 + z18 + z19 + z20 + z21 + z23 + z24 + z25 + z26 + z27 + z28 + z29 + z30 + z31 + z32 + z33 + z34 + z35 + z36 + z37 + z39;
z43 = pow(z42,2.);
z44 = z41 + z43;
z45 = sqrt(0.i+z44);
z46 = z17 + z18 + z19 + z20 + z21 + z23 + z24 + z25 + z26 + z27 + z28 + z29 + z30 + z31 + z32 + z33 + z34 + z35 + z36 + z37 + z39 + z45;
z47 = pow(z46,-1./3.);
z48 = pow(z46,1./3.);
z49 = 1/a31;
z50 = -a11;
z51 = -a22;
z52 = -a33;
z53 = z50 + z51 + z52;
z54 = 0.333333333333333*z53;
z55 = 0.419973683298291*z15*z47;
z56 = -0.2645668419947*z48;
z57 = a33 + z54 + z55 + z56;
z58 = -(a21*a32);
z59 = a22 + z54 + z55 + z56;
z60 = a31*z59;
z61 = z58 + z60;
z62 = 1/z61;
z63 = a23*a31;
z64 = -(a21*z57);
z65 = z63 + z64;
z66 = (-0.20998684164914555-0.3637078786572405i)*z15*z47;
z67 = (0.13228342099734997-0.22912160616643376i)*z48;
z68 = a33 + z54 + z66 + z67;
z69 = a22 + z54 + z66 + z67;
z70 = a31*z69;
z71 = z58 + z70;
z72 = 1/z71;
z73 = -(a21*z68);
z74 = z63 + z73;
z75 = (-0.20998684164914555+0.3637078786572405i)*z15*z47;
z76 = (0.13228342099734997+0.22912160616643376i)*z48;
z77 = a33 + z54 + z75 + z76;
z78 = a22 + z54 + z75 + z76;
z79 = a31*z78;
z80 = z58 + z79;
z81 = 1/z80;
z82 = -(a21*z77);
z83 = z63 + z82;
out(..,1)= z2 - 0.419973683298291*z15*z47 + 0.2645668419947*z48;
out(..,2)= z2 + (0.20998684164914555+0.3637078786572405i)*z15*z47 - (0.13228342099734997-0.22912160616643376i)*z48;
out(..,3)= z2 + (0.20998684164914555-0.3637078786572405i)*z15*z47 - (0.13228342099734997+0.22912160616643376i)*z48;
out(..,4)= -(z49*z57) + a32*z49*z62*z65;
out(..,5)= -(z62*z65);
out(..,6)= 1.;
out(..,7)= -(z49*z68) + a32*z49*z72*z74;
out(..,8)= -(z72*z74);
out(..,9)= 1.;
out(..,10)= -(z49*z77) + a32*z49*z81*z83;
out(..,11)= -(z81*z83);
out(..,12)= 1.;
return out; }


func pad_periodic(x)
/* DOCUMENT 
     pad a cube with periodic value;
   SEE ALSO:
 */
{
  d=dimsof(x)(2:);
  y=array(0.,[3,d(1)+1,d(2)+1,d(3)+1]);
  y(:d(1),:d(2),:d(3))=x;
  y(d(1)+1,:d(2),:d(3))=x(1,,);
  y(:d(1),d(2)+1,:d(3))=x(,1,);
  y(:d(1),:d(2),d(3)+1)=x(,,1);
  y(d(1)+1,d(2)+1,:d(3))=x(1,1,);
  y(d(1)+1,:d(2),d(3)+1)=x(1,,1);
  y(:d(1),d(2)+1,d(3)+1)=x(,1,1);
  y(d(1)+1,d(2)+1,d(3)+1)=x(1,1,1);
  return y;
}


__rgb=[[0,0,0],[48,48,48],[88,88,88],[128,128,128],[160,160,160],[195,195,195],[220,
220,220],[255,255,255],[64,0,0],[128,0,0],[192,0,0],[255,0,0],[255,192,192],[0,
64,0],[0,128,0],[0,192,0],[0,255,0],[192,255,192],[0,0,128],[0,0,192],[0,0,
255],[192,192,255],[64,64,0],[128,128,0],[192,192,0],[255,255,0],[255,255,192],
[0,64,64],[0,128,128],[0,192,192],[0,255,255],[192,255,255],[64,0,64],[128,0,
128],[192,0,192],[255,0,255],[255,192,255],[192,88,0],[255,128,0],[255,168,88],
      [253,220,168]];
