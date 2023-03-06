#include "Bastien/string_utils.i"
#include "Bastien/list_utils.i"
#include "Bastien/plot_utils.i"
#include "Bastien/system_utils.i"
#include "Bastien/ascii_utils.i"
#include "Bastien/fits_utils.i"
#include "Bastien/pdb_utils.i"
#include "Bastien/sky_utils.i"
#include "Bastien/graph_utils.i"
#include "Bastien/struct_utils.i"


func INCLUDE(filename)
{
  include,filename,1;
  listFunc(filename);
}

func pos2Coord(pos,dims)
{
  if(dims(1)!=2) error,"Work only on 2D array...!";
  pos--;
  i=(pos%dims(2))+1;
  j=(pos/dims(2))+1;
  return [i,j];
}


func coord2Pos(coord,dims)
{
  if(dims(1)!=2) error,"Work only on 2D array...!";
  return coord(..,1)+dims(2)*(coord(..,2)-1);
}


func mkdirhier(dirname)
/* DOCUMENT mkdirhier(dirname)

  The mkdirhier command creates the specified directories. Unlike mkdir if
  any of the parent directories of the specified directory do not exist, it
  creates them as well.

   SEE ALSO: mkdir
 */
{
  local nbdir,flag;

  nbdir=numberof(dirname);
  for(i=1;i<=nbdir;i++)
    {
      if(strpart(dirname(i),0:0)=="/")
        curdir=dirname(i);
      else curdir=dirname(i)+"/";

      curdir="/"+curdir;
      tmp=*pointer(curdir);

      idx=where(tmp=='/');
      j=numberof(idx);
      while((--j>1)&&(typeof(lsdir(strpart(curdir,2:idx(j)-1)))=="long"));

      while(j<numberof(idx)) mkdir,strpart(curdir,2:idx(++j)-1);

    }
  return dirname;
}



func eq_copy(&x,y)
  /* DOCUMENT  eq_copy,&x,y
     Do x=y and return nothing.
     The only use is to avoid automatic local definition
 */
{
  x=y;
  return;
}

func last_ref(&x)
/* DOCUMENT last_ref(x)
     returns X, destroying X in the process.
   SEE ALSO: eq_nocopy
 */
{
  local y;
  eq_nocopy, y, x;
  x= [];
  return y;
}


func drawline(a,b,c,color=,type=,width=)
/* DOCUMENT drawline(a,b)
         or drawline(a,b,c)

	    draw line y=ax+b
         or      line ax+by+c=0
   KEYWORDS: color,type,width
*/
{
  if(is_void(c))
    {
      c=b;
      b=a*0.-1.;
    }

  if(is_array(dimsof(a,b,c)))
    {
      a=a+0*c+0*b;
      b=b+0*a+0*c;
      c=c+0*a+0*b;
    }
  else
    error,"Dimension problem !!";

  indx=where((b!=0)+(a!=0));
  if(numberof(indx))
    {
      aa=a(indx);
      bb=b(indx);
      cc=c(indx);
    }
  else
    {
      error,"can't plot equation c=0 !";
    }
  indxb=(bb==0);
  indxa=(aa==0);
  indxab=(1-indxa)*(1-indxb);
  aa+=indxa;
  bb+=indxb;

  rl=limits();
  x0=(indxab+indxa)*rl(1)+indxb*(-cc/aa);
  x1=(indxab+indxa)*rl(2)+indxb*(-cc/aa);
  y0=indxab*(-(aa*x0+cc)/bb)+indxb*rl(3)+indxa*(-cc/bb);
  y1=indxab*(-(aa*x1+cc)/bb)+indxb*rl(4)+indxa*(-cc/bb);

  pldj,x0,y0,x1,y1,color=color,type=type,width=width;
  return;
}

func markvert(maxline=,color=,type=)
/* DOCUMENT markvert(maxline=,color=,type=)
   Use Left button of the mouse to mark verticale line (other button
   will quit).

   KEYWORDS maxline Set the maximal number of vertical
                    lines to plot. By default, no
                    limit.
*/
{
  again=1;
  if(is_void(maxline)||maxline<1)
    {
      maxline=-1;
    }
  nb=0;
  write,"Mark a point of each line.";
  a=[];
  while(again&&nb!=maxline)
    {
      nb++;
      r=mouse(1,1,"");
      again=r(10)==1;
      if (again)
	{
	  write,format=" %3i    %g\n",nb,r(1);
	  drawline,1,0,-r(1),color=color,type=type;
	  grow,a,r(1);
	}
    }
  return a;
}


func plgmk(y,x,color=,marker=,msize=)
/* DOCUMENT plgmk(y,x,color=,marker=,msize=)
   do plg and plmk on y and x
*/
{
  if(is_void(marker)) marker=2;
  if(is_void(msize)) msize=0.5;
  plg,y,x,color=color;
  plmk,y,x,marker=marker,msize=msize,color=color;
  return;
}

func plg2(y,x,m,color=,msize=)
/* DOCUMENT plg2(y,x,m,color=,msize=)
     do plg,y,x,type=0,marker=m,color=color,msize=msize
 */
{
  if(is_void(x)) x=span(1,numberof(y),numberof(y));
  if(is_void(m)) m=1;

  plg,y,x,type=0,marker=m,color=color,msize=msize;
  return;
}

func mplg(y,..,color=,which=,over=)
{
  local x,yy;
  local r,g,b,nbc;
  local offset;

  palette,r,g,b,query=1;
  nbc=numberof(r);

  if(is_void(which)) which=2;
  if(is_void(over)) over=0;
  if(over) over=1;

  yy=transpose(y,[1,which])
  ncol=dimsof(yy)(2);
  npt=numberof(yy)/ncol;
  
  if(more_args()) x=next_arg();
  else x=span(1,npt,npt);
  
  if(is_void(color)) c=int(span(1,nbc,ncol));
  else c=color;
  
  if(numberof(c)!=ncol||numberof(x)!=npt) exit,"dimension problem!";

  offset=0;
  for(i=1;i<=ncol;i++)
    {
      offset-=yy(i,*)(min);
      plg,yy(i,*)+offset,x,color=c(i);
      offset+=over*(yy(i,*)(max));
    }
  return;
}



func rancar(n)
{
  indx=indgen(1:(n*n));
  yy=[];
  while(ni=numberof(indx))
    {
      p=int((random()*ni)+1);
      grow,yy,indx(p);
      indx=indx(where(indgen(1:ni)!=p));
    }
  return yy;
}

func deletebad(good_indx,bad_indx,growp=)
{
  if (is_void(growp)) growp=0;
  if(numberof(bindx))
    {
      indx=(bad_indx(,-)+indgen(-growp:growp)(-,))(*);
      indx=indx(where(indx>0));
      difindx=where2((good_indx(,-)-indx(-,))==0);
      if(numberof(difindx)) difindx=difindx(1,);
      grow,difindx,numberof(good_indx)+1;
      indx=where(!histogram(difindx));
      if(numberof(indx)) return good_indx(indx);
      else return [];
    }
  else
    return good_indx;
}

func same(a,b)
{
  na=numberof(a);
  nb=numberof(b);
  if(na==nb) return allof(a==b);
  return 0;
}

func complete(&a,n,fill)
/* DOCUMENT  complete,&a,n,fill
         or  complete(a,n,fill)
     complete a vector a with value fill (default 0)
     until its size is equal to n.

   SEE ALSO: completep2
*/
{
  if(is_void(fill)) fill=0.;
  na=numberof(a);
  if(na!=n)
    {
      if(am_subroutine()) grow,a,array(a(1)*0.+fill,n-na);
      else return grow(a,array(a(1)*0.+fill,n-na));
    }
  return a;
}

func completep2(&a,fill)
/* DOCUMENT completep2,&a,fill
         or completep2(a,fill)
    complete vector a with value fill (default 0) until its size
    is equal to a power of 2.
 
    SEE ALSO: complete
 */
{
  na=numberof(a);
  if(na==0) na=1;
  pa=2^int(log(na)/log(2));
  if(pa<na) pa<<=1;
  if(am_subroutine()) complete,a,pa,fill;
  else return complete(a,pa,fill);
}



func fact(n)
{
  local m;
  local i,mask;

  if(dimsof(n)(1)==0) m=[n];
  else m=n;
  
  i=array(1.,dimsof(m));
  while(numberof((mask=where(m>0))))
    {
      i(mask)*=m(mask);
      m(mask)--;
    }
  return i;
}

func displayFont(typ=,font=,hexa=)
{
  if(is_void(typ)) typ=0;
  if(is_void(hexa)) hexa=0;
  fma;
  limits,0.6,11,0,25;
  for(i=32;i<256;i++)
    {
      pts1=char([i,0x00]);
      pts2=char([0x21,i,0x00]);
      
      if(hexa) st1=swrite(format="%4X",i);
      else st1=swrite(format="%03o",i);
      if(typ) st2=string(&pts2);
      else st2=string(&pts1);
      px=(i-32)/24+1;
      py=24-(i-32)%24;
      plt,st1,px,py,tosys=1;
      plt,st2,px+0.5,py,tosys=1,font=font;
    }
}

func sort2(f)
{
  slt=sort(f(*));
  n1=dimsof(f)(2);
  n2=dimsof(f)(3);
  rslt=array(long,2,numberof(slt));
  rslt(,1)=(slt-1)%n1+1;
  rslt(,2)=(slt-1)/n1+1;
  return rslt;
}

  
func ellp(x,y,a,b,rho)
{
  local a2,b2,f;
  a2=a*a;
  b2=b*b;
  f=rho/a/b;
  
  return x(,-)^2/a2+y(-,)^2/b2-f*x(,-)*y(-,);
}

func sigprob(ksig)
/* DOCUMENT  sigprob(ksig)

   ksig : number of sigma of the detection below the mean.
   return the probability that the detection is real for a gaussian
   distribution.

   SEE ALSO: probsig
 */
{
  local x,y;

  x=span(-ksig,ksig,1001);
  y=1/sqrt(pi*2)*exp(-x^2/2);
  return y(sum)*x(1:2)(dif);
}

func sigqua(ksig)
{
  local x,y;

  nb=10001;
  x=span(-20,20,nb);
  dx=x(1:2)(dif)(1);
  y=1/sqrt(pi*2)*exp(-x^2/2);
  y=integ(y,x,x);
  return y(abs(x(-,)-ksig(,-))(,mnx));
}

func quasig(q)
{
  local x,y;
  nb=10001;
  x=span(-20,20,nb);
  dx=x(1:2)(dif)(1);
  y=1/sqrt(pi*2)*exp(-x^2/2);
  y=integ(y,x,x);
  return x(abs(y(-,)-q(,-))(,mnx));
}

func sigprobc(ksig)
{
  local x,y;

  nb=1001;
  x=span(-ksig,ksig,nb);
  y=1./x(1:2)(dif)/nb-1/sqrt(pi*2)*exp(-x^2/2);
  return (y)(sum)*x(1:2)(dif);
}

func probsig(prob)
{
  nb=10001;
  x=span(-20,20,nb);
  dx=x(1:2)(dif)(1);
  y=1/sqrt(pi*2)*exp(-x^2/2);
  sy=y(sort(y));
  csy=sy(0::-1)(cum)(0:2:-1)*dx;
  plg,csy,sy;
  
  id=abs(csy-prob)(mnx);
  id=where(y>=sy(id));
  if(numberof(id)==0) return 0;
  else return (x(id(0))-x(id(1)))*0.5;
}


func probsigall(ldist,lx,k)
{
  mn=(ldist*lx)(sum)/ldist(sum);
  sg=sqrt((ldist*lx*lx)(sum)/ldist(sum)-mn*mn);
  idx=where(abs(lx-mn)<(k*sg));
  write,format="mean : %f\nsigm : %f\n",mn,sg;
  return 100*ldist(idx)(sum)/ldist(sum);
}


func compsnc(nflux,nsigm,pixc,&nbpx,&sng,width=,draw=)
{
  local signal,noise,rslt;
  local idxs,idxi,nidx,nidxs,nidxi;
  local nsig2;

  draw=is_void(draw)?1:draw;
  width=is_void(width)?100:int(width);
  rslt=[];sng=[];
  if(numberof(pixc)==0) pixc=[pixc];

  for(i=1;i<=numberof(pixc);i++)
    {
      infidx=max(pixc(i)-width,1);
      supidx=min(pixc(i)+width,numberof(nflux));
      nidxi=numberof((idxi=indgen(pixc(i):infidx:-1)));
      nidxs=numberof((idxs=indgen(pixc(i):supidx)));
      
      if(nidxi>nidxs)
        {
          signal=1-nflux(idxi);
          signal(1:nidxs)+=1-nflux(idxs);
          
          noise=nsigm(idxi)^2;
          noise(1:nidxs)+=nsigm(idxs)^2;
        }
      else
        {
          signal=1-nflux(idxs);
          signal(1:nidxi)+=1-nflux(idxi);
          
          noise=nsigm(idxs)^2;
          noise(1:nidxi)+=nsigm(idxi)^2;
        }
      signal(1)-=1-nflux(pixc(i));
      noise(1)-=nsigm(pixc(i))^2;
      signal=signal(cum)(2:);
      noise=noise(cum)(2:);
      noise=sqrt(noise);
      if(draw)
        {
          plg,signal,color=-7;
          plg,noise,color=-5;
          plg,signal/noise,color=-5,type=2;
        }
      grow,sng,signal((signal/noise)(mxx));
      grow,rslt,(signal/noise)(max);
      grow,nbpx,(signal/noise)(mxx)
    }
  return rslt;
}

  
      
func plsn(lines,regions,nsigm,draw=)
{
  local mx,yc;
  local invf,sg2;

  if(is_void(draw)) draw=1;
  rslt=[];
  for(i=1;i<=numberof(lines);i++)
    {
      yc=CompConvFlux(lines(i),regions,755);
      mx=yc(mnx);
      grow,rslt,compsnc(yc,nsigm,mx,draw=draw);
    }
  return rslt;
}

func plsn2(lines,twave,nsigm,&nbpx,draw=)
{
  local mx,yc;
  local invf,sg2;

  nbpx=[];
  rg=[];
  AddRegion,rg,twave,twave(1),twave(0),twave,twave,1;
  if(is_void(draw)) draw=1;
  rslt=[];
  for(i=1;i<=numberof(lines);i++)
    {
      yc=CompConvFlux(lines(i),rg,755);
      mx=yc(mnx);
      grow,rslt,compsnc(yc,nsigm,mx,nbpx,draw=draw);
    }
  return rslt;
}

func echant(wave,flux,cont,sigm,nbscale)
{
  twave=array(0.,numberof(wave)*nbscale);
  hwave=wave(1)+(wave(1)-wave(2))/2.;
  grow,hwave,wave,2*wave(0)-wave(-1);
  hwave=hwave(zcen);
  dwave=hwave(dif)/nbscale;
  for(i=1;i<=nbscale;i++)
    twave(i::nbscale)=hwave(1:-1)+(i-1)*dwave;
  grow,twave,hwave(0);
  twave=twave(zcen);
  tflux=twave;
  tcont=twave;
  tsigm=twave;
  for(i=1;i<=nbscale;i++)
    {
      tflux(i:0:nbscale)=flux;
      tcont(i:0:nbscale)=cont;
      tsigm(i:0:nbscale)=sigm*sqrt(nbscale);
    }
  return [twave,tflux,tcont,tsigm];
}


func cldnsum(cldn)
{
  local mcldn;
  mcldn=min(cldn(*));
  return mcldn+log10((10^(cldn-mcldn))(*)(sum));
}

    
func blockmat(a11,a12,a21,a22)
{
  local n11,n12,n21,n22,n;
  local rslt;

  n11=dimsof(a11)(2:);
  n12=dimsof(a12)(2:);
  n21=dimsof(a21)(2:);
  n22=dimsof(a22)(2:);
  if(n12(1)!=n11(1)||n12(2)!=n22(2))
    error,"Bad dimension for 2nd argument !";
  if(n21(1)!=n22(1)||n21(2)!=n11(2))
    error,"Bad dimension for 3rd argument !";

  n=2;
  grow,n,n11+n22;
  
  rslt=array(a11(1)+a12(1)+a21(1)+a22(1),n);
  rslt(1:n11(1),1:n11(2))=a11;
  rslt(1:n11(1),-n12(2)+1:0)=a12;
  rslt(-n21(1)+1:0,1:n21(2))=a21;
  rslt(-n22(1)+1:0,-n22(2)+1:0)=a22;
  return rslt;
}



func myselect(y,x,&path,free=,disp=)
/* DOCUMENT myselect(y,x,&path,free=,disp=)
     given a set of points [y,x] return the index of the point in path.
     If path is void, it is asked with the mouse.
     KEYWORDS: free : if set use a free polygone
                      else use rectangle that contain all the entered points.
   SEE ALSO:
 */
{
  local ptx,pty;
  local fin,r;

  if(is_void(disp)) disp=1;
  if(is_void(free)) free=1;
  if(is_void(   x)) x=span(1,numberof(y),numberof(y));
  if(is_void(path))
    {
      fin=0;
      ptx=pty=[];
      while(!fin)
        {
          r=mouse(,0);
          fin=r(-1)==3;
          if(fin) continue;
          grow,ptx,r(1);
          grow,pty,r(2);
          plg,pty(0),ptx(0),type=0,color=-5,marker='\1';
          if(numberof(ptx)>1)
            plg,pty(-1:0),ptx(-1:0),color=-5;
        }
      nptx=numberof(ptx);
      if(free&&nptx>2) plg,pty([1,nptx]),ptx([1,nptx]),color=-5;
      path=[ptx,pty];
    }
  else
    {
      ptx=path(,1);
      pty=path(,2);
      nptx=numberof(ptx);
    }
  if(!free&&nptx<2)
    error,"Need almost two points to define a rectangle !!";
  if(free&&nptx<3)
    error,"Need almost three points to define a free area !!";
  r=2*pi;
  if(free)
    {
      flux=array(0.,numberof(x));
      for(i=0;i<numberof(ptx);i++)
        {
          diffa=atan(pty(i+1)-y,ptx(i+1)-x)-atan(pty(i)-y,ptx(i)-x);
          diffa=(diffa+r)%r; // here r=2*PI
          diffb=r-diffa;
          flux+=merge2(diffa,-diffb,diffa<diffb);
          //      return where(abs(flux)>=pi);
        }
      flux=where(abs(flux)>1);
    }
  else
    {
      xmin=min(ptx);
      xmax=max(ptx);
      ymin=min(pty);
      ymax=max(pty);
      flux=where((x>=xmin)*(x<=xmax)*(y>=ymin)*(y<=ymax));
    }
  if(disp&is_array(flux))
    plg,y(flux),x(flux),color=-7,type=0,marker='\2';
  return flux;
}


func clusterposi(h,s=)
/* DOCUMENT Obsolete use find_cluster !
     
 */
{
  error,"Obsolete: use find_cluster(h,size=s). Beware result is now r(,1) r(,2)!"
}
func find_cluster(logical_array,size=)
/* DOCUMENT  find_cluster(logical_array)
     logical_array must be a 1D array of logical value (0 or 1,
     mainly logical operation on two array of 1D)
     return 2D array r so that r(,1) and r(,2) give start and end
     respectively of following true regions 
 */
{
  local dh,infi,supi;
  local nlga,rlga;

  nlogical_array=numberof(logical_array);
  
  if(noneof(logical_array)) return [];
  if(allof( logical_array)) return [[1],[nlogical_array]];
  
  dh=logical_array(dif);
  supi=where(dh==-1);
  infi=where(dh==1);
  if(!numberof(infi)) infi=[1];
  else infi++;
  if(!numberof(supi)) supi=[nlogical_array];

  if(infi(1)>supi(1)) {grow,infi,[1n];infi=roll(infi,1n);}
  if(supi(0)<infi(0)) grow,supi,nlogical_array;

  if(!is_void(size))
    {
      idx=where((supi-infi+1)>=size);
      if(!is_array(idx)) return [];
      infi=infi(idx);
      supi=supi(idx);
    }
  return long([infi,supi]);
}

func make_idx(nbpts,range_size,over_size)
/* DOCUMENT rslt=make_idx(nbpts,range_size,over_size)

   return an 2D array rslt of  index so that the 1D array rslt(,i) are
   ranges  of   index  of  size   range_size  that  cover   the  range
   [1...nbpts]. Moreover  rslt(,i) overlaps rslt(,i+1)  over over_size
   pixels.

   EXAMPLE:  rslt=make_idx(52,10,3)
             info,rslt  // return array(long,10,70)
             rslt(,1)   // return [ 1, 2, 3, 4, 5, 6, 7, 8, 9,10]
             rslt(,2)   // return [ 8, 9,10,11,12,13,14,15,16,17]
             ...
             rslt(,7)   // return [43,44,45,46,47,48,49,50,51,52]
*/
{
  local nrg,n;

  range_size=is_void(range_size)?(nbpts/50):long(range_size);
  over_size=is_void(range_size)?0:min(max(0,long(over_size)),range_size-1);

  over_size=range_size-over_size;
  
  n=long((nbpts-range_size)/over_size)+1;
  nrg=(indgen(1:n)-1)*over_size;
  if((nrg(0)+range_size)<nbpts) grow,nrg,nbpts-range_size;
    
  return (indgen(1:range_size))(,-)+nrg(-,);
}

func merge_result(input,idx,brd=)
{
  local poids,fcont;
  if(is_void(brd)) brd=0;
  
  fcont=array(0.,max(idx(*)));
  poids=fcont;
  
  fcont(idx(1:-brd,1))+=input(1:-brd,1);
  poids(idx(1:-brd,1))++;  

  for(i=2;i<dimsof(input)(0);i++)
    {
      fcont(idx(brd+1:-brd,i))+=input(brd+1:-brd,i);
      poids(idx(brd+1:-brd,i))++;
    }

  i=dimsof(input)(0)
  fcont(idx(brd+1:0,i))+=input(brd+1:0,i);
  poids(idx(brd+1:0,i))++;  


  return fcont/poids;
}

func mederr(x,pourc=)
{
  local hcut,lcut,med,err;
  if(is_void(pourc)) pourc=90.;
  pourc=is_void(pourc)?0.9:pourc/100.;

  idx=sort(x);
  hnidx=(nidx=numberof(idx))>>1;
  //compute mediane
  if(nidx&1) med=x(idx(hnidx+1));
  else       med=avg(x(idx(hnidx:hnidx+1)));

  pourc=0.5*(1-pourc);
  lcut=long(max(1,1+nidx*pourc));
  hcut=long(min(nidx,1+nidx-lcut));
  return [med,x(idx(lcut)),x(idx(hcut))];
}
  
  
func glue(&head,tail)
{
  local dh,dt,rslt;
  local ams;

  ams=am_subroutine();
  
  dh=dimsof(head);
  dt=dimsof(tail);
  if(!numberof(dh))
    {
      if(ams) head=tail;return;
      return tail;
    }
  if(dh(1)==1)
    {
      if(dt(1)!=1) error,"later arguments not conformable with 1st in glue";
      if(ams) {grow,head,tail;return;}
      return grow(head,tail);
    }


  if((dt(1)!=dh(1))||(abs(dh(3:)-dt(3:))(max)))
    error,"later arguments not conformable with 1st in glue";
  dms=dh;
  dms(2)+=dt(2);
  rslt=array(head(1),dms);
  rslt(      1:dh(2),*)=head(,*);
  rslt(dh(2)+1:     ,*)=tail(,*);
  if(ams) {head=rslt;rslt=[];return;}
  return rslt;
}

func getMIndex(data,&idmin,&idmax)
{
  r1=mouse();
  drawline,-1,0,r1(1),color=-5;
  r2=mouse();
  drawline,-1,0,r2(1),color=-5;
  
  return getIndex(data,min(r1(1),r2(1)),max(r1(1),r2(1)),idmin,idmax);
}


func getIndex(range,xmin,xmax,&idxmin,&idxmax)
{
  idxmin=abs(range-xmin)(mnx);
  idxmax=abs(range-xmax)(mnx);
  return idxmin:idxmax;
}


func sbit2int(str,lowb=,highb=)
{
  local rslt;

  if(is_void(lowb)) lowb='.';
  if(is_void(highb)) highb='#';

  rslt=array(long,dimsof(str));
  for(i=1;i<=numberof(str);i++)
    {
      ptr=*pointer(str(i));
      ptr(1:-1)=(ptr(1:-1)-lowb)/(highb-lowb)+'0';
      rslt(i)=str2int(string(&ptr));
    }
  return rslt;
}

func int2sbit(nb,lowb=,highb=,nbdg=)
{
  local rslt;

  if(is_void(lowb)) lowb='.';
  if(is_void(highb)) highb='#';
  if(is_void( nbdg)) nbdg=max(1+int(log10(nb)));
  
  rslt=array(string,dimsof(nb));
  frmt=swrite(format="%%0%dd",nbdg);
  tmp=swrite(format=frmt,nb);
  for(i=1;i<=numberof(tmp);i++)
    {
      ptr=*pointer(tmp(i));
      ptr(1:-1)=(ptr(1:-1)-'0')*(highb-lowb)+lowb;
      rslt(i)=string(&ptr);
    }
  return rslt;
}

func findOverlap(wmn,wmx,lstwmn,lstwmx,&liminf,&limsup)
{
  idx=where((wmn<=lstwmx)&(wmx>=lstwmn));

  liminf=limsup=[];
  if(!numberof(idx)) return [];

  nidx=numberof(idx);
  area=array(0.,nidx);
  liminf=max(array(wmn,nidx),lstwmn(idx));
  limsup=min(array(wmx,nidx),lstwmx(idx));
  return idx;
  
}



func wgetIdxInIncreasingTable(x,val,range=)
{
  if(is_func(getIdxInIncreasingTable))
    rslt=getIdxInIncreasingTable(x,val,range=range);
  else
    {
      rslt=array(long,dimsof(val));
      for(i=1;i<=numberof(val);i++)
        rslt(i)=abs(x-val(i))(mnx);
      if(range) return rslt(1):rslt(2);
    }
  return rslt;
}

