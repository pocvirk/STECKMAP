func mvar(perc)
{
  return mvax(1,perc);
}

func mval(perc)
{
  return mvax(-1,perc);
}

func mvax(dir,perc)
{
  local a,b,c,d;
  get_style,a,b,c,d;
  nbsys=numberof(b);
  sys=plsys();
  for(i=1;i<=nbsys;i++)
    {
      plsys,i;
      mvx,dir,perc;
    }
  plsys,sys;
}

func mvr(perc)
{
  return mvx(1,perc);
}

func mvl(perc)
{
  return mvx(-1,perc);
}

func mvx(dir,perc)
{
  if(is_void(perc)) perc=0.5;
  ll=limits();
  islogx=(int(ll(5))&128);
  xmin=ll(1);
  xmax=ll(2);
  if(islogx)
    {
      xmin=log10(xmin);
      xmax=log10(xmax);
    }
  xmove=(xmax-xmin)*perc*dir;
  xmin+=xmove;
  xmax+=xmove;
  if(islogx)
    {
      xmin=10^xmin;
      xmax=10^xmax;
    }
  limits,xmin,xmax;
}


func zoom(fct,x=,y=)
{
  x=is_void(x)?1:(x?1:0);
  y=is_void(y)?1:(y?1:0);

  ll=limits();
  dx=0.5*(ll(2)-ll(1));
  dy=0.5*(ll(4)-ll(3));
  cx=0.5*(ll(2)+ll(1));
  cy=0.5*(ll(4)+ll(3));
  fct-=1;
  dx+=fct*x*dx;
  dy+=fct*y*dy;
  limits,cx-dx,cx+dx,cy-dy,cy+dy;
}



func plmkmc(x,y,color=,marker=,width=,msize=)
{
  if(is_void(marker)) {
    marker= (_plmk_count-1)%7 + 1;
    _plmk_count++;
  }
  
  typec=typeof(color);
  if(typec=="string") {
    ref_lst=color;
    ref_uni=suniq(color);
  }
  else if((typec=="long")||(typec=="int")) {
    ref_lst=color;
    ref_uni=where(histogram(color));
  }
  else if(typec=="char") {
    dms=dimsof(color);
    ref_lst=array(string,dms(0));
    for(i=1;i<=dms(0);i++) ref_lst(i)=string(&color(,i));
    ref_uni=suniq(ref_lst);
  }

  for(i=1;i<=numberof(ref_uni);i++)
    {
      idx=where(ref_lst==ref_uni(i));
      plmk,x(idx),y(idx),color=color(..,idx(1)),marker=marker,width=width,msize=msize;
    }
}

func get_inner_box(sys)
{
  local a,b,c,d;

  get_style,a,b,c,d;
  b=b(sys);

  box=b.viewport;
  box(1:2)+=[1,-1]*b.ticks.horiz.tickLen(max);
  box(1:2)+=[1,-1]*b.ticks.horiz.tickLen(max);

  box(3:4)+=[1,-1]*b.ticks.vert.tickLen(max);
  box(3:4)+=[1,-1]*b.ticks.vert.tickLen(max);
  return box;
  
}

func plTitle(title,offset=,font=,height=)
{
  local a,b,c,d;
  if(is_void(font  )) font=pltitle_font;
  if(is_void(height)) height=pltitle_height;
  if(is_void(offset)) offset=0.0;
  get_style,a,b,c,d;

  xmin=b.viewport(1,)(min);
  xmax=b.viewport(2,)(max);
  ypos=b.viewport(4,max)+offset;
  xpos=0.5*(xmin+xmax);

  plt,title,xpos,ypos,justify="CB",
    font=font,height=height;
}

func colorBar(sys,levels,cmin=,cmax=,flag=,dticks=,offset=)
{
  local a,b,c,d;
  local red,green,blue;
  if((numberof(b)>sys)||(sys==0)) error,"Bad system number !";
  if(is_void(flag  )) flag="r";
  flag=strlower(flag);

  verti_bar=1-(horiz_bar=(flag=="b")||(flag=="t"));
  
  width=0.02;
  get_style,a,b,c,d;
  vwp=b(sys).viewport;

  if(is_void(offset)) {
    if(verti_bar) offset=max(0.01,b.ticks.vert.tickLen(max));
    if(horiz_bar) offset=max(0.01,b.ticks.horiz.tickLen(max));
  }

  palette,red,green,blue,query=1;
  ncolor=numberof(red);
  
  xmin=vwp(1);xmax=vwp(2);nx=1;
  ymin=vwp(3);ymax=vwp(4);ny=1;

  if(flag=="t") {ymin=ymax+offset;ymax=ymin+width;nx=ncolor;}
  if(flag=="b") {ymax=ymin-offset;ymin=ymax-width;nx=ncolor;}
  if(flag=="r") {xmin=xmax+offset;xmax=xmin+width;ny=ncolor;}
  if(flag=="l") {xmax=xmin-offset;xmin=xmax-width;ny=ncolor;}

  if(dticks)
    {
      tmp=b(sys);
      if(flag=="b") tmp.ticks.horiz.flags=tmp.ticks.horiz.flags&~0x01;
      if(flag=="t") tmp.ticks.horiz.flags=tmp.ticks.horiz.flags&~0x02;
      if(flag=="l") tmp.ticks.vert.flags=tmp.ticks.vert.flags&~0x01;
      if(flag=="r") tmp.ticks.vert.flags=tmp.ticks.vert.flags&~0x02;
      b(sys)=tmp;
      set_style,a,b,c,d;
      redraw;
    }
  
  clr=[span(0,1,ncolor)];
  if((flag=="l")||(flag=="r")) clr=transpose(clr);
  
  old_sys=plsys(0);
  X=array(span(xmin,xmax,nx+1),ny+1);
  Y=transpose(array(span(ymin,ymax,ny+1),nx+1));
  plf,clr,Y,X;

  if(is_void(levels)) levels=spann(cmin,cmax,11);
  idx=where((levels>=cmin)&(levels<=cmax));
  if(is_void(idx))
    {
      plsys,old_sys;
      return;
    }
  levels=levels(idx);
  nlev=numberof(levels);
  text_offset=0.003;
  if(horiz_bar)
    {
      tymin=array(ymin,nlev);tymax=array(ymax,nlev);
      txmin=txmax=interp([xmin,xmax],[cmin,cmax],levels);
      if(flag=="b")
        {
          tymin-=text_offset;
          txt_px=txmin;
          txt_py=tymin-text_offset;
          just="CT";
        }
      else
        {
          tymax+=text_offset;
          txt_px=tyxmin;
          txt_py=tymax+text_offset;
          just="CB";
        }
    }
  else
    {
      txmin=array(xmin,nlev);txmax=array(xmax,nlev);
      tymin=tymax=interp([ymin,ymax],[cmin,cmax],levels);
      if(flag=="r")
        {
          txmax+=text_offset;
          txt_px=txmax+text_offset;
          txt_py=tymin;
          just="LH";
        }
      else
        {
          txmin-=text_offset;
          txt_px=txmin-text_offset;
          txt_py=tymin;
          just="RH";
        }
    }
  pldj,txmin,tymin,txmax,tymax;
  d=int(max(0,-floor(log10(levels(dif)(avg)))));
  fmt=swrite(format="%%.%df",d);

  labels=swrite(format=fmt,float(levels));
  zerol=swrite(format=fmt,-0.);
  idx=where(labels==zerol);
  if(is_array(idx)) labels(idx)=strpart(zerol,2:);
  plt1,labels,txt_px,txt_py,justify=just,tosys=0,height=10;
  plsys,old_sys;
}








