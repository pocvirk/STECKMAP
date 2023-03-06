#include "style.i"

func copyWin(win,wout,lmt=,clear=)
/* DOCUMENT copyWin,window_source,window_target
     Copy the style and the graphic elements of window window_source to
     window window_target
   KEYWORDS: lmt=   if set, copy also the limits (the default)
             clear= if set erase the target window before copying (the default)
   SEE ALSO:replotAll,replot1Sys
 */
{
  local a,b,c,d;

  if(is_void(clear)) clear=1;
  old_win=current_window();
  window,win;
  get_style,a,b,c,d;
  window,wout;
  set_style,a,b,c,d;
  if(clear) fma;
  replotAll,win,wout,lmt=lmt;
  if(old_win<0) old_win=wout;
  window,old_win;
}


func replotAll(win,wout,lmt=)
/* DOCUMENT replotAll(window_source,window_target)
     Copy graphical elements of all systems (but not the style)
     of window window_source to window_target
   KEYWORDS: lmt : if set, copy also the limits (the default)
   SEE ALSO: copyWin,replot1Sys
 */
{
  old_win=current_window();
  nbsin =getNbSys(win);
  nbsout=getNbSys(wout);
  if(nbsin!=nbsout) error,"The two windows have not the same number of systems !!";
  for(i=0;i<=nbsin;i++)
    replot1Sys,win,i,wout,i,lmt=lmt;
  window,wout;
  redraw;
  if(old_win<0) old_win=wout;
  window,old_win;
}

func replot1Sys(win,sin,wout,sout,lmt=)
/* DOCUMENT replot1Sys(window_source,sys_source,winow_target,system_target)
   Copy all the graphical elements of the system sys_source
   of window_source to the system sys_target of window window_target

   KEYWORDS: lmt : if set, copy also the limits (the default)

   SEE ALSO:copyWin,replotAll
 */
{
  if(is_void(lmt)) lmt=1;
  window,win;
  old_sin=plsys(sin);
  lmtin=limits();
  nbobj=numberof(plq());
  window,wout;
  old_sout=plsys(sout);
  if(lmt) limits,lmtin;
  for(i=1;i<=nbobj;i++)
    {
      window,win;
      plsys,sin;
      prop=plq(i);
      window,wout;
      plsys,sout;
      replot,prop;
    }
  plsys,old_sout;
  window,win;
  plsys,old_sin;
}

func savePlot(filename,win)
/* DOCUMENT savePlot(filename,win)
   Save Yorick plot of window win in file 'filename'.
   The plot can be reload with loadPlot.

   EXAMPLE:
   window,0;
   pli,random(100,100);
   pltitle,"Random array";
   limits,20,60,10,70;
   plt,"Zoom In",40,40,tosys=1,color="yellow",height=18;
   savePlot,"rand_array.gdb",0;
   loadPlot,"rand_array.gdb",1;
   
   SEE ALSO: loadPlot,copyWin
 */
{
  local a,b,c,d,x,y,z;
  local x0,x1,y0,y1,txt;
  local ireg;
  local p1,p2,p3,p4,p5;

  fstrm=createb(filename);
  old_win=current_window();
  if(old_win>=0) old_sys=plsys();
  window,win;
  get_style,a,b,c,d;
  saveVariable,fstrm,"getstyle_p1",a;
  saveVariable,fstrm,"getstyle_p2",b;
  saveVariable,fstrm,"getstyle_p3",c;
  saveVariable,fstrm,"getstyle_p4",d;

  nbsys=getNbSys(win);
  for(i=0;i<=nbsys;i++)
    {
      plsys,i;
      lmt=limits();
      nbobj=numberof(plq());
      saveVariable,fstrm,swrite(format="system_%d",i),i;
      saveVariable,fstrm,swrite(format="limits_%d",i),lmt;
      for(j=1;j<=nbobj;j++)
        {
          prop=plq(j);
          decompProp,prop,p1,p2,p3,p4,p5;
          saveVariable,fstrm,swrite(format="prop1_%d_%d",i,j),(is_void(p1)?"dummy":p1);
          saveVariable,fstrm,swrite(format="prop2_%d_%d",i,j),(is_void(p2)?"dummy":p2);
          saveVariable,fstrm,swrite(format="prop3_%d_%d",i,j),(is_void(p3)?"dummy":p3);
          saveVariable,fstrm,swrite(format="prop4_%d_%d",i,j),(is_void(p4)?"dummy":p4);

          rslt=reshapeProp(prop);
          saveVariable,fstrm,swrite(format="prop5_%d_%d",i,j),(is_void(rslt)?"dummy":rslt);
        }
    }
  close,fstrm;
  if(old_win>=0)
    {
      window,old_win;
      plsys,old_sys;
    }
}

func loadPlot(filename,wout,clear=,lmt=)
/* DOCUMENT loadPlot(filename,wout)
   Load Yorick plot form file 'filename' in window wout
   The plot have to be saved with savePlot.

   EXAMPLE:
   window,0;
   pli,random(100,100);
   pltitle,"Random array";
   limits,20,60,10,70;
   plt,"Zoom In",40,40,tosys=1,color="yellow",height=18;
   savePlot,"rand_array.gdb",0;
   loadPlot,"rand_array.gdb",1;

   KEYWORDS: lmt=   if set (default), restore also the
                    imits
             clear= if set (default) erase the window
                    before loading 
   
   SEE ALSO: savePlot,copyWin
 */
{
  if(is_void(clear)) clear=1;
  if(is_void(  lmt)) lmt=1;
  fstrm=openb(filename);
  old_win=current_window();
  if(old_win>=0) old_sys=plsys();

  window,wout;
  set_style,fstrm.getstyle_p1,fstrm.getstyle_p2,fstrm.getstyle_p3,fstrm.getstyle_p4;
  if(clear) fma;
  names=*get_vars(fstrm)(1);
  nnames=numberof(names);
  idx=4;
  while(++idx<=nnames)
    {
      if(strmatch(names(idx),"system_"))
        {
          plsys,get_member(fstrm,names(idx));
          limits;
          continue;
        }
      if(strmatch(names(idx),"limits_"))
        {
          if(lmt) limits,get_member(fstrm,names(idx));
          continue;
        }

      if(strmatch(names(idx),"prop1_"))
        {
          p1  =get_member(fstrm,names(idx));
          p2  =get_member(fstrm,names(++idx));
          p3  =get_member(fstrm,names(++idx));
          p4  =get_member(fstrm,names(++idx));
          rslt=get_member(fstrm,names(++idx));
          
          replot,p1,p2,p3,p4,rslt;
          continue;
        }
      write,format="[WARNING] Unknown variable flag %s !!\n",names(idx);
    }
  close,fstrm;
  redraw;
  if(old_win>=0)
    {
      window,old_win;
      plsys,old_sys;
    }
}

func reshapeProp(prop)
/* DOCUMENT reshapeProp(prop)
     Return an array of pointers that point to the data useful for
     drawing the graphical element. The size of the array depend
     of the graphical type, and the order of the pointers is the
     same than the one of the plq() fifth element.
 */
{
  aig=(*prop(1))(1);
  if(aig==1)      return reshape_plg(prop);
  else if(aig==2) return reshape_pldj(prop);
  else if(aig==3) return reshape_plt(prop);
  else if(aig==4) return reshape_plm(prop);
  else if(aig==5) return reshape_plf(prop);
  else if(aig==6) return reshape_plv(prop);
  else if(aig==7) return reshape_plc(prop);
  else if(aig==8) return reshape_pli(prop);
  else if(aig==9) return reshape_plfp(prop);
  else if(aig!=0) error,"Unknown graphical element !!";
}


func replot(p1,p2,p3,p4,rslt)
/* DOCUMENT replot(prop)
            replot(p1,p2,p3,p4,rslt)
   Replot a graphical element with its plq properties.
   prop is the array return by plq(i). then replot(plq(i))
   will replot the graphical element i.
     
   SEE ALSO: copyWin,replotAll,replot1Sys
 */
{
  if(is_void(p2)||is_void(p3)||is_void(p4)||is_void(rslt))
    {
      aig=(*prop(1))(1);
      p1=*prop(1);
      p2=*prop(2);
      p3=*prop(3);
      p4=*prop(4);
      rslt=reshapeProp(prop);
    }
  else
    aig=p1(1);
  if(aig==1)      replot_plg ,p1,p2,p3,p4,rslt;
  else if(aig==2) replot_pldj,p1,p2,p3,p4,rslt;
  else if(aig==3) replot_plt ,p1,p2,p3,p4,rslt;
  else if(aig==4) replot_plm ,p1,p2,p3,p4,rslt;
  else if(aig==5) replot_plf ,p1,p2,p3,p4,rslt;
  else if(aig==6) replot_plv ,p1,p2,p3,p4,rslt;
  else if(aig==7) replot_plc ,p1,p2,p3,p4,rslt;
  else if(aig==8) replot_pli ,p1,p2,p3,p4,rslt;
  else if(aig==9) replot_plfp,p1,p2,p3,p4,rslt;
  else if(aig!=0) error,"Unknown graphical element !!";
}

func reshape_plg(prop)
/* DOCUMENT reshape_plg(prop)
   Return [&x,&y], where x and y are the array used for the plot (plg,y,x).
     
   SEE ALSO: reshapeProp
 */
{
  local y,x;
  local p1,p2,p3,p4,p5;
  decompProp,prop,p1,p2,p3,p4,p5;
  if(p1(1)!=1) error,"Bad properties !!";
  adx=p5(2);reshape,x,adx,double,[1,p5(1)];
  ady=p5(3);reshape,y,ady,double,[1,p5(1)];
  return [&x,&y];
}

func replot_plg(p1,p2,p3,p4,rslt)
/* DOCUMENT replot_plg(prop)
   Same as replot but only for plg graphical elements
   SEE ALSO: replot
 */
{
  if(p1(1)!=1) error,"Bad properties !!";

  plg,*rslt(2),*rslt(1),hide=p1(2),legend=p2,
    color= getColor(p3(1)),type=p3(2),marks= p3(3),
    mcolor=getColor(p3(4)),marker=p3(5),rays=p3(6),
    closed=p3(7),smooth=p3(8),width= p4(1),msize=p4(2),
    mspace=p4(3),rspace=p4(4),rphase=p4(5),arrowl=p4(6),arroww=p4(7);
}


func reshape_pldj(prop)
{
  local x0,y0,x1,y1;
  local p1,p2,p3,p4,p5;
  decompProp,prop,p1,p2,p3,p4,p5;
  if(p1(1)!=2) error,"Bad properties !!";
  adx0=p5(2);reshape,x0,adx0,double,[1,p5(1)];
  ady0=p5(3);reshape,y0,ady0,double,[1,p5(1)];
  adx1=p5(4);reshape,x1,adx1,double,[1,p5(1)];
  ady1=p5(5);reshape,y1,ady1,double,[1,p5(1)];
  return [&x0,&y0,&x1,&y1];
}

func replot_pldj(p1,p2,p3,p4,rslt)
/* DOCUMENT replot_pldj(prop)
   Same as replot but only for pldj graphical elements
   SEE ALSO: replot
 */
{
  if(p1(1)!=2) error,"Bad properties !!";
  pldj,*rslt(1),*rslt(2),*rslt(3),*rslt(4),hide=p1(2),legend=p2,
    color= getColor(p3(1)),type=p3(2),
    width= p4(1);
}

func reshape_plt(prop)
{
  local txt;
  local p1,p2,p3,p4,p5;
  decompProp,prop,p1,p2,p3,p4,p5;
  if(p1(1)!=3) error,"Bad properties !!";
  adt=p5(2);reshape,txt,adt,char,[1,p5(1)];
  return [&txt];
}


func replot_plt(p1,p2,p3,p4,rslt)
/* DOCUMENT replot_plt(prop)
   Same as replot but only for plt graphical elements
   SEE ALSO: replot
 */
{
  if(p1(1)!=3) error,"Bad properties !!";
  plt,string(rslt(1)),p4(2),p4(3),hide=p1(2),legend=p2,tosys=plsys(),
    color= getColor(p3(1)),font=p3(2),orient= p3(3),justify=p3(4),opaque=p3(5),
    height= p4(1);
}

func reshape_plm(prop)
{
  local p1,p2,p3,p4,p5;
  local x,y,ireg;
  decompProp,prop,p1,p2,p3,p4,p5;
  if(p1(1)!=4) error,"Bad properties !!";
  adx=p5(3);reshape,x,adx,double,[2,p5(1),p5(2)];
  ady=p5(4);reshape,y,ady,double,[2,p5(1),p5(2)];
  adi=p5(5);reshape,ireg,adi,int,[2,p5(1),p5(2)];
  return [&x,&y,&ireg];
}

func replot_plm(p1,p2,p3,p4,rlst)
/* DOCUMENT replot_plm(prop)
   Same as replot but only for plm graphical elements
   SEE ALSO: replot
 */
{
  if(p1(1)!=4) error,"Bad properties !!";
  plm,*rslt(2),*rslt(1),*rslt(3),legend=p2,hide=p1(2),
    color=getColor(p3(1)),type=p3(2),region=p3(3),
    boundary=p3(4),inhibit=p3(5),width=p4(1);
}

func reshape_plf(prop)
{
  local p1,p2,p3,p4,p5;
  local z,y,x,ireg;
  decompProp,prop,p1,p2,p3,p4,p5;
  if(p1(1)!=5) error,"Bad properties !!";
  adx=p5(3);reshape,x,adx,double,[2,p5(1),p5(2)];
  ady=p5(4);reshape,y,ady,double,[2,p5(1),p5(2)];
  adi=p5(5);reshape,ireg,adi,double,[2,p5(1),p5(2)];
  adc=p5(6);
  if(p3(4))
    reshape,z,adc,char,[3,3,p5(1)-1,p5(2)-1];
  else
    reshape,z,adc,char,[2,p5(1)-1,p5(2)-1];
  return [&x,&y,&ireg,&z];
}


func replot_plf(p1,p2,p3,p4,rslt)
/* DOCUMENT replot_plf(prop)
   Same as replot but only for plf graphical elements
   SEE ALSO: replot
 */
{
  if(p1(1)!=5) error,"Bad properties !!";
  plf,*rslt(4),*rslt(2),*rslt(1),legend=p2,hide=p1(2),
    region=p3(1),edges=p3(2),ecolor=getColor(p3(3)),
    ewidth=p4(1);
}

func reshape_plv(prop)
{
  local p1,p2,p3,p4,p5;
  local x,y,ireg,vx,vy;
  decompProp,prop,p1,p2,p3,p4,p5;
  if(p1(1)!=6) error,"Bad properties !!";
  adx =p5(3);reshape,x,adx,double,[2,p5(1),p5(2)];
  ady =p5(4);reshape,y,ady,double,[2,p5(1),p5(2)];
  adi =p5(5);reshape,ireg,adi,int,[2,p5(1),p5(2)];
  advx=p5(6);reshape,vx,advx,double,[2,p5(1),p5(2)];
  advy=p5(7);reshape,vy,advy,double,[2,p5(1),p5(2)];
  return [&x,&y,&ireg,&vx,&vy];
}

func replot_plv(p1,p2,p3,p4,rslt)
/* DOCUMENT replot_plv(prop)
   Same as replot but only for plv graphical elements
   SEE ALSO: replot
 */
{
  if(p1(1)!=6) error,"Bad properties !!";
  plv,*rslt(5),*rslt(4),*rslt(2),*rslt(1),*rslt(3),hide=p1(2),legend=p2,
    region=p3(1),color=getColor(p3(2)),hollow=p3(3),
    width=p4(1),aspect=p4(2),scale=p4(3);
}

func reshape_plc(prop)
{
  local p1,p2,p3,p4,p5;
  local x,y,ireg,z,tria,levels;
  decompProp,prop,p1,p2,p3,p4,p5;
  if(p1(1)!=7) error,"Bad properties !!";
  adx=p5(3);reshape,x,adx,double,[2,p5(1),p5(2)];
  ady=p5(4);reshape,y,ady,double,[2,p5(1),p5(2)];
  adi=p5(5);reshape,ireg,adi,int,[2,p5(1),p5(2)];
  adz=p5(6);reshape,z,adz,double,[2,p5(1),p5(2)];
  adt=p5(7);reshape,tria,adt,short,[2,p5(1),p5(2)];
  adl=p5(9);reshape,levels,adl,double,[1,p5(8)];
  return [&x,&y,&ireg,&z,&tria,&levels];
}


func replot_plc(p1,p2,p3,p4,rslt)
/* DOCUMENT replot_plc(prop)
   Same as replot but only for plc graphical elements
   SEE ALSO: replot
 */
{
  if(p1(1)!=7) error,"Bad properties !!";
  plc,*rslt(4),*rslt(2),*rslt(1),*rslt(3),levs=*rslt(6),legend=p2,hide=p1(2),
    region=p3(1),color=getColor(p3(2)),type=p3(3),
    marks=p3(4),mcolor=p3(5),marker=p3(6),smooth=p3(7),
    width=p4(1),msize=p4(2),mspace=p4(3),mphase=p4(4),
    triangle=*rslt(5);
}

func reshape_pli(prop)
{
  local p1,p2,p3,p4,p5;
  local colors;
  decompProp,prop,p1,p2,p3,p4,p5;
  if(p1(1)!=8) error,"Bad properties !!";
  adc=p5(3);reshape,colors,adc,char,[2,p5(1),p5(2)];
  return [&colors];
}

func replot_pli(p1,p2,p3,p4,rslt)
/* DOCUMENT replot_pli(prop)
   Same as replot but only for pli graphical elements
   SEE ALSO: replot
 */
{
  if(p1(1)!=8) error,"Bad properties !!";
  pli,*rslt(1),p4(1),p4(2),p4(3),p4(4),legend=p2,hide=p1(2);
}

func reshape_plfp(prop)
{
  local p1,p2,p3,p4,p5;
  local pn,x,y,z;
  decompProp,prop,p1,p2,p3,p4,p5;
  if(p1(1)!=9) error,"Bad properties !!";
  adpn=p5(5);reshape,pn,adpn,long,[1,p5(1)];
  adx =p5(3);reshape,x,adx,double,[1,sum(pn)];
  ady =p5(2);reshape,y,ady,double,[1,sum(pn)];
  adz =p5(4);
  if(p5(1)&&p5(4))
    if(p3(3)) //RGB mode
      reshape,z,adz,char,[2,3,p5(1)];
    else
      reshape,z,adz,char,[1,p5(1)];
  else
    z=[];
  
  if(is_void(z)) return [&x,&y,&pn];
  return [&x,&y,&z,&pn];
}

func replot_plfp(p1,p2,p3,p4,rslt)
/* DOCUMENT replot_plfp(prop)
   Same as replot but only for plfp graphical elements
   SEE ALSO: replot
 */
{
  if(p1(1)!=9) error,"Bad properties !!";
  if(numberof(rslt)==3)
    plfp,,*rslt(1),*rslt(2),*rslt(3),hide=p1(2),legend=p2,
      edges=p3(1),ecolor=getColor(p3(2)),ewidth=p4(1);
  else
    plfp,*rslt(3),*rslt(1),*rslt(2),*rslt(4),hide=p1(2),legend=p2,
      edges=p3(1),ecolor=getColor(p3(2)),ewidth=p4(1);
}

func decompProp(prop,&prop1,&prop2,&prop3,&prop4,&prop5)
/* DOCUMENT  decompProp(prop,&p1,&p2,&p3,&p4,&p5)
   Put in propi the array *prop(i);
 */
{
  prop1=*prop(1);
  prop2=*prop(2);
  prop3=*prop(3);
  prop4=*prop(4);
  prop5=*prop(5);
}

func getNbSys(win)
/* DOCUMENT getNbSys(win)
   Return the numberof og system in the window 'win'. If win is ommited, return
   the numberof of system in the current window.
   This routine does not change the selected window.
   SEE ALSO:
 */
{
  local a,b,c,d;
  old_win=current_window();
  window,win;
  get_style,a,b,c,d;
  window,old_win;
  return numberof(b);
}

func getColor(color)
/* DOCUMENT getColor(color)
   Check if color is a rgb encoded color. If so return the array
   [red,green,blue], else just return color.
 */
{
  if(!(color&0xF000000)) return color;
  red  =(color&0x00000FF);
  green=(color&0x000FF00)>>8;
  blue =(color&0x0FF0000)>>16;
  return char([red,green,blue]);
}

func saveVariable(fstrm,name,value)
/* DOCUMENT saveVariable(fstrm,name,value)
   Save the array 'value' in the binary file stream 'fstrm' under the name
   'name'.
     
   SEE ALSO:
 */
{
  add_variable,fstrm,-1,name,structof(value),dimsof(value);
  get_member(fstrm,name)=value;
}

