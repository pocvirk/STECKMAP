#include "Bastien/utils.i"

if(is_void(__SetNameOfColors))
  __SetNameOfColors=0;

struct colorStruct {
  string name;
  char rgb(3);
}

_color_file_="/usr/lib/X11/rgb.txt";


func readRGB(filename)
{
  lines=sreadAscii(filename,cmt="!",verb=0);
  flag=array(0n,numberof(lines));
  for(i=1;i<=numberof(lines);i++)
    {
      words=split2words(lines(i));
      if(numberof(words)==4) flag(i)=1;
    }
  idx=where(flag);
  lines=lines(idx);
  a=b=c=array(0n,numberof(idx));
  s=    array(string,numberof(idx));
  sread,format="%d %d %d %s",lines,a,b,c,s;
  rslt=array(colorStruct,numberof(idx));
  rslt.name=s;
  rslt.rgb=transpose(char([a,b,c]));
    
  return rslt;
}

colorList=readRGB(_color_file_);

if(__SetNameOfColors&&is_void(__grey60))
  for(i=1;i<=numberof(colorList);i++)
    symbol_set,"__"+colorList(i).name,colorList(i).rgb;


func dispColor(lumi,type=)
{
  if(is_void(type)) type=1;

  nbc=numberof(colorList);
  xc=[0,0,1,1];
  yc=[0,1,1,0];
  if(type==1)
    {
      n=int(sqrt(nbc));
      if(n*n<nbc) ++n;
      k=0;
      for(i=1;i<=n;i++)
        for(j=1;j<=n;j++)
          {
            k++;
            if(k>nbc) break;
            plfp,[colorList(k).rgb],xc+i,yc+j,[4];
          }
    }
  else if (type==2)
    {
      for(i=1;i<=nbc;i++)
        {
          crgb=colorList(i).rgb;
          plfp,[crgb],xc+crgb(1)+crgb(3)*256,yc+crgb(2),[4];
        }
    }
  else if (type==3)
    {
      phi=span(0,255,256);
      the=span(0,255,256);
      r=char(255*lumi*sin(pi*the*0.5/255)(,-)*cos(pi*phi*0.5/255)(-,));
      g=char(255*lumi*sin(pi*the*0.5/255)(,-)*sin(pi*phi*0.5/255)(-,));
      b=char(255*lumi*cos(pi*the*0.5/255));
      xc=[0,0,1,1];
      yc=[0,1,1,0];
      for(i=1;i<=256;i++)
        for(j=1;j<=256;j++)
          plfp,[char([r(j,i),g(j,i),b(j)])],xc+the(i),yc+phi(j),[4];
      return transpose([r,g,b]);
    }
}
