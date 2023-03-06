#include "STECKMAP/Pierre/POP/sfit.i"
xnpal(0.);
ncolors=240;
x=span(0.,1.,ncolors);
c=char(3);
allc=array(char,[2,3,ncolors]);
for(i=1;i<=numberof(x);i++){allc(,i)=char(255.*xnpal(x(i)));};

f=open("~/Yorick/Gist/xn.gp","w");
write,f,"# Gist xn palette (for x neutral)\n",format="%s";
write,f,"# blabla\n",format="%s";
write,f,"\n",format="%s";
write,f,"ncolors=",format="%s";
write,f,ncolors;
write,f,"\n";
write,f,"\n";
write,f,"#ntsc blabla\n",format="%s";
write,f,"ntsc= 1\n",format="%s";
write,f,"\n",format="%s";
write,f,"#  r  g  b\n",format="%s";
for(i=1;i<=ncolors;i++){
  write,f,int(255.*xnpal(x(i)));
 };
close,f;

ws,1;
palette,"xn.gp";
plk,random(10,10);
