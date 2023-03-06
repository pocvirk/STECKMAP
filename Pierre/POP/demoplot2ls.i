#include "plot2ls.i"
//#include "Pierre/POP/sfit.i"
// faire en sorte que ws, plh, random soient definis 

ws;
plh,random(100),type=1;
plh,random(100),type=2;
plh,random(100),type=3;
           
ylims;  // fait un peu de place en haut
// si ca ne suffit pas, bidouille avec range,y1,y2
typ=[1,2,3];
s=["machin","truc","chose"];
mlabel,s,[10.,1.06],[30.,1.4],t=typ;

