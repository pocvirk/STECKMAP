#include "STECKMAP/Pierre/POP/sfit.i"
#include "ascii2ls.i"

func write_psfrag(filename){
  /* DOCUMENT
     writes a tex file for psfrag containing a number of default tags
     including the ps figure filename
  */

  //  f= open("~/STECKMAP/Pierre/POP/tex/psfragcrap.tex","w");
  //  a=stextread("~/Yorick/STECKMAP/Pierre/POP/tex/psftest.tex");
  //  a=stextread(pidir+"/tex/psftest.tex",0);
  a=exec("cat "+pidir+"/tex/psftest.tex");
  a=strreplace(a,"SFR.ps",filename);
  f=open("psfragcrap.tex","w");
  write,f,a+"\n";
  close,f;
  return a;
};
  
func replace_tags(filename,targetps,nodisp=){
  /*DOCUMENT
    replaces the tags in filename with the default tags defined in psftest.tex
    then run latex on the temp file
    then dvips the dvi
    then displays it if nodisp!=1
    notice you need to have a ps file first!
  */

  if(is_void(targetps)) targetps = "psfragcrap.ps";
  bar=write_psfrag(filename);
  foo=exec("latex psfragcrap.tex");
  //bar=exec("dvips -o "+filename+ "psfragcrap.dvi"); // non cant do that cause the dvi itself invokes the ps. need to create a temp output
  bar=exec("dvips -o "+targetps+ " psfragcrap.dvi");
  if (nodisp!=1) foo=exec("gv "+targetps);
};
    
  

  

  
    
