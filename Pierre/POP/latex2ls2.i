func WSLpsnup(file,nodisp=){
  /*DOCUMENT
    handles WSL-like plots: creates a ps with 1 WSL plots per page, to hepl browsing through, using latex
  */
  
  if(is_void(nodisp)) nodisp=0;
  
  f=open("crap.tex","w");
  write,f,"\\documentclass{article}";
  write,f,"\\usepackage[T1]{fontenc}";
  write,f,"\\usepackage[latin1]{inputenc}";
  write,f,"\\usepackage{graphics}";
  write,f,"\\textwidth=16.5cm";
  write,f,"\\textheight=27.5cm";
  write,f,"\\topmargin=-3.cm";
  write,f,"\\evensidemargin=-2.5cm";
  write,f,"\\oddsidemargin=-2.5cm";
  write,f,"\\begin{document}";
  for(i=1;i<=numberof(file);i++){
    write,f,"\\begin{figure}";
    write,f,"\\includegraphics{"+file(i)+"}";
    write,f,"\\end{figure}";
  };
  write,f,"\\end{document}";
  close,f;
  exec("latex crap.tex");
  exec("dvips -o crap.ps crap.dvi");
  if(nodisp==0) exec("gv -landscape crap.ps");
  
};
  
