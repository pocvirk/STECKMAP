// SOME ADDITIONAL PATHS
//usdir=exec("echo $HOME")(1);
//yodir=usdir+"/Yorick/";
//sdir=yodir+"STECKMAP/";
#include "STECKMAP/Eric/system.i"
sdir="~/Yorick/STECKMAP/";
ericdir=sdir+"Eric";
chrisdir=sdir+"Chris";
bastdir=sdir+"Bastien";
pidir=sdir+"Pierre/POP";

set_path,".:~/yorick:~/Yorick:"+Y_SITE+"i:"+Y_SITE+"contrib:"+Y_SITE+"i0:"+Y_HOME+"lib:"+ericdir+":"+bastdir+":"+chrisdir+":"+pidir;

//include, "utls.i",1;

usdir=exec("echo $HOME")(1);
yodir=usdir+"/Yorick/";
basisdir=usdir+"/Yorick/STECKMAP/Pierre/POP/BASE/";
examplesdir=usdir+"/Yorick/STECKMAP/Pierre/POP/EXAMPLES/";
//GISTPATH=GISTPATH+":"+usdir+"/Yorick/Gist/";


pldefault, style=usdir+"/Yorick/STECKMAP/Gist/bboxed.gs", marks=0, width=1, palette=usdir+"/Yorick/STECKMAP/Gist/ridl-03.gp",legends=0;
if((!open("~/.init.i", "", 1))==0) include, "~/.init.i";

if(is_void(builtin_help)) builtin_help=help;
func help(fct)
{
  info,fct;
  if((typeof(fct)=="function")||(typeof(fct)=="builtin")) builtin_help,fct;
  else if (numberof(fct)==1) fct;
}

//set_path,"~/yorick:~/Yorick:"+Y_SITE+"i:"+Y_SITE+"contrib:"+Y_SITE+"i0:"+Y_HOME+"/lib:"+ericdir+":"+bastdir+":"+chrisdir+":"+pidir;


