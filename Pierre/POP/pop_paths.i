// SOME ADDITIONAL PATHS
//usdir=exec("echo $HOME")(1);
//yodir=usdir+"/Yorick/";
//sdir=yodir+"STECKMAP/";
//#include "STECKMAP/Eric/system.i"

// first of all get the STECKMAPDIR from environment variable!!
$echo $STECKMAPROOTDIR > ~/yotmpfile
f=open("~/yotmpfile");
steckmaprootdir=rdline(f);
close,f;

steckmapdir=steckmaprootdir+"STECKMAP/";
//sdir="~/Yorick/STECKMAP/";
sdir=steckmapdir;
ericdir=sdir+"Eric";
chrisdir=sdir+"Chris";
bastdir=sdir+"Bastien";
pidir=sdir+"Pierre/POP";

set_path,".:"+Y_SITE+"i:"+Y_SITE+"contrib:"+Y_SITE+"i0:"+Y_HOME+"lib:"+ericdir+":"+bastdir+":"+chrisdir+":"+pidir+":"+steckmaprootdir;

//include, "utls.i",1;

//usdir=exec("echo $HOME")(1);
//yodir=usdir+"/Yorick/";   // yodir should not be used any more
basisdir=steckmapdir+"Pierre/POP/BASE/";
examplesdir=steckmapdir+"Pierre/POP/EXAMPLES/";
//GISTPATH=GISTPATH+":"+usdir+"/Yorick/Gist/";


pldefault, style=steckmapdir+"Gist/bboxed.gs", marks=0, width=1, palette=steckmapdir+"Gist/ridl-03.gp",legends=0;
//pldefault, style=libFunc+"/Gist/bboxed.gs", marks=0, width=1, palette=libFunc+"/Gist/ridl-03.gp",legends=0;

//if((!open("~/.init.i", "", 1))==0) include, "~/.init.i";

if(is_void(builtin_help)) builtin_help=help;
func help(fct)
{
  info,fct;
  if((typeof(fct)=="function")||(typeof(fct)=="builtin")) builtin_help,fct;
  else if (numberof(fct)==1) fct;
}

//set_path,"~/yorick:~/Yorick:"+Y_SITE+"i:"+Y_SITE+"contrib:"+Y_SITE+"i0:"+Y_HOME+"/lib:"+ericdir+":"+bastdir+":"+chrisdir+":"+pidir;


