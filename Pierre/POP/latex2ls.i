func gencomarc(filename,n=){
  /* DOCUMENT
     generates a concatenated filelist containing all the paths to the figures used in the latex file by includegraphics incantation to feed a tar command
     the \includegraphics incantations are counted from the n-th line (default n=20)
     writes the list in crap.lst
  */
  if(is_void(n)) n=20;
  nl=wordcount(filename)(2);
  f=open(filename);
  a=rdline(f,nl);
  a=a(n:);
  ni=where((strmatch(a,"includegraphics")==1)&(strmatch(a,"%")!=1));
  res=[];
  
  for(i=1;i<=numberof(ni);i++){
  grow,res," "+split2words(split2words(strreplace(a(ni(i)),"includegraphics",""),sep="{")(0),sep="}")(1)+"*";
  };

  f=open("crap.lst","w");
  write,f,strglue(res,sep=" ");
  close,f;
  
  return res;

};

func truncarray(a,n){
  /* DOCUMENT
     used to truncate arrays to the n-th decimal number after point.
     useful for printing arrays nicely
  */

  return a-((10.^(-n))*((10.^(n))*a-int(a*(10.^(n)))));
};

func prtable(a){
  /* DOCUMENT
     does pr1 for a whole table, unlike pr1
  */

  b=array(string,dimsof(a));
  na=numberof(a);
  
  for(i=1;i<=na;i++){b(i)=pr1(a(i));};
  return b;
};
  

func textable(a,h=,f=){
  /* DOCUMENT
     produces a latex formatted table of 2d string array a
     with header h and footer h
  */

  nl=dimsof(a)(0);
  nc=dimsof(a)(-1);

  betab="\\begin{tabular}{c";
  for(i=1;i<=nc-1;i++){betab+="c";};
  betab+="} \n";
  if(!is_void(h)) grow,betab,h;
  for(j=1;j<=nl;j++){
    ltab="";
    for(i=1;i<=nc;i++){
      ltab+=pr1(a(i,j));
      if(i==nc) continue;
      ltab+=" & ";
    };
    grow,betab,ltab+" \\\\ \n";
  };
  if(!is_void(f)) grow,betab,f;
  grow,betab,"\\end{tabular}";
  res=betab;
  return res;
}



