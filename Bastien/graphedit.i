
func putt(text,cosys=,legend=,hide=,color=,font=,height=,opaque=,orient=,justify=)
/* DOCUMENT putt,text

     same as plt but uses mouse to put text

   KEYWORDS: cosys=  selects coordinate system.
                     default: -1

   SEE ALSO: movet, textinfo, textedit, textselect,
             mousefindt, textfindt, put_indx, savet.
*/
{
  if(is_void(justify)) justify="CH";
  if(is_void(legend)) legend=text;
  if(is_void(cosys)) cosys=-1;
  r=mouse(cosys,0,"Click at the position to put text :'"+text+"'");

  cosysbck=plsys(int(r(-2)));
  plt,text,r(1),r(2),tosys=int(r(-2)),justify=justify,legend=legend,hide=hide,color=color,font=font,height=height,opaque=opaque,orient=orient;
  plsys,cosysbck;
}

func movet(..,cases=,cosys=)
/* DOCUMENT    movet
            or movet,text

    selects a text with mouse (first form) or with
    value (second form) and moves it.

   KEYWORDS: cases=  with second form
                     1 case sensitive
                     0 not case sensitive (default)
             cosys=  selects coordinate system.
                     default: -1

   SEE ALSO: putt, textinfo, textedit, textselect,
             mousefindt, textfindt, put_indx, savet.
*/
{
  local mode,indx,txtadr;
  local str2find,txtp,txts;

  mode=0;
  if(more_args())
    {
      mode=1;
      str2find=next_arg();
    }
  if(is_void(cases)) cases=0;
  if(is_void(cosys)) cosys=-1;
  cosysbck=plsys();

  if(mode==0) indx=textselect(mousefindt(cosys=cosys));
  if(mode==1) indx=textselect(textfindt(str2find,cases=cases,cosys=cosys));

  indx=indx(*);
  plsys,indx(2);
  par=plq(indx(1));
  if ((*par(1))(2)==1)
    {
      write,"!! Text hide. Show it !!";
      pledit,indx(1),hide=0;
    }

  postx=(*par(4))(2);
  posty=(*par(4))(3);
  lmts=limits();
  if (postx<lmts(1)||postx>lmts(2)||posty<lmts(3)||posty>lmts(4))
    { write,"!! Text out of the window !!"; }

  colorbck=(*par(3))(1);
  colornew=-5-2*(colorbck==-5);
  pledit,indx(1),color=colornew;

  txtadr=(*par(5))(2);
  reshape,txtp,txtadr,char,(*par(5))(1);
  txts=string(&txtp);
  write,format=" - Text selected : '%s'\n",txts;
  
  r=mouse(-1,0," ## Click at the new position ##");
  while(r(-1)==1)
    {
      par=plq(indx(1));
      postx=(*par(4))(2);
      posty=(*par(4))(3);
      convertsys,r(1),r(2),int(r(-2)),dx,dy,indx(2);
      dx=dx-postx;
      dy=dy-posty;
      pledit,indx(1),dx=dx,dy=dy,justify="CH";
      r=mouse(-1,0," ## Click at the new position ##");
    }
  pledit,indx(1),color=colorbck;
  plsys,cosysbck;
  return;
}

func textinfo(indx)
/* DOCUMENT textinfo,indx

     Displays informations about text of index indx.
     indx: index of the text (result of mousefindt or textfindt)

   KEYWORD: cosys=  selects coordinate system.
                    default: 1

   SEE ALSO: putt, movett, textedit, textselect,
             mousefindt, textfindt, put_indx, savet.
*/
{
  local txtadr,txtp,txts,par,vindx,cosysbck;

  cosysbck=plsys();
  vindx=indx(*);

  plsys,vindx(2);
  par=plq(vindx(1));
  txtadr=(*par(5))(2);
  reshape,txtp,txtadr,char,[1,(*par(5))(1)];
  txts=string(&txtp);
  noyes=["No","Yes"];
  jsth=["N","L","C","R"];
  jstv=["N","T","C","H","A","B"];
  ttmp=*par(3);
  write,"-------------";
  write,format=" [1] Text value  : '%s'\n",txts;
  write,format=" [2] hide        : %s\n",noyes((*par(1))(2)+1);
  write,format=" [3] color       : %i\n",ttmp(1);
  write,format=" [4] font        : %i\n",ttmp(2);
  write,format=" [5] orient      : %i\n",ttmp(3);
  write,format=" [6] justify     : \"%1s%1s\"\n",jsth(ttmp(4)%4+1),jstv(int(ttmp(4)/4)+1);
  write,format=" [7] opaque      : %i\n",ttmp(5);
  write,format=" [8] height      : %-6.2f\n",(*par(4))(1);
  write,format=" [9] x,y         : %g,%g\n",(*par(4))(2),(*par(4))(3);
  write,"-------------";
  plsys,cosysbck;
  return;
}

func textedit(indx,cases=,cosys=)
/* DOCUMENT textedit,indx
        or  textedit,txt
        or  textedit

     - first form,edits text of index 'indx' (result of
       mousefindt or textfindt).
     - second form, same as:
         textedit,textfindt(txt)
     - third form, same as:
         textedit,mousefindt()

   KEYWORDS: cases=  1 case sensitive
                     0 not case sensitive (default)
             cosys=  selects coordinate system.
                     default: 1

   SEE ALSO: putt, movett, textinfo, textselect,
             mousefindt, textfindt, put_indx, savet.
*/
{
  local txtadr,txtp,txts,chx,notend;
  local cx,cy,nx,ny;
  local ccol;

  if(is_void(indx)) indx=mousefindt(cosys=cosys);
  if(typeof(indx)=="string") indx=textfindt(indx,cases=cases,cosys=cosys);

  vindx=textselect(indx)(*);
  notend=1;
  chx=0n;
  vindx;
  sysbck=plsys(vindx(2));
  while(notend)
    {
      textinfo,vindx;
      read,prompt="Your choose ([0] quit) : ",chx;
      if(!numberof(chx)) chx=-1;
      chx=int(chx);
      if(chx==0) notend=0;
      if (chx<=0||chx>9) continue;

      par=plq(vindx(1));

      if(chx==1)
	{
	  txtadr=(*par(5))(2);
	  ccol="";
	  read,prompt="New text : ",ccol;
	  txtp=(*pointer(ccol));
	  nc2=numberof(txtp);
	  reshape,txts,txtadr,char,[1,nc2];
	  (*par(5))(1)=nc2-1;
	  txts=txtp;
	  redraw;
	  continue;
	}
      if(chx==2)
	{
	  pledit,vindx(1),hide=(1-(*par(1))(2));
	  continue;
	}
      if(chx==3)
	{
	  ccol=0n;
	  read,prompt="New color index : ",ccol;
	  if(numberof(ccol)) pledit,vindx(1),color=ccol;
	  continue;
	}
      if(chx==4)
	{
	  ccol=0n;
	  read,prompt="New font index : ",ccol;
	  if(numberof(ccol)) pledit,vindx(1),font=ccol;
	}
      if(chx==5)
	{
	  ccol=0n;
	  read,prompt="New orientation : ",ccol;
	  if(numberof(ccol)) pledit,vindx(1),orient=ccol;
	}
      if(chx==6)
	{
	  ccol="";
	  read,prompt="New justification : ",ccol;
	  if(numberof(ccol)) pledit,vindx(1),justify=ccol;
	}
      if(chx==7)
	{
	  pledit,vindx(1),opaque=(1-(*par(3))(5));
	  continue;
	}
      if(chx==8)
	{
	  ccol=0n;
	  read,prompt="New height : ",ccol;
	  if(numberof(ccol)) pledit,vindx(1),height=ccol;
	  continue;
	}
      if(chx==9)
	{
	  cx=(*par(4))(2);
	  cy=(*par(4))(3);
	  nx=cx;
	  ny=cy;
	  ccol="";
	  read,prompt="Use the mouse (Y/N) : ",ccol;
	  if(ccol=="Y"||ccol=="y")
	    {
              while((r=mouse(-1,0," ## Click at the new position ##"))(-1)!=2)
                {
                  convertsys,r(1),r(2),int(r(-2)),nx,ny,vindx(2);
                  r(1)=nx;
                  r(2)=ny;
                  nx=r(1)-cx;
                  ny=r(2)-cy;
                  cx=r(1);
                  cy=r(2);
                  pledit,vindx(1),dx=nx,dy=ny;
                }
	    }
	  else
	    {
	      read,prompt="New x position : ",format="%f",nx;
	      read,prompt="New y position : ",format="%f",ny;
              nx-=cx;
              ny-=cy;
              pledit,vindx(1),dx=nx,dy=ny;
	    }
	  continue;
	}
    }
  plsys,cosysbck;
}


func textselect(indx)
/* DOCUMENT  textselect(indx)

     Selects a text in the list of texts.
     indx is an array of index of text element (can
     be the result of mousefindt or textfindt)

   KEYWORDS: cosys=  selects coordinate system.
                    default: -1

   SEE ALSO: putt, movett, textinfo, textedit,
             mousefindt, textfindt, put_indx, savet.
*/
{
  totindx=dimsof(indx)(3);
  if (totindx==0) exit,"!! Not text to select !!";
  if (totindx==1) return indx(,1);
  
  
  cosysbck=plsys();
  curindx=1;
  notend=1;
  while(notend)
    {
      plsys,indx(2,curindx);
      par=plq(indx(1,curindx));
      colorbck=(*par(3))(1);
      colornew=-5-2*(colorbck==-5);
      pledit,indx(1,curindx),color=colornew;
      r=mouse(-1,0,"Buttons: | Previous | Select | Next |");
      pledit,indx(1,curindx),color=colorbck;
      if (r(-1)==1) curindx=curindx%totindx+1;
      if (r(-1)==2) notend=0;
      if (r(-1)==3) curindx=(curindx-2+totindx)%totindx+1;
    }
  plsys,cosysbck;
  return indx(,curindx);
}

func mousefindt(dummy,cosys=)
/* DOCUMENT mousefindt(&cosys)

     returns index of the text nearest the click
     of the mouse.

   SEE ALSO: putt, movett, textinfo, textedit, textselect,
             textfindt, put_indx, savet.
*/
{
  if(is_void(cosys)) cosys=-1;

  r=mouse(cosys,0," ## Click at the center of the text ##");
  if(cosys==-1) cosys=int(r(-2));

  cosysbck=plsys(cosys);
  nel=numberof(plq());
  if(nel==0) exit,"!! Nothing is drawn !!";
  
  posms=r(1:2);
  indxt=[];
  distt=[];
  nnt=0;
  rl=limits();
  rfact=[abs(rl(2)-rl(1)),abs(rl(4)-rl(3))];
  for(i=1;i<=nel;i++)
    {
      par=plq(i);
      if ((*par(1))(1)!=3||(*par(1))(2)!=0) continue;
      nnt++;
      postxt=[(*par(4))(2),(*par(4))(3)];
      grow,distt,(((posms-postxt)/rfact)^2)(sum);
      grow,indxt,[i,cosys];
    }
  plsys,cosysbck;
  indxt=reform(indxt,[2,2,numberof(indxt)/2]);
  if(nnt==0) exit,"!! No text drawn !!";
  return indxt(,where(distt==distt(min)));
}

func textfindt(txt,cosys=,cases=)
/* DOCUMENT textfindt(txt)

     returns index of the text(s) which value match
     with 'txt'.

   KEYWORDS: cases=   1 case sensitive
                      0 not case sensitive (default)

   SEE ALSO: putt, movett, textinfo, textedit, textselect,
             mousefindt, put_indx, savet.
*/
{
  local txtp,txtadr,par,nel,cosysbck,str2find,indx,txts;
  local sys_d,sys_f;
  
  if(is_void(cases)) cases=0;
  if(is_void(cosys)||cosys<0)
    {
      sys_d=1;
      sys_f=GetNbSys();
    }
  else
    {
      sys_d=cosys;
      sys_f=cosys;
    }
      
  if(typeof(txt)!="string") exit,"!! Argument must be a string !!";
  if(!cases) str2find=strtolower(txt);

  cosysbck=plsys();
  nelt=0;
  indx=[];
  for(sys=sys_d;sys<=sys_f;sys++)
    {
      plsys,sys;
      nelt+=(nel=numberof(plq()));

      for(i=1;i<=nel;i++)
        {
          par=plq(i);
          if ((*par(1))(1)!=3) continue;
          txtadr=(*par(5))(2);
          reshape,txtp,txtadr,char,(*par(5))(1);
          txts=string(&txtp);
          if(!cases) txts=strtolower(txts);
          if(txts==str2find) grow,indx,[i,sys];
        }
    }
  plsys,cosysbck;
  if(nelt==0) exit,"!! Nothing is drawn !!";
  if (!numberof(indx)) exit,"!! Text '"+txt+"' not found !!";
  indx=reform(indx,[2,2,numberof(indx)/2]);
  return indx;
}

func put_indx(dummy,cosys=)
/* DOCUMENT put_indx()

   returns index of texts in the graphic.
   The result is an array of dimension
   [2,2,nt] where nt is the number of texte
   drawn.
   result(2,) give the system where the text
              is plotted.
   result(1,) give index in plq command
              (negative value for hidden text)
   
   KEYWORDS: cosys=  selects coordinate system.
                    default: 1
     
   SEE ALSO: putt, movett, textinfo, textedit, textselect,
             mousefindt, textfindt, savet.
 */
{
  local sysbck,nel,nelt;
  local indxt,par;
  local sys,sys_d,sys_f;
  
  if(is_void(cosys))
    {
      sys_d=1;
      sys_f=GetNbSys();
    }
  else
      sys_d=sys_f=cosys;

  sysbck=plsys();
  nelt=0;
  indxt=[];
  for(sys=sys_d;sys<=sys_f;sys++)
    {
      plsys,sys;
      nelt+=(nel=numberof(plq()));
      for(i=1;i<=nel;i++)
        {
          par=plq(i);
          if ((*par(1))(1)!=3) continue;
          if ((*par(1))(2)==0) grow,indxt,[i,sys];
          else grow,indxt,[-i,sys];
        }
    }
  plsys,sysbck;
  if(nelt==0) exit,"!! Nothing is drawn !!";
  if (!numberof(indxt)) exit,"!! No Text drawn !!";
  
  indxt=reform(indxt,[2,2,numberof(indxt)/2]);
  return indxt;
}

func textarray(dummy,cosys=)
{
  local i,idx,rslt;
  idx=abs(put_indx(cosys=cosys));
  rslt=array(string,numberof(idx)/2);
  for(i=1;i<=numberof(idx)/2;i++)
    {
      plsys,idx(2,i);
      par=plq(idx(1,i));
      par5=(*par(5));
      adrtxt=par5(2);
      reshape,txtp,adrtxt,char,[1,par5(1)];
      rslt(i)=string(&txtp);
    }
  return rslt;
}

func tposxarray(dummy,cosys=)
{
  local i,idx,rslt;
  idx=abs(put_indx(cosys=cosys));
  rslt=array(0.,numberof(idx)/2);
  for(i=1;i<=numberof(idx)/2;i++)
    {
      plsys,idx(2,i);
      par=plq(idx(1,i));
      rslt(i)=(*par(4))(2);
    }
  return rslt;
}

func tposyarray(dummy,cosys=)
{
  local i,idx,rslt;
  idx=abs(put_indx(cosys=cosys));
  rslt=array(0.,numberof(idx)/2);
  for(i=1;i<=numberof(idx)/2;i++)
    {
      plsys,idx(2,i);
      par=plq(idx(1,i));
      rslt(i)=(*par(4))(3);
    }
  return rslt;
}

func sysarray(dummy,cosys=)
{
  return put_indx(cosys=cosys)(2,);
}

func hidearray(dummy,cosys=)
{
  return put_indx(cosys=cosys)(1,)<0;
}

    

func savet(cmd_nb,cosys=,form=,file=)
/* DOCUMENT savet()
         or savet,[0,sys]
         or savet,indx

         Make a string of plt commands.
         First  form, for all plt commands.
         Second form, for the last plt command.
         Third  form, for the listed plt commands.

   KEYWORDS: cosys=   selects coordinate system.
                      default: -1

             form=0/1 without or with a newline after
                      each commands.

             file="file.i" save commands in file "file.i"
     
   SEE ALSO: putt, movett, textinfo, textedit, textselect,
             mousefindt, textfindt, put_indx.
 */
{
  local par,nel;
  local cmd,indxt;
  local cmds,txtp,txts;
  local sys,sys_d,sys_f,old_sys;
  
  if(is_void(form)) form=1;
  indxt=abs(put_indx(cosys=cosys));
  nel=numberof(indxt);
  if(nel==0) error,"No text is drawn !!";

  cmd=cmd_nb;
  if(is_void(cmd_nb)) cmd=indxt;
  if(cmd(1,1)==0) cmd=indxt(,0);
  
  tmp=where2((indxt(,,-)==cmd(,-,))(sum,,)==2)(1,);
  if(numberof(tmp)==0) error,"Bad value for text index !!";
  tmp=indxt(,(tmp-1)%nel+1);
  lst=sort(tmp(2,));
  tmp=tmp(,lst);
  nel=dimsof(tmp)(3);

  fmts="plsys,%i;";
  fmt ="plt,\"%s\",%f,%f,tosys=%i,justify=%i,legend=\"%s\"";
  fmt+=",hide=%i,color=%i,font=%i,height=%f,opaque=%i,orient=%i;";
  if(form) {fmt+="\n";fmts+="\n";}
  cmds="";
  sysbck=plsys();
  old_sys=tmp(2,1)-1;
  for(i=1;i<=nel;i++)
    {
      if(tmp(2,i)!=old_sys)
        {
          cmds+=swrite(format=fmts,tmp(2,i));
          plsys,tmp(2,i);
          old_sys=tmp(2,i);
        }
      par=plq(tmp(1,i));
      par1=(*par(1));par2=(*par(2));par3=(*par(3));par4=(*par(4));par5=(*par(5));
      adrtxt=par5(2);
      reshape,txtp,adrtxt,char,[1,par5(1)];
      txts=string(&txtp);
      hide    = par1(2);
      legend  = par2;
      color   = par3(1);
      font    = par3(2);
      orient  = par3(3);
      justify = par3(4);
      opaque  = par3(5);
      height  = par4(1);
      x       = par4(2);
      y       = par4(3);

      cmds+=swrite(format=fmt,txts,x,y,tmp(2,i),justify,legend,hide,color,font,height,opaque,orient);
    }
  if (file)
    {
      ff=open(file,"w");
      write,ff,cmds;
      close,ff;
      return;
    }
  return cmds;
}


func convertsys(xi,yi,sysin,&xo,&yo,sysout)
{
  local lmi,vwi,lmo,vwo;
  local sysbck;
  local xyt;

  sysbck=plsys(sysin);
  lmi=limits()(1:4);
  vwi=[0.,1.,0.,1.];
  if(sysin) vwi=viewport();

  plsys,sysout;
  lmo=limits()(1:4);
  vwo=[0.,1.,0.,1.];
  if (sysout) vwo=viewport();

  xyt=([xi,yi]-lmi(1::2))*(vwi(2::2)-vwi(1::2))/(lmi(2::2)-lmi(1::2))+vwi(1::2);
  xyt=(xyt-vwo(1::2))*(lmo(2::2)-lmo(1::2))/(vwo(2::2)-vwo(1::2))+lmo(1::2);
  xo=xyt(1);
  yo=xyt(2);
  return xyt;
}

func convertSysTo0(xi,yi,sysin,&xo,&yo)
{
  local lmi,vwi;
  local a,b,c,d;

  old_sys=plsys(sysin);
  vwi=viewport();
  lmi=limits();
  plsys,old_sys;
  xo=(xi-lmi(1))/(lmi(2)-lmi(1))*(vwi(2)-vwi(1))+vwi(1);
  yo=(yi-lmi(3))/(lmi(4)-lmi(3))*(vwi(4)-vwi(3))+vwi(3);
  return [xo,yo];
}
