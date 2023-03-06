func myxytitles(xtitle, ytitle, adjust,color=,invy=)
/* DOCUMENT xytitles, xtitle, ytitle
       -or- xytitles, xtitle, ytitle, [deltax,deltay]
     Plot XTITLE horizontally under the viewport and YTITLE vertically
     to the left of the viewport.  If the tick numbers interfere with
     the labels, you can specify the [DELTAX,DELTAY] in NDC units to
     displace the labels.  (Especially for the y title, the adjustment
     may depend on how many digits the numbers on your scale actually
     have.)  Note that DELTAX moves YTITLE and DELTAY moves XTITLE.
     WARNING: There is no easy way to ensure that this type of title
              will not interfere with the tick numbering.  Interference
              may make the numbers or the title or both illegible.
   SEE ALSO: plt, pltitle
 */
{
  invy=(is_void(invy)?0:!(invy==0));
  invf=1-2*invy;
  if (is_void(adjust)) adjust= [0.,0.];
  port= viewport();
  if (xtitle && strlen(xtitle))
    plt, xtitle, port(zcen:1:2)(1), port(3)-0.050+adjust(2),
      font=pltitle_font, justify="CT", height=pltitle_height,color=color;
  if (ytitle && strlen(ytitle))
    plt, ytitle, port(1+invy)-0.050*invf+adjust(1), port(zcen:3:4)(1),
      font=pltitle_font, justify="CB", height=pltitle_height, orient=1,color=color;
}

func plid(z,cmin=,cmax=)
{
  dms=dimsof(z);
  pli,z,0.5,0.5,dms(2)+0.5,dms(3)+0.5,cmin=cmin,cmax=cmax;
}
