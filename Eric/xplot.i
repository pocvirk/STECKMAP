/*
 * xplot.i --
 *
 *	Extended plot routines for Yorick.
 *
 * Copyright (c) 1996-2003, Eric THIEBAUT.
 *
 * History:
 *	$Id$
 *	$Log$
 */

/*---------------------------------------------------------------------------*/
/* MOUSE ROUTINES */

func xmouse_point(nil,win=,prompt=,system=,forever=)
/* DOCUMENT xmouse_point()
       -or- xmouse_point
     Interactively  define  a  point  as with  xmouse("point").   The  same
     keywords as xmouse (which see) can be specified.
     
   SEE  ALSO: xmouse. */
{
  m = __xmouse(0, "click mouse to choose a position");
  if (am_subroutine()) {
    if (is_void(m)) write, "aborted";
    else write, format="position is (%g,%g)\n", m(1), m(2);
  } else if (! is_void(m)) return m(1:2);
}

func xmouse_box(nil,win=,prompt=,system=,forever=)
/* DOCUMENT xmouse_box()
       -or- xmouse_box
     Interactively  define a  rectangular box  as with  xmouse("box").  The
     same keywords as xmouse (which see) can be specified.
     
   SEE  ALSO: xmouse. */
{
  m = __xmouse(1, "click and drag mouse to select a region");
  if (am_subroutine()) {
    if (is_void(m)) write, "aborted";
    else write, format="region is [%g : %g]x[%g : %g]\n",
           min(m(1), m(3)), max(m(1), m(3)),
           min(m(2), m(4)), max(m(2), m(4));
  } else if (! is_void(m)) return [min(m(1), m(3)), max(m(1), m(3)),
                                   min(m(2), m(4)), max(m(2), m(4))];
}

func xmouse_line(nil,win=,prompt=,system=,forever=)
/* DOCUMENT xmouse_line()
       -or- xmouse_line
     Interactively define a line as with xmouse("line").  The same keywords
     as xmouse (which see) can be specified.
     
   SEE  ALSO: xmouse. */
{
  m = __xmouse(2, "click and drag mouse to choose a line");
  if (am_subroutine()) {
    if (is_void(m)) write, "aborted";
    else write, format="line endpoints are: (%g,%g) and (%g,%g)\n",
           m(1), m(2), m(3), m(4);
  } else if (! is_void(m)) return m(1:4);
}

func xmouse_length(nil,win=,prompt=,system=,forever=)
/* DOCUMENT xmouse_length()
       -or- xmouse_length
     Interactively  measure a  length as  with xmouse("length").   The same
     keywords as xmouse (which see) can be specified.
     
   SEE  ALSO: xmouse. */
{
  m = __xmouse(2, "click and drag mouse to measure a distance");
  if (am_subroutine()) {
    if (is_void(m)) write, "aborted";
    else write, format="distance from (%g,%g) to (%g,%g) = %g\n",
           m(1), m(2), m(3), m(4), abs(m(1) - m(3), m(2) - m(4));
  } else if (! is_void(m)) return abs(m(1) - m(3), m(2) - m(4));
}

func xmouse(type,win=,prompt=,system=,forever=)
/* DOCUMENT xmouse()
       -or- xmouse(type)
       -or- xmouse
       -or- xmouse, type
     Asks the  user to click into  a graphic window to  indicate a position
     (the default  if TYPE is not specified),  or to define a  segment or a
     rectangular box,  or to measure  a distance.  The possible  values for
     TYPE are:

         TYPE          RESULT
       ---------   ----------------
       0 "point"   [X,Y]
       1 "box"     [XMIN, XMAX, YMIN, YMAX]
       2 "line"    [X0, Y0, X1, Y1]
       3 "length"  LENGHT

     When called as a function, the coordinates of the result get returned;
     when called as a subroutine, the position is printed out.  If the user
     cancel the operation (see mouse function) or click into another window
     than the target one, no position is selected (nil get returned).  This
     behaviour can be  changed by setting keyword FOREVER  to true in which
     case the function loop until a valid selection is done.

     Keyword WIN can be used to  specify the target graphic window which is
     by default the curent graphic window.

     Keyword PROMPT can be used to specify a prompt string different from
     the default one.

     Keyword SYSTEM  can be used to  specify the coordinate  system to use,
     the  default being  to use  the coordinate  system that  is  under the
     mouse.   The   returned  coordinates/lenght  are  in   units  of  that
     coordinate system.

   SEE ALSO: mouse, xmouse_point, xmouse_box, xmouse_line, xmouse_length. */
{
  if (is_void(type)) {
    op = xmouse_point;
  } else if (is_array(type) && ! dimsof(type)(1)) {
    if ((s = structof(type)) == string) {
      op = symbol_def("xmouse_"+type);
      if (is_func(op) != 1) error, "unrecognized type name: \""+type+"\"";
    } else if (s == long || s == int || s == char || s == short) {
      /**/ if (type == 0) op = xmouse_point;
      else if (type == 1) op = xmouse_box;
      else if (type == 2) op = xmouse_line;
      else if (type == 3) op = xmouse_length;
      else error, "bad mouse selection type";
    } else {
      error, "bad mouse selection type must be a scalar integer of string";
    }
  } else {
    error, "bad mouse selection type";
  }
  if (am_subroutine()) op, win=win, system=system, prompt=prompt, forever=forever;
  else return op(win=win, system=system, prompt=prompt, forever=forever);
}

func __xmouse(style, default_prompt)
/* DOCUMENT __xmouse()
     Private function used by xmouse_* routines.
     
   SEE  ALSO: xmouse. */
{
  extern win, prompt, sustem, forever;
  if (! is_void(win)) window, win;
  if (is_void(system)) system = -1;
  if (is_void(prompt)) prompt = default_prompt;
  for (;;) {
    m = mouse(system, style, prompt);
    if (! forever || (! is_void(m) && m(10))) return m;
  }
}

func xmouse_test
{
  t=xmouse_point();
  plp,t(2),t(1),symbol=8,color="red",fill=1;
  
  t=xmouse_line();
  pldj,t(1),t(2),t(3),t(4),color="green";
  
  t=xmouse_box();
  plbox,t,color="green";
  
  xmouse_length;
}

/*---------------------------------------------------------------------------*/

func copy_limits(src=, dst=)
/* DOCUMENT copy_limits, src=..., dst=...;
     Copy plot limits from source graphic window SRC to destination graphic
     window DST.  SRC and DST are  both specified by keyword and default to
     the current graphic window if not set.

   SEE ALSO: limits, window. */
{
  if ((cur = current_window()) < 0) return; /* no current window */
  if (is_void(dst)) dst = cur;
  if (is_void(src)) src = cur;
  if (dst == src) return;
  window, src;
  lim = limits();
  window, dst;
  limits, lim;
  window, cur;
}

/*---------------------------------------------------------------------------*/

func plbox(xmin,xmax,ymin,ymax,color=,width=,type=)
/* DOCUMENT plbox, xmin, xmax, ymin, ymax;
       -or- plbox, [xmin, xmax, ymin, ymax];
     Plots a rectangular  box onto the current graphic  device.  As for plg
     (which see), keywords COLOR, WIDTH and TYPE can be specified.
     
   SEE ALSO: plg, pl_get_color. */
{
  if (is_void(xmax) && numberof(xmin)==4) {
    ymax = xmin(4);    
    ymin = xmin(3);
    xmax = xmin(2);
    xmin = xmin(1);
  }
  plg, [ymin,ymin,ymax,ymax], [xmin,xmax,xmax,xmin], closed=1,
    color=pl_get_color(color), width=width, type=type;
}

#if 0

func xplt(text, win=,prompt=,system=,legend=,
          type=,width=,marker=,marks=,
          justify=,color=,font=,height=,orient=,opaque=,
          anchor=)
/* DOCUMENT

     Keyword SYSTEM can be used to specify the coordinate system for the
     anchor point.  If SYSTEM>=0, that coordinate system is used; if
     SYSTEM<0, the coordinate system is the one under the mouse when the
     button was pressed.  By default, SYSTEM=0 (NDC coordinate system).
*/
{
  if (! is_void(win)) window, win;
  if (is_void(system)) system = 0;
  if (is_void(prompt)) prompt = "Click with mouse where you want to anchor the text."
  m = mouse(system, 0, prompt);
  if (is_void(m)) return;
  x = m(1);
  y = m(2);
  system = long(m(9));

  /* set coordinate system for the subsequent plots and save previous one */
  previous_system = plsys(system);

  /* Text justification is "HV" where:
      H (horizontal) is one of: N (normal = left), L (left), C (center),
                                or R (right);
      V (vertical) is one of: N (normal = baseline), T (top), C (capline),
                              H (half), A (baseline), or B (bottom).
  */
  if (is_void(anchor) || anchor == "w" || anchor == "west" || anchor == "left") {
    justify = "LH";
  } else if (anchor == "c" || anchor == "center") {
    justify = "CH";
  } else if (anchor == "n" || anchor == "north" || anchor == "top") {
    justify = "CT";
  } else if (anchor == "e" || anchor == "east" || anchor == "right") {
    justify = "RH";
  } else if (anchor == "s" || anchor == "south" || anchor == "bottom") {
    justify = "CB";
  } else if (anchor == "nw") {
    justify = "LT";
  } else if (anchor == "ne") {
    justify = "RT";
  } else if (anchor == "sw") {
    justify = "LB";
  } else if (anchor == "se") {
    justify = "RB";
  } else {
    error, "bad value for keyword ANCHOR";
  }

  /* plt keywords: legend, hide, color, font, height, opaque, orient, justify */
  plt, text, x, y, tosys=1, justify=justify,legend=legend,
    color=color, font=font, height=height, orient=orient, opaque=opaque;

  /* plg keywords: legend, hide,
              type, width, color, closed, smooth,
              marks, marker, mspace, mphase,
              rays, arrowl, arroww, rspace, rphase. */
  //plg, [y0, y0], [x0, x1], color=color, width=width, marker=marker, marks=marks;

  /* restore coordinate system */
  plsys, previous_system;
}

#endif

/*---------------------------------------------------------------------------*/

func plp(y, x, dx=, xlo=, xhi=, dy=, ylo=, yhi=, size=, symbol=, ticks=,
	 legend=, type=, width=, color=, fill=)
/* DOCUMENT plp, y, x;
       -or- plp, y, x, dx=sigma_x, dy=sigma_y;
     Plots points (X,Y) with symbols and/or  error bars.  X, and Y may have
     any dimensionality, but  must have the same number  of elements.  If X
     is nil, it defaults to indgen(numberof(Y)).

     Keyword SYMBOL may be used to choose the shape of each symbol (see
     pl_get_symbol).
     
     Keyword  SIZE may be used to  change the size of the  symbols and tick
     marks (SIZE acts as a multiplier, default value is 1.0).

     If value of  keyword FILL is true (non-nil  and non-zero), symbols are
     filled with COLOR (default is to draw open symbols).

     Keywords XLO, XHI, YLO, and/or YHI  can be used to indicate the bounds
     of the optional  error bars (default is to draw  no error bars).  Only
     specified bounds get plotted as  error bars. If value of keyword TICKS
     is true  (non-nil and non-zero), ticks  get drawn at  the endpoints of
     the error bars.   Alternatively, keywords DX and/or DY  can be used to
     plot  error bars  as segments  from XLO=X-DX  to XHI=X+DX  and/or from
     YLO=Y-DY to  YHI=Y+DY.  If keyword  DX (respectively DY) is  used, any
     value of XLO and XHI (respectively YLO and YHI) is ignored.

     The other keywords are the same as for pldj (TYPE is only used to draw
     error bars):
   KEYWORDS: legend, type, width, color.

   SEE ALSO: pl_get_symbol, pl_get_color,
             pldj, plg, plm, plc, plv, plf, pli, plt, pldj, plfp, plmk,
             limits, logxy, range, fma, hcp. */
{
  /* NDC units for symbols/ticks (one pixel = 0.00125268 NDC at 75 DPI) */
  u0 = 0.0;       // zero
  u1 = 0.0100214; // radius of about 8 pixels at 75 DPI
  if (! is_void(size)) u1 *= size;

  /* parse color */
  color = pl_get_color(color);

  /* default X */ 
  if (is_void(x)) (x = array(double, dimsof(y)))(*) = indgen(numberof(y));
  
  /* error bars */
  if (is_void(dx)) {
    err = (! is_void(xlo)) + 2*(! is_void(xhi));
  } else {
    xlo = x - dx;
    xhi = x + dx;
    err = 3;
  }
  if (err) {
    pldj, (is_void(xlo) ? x : xlo), y, (is_void(xhi) ? x : xhi), y,
      type=type, width=width, color=color;
    if (ticks) {
      xm = [ u0, u0];
      ym = [-u1, u1];
      if      (err == 1) __plp,   y,       xlo;
      else if (err == 2) __plp,   y,       xhi;
      else               __plp, [y, y], [xlo, xhi];
    }
    xhi = xlo = [];
  }
  if (is_void(dy)) {
    err = (! is_void(ylo)) + 2*(! is_void(yhi));
  } else {
    ylo = y - dy;
    yhi = y + dy;
    err = 3;
  }
  if (err) {
    pldj, x, (is_void(ylo) ? y : ylo), x, (is_void(yhi) ? y : yhi),
      type=type, width=width, color=color;
    if (ticks) {
      xm = [-u1, u1];
      ym = [ u0, u0];
      if      (err == 1) __plp,    ylo,       x;
      else if (err == 2) __plp,    yhi,       x;
      else               __plp, [ylo, yhi], [x, x];
    }
    yhi = ylo = [];
  }

  /* symbols */
  symbol = pl_get_symbol(symbol);
  if (! symbol) return;
  if (symbol == 1) {
    /* square */
    u2 = u1*sqrt(0.5);
    xm = [-u2, u2, u2,-u2];
    ym = [ u2, u2,-u2,-u2];
  } else if (symbol == 2) {
    /* + cross */
    xm = [-u1, u1, u0, u0, u0, u0];
    ym = [ u0, u0, u0, u1,-u1, u0];
    fill = 0;
  } else if (symbol == 3) {
    /* triangle */
    u2 = u1*0.5;
    u3 = u1*sqrt(0.75);
    xm = [u0, u3,-u3];
    ym = [u1,-u2,-u2];
  } else if (symbol == 4) {
    /* hexagon */
    u2 = u1*0.5;
    u3 = u1*sqrt(0.75);
    xm = [ u1, u2,-u2,-u1,-u2, u2];
    ym = [ u0, u3, u3, u0,-u3,-u3];
  } else if (symbol == 5) {
    /* diamond */
    xm = [u1, u0,-u1, u0];
    ym = [u0, u1, u0,-u1];
  } else if (symbol == 6) {
    /* x cross (rotated 45 degrees) */
    u2 = u1*sqrt(0.5);
    xm = [u2,-u2, u0, u2,-u2, u0];
    ym = [u2,-u2, u0,-u2, u2, u0];
    fill = 0;
  } else if (symbol == 7) {
    /* triangle (upside down) */
    u2 = u1*0.5;
    u3 = u1*sqrt(0.75);
    xm = [ u0, u3,-u3];
    ym = [-u1, u2, u2];
  } else if (symbol == 8) {
    /* 5 branch star
     *   C18 = cos(18*ONE_DEGREE)
     *   S18 = sin(18*ONE_DEGREE)
     *   C54 = cos(54*ONE_DEGREE)
     *   S54 = sin(54*ONE_DEGREE)
     */
    u2 = 0.224514*u1; // C54*S18/S54
    u3 = 0.309017*u1; // S18
    u4 = 0.951057*u1; // C18
    u5 = 0.363271*u1; // C18*S18/S54
    u6 = 0.118034*u1; // S18*S18/S54
    u7 = 0.587785*u1; // C54
    u8 = 0.809017*u1; // S54
    u9 = 0.381966*u1; // S18/S54
    xm = [ u0, u2, u4, u5, u7, u0,-u7,-u5,-u4,-u2];
    ym = [ u1, u3, u3,-u6,-u8,-u9,-u8,-u6, u3, u3];
  } else {
    /* N-side polygon in unit circle */
    PI = 3.141592653589793238462643383279503;
    a = (2.0*PI/symbol)*indgen(0:symbol-1);
    xm = u1*cos(a);
    ym = u1*sin(a);
  }
  __plp, y, x;
}

func __plp(y, x)
/* DOCUMENT __plp, x, y;
     Private routine used by plp. */
{
  extern xm, ym, color, fill, legend, width;
  local z;
  n = array(1, 1 + numberof(y));
  n(1) = numberof(ym);
  if (fill && n(1) > 2) {
    if (numberof(color) == 3) {
      z = array(char, 3, numberof(n));
      z(,) = color;
    } else {
      z = array(char(color), numberof(n));
    }
  }
  plfp, z, grow(ym,y(*)), grow(xm,x(*)), n,
    legend=legend, edges=1, ewidth=width, ecolor=color;
}

func pl_get_symbol(symbol)
/* DOCUMENT pl_get_symbol(symbol)
     Get symbol value as an integer, SYMBOL must be a scalar and may be
     either an integer, a character or a string:

       INT CHAR  STRING                 DESCRIPTION
       ----------------------------------------------------------------------
        0                               nothing (just draw error bars if any)
        1   #  "square"                 a square
        2   +  "plus"                   a plus sign
        3   ^  "triangle" "uptriangle"  a triangle
        4   o  "circle"                 a circle (actually an hexagon)
        5   @  "diamond"                a square rotated by 45 degrees 
        6   x  "cross"                  an X-cross    <- this is the default
        7   v  "downtriangle"           an upside down triangle
        8   *  "star"                   a star
       >=9                              a polygon with SYMBOL sides
       ----------------------------------------------------------------------

     The one-character symbol may given as loweer/upper case and as a
     string or a char; e.g. 'v', 'V', "v" and "V" all stand for an upside
     down triangle.

  SEE ALSO: plp, pl_get_color. */
{
  if (is_void(symbol)) return 6;
  if (! is_array(symbol) || dimsof(symbol)(1))
    error, "symbol must be a scalar";
  if ((s = structof(symbol)) == string) {
    len = strlen(symbol);
    if (len >= 1) {
      c = *pointer(symbol);
      if (len == 1) {
        if (c=='#') return 1;
        if (c=='+') return 2;
        if (c=='^') return 3;
        if (c=='o' || c =='O') return 4;
        if (c=='@') return 5;
        if (c=='x' || c=='X') return 6;
        if (c=='v' || c=='V') return 7;
        if (c=='*') return 8;
      } else {
        if (c=='s') {
          if (symbol=="square") return 1;
          if (symbol == "star") return 8;
        } else if (c=='p') {
          if (symbol=="plus") return 2;
        } else if (c=='t') {
          if (symbol == "triangle") return 3;
        } else if (c=='c') {
          if (symbol == "circle") return 4;
          if (symbol == "cross") return 6;
        } else if (c=='d') {
          if (symbol == "diamond") return 5;
          if (symbol == "downtriangle") return 7;
        } else if (c=='u') {
          if (symbol=="uptriangle") return 3;
        }
      }
    }
    error, "unrecognized symbol name: \""+symbol+"\"";
  }
  if (s == char) {
    if (c=='#') return 1;
    if (c=='+') return 2;
    if (c=='^') return 3;
    if (c=='o' || c =='O') return 4;
    if (c=='@') return 5;
    if (c=='x' || c=='X') return 6;
    if (c=='v' || c=='V') return 7;
    if (c=='*') return 8;
  } else if (s!=long && s!=int && s!=short) {
    error, "bad data type for symbol value";
  }
  if (symbol < 0) error, "bad symbol value";
  return long(symbol);
}

/*---------------------------------------------------------------------------*/
/* COLORS */

local __pl_color_list;
__pl_color_list = ["bg","fg","black","white","red","green","blue","cyan",
                   "magenta","yellow","grayd","grayc","grayb","graya",
                   "xor","extra"];
/* DOCUMENT __pl_color_list - private list of color names */

func pl_get_color(color)
/* DOCUMENT pl_get_color(color)
     Get color  value as an  integer or a  triplet [r,g,b]; COLOR can  be a
     scalar sting,  integer or triplet [r,g,b].  The  following color names
     and equivalent integer values are supported:
     
        NAME      INTEGER         NAME      INTEGER
       -------------------       -------------------
       "bg"       255   -1       "magenta"  247   -9
       "fg"       254   -2       "yellow"   246  -10
       "black"    253   -3       "grayd"    245  -11
       "white"    252   -4       "grayc"    244  -12
       "red"      251   -5       "grayb"    243  -13
       "green"    250   -6       "graya"    242  -14
       "blue"     249   -7       "xor"      241  -15
       "cyan"     248   -8       "extra"    240  -16
       -------------------       -------------------

     If COLOR is nil, the returned value is 254 ("fg" color).

  SEE ALSO: color, pl_get_symbol. */
{
  if (is_void(color)) return 254; /* fg */
  if ((s = structof(color)) == string) {
    if (dimsof(color)(1)) error, "color name must be a scalar";
    n = where(color == __pl_color_list);
    if (numberof(n)!=1) error, "unrecognized color name: \""+color+"\"";
    return 256 - n(1);
  } else if (s == long || s == int || s == char || s == short) {
    if (! (ndims = dimsof(color)(1)) ||
        ndims == 1 && numberof(color) == 3) return color;
    error, "color must be a scalar integer or a [r,g,b] triplet";
  }
  error, "bad data type for color value";
}

/*---------------------------------------------------------------------------*/
