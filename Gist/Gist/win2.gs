# Gist drawing style
# $Id$
# $Log$

# Two viewports on a portrait page.
# Legends: none

# See the header file gist.h for a more complete description of the
# meanings of the various keywords set in this file; they coincide with the
# values of the corresponding C data structure members defined in gist.h.
# Here is a brief description:


# Conversion table:
#          100 dpi     75 dpi    point      cm        inch      NDC
#         ---------  --------- --------  --------- ---------   ----------
# 100 dpi    1         0.75      0.7227   0.0254    0.01       0.00093951
#  75 dpi    1.3333    1         0.9636   0.033867  0.013333   0.001253
#   point    1.3837    1.03778   1        0.035146  0.013837   0.0013
#      cm   39.3701   29.5276   28.4528   1         0.393701   0.036989
#    inch  100        75        72.27     2.54      1          0.093951
#     NDC 1064.38    798.288   769.231   27.0354   10.6438     1

# 

# Units and coordinate systems:
#   viewport, tickOff, labelOff, tickLen, xOver, yOver, height (of text)
#   legends.(x,y,dx,dy)
#     Coordinates are in Gist's "NDC" coordinate system.  In this system,
#     0.0013 unit is 1.000 point, and there are 72.27 points per inch:
#       pts=    1      6      8     10     12     14     16     18     24
#       NDC= 0.0013 0.0078 0.0104 0.0130 0.0156 0.0182 0.0208 0.0234 0.0312
#     The lower left hand corner of the sheet of paper is at (0,0).
#     For landscape=0, the 8.5 inch edge of the paper is horizontal;
#     for landscape=1, the 11 inch edge of the paper is horizontal.
#   width
#     Line width is measured in relative units, with 1.0 being 1/2 point.

# Ticks flags (add together the ones you want):
#   0x001  Draw ticks on bottom or left edge of viewport
#   0x002  Draw ticks on top or right edge of viewport
#   0x004  Draw ticks centered on origin in middle of viewport
#   0x008  Ticks project inward into viewport
#   0x010  Ticks project outward away from viewport (0x18 for both)
#   0x020  Draw tick label numbers on bottom or left edge of viewport
#   0x040  Draw tick label numbers on top or right edge of viewport
#   0x080  Draw all grid lines down to gridLevel
#   0x100  Draw single grid line at origin

# Line types:
#   solid        1
#   dash         2
#   dot          3
#   dash-dot     4
#   dash-dot-dot 5

# Font numbers:
#   Courier    0x00
#   Times      0x04
#   Helvetica  0x08
#   Symbol     0x0c
#   Schoolbook 0x10
# Add 0x01 for bold, 0x02 for italic

#    The top-level keywords are:
# 
#      landscape   - 0 or 1 to get portrait or landscape mode
#      default     - to change the default coordinate system
# 		     parameters
#      system      - to specify a particular coordinate system
# 		     The first system will be system number 0, the
# 		     second, system number 1, and so on.  Usually,
# 		     you should specify a default= which includes all
# 		     of the properties shared by all your coordinate
# 		     systems, then each system= keyword should
# 		     specify only differences between the default
# 		     and the particular system.  See work2.gs.
#      legends     - to specify the layout of curve legends
#      clegends    - to specify the layout of contour legends
# 
# The keywords for a system (or default) are:
# 
#      legend      - 0 or a quoted string, this will be the name of
# 		     the coordinate system
#      viewport    - { xmin, xmax, ymin, ymax } in NDC
#      ticks       - to specify the layout of tick marks and their labels
# 
# The keywords for ticks are:
# 
#      horiz       - to specify the layout of the horizontal ticks
#      vert        - to specify the layout of the vertical ticks
#      frame       - 1 or 0 to draw or not draw a box around the viewport
#      frameStyle  - line style to use in drawing frame box
# 
# The keywords for horiz or vert are the same as the GaTickStyle data
# structure:
# 
#      nMajor      - (floating point) maximum number of major ticks
#      nMinor      - (floating point) maximum number of minor ticks
#      logAdjMajor - factor to increase nMajor by for log axes
#      logAdjMinor - factor to increase nMinor by for log axes
#      nDigits     - number of digits after decimal point to allow
# 		     before switching to oveflow format
#      gridLevel   - level in tick heirarchy at which to stop full
# 		     grid lines (level 0 is major ticks, 1 next smaller,
# 		     and so on)
#      flags       - are flags as in the GaTickStyle data structure
# 		     (see section 2G above, TICK_L etc.)
#      tickOff     - NDC offset of ticks from viewport edge (+ out - in)
#      labelOff    - NDC offset of labels from viewport edge (+ out - in)
#      tickLen     - { len0, len1, len2, len3, len4 }
# 		     lengths of the five levels of ticks in NDC
#      tickStyle   - to specify linestyle for ticks
#      gridStyle   - to specify linestyle for grid lines
#      textStyle   - to specify text style for tick labels
#      xOver       - NDC x coordinate of overflow text (after nDigits)
#      yOver       - NDC y coordinate of overflow text (after nDigits)
# 
# The frameStyle, tickStyle, and gridStyle keywords are:
# 
#      color       - -1 to -10, line color (-2 is foreground)
#      type        - 1 to 5, 1 for solid, others as in L_DASH, etc.
#      width       - relative line width, 1.0 is 0.5*ONE_POINT
# 
# The textStyle keywords are:
# 
#      color      - -1 to -10, line color (-2 is foreground)
#      font       - 0 to 19 (0 is Courier, 4 Times-Roman, 8 Helvetica,
# 		    16 New Century Schoolbok)
#      height     - font size in NDC (ONE_POINT=0.0013 NDC units)
#      path       - 0 for T_RIGHT, 3 for T_DOWN
#      alignH     - 0 to 3 to specify horizontal alignment
# 		    0 normal, 1 left, 2 center, 3 right
#      alignV     - 0 to 5 to specify vertical alignment
# 		    0 normal, 1 top, 2 cap, 3 half, 4 base, 5 bottom
#      opaque     - 0 or 1 for transparent or opaque text
#      prec, expand, spacing, upX, upY - (unused)
#
 
# This actually repeats the default values in gread.c

# Paper orientation:
#   Paper size is 8.5 x 11 inch = 0.8 x 1.0 NDC
#   For landscape=0, the 8.5 inch edge of the paper is horizontal;
#   for landscape=1, the 11 inch edge of the paper is horizontal.
landscape= 0

# This viewport allow for a bottom horizontal colorbar.
default = {
  legend= 0,

  ticks= {

    horiz= {
      nMajor= 7.5,  nMinor= 50.0,  logAdjMajor= 1.2,  logAdjMinor= 1.2,
      nDigits= 3,  gridLevel= 1,  flags= 0x000,
      tickOff= 0.0091,  labelOff= 0.0104,
      tickLen= { 0.0065, 0.0039, 0.0013, 0.0000, 0.0000 },
      tickStyle= { color= -2,  type= 1,  width= 1.0 },
      gridStyle= { color= -2,  type= 3,  width= 1.0 },
      textStyle= { color= -2,  font= 0x08,  prec= 2,  height= 0.0156,
        expand= 1.0,  spacing= 0.0,  upX= 0.0,  upY= 1.0,
        path= 0,  alignH= 0,  alignV= 0,  opaque= 0 },
      xOver= 0.395,  yOver= 0.370 },

    vert= {
      nMajor= 7.5,  nMinor= 50.0,  logAdjMajor= 1.2,  logAdjMinor= 1.2,
      nDigits= 4,  gridLevel= 1,  flags= 0x000,
      tickOff= 0.0091,  labelOff= 0.0117,
      tickLen= { 0.0065, 0.0039, 0.0013, 0.0000, 0.0000 },
      tickStyle= { color= -2,  type= 1,  width= 1.0 },
      gridStyle= { color= -2,  type= 3,  width= 1.0 },
      textStyle= { color= -2,  font= 0x08,  prec= 2,  height= 0.0156,
        expand= 1.0,  spacing= 0.0,  upX= 0.0,  upY= 1.0,
        path= 0,  alignH= 0,  alignV= 0,  opaque= 0 },
      xOver= 0.150,  yOver= 0.370 },

    frame= 1,
    frameStyle= { color= -2,  type= 1,  width= 1.0 }
  }
}

system= {
  legend= "Model vs. Data",
  viewport= { 0.1703, 0.6604, 0.5785, 0.8606 },
  ticks= {
    horiz= { flags= 0x04B }, vert= { flags= 0x02B }
    frameStyle= { color= -2,  type= 1,  width= 1.0 }
  }
}

system= {
  legend= "Residuals",
  viewport= { 0.1703, 0.6604, 0.4095, 0.5603 },
  ticks= {
    horiz= { flags= 0x02B }, vert= { flags= 0x12B }
    frameStyle= { color= -2,  type= 1,  width= 1.0 }
  }
}


# The legends and clegends keywords are:
# 
#      x          - NDC x coordinate for legend box
#      y          - NDC y coordinate for legend box
#      dx         - NDC dx to second legend box (none if dx and dy both 0.0)
#      dy         - NDC dy to second legend box
# 		    (Use dy==0.0, dx>0 to put the legends in 2 columns)
#      textStyle  - to specify textstyle for legends.  Courier works best.
#      nchars     - maximum number of characters per line
#      nlines     - maximum number of lines (use nlines==0 to get no legends)
#      nwrap      - maximum number of lines before truncating long legend
# 
# No legends (nlines= 0):
legends= { nlines= 0 }
#
clegends= { nlines= 0 }

