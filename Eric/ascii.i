/*
 * ascii.i --
 *
 *	Reading ascii data from a file.
 *
 * $Id: ascii.i,v 1.3 1996/03/11 15:34:14 eric Exp $
 *
 * Copyright (c) 1996, Eric THIEBAUT (thiebaut@obs.univ-lyon1.fr, Centre de
 * Recherche Astrophysique de Lyon, 9 avenue Charles  Andre,  F-69561 Saint
 * Genis Laval Cedex).
 *
 * This program is free software; you can redistribute it and/or  modify it
 * under the terms of the GNU General Public License  as  published  by the
 * Free Software Foundation; either version 2 of the License,  or  (at your
 * option) any later version.
 *
 * This program is distributed in the hope  that  it  will  be  useful, but
 * WITHOUT  ANY   WARRANTY;   without   even   the   implied   warranty  of
 * MERCHANTABILITY or  FITNESS  FOR  A  PARTICULAR  PURPOSE.   See  the GNU
 * General Public License for more details (to receive a  copy  of  the GNU
 * General Public License, write to the Free Software Foundation, Inc., 675
 * Mass Ave, Cambridge, MA 02139, USA).
 *
 * ------------------------------------------------------------------------
 * History:
 *	02/20/96: release 1.1
 *	02/20/96 by Eric THIEBAUT: improve performances of asciiRead()
 *	    execution time almost divided by 2, "string.i" no more needed.
 *      02/21/96 by Dave Munro: improved performance by a factor of 250 or so.
 */

func asciiRead(file)
/* DOCUMENT data= asciiRead(name)
            data= asciiRead(file)
	read ascii numeric data in columns from file NAME, or the already
	open file FILE.
	The result is a NCOLUMNS-by-NLINES array of doubles.

	Data are read as double values arranged in columns separated
	by any number of spaces or tabs.  Comments starting with a "#"
	or any other character which is not part of a number are ignored
	up to the end-of-line.  Blank lines are ignored.
	The first non-blank/commented line gives the number of values per
	column, for subsequent lines.  Subsequent lines must have the
	same number of columns -- blanks in columns are not permitted,
	use 0.0 instead.  However, minimal error checking is performed,
	and if the data is not really in columns, asciiRead can silently
	fail to interpret your file as you would scanning it by eye.

	The read operation will be much faster if the number of commented
	lines is relatively small.  Blank lines cost nothing, while a line
	containing just a "#" is expensive.

   SEE ALSO: read
*/
{
  /* open the file if it's not already open */
  if (structof(file)==string) file= open(file);

  /* read lines one at a time until the "model" line which
   * determines the number of columns is discovered
   * assume the number of columns is less than 128 */
  x= array(0.0, 128);
  ncols= 0;
  while ((line= rdline(file))) {
    ncols= sread(line, x);
    if (ncols) break;          /* got a line with numbers */
  }
  if (!ncols) return [];

  nrows= 1;
  list= _lst([x(1:ncols)]);
  x= array(0.0, ncols, 10000/ncols + 1);
  for(;;) {
    /* try to grab at least 10000 numbers from the file
     * blank lines will be skipped, but any comments will
     * interrupt the read */
    n= read(file, x);
    if (!n) {
      /* if didn't get any, drop back to reading comments one
       * line at a time until we get some more numbers */
      while ((line= rdline(file))) {
	n= sread(line, x);
	if (n) break;
      }
      if (!line) break;    /* rdline detected end-of-file, n==0 too */
    }
    if (n%ncols) error, "data is not in columns";
    n/= ncols;

    /* grow the list the fast way, adding new values to its head
     * (adding to the tail would make growth an n^2 proposition,
     *  as would using the grow function) */
    list= _cat(x(,1:n), list);
    nrows+= n;
  }

  /* pop chunks off list and reassemble result */
  x= array(0.0, ncols, nrows);
  for (i=nrows ; list ; list=_cdr(list)) {
    n= numberof(_car(list))/ncols;
    x(,i-n+1:i)= _car(list);
    i-= n;
  }

  return x;
}
