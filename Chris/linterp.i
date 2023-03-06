//----------------------------------------------------------------------
// @(#) linterp.i: (bi)linear interpolation for Yorick.
//----------------------------------------------------------------------
// Copyright (c) 1995 Eric THIEBAUT.
//	$Id: linterp.i,v 1.2 1995/12/18 11:32:06 eric Exp $
//----------------------------------------------------------------------

func LInterp(a, xp, yp, zp, x=, y=, z=, dbl=)
/* DOCUMENT LInterp -- compute (bi/tri)linear interpolation of an array.
|
| SYNOPSIS:
|   LInterp(a, xp [, yp ][, zp ] [, x=, y=,z=, dbl= ] )
|
| ARGUMENTS:
|   A   = Array to interpolate (1-D, 2-D or 3-D array).
|   XP  = ``X prime'': abscissa where the function is to be interpolated
|           (must be specified for 1-D and 2-D and 3-D input array).
|   YP  = ``Y prime'': ordinate where the function is to be interpolated
|           (must be specified for 2-D & 3-D input array only).
|   ZP  = ``Z prime'': alt-ordinate where the function is to be interpolated
|           (must be specified for 3-D input array only).
|   X   = Abscissa where the function is sampled (optional, default
|           is to use the column index).  The sampling abscissa must
|           be in increasing order and with a constant step size.
|   Y   = Ordinate where the function is sampled (optional, default
|           is to use the row index).    The sampling ordinate must
|           be in increasing order and with a constant step size.
|   Z   = altOrdinate where the function is sampled (optional, default
|           is to use the row index).    The sampling ordinate must
|           be in increasing order and with a constant step size.
|   DBL = Flag to force computation in double precision.
|
| DESCRIPTION:
|   If input array is 1-D, the ouput is the linear interpolation for
|   coordinate XP. E.g.:
|       XP = span(-1,40,numberof(A))
|       B = LInterp(A, XP)
|   Then B(i) is the linear interpolated approximation of A at XP(i).
|
|   If input array is 2-D, the ouput is the bilinear interpolation for
|   coordinate (XP, YP). E.g.:
|       XP = span(-1,40,100)(,-:1:50)
|       YP = span(1,40,50)(-:1:100,)
|       B  = LInterp(A, XP, YP)
|   Then B(i) is the bilinear interpolated approximation of A at (XP(i),YP(i)).
|
|   If input array is 3-D, the ouput is the bilinear interpolation for
|   coordinate (XP, YP). E.g.:
|       XP = span(-1,40,100)(,-:1:50)(,,-:1:40)
|       YP = span(1,40,50)(-:1:100,)(,,-:1:40)
|       ZP = span(1,40,40)(,-:1:40)(,,-:1:100)
|       B  = LInterp(A, XP, YP,ZP)
|   Then B(i) is the bilinear interpolated approximation of A at
|    (XP(i),YP(i),ZP(i)).
|
| NOTES:
|   XP,YP and ZP may have any geometry (but the geometry of XP,YP and ZP must be
|   the same); the output geometry will be that geometry.
|
|   If XP (or YP or ZP) runs outside the sampling limits, the extrapolated value
|   of A is computed assuming A remains constant from the boundary of its
|   support to the infinity in a direction parallel to one of the axis (see
|   LInterpC(eric)).
|
|   The sampling coordinates (if specified) *MUST* be in increasing order
|   and with a constant step size.
|
|
| SEE ALSO:
|    LInterpC(eric).
|
| HISTORY:
|    July 23, 1995 by Eric THIEBAUT.
|    September 8, 1995 by Eric THIEBAUT: avoid reference to non-local variables
|                                        (i.e., `i', `j', `u', `v').
*/
{
  local nx, ny, nz, i, j, k, u, v, w;

  d = dimsof(a);
  n = d(1);
  coords = is_void(xp) + 2*is_void(yp) + 4*is_void(zp);
  if (n == 1) {
    nx = d(2);
    if (coords != 6)
      error, "XP and only XP must be specified for linear interpolation";
    if (is_void(x)) x = indgen(nx);
    else if (numberof(x) != nx || dimsof(x)(1) > 1) error, "bad size for X";
    LInterpC, i, u, xp, x, dbl=dbl;
    return (1 - u)*a(i) + u*a(i + 1);
  } else if (n == 2) {
    nx = d(2);
    ny = d(3);
    if (coords != 4)
      error, "XP and YP (only) must be specified for linear interpolation";
    if (is_void(x)) x = indgen(nx);
    else if (numberof(x) != nx || dimsof(x)(1) > 1) error, "bad size for X";
    if (is_void(y)) y = indgen(ny);
    else if (numberof(y) != ny || dimsof(y)(1) > 1) error, "bad size for Y";
    LInterpC, i, u, xp, x, dbl=dbl;
    LInterpC, j, v, yp, y, dbl=dbl;
    k = i + nx*(j - 1);
    j = [];
    w = 1 - v;
    return (a(k)*w + a(k+nx)*v)*(1 - u) + (a(k+1)*w + a(k+(nx+1))*v)*u;
  } else if (n == 3) {
    nx = d(2);
    ny = d(3);
    nz = d(3);
    if (coords != 0)
      error, "XP, YP and ZP must be specified for linear interpolation";
    if (is_void(x)) x = indgen(nx);
    else if (numberof(x) != nx || dimsof(x)(1) > 1) error, "bad size for X";
    if (is_void(y)) y = indgen(ny);
    else if (numberof(y) != ny || dimsof(y)(1) > 1) error, "bad size for Y";
    if (is_void(z)) z = indgen(nz);
    else if (numberof(z) != nz || dimsof(z)(1) > 1) error, "bad size for Z";
    LInterpC, i, u, xp, x, dbl=dbl;
    LInterpC, j, v, yp, y, dbl=dbl;
    LInterpC, k, w, zp, z, dbl=dbl;
#if 0
    /* 3D interpolation using 3D indices (not supposed to be working...): */
    return ((1-u)*(1-v)*(1-w)*a(i,j,k) +
            u*(1-v)*(1-w)*a(i+1,j,k) +
            (1-u)*v*(1-w)*a(i,j+1,k) +
            (1-u)*(1-v)*w*a(i,j,k+1) +
            (1-u)*v*w*a(i,j+1,k+1) +
            u*(1-v)*w*a(i+1,j,k+1) +
            u*v*(1-w)*a(i+1,j+1,k) +
            u*qv*qw*a(i+1,j+1,k+1));
    /* 3D interpolation using 1D index (with parenthesis to speed up
       computation of indices): */
    i += nx*(j - 1 + ny*(k - 1));
    j = k = [];
    return ((1-u)*(1-v)*(1-w)*a(i) +
            u*(1-v)*(1-w)*a(i+1) +
            (1-u)*v*(1-w)*a(i+nx) +
            (1-u)*(1-v)*w*a(i+nx*ny) +
            (1-u)*v*w*a(i+nx*(1+ny)) +
            u*(1-v)*w*a(i+(1+nx*ny)) +
            u*v*(1-w)*a(i+(1+nx)) +
            u*v*w*a(i+(1+nx*(1+ny))));
#endif
    /* 3D interpolation using 1D index and factorization: */
    i += nx*(j - 1 + ny*(k - 1));
    j = k = [];
    vp = 1 - v;
    wp = 1 - w;
    return ((1-u)*(vp*(wp*a(i)        + w*a(i+nx*ny)) +
                   v *(wp*a(i+nx)     + w*a(i+nx*(1+ny)))) +
            u    *(vp*(wp*a(i+1)      + w*a(i+(1+nx*ny))) +
                   v *(wp*a(i+(1+nx)) + w*a(i+(1+nx*(1+ny))))));
  } else {
    error, "array A must be 1-D or 2-D or 3-D";
  }
}

//----------------------------------------------------------------------

func LInterpC(&i, &u, xp, x, dbl=)
/* DOCUMENT LInterpC -- compute index and weights for linear interpolation.
|
| SYNOPSIS:
|   LInterpC, &i, &u, xp, x [, dbl=]
|
| INPUT:
|   XP  = ``X prime'': coordinates where the function will be evaluated.
|   X   = coordinates where the function is sampled (must be a vector of,
|           at least 2 elements).
|   DBL = flag to force computation in double precision.
|
| OUTPUT:
|   I  = index of lower bounds, INDEX is guaranted to be in the range
|            0, ..., N_Elements(X)-2L
|   U  = weights for the interpolation:
|
| DESCRIPTION:
|   The linearly interpolated value `Fi' of `F' at `xnew' can be computed by:
|       LInterpC, i, u, xp, x
|       Fi = (1-u) * F(i) + u * F(i+1)
|
| WARNING:
|   It is assumed that:
|    1. The sampling coordinate is in increasing order and with a constant
|       step size.
|    2. If XP is less than min(X), than the interpolated value will be
|       F(1).
|    3. If XP is greater than max(X), than the interpolated value will be
|       F(numberof(X)).
|       
| SEE ALSO:
|    LInterp(eric), digitize(Yorick).
|
| HISTORY:
|    July 20, 1995 by Eric THIEBAUT.
*/
{
  if (numberof(x) < 2 || dimsof(x)(1) != 1) {
    error, "X must be a vector of, at least, 2 elements";
  }
  if (numberof(xp) < 2) error, "Bad size for XP";
  xmin = x(1);
  xmax = x(0);
  floating_point = (dbl || structof(xp) == double ? double : float);
  u = array(floating_point, dimsof(xp));
  step = (xmax - xmin)/floating_point(numberof(x)-1);
  if (step <= 0) error, "Elements of X must be in strictly increasing order";
  i = 1 + long(floor((xp - xmin)/step));
  if (is_array((w = where((xp > xmin)&(xp < xmax))))) {
    u(w) = (xp(w) - x(i(w)))/step;
  }
  if (is_array((w = where(xp < xmin)))) {
    u(w) = 0;
    i(w) = 1;
  }
  if (is_array((w = where(xp >= xmax)))) {
    u(w) = 1;
    i(w) = numberof(x) - 1;
  }
}
