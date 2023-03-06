/*
 * regul.i --
 *
 *	Regularization and constraints for optimization in Yorick.
 *	Provides routines:
 *	  - regul:            regularization term.
 *	  - entropy:          regularization term based on entropy.
 *	  - negentropy:       opositive of entropy.
 *	  - ln:               approximation of natural logarithm function
 *	                      used in entropy to allow negative values.
 *	  - smooth:           smooth an array
 *	  - smoothness:       measure the smoothness af an array
 *	  - finiteDifference: compute finite difference for smoothness()
 *
 *
 * Copyright  (c) 1996,  Eric THIEBAUT (thiebaut@obs.univ-lyon1.fr,  Centre
 * de  Recherche  Astrophysique  de  Lyon,  9 avenue Charles Andre, F-69561
 * Saint Genis Laval Cedex).
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
 */

/*
 * History:
 *
 * 	$Id: regul.i,v 1.6 1996/05/16 08:17:31 eric Exp eric $
 *	$Log$
 *	1996/04/16 (Eric THIEBAUT): added negentropy function.
 *	1996/04/17 (Eric THIEBAUT): added smooth function.
 *	1996/04/17 (Eric THIEBAUT): removed positivity reparametrization
 *	  functions (now in "Eric/positivity").
 *	1996/05/16 (Eric THIEBAUT): added ln(), entropy() and regul()
 *	  functions, also redefine negentropy() to allow negative values.
 */

/*---------------------------------------------------------------------------*/
func regul(x, &dx, deriv=, level=, normalize=, zero=)
/* DOCUMENT q = regul(x)
            q = regul(x, dx, deriv=1)
     regularization penalizing function for inversion of ill-conditionned
     problems.

     If LEVEL is positive or nil, returns negentropy(X, ...); otherwise,
     return returns smoothness(X, level=-LEVEL).
     
     See smoothness() and entropy() for the meaning of keywords.
 */
{
  if (level && level < 0) {
    return smoothness(x, dx, deriv=deriv, normalize=normalize,
		      zero=zero, weight=weight, order=-level);
  }
  return negentropy(x, dx, deriv=deriv, normalize=normalize, level=level);
}

/*---------------------------------------------------------------------------*/
local _ln_eps;
/* DOCUMENT local _ln_eps;
     threshold for approximation of natural logarithm; should be
     small and *MUST* be > 0; see ln().
*/
_ln_eps = 1e-100;

/*---------------------------------------------------------------------------*/
func ln(x, deriv=)
/* DOCUMENT y= ln(x)
     return approximation of the natural logarithm of X.  The approximation
     consists in a 2nd order polynomial extrapolation of log() for values
     less than _ln_eps:
       ln(X) = log(X)                                     where X >= EPS
               log(EPS)-1.5 + 2.0*X/EPS - 0.5*(X/EPS)^2   elsewhere
     where EPS=_ln_eps.  The coefficients of the polynomial are such as
     to preserve the continuity of ln() as well as its 1st and 2nd
     derivatives.

     This approximation is intended to allow the computation of the entropy
     even when the input has negative values.  Although zero and negative
     values are still very penalizing.

     Keyword DERIV may be set to 1 or 2 to return the 1st or 2nd derivatives
     of ln() respectively.
     
   SEE ALSO: _ln_eps, entropy.
*/
{
  extern _ln_eps;
  if (!deriv) {
    if (numberof((i = where(x < _ln_eps)))) {
      y = array(double, dimsof(x));
      xi = x(i);
      y(i) = log(_ln_eps) - 1.5 +
	0.5 / _ln_eps / _ln_eps * (4.0 * _ln_eps - xi) * xi;
      if (numberof((i = where(x >= _ln_eps))))
	y(i) = log(x(i));
      return y;
    }
    return log(x);
  } else if (deriv == 1) {
    if (numberof((i = where(x < _ln_eps)))) {
      y = array(double, dimsof(x));
      y(i) = 2.0 / _ln_eps - x(i) / _ln_eps^2;
      if (numberof((i = where(x >= _ln_eps))))
	y(i) = 1.0 / x(i);
      return y;
    }
    return 1.0 / x;
  } else if (deriv == 2) {
    if (numberof((i = where(x < _ln_eps)))) {
      y = array(double, dimsof(x));
      y(i) = -1.0 / _ln_eps^2;
      if (numberof((i = where(x >= _ln_eps)))) {
	xi = x(i);
	y(i) = -1.0 / xi / xi;
      }
      return y;
    }
    return -1.0 / x / x;
  }
  error, "DERIV keyword must be one of nil, 0, 1 or 2";
}

/*---------------------------------------------------------------------------*/
func entropy(x, &dx, deriv=, normalize=, level=)
/* DOCUMENT s= entropy(x, level=level_of_smoothness)
            s= entropy(x, dx, deriv=1, level=level_of_smoothness)     
     return the entropy of array X:
       S = sum(X - PRIOR - X * (ln(X) - ln(PRIOR))
     where PRIOR is (uniform prior):
       PRIOR = 1 / OMEGA
     where OMEGA is the number of non-zero values in X, or (floating smooth
     prior):
       PRIOR = smooth(X, LEVEL)

     X may have zero and negative values (see ln()).

     If keyword NORMALIZE is non-nil and non-zero, return the entropy (and
     possibly gradient) of the normalized array:
       XP = X / sqrt(sum(X^2)).
     With NORMALIZE=1, the entropy is only sensitive to the morphology of X.

     Keyword LEVEL set the level of smoothness for the floating  prior: the
     floating prior is smooth(X, LEVEL).  If  LEVEL  is  zero  or  nil (the
     default) a uniform prior constraint is used.

     If  keyword  DERIV  is  non-nil  and  non-zero,  the  gradient  of the
     neg-entropy with respect to X is computed and stored into DX.

   SEE ALSO: negentropy, smooth, ln, optimBrent, optimConjGrad.
*/
{
  /*
   * Entropy is:
   *   S = sum(X - P - X * log(X / P))
   * which can be approximated by:
   *   S # sum(X - P - X * ln(X) + X * ln(P))
   * where ln() is log() for positive values but allows negative
   * values.  The floating prior is given by:
   *   P = R(,+) * X(+) = smooth(X)
   * The gradient is:
   *   dS_dX = 1 - R(sum,) - ln(X) + ln(P) - X * lnp(X)
   *           + R(+,) * (X * lnp(P))(+)
   * where lnp() is the derivative of ln().  Since R is symmetric:
   *   R(+,) * (X * lnp(P))(+) = smooth(X * lnp(P)).
   * and since  R(sum,) = 1, the gradient becomes:
   *   dS_dX = - ln(X) + ln(P) - X * lnp(X) + smooth(X * lnp(P))
   *
   * If the input array is normalized:
   *   S = S(XP) = sum(XP * (ln(PP) - ln(XP)))
   * with:
   *   XP = X / sqrt(sum(X^2))   and  PP = smoth(XP)
   * the gradient becomes:
   *   dS_dX = (dS_dXP - XP * sum(XP * dS_dXP)) / sqrt(sum(X^2)) 
   * where:
   *   dS_dXP = - ln(XP) + ln(PP) - XP * lnp(XP) + smooth(XP * lnp(PP)).
   */
  
  if (normalize) {
    if ((s = sqrt(sum(x * x))) == 0.0)
      error, "input array is zero everywhere";
    x /= s;
  }

  if (level) {

    /*
     * Floating smooth prior
     */
    
    p = smooth(x, level);
    dl = ln(p) - ln(x);
    if (deriv) {
      dx = dl + smooth(x * ln(p, deriv=1), level) - x * ln(x, deriv=1);
      if (normalize)
	dx = (dx - x * sum(x * dx)) / s;
    }
    return (normalize ? sum(x * dl) : sum(x) - sum(p) + sum(x * dl));
  }

  
  /*
   * Uniform prior
   */

  n = double(numberof(where(x)));
  p = sum(x) / n;
  dl = ln(p) - ln(x);
  if (deriv) {
    dx = dl + ln(p, deriv=1) / n * x - x * ln(x, deriv=1);
    if (normalize)
      dx = (dx - x * sum(x * dx)) / s;
  }
  return sum(x * dl);
}

/*--------------------------------------------------------------------------*/
func negentropy(x, &dx, deriv=, normalize=, level=)
/* DOCUMENT q= negentropy(x, level=level_of_smoothness)
            q= negentropy(x, dx, deriv=1, level=level_of_smoothness)     
     return the opposite of the entropy of array X:
       Q = -entropy(X, ...).

   SEE ALSO: entropy.
*/
{
  if (normalize) {
    if ((s = sqrt(sum(x * x))) == 0.0)
      error, "input array is zero everywhere";
    x /= s;
  }

  if (level) {

    /*
     * Floating smooth prior
     */
    
    p = smooth(x, level);
    dl = ln(x) - ln(p);
    if (deriv) {
      dx = dl - smooth(x * ln(p, deriv=1), level) + x * ln(x, deriv=1);
      if (normalize)
	dx = (dx - x * sum(x * dx)) / s;
    }
    return (normalize ? sum(x * dl) : sum(p) - sum(x) + sum(x * dl));
  }

  
  /*
   * Uniform prior
   */

  n = double(numberof(where(x)));
  p = sum(x) / n;
  dl = ln(x) - ln(p);
  if (deriv) {
    dx = dl - ln(p, deriv=1) / n * x + x * ln(x, deriv=1);
    if (normalize)
      dx = (dx - x * sum(x * dx)) / s;
  }
  return sum(x * dl);
}

/*--------------------------------------------------------------------------*/
func finiteDifference(x, zero=, order=, transp=)
/* DOCUMENT finiteDifference -- compute finite difference for smoothness()
       dx= finiteDifference(x);
   returns finite difference derivatives of X for its first dimension.
       Keyword ORDER  specifies  order  of  redivatives  (1,  2  or 3,
   default is 1).
       Keyword ZERO can be either nil, 0, 1, 2, or 3 and  indicates if
   left (bit 1 of ZERO is set?) and/or right (bit 2 of  ZERO  is set?)
   limits of X (in its first dimension) should be considered  as being
   zero.
       If keyword TRANSP  is  true,  the  tranposition  of  the finite
   difference matrix is used. This is useful to  compute  the gradient
   of the smooth error and has no effect when ORDER  is  2,  i.e., the
   finite difference matrix is symmetrical.
*/
{
    zero= is_void(zero) ? 0n : (int(zero) & 3);	// evaluates ZERO keyword

    if (is_void(order) || order == 1) {
	if (zero == 0) {
	    if (transp) {
		(dims= dimsof(x))(2)++;
		y= array(0., dims);
		y(2:, ..)   = x;
		y(:-1, ..) -= x;
	    } else {
		return x(dif, ..);
	    }
	} else if (zero == 1) {
	    // left limit is zero
	    y= x;
	    if (transp) {
		y(:-1, ..) -= x(2:, ..);
	    } else {
		y(2:, ..) -= x(:-1, ..);
	    }
	} else if (zero == 2) {
	    // right limit is zero
	    y= -x;
	    if (transp) {
		y(2:, ..) += x(:-1, ..);
	    } else {
		y(:-1, ..) += x(2:, ..);
	    }
	} else {
	    // left and right limits are zero
	    if (transp) {
		return x(:-1, ..) - x(2:, ..);
	    } else {
		(dims= dimsof(x))(2)++;
		y= array(0., dims);
		y(:-1, ..) = x;
		y(2:, ..) -= x;
	    }
	}
    } else if (order == 2) {
	y= (-2.) * x;
	if (!(zero & 1)) y(1, ..) += x(1, ..);
	if (!(zero & 2)) y(0, ..) += x(0, ..);
	y(:-1, ..) += x(2:, ..);
	y(2:, ..)  += x(:-1, ..);
    } else if (order == 3) {
	error, "ORDER=3 not yet implemented: never trust the documentation!";
    } else {
	error, "invalid ORDER";
    }
    return y;
}

/*--------------------------------------------------------------------------*/
func smoothness(x, &df, deriv=, normalize=, order=, weight=, which=, zero=)
/* DOCUMENT smoothness -- measure the smoothness af an array

 PROTOTYPE:
   func smoothness(x, &df, deriv=, order=, weight=, which=, zero=)

 SYNOPSIS:
   error= smoothness(x, ...)

 ARGUMENTS:
   x	    input array (should be float or double).
   
   df       output gradient of smoothness(X) with respect to X.

   deriv=   set it to  non-zero  to  request  the  computation  of the
            gradient (optional, default is `[]': the  gradient  is not
            computed).

   normalize=
            (optional input flag) compute  the  normalized smoothness,
	    i.e., the smoothness of X/sqrt(sum(X^2)).

   order=   (optional scalar integer input) order of  derivative: 1, 2
	    or 3 are implemented; 1st derivative is the default.

   weight=  (optional   scalar/vector   input)   weighting    of   the
            regularization term (default is 1).  If  it  is  a scalar,
            the weight  is  the  same  for  all  the  dimension  of X;
            otherwise, it must be a vector with as many elements  as X
            has  dimensions,  each  element  set  the  weight  for the
            smoothness of X along the corresponding dimension.

   which=   scalar  indicating for  which  actual dimension  of  X the
            regularization error should  be  computed  (default  is to
            compute  the  error  for  all  dimensions).   If  WHICH is
            specified, WEIGHT and ZERO must be scalar.

   zero=    (optional scalar or vector integer input) 1st (2nd) bit is
            set to assume that the array is extended by zeroes outside
            its left (right) limit.  Default is to make  no assumption
            about the outside values.  If X is multi-dimensional, ZERO
            should be either scalar or vector with as many elements as
            the number of dimensions of X,  in  the  later  case, each
            element flags the corresponding X dimension.

 DESCRIPTION:
   The calling sequence  of  `smoothness'  is  designed  for  use with
   optimisation  of  function  routines  (see  `optimConjGrad').   The
   returned  value  is  the  squared  norm  sum(DX^2)  of  the  finite
   difference derivative of X, where (in the mono-dimensional case):

       DX(i) = X(i) - X(i-1)                           if ORDER = 1
       
       DX(i) = X(i-1) - 2 * X(i) + X(i+1)              if ORDER = 2

       DX(i) = -X(i-1) + 3 * X(i) - X(i+1) + X(i+2)    if ORDER = 3

   The  lower  is  this  returned  value  the  smoother  is  the input
   array. The  returned  value  is  maximum  for  isolated  spikes and
   minimum for a uniform array.  This function may be used  to enforce
   smoothness in an inversion process.

 HISTORY:
   August 31, 1995 by Eric THIEBAUT.
   December 4, 1995 by Eric THIEBAUT: rename  the function and make it
       applicable to any number of dimensions.

 BUGS: None of the dimension of X can be smaller than 2.   For reasons
   of efficiency, the checking of the arguments may be poor...

 SEE ALSO: negentropy, smooth, optimConjGrad, optimBrent.
*/
{
    // normalization?
    if ((normalize= normalize ? 1 : 0)) {
	sx2= sum(x * x);
	if (sx2 <= 0.) {
	    if (deriv) df= array(0., dimsof(x));
	    return 0.;
	}
	x /= sqrt(sx2);
    }
    
    if (is_void(which)) {

	ndims= (dims= dimsof(x))(1);	// get dimension list of X

	// parse ZERO
	if (numberof(zero) == ndims) {
	    zero= int(zero & 3);
	} else if (numberof(zero) == 1) {
	    zero= array(int(zero & 3), ndims);
	} else if (numberof(zero) == 0) {
	    zero= array(0n, ndims);
	} else {
	    error, "bad dimension for keyword ZERO";
	}

	// parse WEIGHT
	if (numberof(weight) == 0) {
	    weight= array(1., ndims);
	} else if (numberof(weight) == 1) {
	    weight= array(double(weight), ndims);
	} else if (numberof(weight) != ndims) {
	    error, "bad dimension for keyword WEIGHT";
	}

	// compute smoothness of 1st dimension
	err= smoothness(x, df, deriv=deriv, order=order,
			weight=weight(1), which=1, zero=zero(1));

	// compute smoothness of other dimensions
	for (n=2; n<=ndims; n++) {
	    err += smoothness(x, tmp, deriv=deriv, order=order,
			      weight=weight(n), which=n, zero=zero(n));
	    if (deriv) {
		df += tmp;
		tmp= [];	// free some memory
	    }
	}

    } else {
    
	if (is_void(weight)) weight= 1.;

	if (which != 1)
	    x= transpose(x, 2-which);	// make the working index the first one
    
	dx= finiteDifference(x, zero=zero, order=order);
    
	if (deriv) {
	    df= (2. * weight) *
		finiteDifference(dx, zero=zero, order=order, transp=1n);
	    if (which != 1)
		df= transpose(df, which);	// replace the working index
	}

	err= weight * sum(dx * dx);	// DX*DX is faster than DX^2 in Yorick

    }

    // Take into account normalization for the gradient
    // and return the smoothness error
    if (deriv && normalize)
	df= (df - x * sum(x * df)) / sqrt(sx2);
    return err;
}
