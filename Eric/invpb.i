/*
 * invpb.i --
 *
 *	Routines to solve inverse-problems in Yorick.
 *
 * Copyright (c) 2000 - 2002, Eric THIEBAUT.
 *
 * Notes:
 *	To prevent conflicts in naming of symbols, all public routines of
 *	this package are prefixed by "ip_" and all private routines are
 *	prefixed by "__ip_".
 *
 * History:
 *	$Id$
 *	$Log$
 */

/*---------------------------------------------------------------------------*/
/* INTRODUCTION */
local ip_intro;
/* DOCUMENT ip_intro -- Introduction to routines for solving large scale
                        inverse problems.

   REGULARIZATION
     Computation of regularization terms and their derivatives is implemented
     via routines that have the same prototype:
       func regul(x, deriv, misc, &dx)       
     where the arguments are:
       X     = Parameters.
       DERIV = 0 => Returns RGL(X).
               1 => Returns RGL(X) and stores gradient of RGL(X) in DX.
               2 => Returns a N-by-N+1 array AB that gives the quadratic
                    approximation of RGL in the subspace defined by the
                    search directions {*DX(1), ..., *DX(N)} around X (see
                    ip_solve_local for the description of AB).
       MISC  = Anything but X that is needed to compute the regularization
               term.
       DX    = Unused (DERIV=0) or output gradient (DERIV=1) or input
               search directions (DERIV=2). 


   BOUNDS
     The parameters can be bounded.

   LOCAL QUADRATIC APPROXIMATION
     see ip_solve_local.
*/

/*---------------------------------------------------------------------------*/
/* 2ND ORDER LOCAL APPROXIMATION */

func ip_solve_local(ab, lp, dx)
/* DOCUMENT    w = ip_solve_local(ab, lp);
       -or- xstp = ip_solve_local(ab, lp, dx);

     Solve local optimization  problem for Levenberg-Marquardt parameter LP
     given  local  2nd order  Taylor  expansion  AB  as computed  by,  e.g.
     ip_ssq_approx.  If DX is omitted,  the return value W is the (guessed)
     optimum step  for the parameters; otherwise, the  optimization is done
     locally for a subspace DX of search directions and the result is:

       XSTP = W(1) * DX(..,1) + W(2) * DX(..,2) + ...

     or, if DX is an array of pointers:

       XSTP = W(1) * *DX(1) + W(2) * *DX(2) + ...

     Let PSI(P)  be the  penalty function to  optimize with respect  to the
     parameters P and  let N=numberof(P) be the number  of parameters, then
     AB is a N by N+1  array and defines the quadratic approximation of PSI
     in the vicinity of P as follow:

       PSI(P + DP) ~ PSI(P) + sum(DP * (0.5 * (A(+,) * DP(+)) - B))

     where DP is a (small) variation  of P, the vector of current parameter
     values, A=AB(,:-1)  and B=AB(,0).   If the quadratic  approximation is
     obtained by Taylor  expansion, A is the Hessian of PSI  and B is minus
     the gradient of PSI.  A=AB(,:-1) must therefore be symmetric.

     Levenberg-Marquardt parameter  LP is used  to control the size  of the
     trust  region,   i.e.   the  region  around  P   where  the  quadratic
     approximation is  assumed to be  valid.  LP must be  non-negative, the
     larger is LP  the smaller will be the  step returned by ip_solve_local
     (i.e.   the   smaller  is  the   trust  region  where   the  quadratic
     approximation  is  valid).  If  LP  is  omitted  or with  LP=0.0,  the
     quadratic approximation is assumed to be exact everywhere (i.e. PSI is
     quadratic with respect to P).

   SEE ALSO LUsolve, ip_ssq_approx, ip_expected_change. */
{
  a = ab(,:-1);
  b = ab(,0);
  n = numberof(b);
  if (lp) a(1:n*n:n+1) *= 1.0 + lp;
  w = (n==1 ? b/a : LUsolve(a, b));
  if (is_void(dx)) return w;
  if (structof(dx) == pointer) {
    if (numberof(dx) != n) error, "bad number of directions in DX";
    xstp = w(1)*(*dx(1));
    for (i=2 ; i<=n ; ++i) xstp += w(i)*(*dx(i));
  } else {
    if (dimsof(dx)(0) != n) error, "bad dimensions for DX";
    xstp = w(1)*dx(..,1);
    for (i=2 ; i<=n ; ++i) xstp += w(i)*dx(..,i);
  }
  return xstp;
}

func ip_expected_change(ab, w)
/* DOCUMENT ip_expected_change(ab);
       -or- ip_expected_change(ab, w);
     Compute expected  change in penalty function from  its local quadratic
     approximation AB (see  ip_solve_local).  W is the step  as returned by
     ip_solve_local, if W is omitted, ip_solve_local is called to get it.

   SEE ALSO ip_solve_local. */
{
  if (is_void(w)) w = ip_solve_local(ab);
  return sum((0.5*(ab(,:-1)(+,)*w(+)) - ab(,0))*w);
}

func ip_ssq_approx(res, jac, wght)
/* DOCUMENT ab = ip_ssq_approx(res, jac);
       -or- ab = ip_ssq_approx(res, jac, wght);

     Build local quadratic  approximation of a penalty function  PSI of the
     form:

       PSI = sum(WGHT * RES^2)

     where  RES are  the so-called  "residuals" assumed  to depend  on some
     parameters PARAM  (unknown and irrelevant  for this routine)  and WGHT
     are optional weights (no WGHT is the same as with WGHT=1).  JAC is the
     Jacobian of the residuals with respect to the parameters:

       *JAC(i)   = d RES / d PARAM(i)    if JAC is an array of pointers
       JAC(..,i) = d RES / d PARAM(i)    otherwise

     Such a  penalty function is common  in, e.g., model  fitting where the
     residuals are  the difference  between some data  and a model  and the
     weights are the inverse of the data variance.

     This routine returns  a N by N+1 array AB  which defines the quadratic
     approximation  of  PSI  (see  ip_solve_local).  The  approximation  is
     obtained by 2nd  order Taylor expansion and is intended  to be used to
     solve   the  non-linear  least   squares  optimization   problem  (see
     ip_solve_local).  The  term involving the second  order derivatives of
     the  residuals with  respect to  the  parameters is  neglected in  the
     approximation.

   SEE ALSO ip_solve_local. */
{
  local ws;
  if (structof(jac)==pointer) {
    n = numberof(jac);
    ab = array(double, n, n+1);
    if (is_void(wght)) {
      for (i=1 ; i<=n ; ++i) {
        eq_nocopy, ws, *jac(i);  
	ab(i,0) = -sum(ws*res);
	for (j=1 ; j<=i ; ++j) ab(j,i) = ab(i,j) = sum(ws*(*jac(j)));
      }
    } else {
      for (i=1 ; i<=n ; ++i) {
	ws = wght*(*jac(i));
	ab(i,0) = -sum(ws*res);
	for (j=1 ; j<=i ; ++j) ab(j,i) = ab(i,j) = sum(ws*(*jac(j)));
      }
    }
  } else {
    n = dimsof(jac)(0);
    ab = array(double, n, n+1);
    for (i=1 ; i<=n ; ++i) {
      ws = is_void(wght) ? jac(..,i) : wght*jac(..,i);
      ab(i,0)= -sum(ws*res);
      for (j=1 ; j<=i ; ++j) ab(j,i) = ab(i,j) = sum(ws*jac(..,j));
    }
  }
  return ab + ab;  /* same as 2*AB but should be faster at least because
                      no array broadcasting is required to comply with
                      conformation rules */
}

/*---------------------------------------------------------------------------*/
/* REGULARIZATION */

func ip_tikonov_regul(x, deriv, &dx)
/* DOCUMENT ip_tikonov_regul(x)
       -or- ip_tikonov_regul(x, deriv, dx)

     Compute Tikonov regularization value  and/or derivatives.  If DERIV is
     nil or DERIV<=1, the result is the value of the regularization term:

       RGL(X) = sum(X*X)

     If  DERIV=1, then  symbol DX  will be  used to  store the  gradient of
     RGL(X) with respect to X and RGL(X) is returned.
     
     If DERIV=2, then DX must be  an array of pointers of search directions
     and  the   result  is  the   local  quadratic  approximation   of  the
     regularization term (see ip_solve_local).

  SEE ALSO ip_solve_local, ip_ssq_approx. */
{
  if (! deriv) return sum(x*x);
  if (deriv == 1) {
    dx = x + x; /* faster than 2*x */
    if (! am_subroutine()) return sum(x*x);
  } else if (deriv==2) {
    return ip_ssq_approx(x, dx);
  } else error, "bad value for DERIV";
}

func ip_periodic_roughness(x, deriv, &dx, order=, which=)
/* DOCUMENT ip_periodic_roughness(x)
       -or- ip_periodic_roughness(x, deriv, dx)

     Compute the regularization term based on a measure of the roughness of
     X.   This measure  is quadratic  in X  and its  definition  depends on
     ORDER.  If  ORDER=1 (or if ORDER  is not specified, since  this is the
     default), the roughness of X is measured by its gradient along all its
     dimensions:
     
       RGL(X) = sum(ip_periodic_dif(X, 1)^2)
              + sum(ip_periodic_dif(X, 2)^2)
              + ...
              + sum(ip_periodic_dif(X, N)^2)
     
     where N is  the number of dimensions of X.   If ORDER=2, the roughness
     of X is measured by its Laplacian all its dimensions:
     
       RGL(X) = sum((   ip_periodic_dif(ip_periodic_dif(X, 1), 1)
                      + ip_periodic_dif(ip_periodic_dif(X, 2), 2)
                      + ...
                      + ip_periodic_dif(ip_periodic_dif(X, N), N))^2)
                      
     These  definitions  ensure  that  the  roughness  is  isotropic  (i.e.
     invariant by a rotation of the axes).

     If DERIV is nil or DERIV<=1,  the result is RGL(X), the measure of the
     roughness of X.
     
     If DERIV=1, then  symbol DX will be used to store  the gradient of the
     roughness with respect to X.
     
     If DERIV=2, then DX must be  an array of pointers of search directions
     and  the   result  is  the   local  quadratic  approximation   of  the
     regularization term (see ip_ssq_approx and ip_solve_local).

     Keyword WHICH can be used to  specify a list of dimensions along which
     the roughness  is measured.  The  default is to measure  the roughness
     along all dimensions of X.  Negative or zero values in WHICH are taken
     as relative  to the  last dimension of  X.  For  instance, WHICH=[2,0]
     means second and last dimensions of X.

  SEE ALSO ip_periodic_dif. */
{
  if (is_void(deriv)) deriv = 0;
  if (deriv==2) {
    if (structof(dx) != pointer) error, "DX must be an array of pointers";
    ndx = numberof(dx);
  } else {
    ndx = 0;
  }

  /* Get list of dimensions along which the roughness must be measured. */ 
  dims = dimsof(x);
  rank = dims(1);
  if (is_void(which)) {
    which = indgen(rank);
  } else {
    which += rank*(which <= 0);
    if (structof(which) != long) error, "bad data type for WHICH";
    if (min(which) < 1 || max(which) > rank)
      error, "out of range dimension in WHICH";
  }
  nw = numberof(which);

  /* Roughness measured by finite difference gradient. */
  if (is_void(order) || order==1) {
    if (deriv<=1) {
      rgl = 0.0;
      if (deriv) dx = array(0.0, dims);
      for (k=1 ; k<=nw ; ++k) {
        this = which(k);
	res = ip_periodic_dif(x, this);
	rgl += sum(res*res);
	if (deriv) dx += ip_periodic_dif(res, this, 1);
      }
      if (deriv) dx += dx; /* same as DX *= 2 but faster */
      return rgl;
    }
    if (deriv==2) {
      ab = 0.0;
      for (k=1 ; k<=nw ; ++k) {
        this = which(k);
	res = ip_periodic_dif(x, this);
	jac = array(pointer, ndx);
	for (i=1 ; i<=ndx ; ++i) jac(i) = &ip_periodic_dif(*dx(i), this);
	ab += ip_ssq_approx(res, jac);
	jac = [];
      }
      return ab;
    }
  }

  /* Roughness measured by Laplacian. */
  if (order==2) {
    /* Residuals. */
    this = which(1);
    res = ip_periodic_dif(ip_periodic_dif(x, this), this);
    for (k=2 ; k<=nw ; ++k) {
      this = which(k);
      res += ip_periodic_dif(ip_periodic_dif(x, this), this);
    }
    
    /* Penalty (and gradient). */
    if (deriv<=1) {
      if (deriv) {
        this = which(1);
	dx = ip_periodic_dif(ip_periodic_dif(res, this, 1), this, 1);
	for (k=2 ; k<=nw ; ++k) {
          this = which(k);
	  dx += ip_periodic_dif(ip_periodic_dif(res, this, 1), this, 1);
        }
	dx += dx; /* same as DX *= 2 but faster */
      }
      return sum(res*res);
    }

    /* 2nd order approximation. */
    if (deriv==2) {
      jac = array(pointer, ndx);
      for (i=1 ; i<=ndx ; ++i) {
        this = which(1);
	jac_i = ip_periodic_dif(ip_periodic_dif(*dx(i), this), this);
	for (k=2 ; k<=nw ; ++k) {
          this = which(k);
	  jac_i += ip_periodic_dif(ip_periodic_dif(*dx(i), this), this);
        }
	jac(i) = &jac_i;
      }
      return ip_ssq_approx(res, jac);
    }
  }

  /* Roughness measured by Laplacian. */
  if (order==3) {
    /* Residuals. */
    res = ip_periodic_dif(x, which(1), 0, 2);
    for (k=2 ; k<=nw ; ++k) res += ip_periodic_dif(x, which(k), 0, 2);
    
    /* Penalty (and gradient). */
    if (deriv<=1) {
      if (deriv) {
	dx = ip_periodic_dif(res, which(1), 1, 2);
	for (k=2 ; k<=nw ; ++k) dx += ip_periodic_dif(res, which(k), 1, 2);
	dx += dx; /* same as DX *= 2 but faster */
      }
      return sum(res*res);
    }

    /* 2nd order approximation. */
    if (deriv==2) {
      jac = array(pointer, ndx);
      for (i=1 ; i<=ndx ; ++i) {
	jac_i = ip_periodic_dif(*dx(i), which(1), 0, 2);
	for (k=2 ; k<=nw ; ++k)
          jac_i += ip_periodic_dif(*dx(i), which(k), 0, 2);
	jac(i) = &jac_i;
      }
      return ip_ssq_approx(res, jac);
    }
  }

  error, "bad value for ORDER or DERIV";
}

func ip_periodic_dif(x, which, transp, repeat)
/* DOCUMENT ip_periodic_dif(x, which)
       -or- ip_periodic_dif(x, which, transp)
       -or- ip_periodic_dif(x, which, transp, repeat)
       
     Apply  linear operator of  periodic differenciation  to array  X along
     dimension WHICH (default WHICH=1), if TRANSP is true, the transpose of
     the  operator is  applied instead.   The linear  operator  of periodic
     differenciation is (for a vector X of length N):
     
       (D.x)(i) = x(i+1) - x(i)     for i=1..N-1 
                  x(1)   - x(N)     for i=N

     its transpose is:
     
       (D'.x)(i)= x(N)   - x(1)     for i=1
                  x(i-1) - x(i)     for i=2..N

     If REPEAT is specified, the operator D (or its transpose D') is applied
     REPEAT times (at least once).
     
   BUG
     The use of the transpose() builtin routine may not be optimal.

   SEE ALSO ip_periodic_roughness, transpose. */
{
  if (is_void(repeat)) repeat = 1;
  if (is_void(which)) which = 1;
  else if (which!=1) x = transpose(x, [1,which]);
  n = 0;
  for ( ; ; ) {
    y = -x;
    if (transp) {
      y(1,..) += x(0,..);
      y(2:,..) += x(:-1,..);
    } else {
      y(:-1,..) += x(2:,..);
      y(0,..) += x(1,..);
    }
    if (++n >= repeat) return (which==1 ? y : transpose(y, [which,1]));
    x = y;
  }
}

func ip_roughness(x, deriv, &dx, order=, which=)
/* DOCUMENT ip_roughness(x)
       -or- ip_roughness(x, deriv, dx)

     Compute the regularization term based on a measure of the roughness of
     X.   This measure  is quadratic  in X  and its  definition  depends on
     ORDER.  If  ORDER=1 (or if ORDER  is not specified, since  this is the
     default), the roughness of X is measured by its gradient along all its
     dimensions:


===============FIXME===================
       RGL(X) = sum(ip_periodic_dif(X, 1)^2)
              + sum(ip_periodic_dif(X, 2)^2)
              + ...
              + sum(ip_periodic_dif(X, N)^2)
     
     where N is  the number of dimensions of X.   If ORDER=2, the roughness
     of X is measured by its Laplacian all its dimensions:
     
       RGL(X) = sum((   ip_periodic_dif(ip_periodic_dif(X, 1), 1)
                      + ip_periodic_dif(ip_periodic_dif(X, 2), 2)
                      + ...
                      + ip_periodic_dif(ip_periodic_dif(X, N), N))^2)
                      
     These  definitions  ensure  that  the  roughness  is  isotropic  (i.e.
     invariant by a rotation of the axes).

     If DERIV is nil or DERIV<=1,  the result is RGL(X), the measure of the
     roughness of X.
     
     If DERIV=1, then  symbol DX will be used to store  the gradient of the
     roughness with respect to X.
     
     If DERIV=2, then DX must be  an array of pointers of search directions
     and  the   result  is  the   local  quadratic  approximation   of  the
     regularization term (see ip_ssq_approx and ip_solve_local).

     Keyword WHICH can be used to  specify a list of dimensions along which
     the roughness  is measured.  The  default is to measure  the roughness
     along all dimensions of X.  Negative or zero values in WHICH are taken
     as relative  to the  last dimension of  X.  For  instance, WHICH=[2,0]
     means second and last dimensions of X.

  SEE ALSO ip_periodic_dif. */
{
  _dif = __ip_dif_op;
  
  if (is_void(deriv)) deriv = 0;
  if (deriv==2) {
    if (structof(dx) != pointer) error, "DX must be an array of pointers";
    ndx = numberof(dx);
  } else {
    ndx = 0;
  }

  /* Get list of dimensions along which the roughness must be measured. */ 
  dims = dimsof(x);
  rank = dims(1);
  if (is_void(which)) {
    which = indgen(rank);
  } else {
    which += rank*(which <= 0);
    if (structof(which) != long) error, "bad data type for WHICH";
    if (min(which) < 1 || max(which) > rank)
      error, "out of range dimension in WHICH";
  }
  nw = numberof(which);

  /* Roughness measured by finite difference gradient. */
  if (is_void(order) || order==1) {
    if (deriv<=1) {
      rgl = 0.0;
      if (deriv) dx = array(0.0, dims);
      for (k=1 ; k<=nw ; ++k) {
        this = which(k);
	res = _dif(x, this);
	rgl += sum(res*res);
	if (deriv) dx += _dif(res, this, 1);
      }
      if (deriv) dx += dx; /* same as DX *= 2 but faster */
      return rgl;
    }
    if (deriv==2) {
      ab = 0.0;
      for (k=1 ; k<=nw ; ++k) {
        this = which(k);
	res = _dif(x, this);
	jac = array(pointer, ndx);
	for (i=1 ; i<=ndx ; ++i) jac(i) = &_dif(*dx(i), this);
	ab += ip_ssq_approx(res, jac);
	jac = [];
      }
      return ab;
    }
  }

  /* Roughness measured by Laplacian. */
  if (order==2) {
    /* Residuals. */
    this = which(1);
    res = _dif(_dif(x, this), this);
    for (k=2 ; k<=nw ; ++k) {
      this = which(k);
      res += _dif(_dif(x, this), this);
    }
    
    /* Penalty (and gradient). */
    if (deriv<=1) {
      if (deriv) {
        this = which(1);
	dx = _dif(_dif(res, this, 1), this, 1);
	for (k=2 ; k<=nw ; ++k) {
          this = which(k);
	  dx += _dif(_dif(res, this, 1), this, 1);
        }
	dx += dx; /* same as DX *= 2 but faster */
      }
      return sum(res*res);
    }

    /* 2nd order approximation. */
    if (deriv==2) {
      jac = array(pointer, ndx);
      for (i=1 ; i<=ndx ; ++i) {
        this = which(1);
	jac_i = _dif(_dif(*dx(i), this), this);
	for (k=2 ; k<=nw ; ++k) {
          this = which(k);
	  jac_i += _dif(_dif(*dx(i), this), this);
        }
	jac(i) = &jac_i;
      }
      return ip_ssq_approx(res, jac);
    }
  }

  error, "bad value for ORDER or DERIV";
}

func __ip_dif_op(x, which, transp, do_not_reorder)
{
  if (which != 1) x = transpose(x, [1,which]);
  if (transp) {
    dims = dimsof(x);
    ++dims(2);
    y = array(double, dims);
    y(2:0,) = x;
    y(1:-1,) -= x;
  } else {
    y = x(dif,);
  }
  return ((which==1 || do_not_reorder) ? y : transpose(y, [1,which]));
}

/*---------------------------------------------------------------------------*/
/* TUNE REGULARIZATION WEIGHT */

func __ip_tune_func(mu)
/* DOCUMENT __ip_tune_func(mu)
     Private function for ip_tune_regul_weight. */
{
  x = SVsolve(al + mu*ar, bl + mu*br);
  dl = sum((0.5*(al(+,)*x(+)) - bl)*x);
  if (update == 1) return l - aim + dl;
  dr = sum((0.5*(ar(+,)*x(+)) - br)*x);
  return (l + mu*r) - aim + (dl + mu*dr);
}

func ip_tune_regul_weight(l, abl, r, abr, mu, aim=, update=, tol=, rate=)
/* DOCUMENT new_mu = ip_tune_regul_weight(l, abl, r, abr, mu)

     Estimate the next  value of the regularization weight  given the local
     approximations of the likelihood and regularization terms and a target
     constraint.

          L = current value of the likelihood term
        ABL = local approximation of the likelihood term
          R = current value of the regularization term
        ABR = local approximation of the regularization term
         MU = current value of the regularization weight
        AIM = target value, default value is 1.0;
     UPDATE = 1 to update MU so that L tends to AIM
              2 to update MU so that L + MU*R tends to AIM
              otherwise, MU remains constant
        TOL = relative tolerance to adjust MU; default is 1e-5;
       RATE = factor to control the maximum rate of change for MU: the
              returned MU will be such that MU/RATE <= MU <= MU*RATE
              default value is 10, must be > 1.  */
{
  extern SVrank;

  /* The default value for MU attempts to balance the magnitudes of the
     likelihood and the regularization terms. */
  if (is_void(mu) && mu <= 0.0) {
    abs_l = abs(l);
    abs_r = abs(r);
    mu = (abs_r > 0.0 && abs_l > 0.0) ?  abs_l/abs_r : 1.0;
  }
  
  if (update == 1 || update == 2) {
    /* Get optional settings. */
    if (is_void(alpha)) alpha = 0.67; 
    if (is_void(tol))   tol = 1e-5;
    if (is_void(aim))   aim = 1.0;
    if (is_void(rate))  rate = 10.0;
  
    /* Extract, Hessian and (minus) gradient parts from local
       approximations of the likelihood and regularization terms.
       Use a bissection method to estimate the next value of MU. */
    al = abl(,:-1);
    bl = abl(,0);
    ar = abr(,:-1);
    br = abr(,0);
    mu1 = mu/rate;
    mu2 = mu*rate;
    t1 = __ip_tune_func(mu1);
    t2 = __ip_tune_func(mu2);
    if (t1*t2 < 0.0) {
      for (;;) {
        mu3 = (mu1 + mu2)/2.0;
        t3 = __ip_tune_func(mu3);
        if (ip_almost_equal(t3, 0.0, tol) ||
            ip_almost_equal(mu1, mu2, tol)) return mu3;
        if (t1*t3 < 0.0) {
          mu2 = mu3;
          t2 = t3;
        } else {
          mu1 = mu3;
          t1 = t3;
        }
      }
    } else if (abs(t1) < abs(t2)) {
      return mu1;
    } else {
      return mu2;
    }
  }
  return mu; /* MU remains always constant. */
}

/*---------------------------------------------------------------------------*/
/* PARAMETER BOUNDS */

func ip_apply_positivity(x, dx)
{
  min_x = min(x);
  if (is_void(dx)) return min_x >= 0.0 ? x : max(0.0, x);
  if (structof(dx) != pointer) error, "DX must be an array of pointers";
  if (min_x <= 0.0) {
    if (min_x < 0.0) error, "non-positive X";
    low = x <= 0.0;
    ndx = numberof(dx);
    for (i=1 ; i<=ndx ; ++i) {
      out = low & (*dx(i) < 0.0);
      if (anyof(out)) dx(i) = &((!out) * *dx(i));
    }
  }
}

func ip_check_bounds(x, &bnd, &lo, &hi)
/* DOCUMENT ip_check_bounds(x, bnd, lo, hi)
   
     Check/setup  bound settings  for  ip_apply_bounds routine.   X is  the
     array of parameters,  BND specifies the type of bounds,  LO and HI are
     the  lower  and  upper  bounds  respectively.   BND,  HI  and  LO  are
     input/output  symbols and  they  may  be nil  otherwise  they must  be
     conformable  with  X.  Their  contents  is  properly defined/fixed  on
     return so that they can be directly used by ip_apply_bounds.

     If BND is undefined, then  all parameters will be considered as having
     a lower  bound if LO is  given and/or an  upper bound if HI  is given.
     Otherwise, BND is an array of integers conformable with X such that:
     
       BND(i) = 0   no bound
                1   lower bound only:  LO(i) <= X(i)
                2   upper bound only:           X(i) <= HI(i)
                3   boxed value:       LO(i) <= X(i) <= HI(i)
                  
     On return,  if there is  no bound constraints  at all, BND, LO  and HI
     will  be nil.   Otherwise, BND  will be  an array  of int's  with same
     dimension list  as the parameters X and  LO and HI will  be either nil
     (i.e.  there is  no  such bound)  or  array(s) of  double's with  same
     dimension list as X.

     SEE ALSO: ip_apply_bounds. */
{  
  if (! is_void(lo)) {
    if ((s = structof(lo)) != double) {
      if (s != char && s != short && s != int && s != long && s != float)
        error, "bad data type for lower bound";
      lo = double(lo);
    }
    if (! is_array(dimsof(x, lo)))
      error, "lower bound not conformable with array of parameters";
    if (numberof(lo) != numberof(x)) lo *= array(1.0, dimsof(x));
  }
  if (! is_void(hi)) {
    if ((s = structof(hi)) != double) {
      if (s != char && s != short && s != int && s != long && s != float)
        error, "bad data type for upper bound";
      hi = double(hi);
    }
    if (! is_array(dimsof(x, hi)))
      error, "upper bound not conformable with array of parameters";  
    if (numberof(hi) != numberof(x)) hi *= array(1.0, dimsof(x));
  }
  flag = (! is_void(lo)) | 2n*(! is_void(hi));
  if (is_void(bnd)) {
    if (flag) bnd = array(flag, dimsof(x));
  } else {
    if ((s = structof(bnd)) != int) {
      if (s != char && s != short && s != long)
        error, "bad data type for bound case";
      bnd = int(bnd);
    }
    if (! is_array(dimsof(x, bnd)))
      error, "bound case not conformable with array of parameters";
    if (max(bnd) > 3n || min(bnd) < 0n) error, "bad value in bound case array";
    if (numberof(bnd) != numberof(x)) bnd *= array(1n, dimsof(x));
    if (anyof(bnd&1n)) {
      if (is_void(lo)) error, "missing lower bound array";
    } else if (! is_void(lo)) {
      write, "WARNING - unused lower bound";
      lo = [];
    }
    if (anyof(bnd&2n)) {
      if (is_void(hi)) error, "missing upper bound array";
    } else if (! is_void(hi)) {
      write, "WARNING - unused upper bound";
      hi = [];
    }
  }
  if (flag == 3n) {
    i = where(bnd == 3n);
    if (is_array(i) && anyof(lo(i) > hi(i)))
      error, "incompatible lower-upper bounds";
  }
}

func ip_apply_bounds(&x, bnd, lo, hi, dx, op=, data=)
/* DOCUMENT ip_apply_bounds, x, bnd, lo, hi;
       -or- ip_apply_bounds, x, bnd, lo, hi, dx;
     Apply  bounds to  parameters X  or to  list of  search  directions DX.
     Parameters BND,  HI and LO are  the bound settings and  must have been
     prepared  by ip_check_bounds.   If DX  is  not given,  the bounds  are
     applied to  the parameters X (used  as an input/output  symbol in this
     case);  otherwise, DX  must  be an  array  of pointers  to the  search
     directions and the bound constraints are applied to theses directions.
     Optionnaly, the bound constraints may be applied to a linear transform
     of X (see the OP and DATA keywords below).

   KEYWORDS
     DATA - relevant data for transformation OP (see below).

     OP   - transformation routine;  the bounds are applied to  Y=OP(X); OP
            must be a  one-to-one linear  transformation; the  prototype of
            function OP is:
            
              func OP(x, dir, data) {
                if (dir < 0) return ...; // inverse transform
                return ...; // transform
              }

            where  DIR is  the  transformation direction  negative for  the
            inverse, otherwise positive or zero and DATA is anything passed
            by the DATA keyword of the ip_apply_bounds routine.

   SEE ALSO: ip_check_bounds. */
{
  if (is_void(bnd)) return; /* nothing to do */  
  if (is_void(dx)) {
    /* Apply bound constraints to the parameters. */
    if (! is_void(op)) x = op(x, 1, data);
    if (! is_void(lo)) {
      i = where((bnd & 1n) * (x < lo));
      if (is_array(i)) x(i) = lo(i);
    }
    if (! is_void(hi)) {
      i = where((bnd & 2n) * (x > hi));
      if (is_array(i)) x(i) = hi(i);
    }
    if (! is_void(op)) x = op(x, -1, data);
  } else {
    /* Apply bound constraints to the search directions. */
    ndx = numberof(dx);
    at_lo = is_void(lo) ? 0n : (bnd & 1n) * (x <= lo);
    at_hi = is_void(hi) ? 0n : (bnd & 2n) * (x >= hi);
    for (i=1 ; i<=ndx ; ++i) {
      s = is_void(op) ? *dx(i) : op(*dx(i), 1, data);
      out = is_void(hi) ? at_lo * (s < 0.0) : \
           (is_void(lo) ? at_hi * (s > 0.0) : \
                         (at_lo * (s < 0.0)) + (at_hi * (s > 0.0)))
      if (anyof(out))
        dx(i) = &(is_void(op) ? (!out) * s : op((!out) * s, -1, data));
    }
  }
}

/*---------------------------------------------------------------------------*/
/* MISCELLANEOUS */

func ip_enorm(v) { return (is_void(v) ? 0.0 : sqrt(sum(v*v))); }
/* DOCUMENT ip_enorm(v);
     Return Euclidian norm of V: sqrt(sum(V*V)) (or 0.0 if V is nil). */

func ip_almost_equal(a, b, tol)
/* DOCUMENT ip_almost_equal(a, b, tol)
     Return true (non-zero) if relative  difference between A and B is less
     or equal TOL.  A and B must be floating point scalars.

   SEE ALSO ip_relative_difference. */
{ return 2.0*abs(a - b) <= tol*(abs(a) + abs(b)); }

func ip_relative_difference(a, b)
/* DOCUMENT ip_relative_difference(a, b)   
     Return  elementwise relative difference  between arrays  A and  B.  To
     avoid division by  zero, the relative difference is  taken to be equal
     to 0 everywhere A and B have the same value.  A and B must be floating
     point conformable arrays.

   SEE ALSO ip_almost_equal. */
{
  r = abs(a - b);
  s = abs(a) + abs(b);
  if (anyof((zero = !s))) s(where(zero)) = 1.0;
  return (r + r) / s;
}

func ip_default(&symbol, default_value)
/* DOCUMENT ip_default, symbol, default_value;
       -or- ip_default(expression, default_value);
     When called as  a subroutine, store DEFAULT_VALUE into  SYMBOL if this
     latter  is  void.   When  called   as  a  function,  return  value  of
     EXPRESSION, if it is non-void; DEFAULT_VALUE, ortherwise. */
{
  if (! am_subroutine()) return (is_void(symbol) ? default_value : symbol);
  if (is_void(symbol)) symbol = default_value;
}

/*---------------------------------------------------------------------------*/
