/*
 * optim_pack.i -
 *
 *	Optimization routines for Yorick.
 *
 *
 * Routines:
 *	Multidimensional minimization:
 *	  optim_vmlm - limited memory variable metric.
 *	  optim_cgmn - conjugate gradients with optional preconditioning.
 *	  optim_pcgmn - preconditioned conjugate gradients.
 *
 *	Line search:
 *	  optim_csrch - finds a step that satisfies a sufficient decrease
 *	  	condition and a curvature condition.
 *	  optim_cstep - computes a safeguarded step for a line search
 *		procedure.
 *
 *
 * Acknowledgments:
 *	This package is partially based on routines written by others:
 *
 *	 - The line search and variable metric limited memory routines
 *	   are derived from the MINPACK-2 project:
 *           ftp://info.mcs.anl.gov/pub/MINPACK-2/
 *
 *	 - Conjugate gradient is based on CG+ package by Gilbert, J.C. and
 *	   Nocedal, J. [1]: http://www.ece.nwu.edu/~rwaltz/CG+.html
 *
 *	The original routines were converted in Yorick and _improved_ in
 *	many aspects:
 *
 *	 1. The multidimensional minimization routines automatically
 *	    resume the search whenever the computed search direction is
 *	    not a descent and should therefore be more robust with respect
 *	    to non-linear problems.
 *
 *	 2. Preconditioning can be used in conjugate gradient method.
 *
 *       3. The multidimensional minimization routines can account for
 *	    simple bound constraints on the parameters.
 *
 *
 * References:
 *	[1] Gilbert, J.C. and Nocedal, J. (1992). "Global Convergence
 *	    Properties of Conjugate Gradient Methods", SIAM Journal on
 *	    Optimization, Vol. 2, pp. 21-42.
 *	[2] Shewchuk, J.R. (1994). "An Introduction to the Conjugate
 *	    Gradient Method Without the Agonizing Pain",
 *	    ftp://warp.cs.cmu.edu/quake-papers/painless-conjugate-gradient.ps
 *
 * History:
 *	January-November 2001.
 *	Observatoire de Lyon (France).
 *	Eric Thiébaut.
 *	---------------------------------------------------------------------
 *	$Id: optim_pack.i,v 1.16 2002/12/19 17:00:39 eric Exp $
 *	$Log: optim_pack.i,v $
 *	Revision 1.16  2002/12/19 17:00:39  eric
 *	 - optim_vmlm: fixed returned values in X, F and G when line search
 *	               returns with a warning.
 *
 *	Revision 1.15  2002/12/16 08:03:43  eric
 *	 - changes internals of optim_driver to make loop more general
 *
 *	Revision 1.14  2002/10/11 16:01:16  eric
 *	 - Fixed a invalid call to optim_cgmn_setup in optim_driver
 *	   (thanks Pierre OCVIRK).
 *
 *	Revision 1.13  2002/07/11 15:37:50  eric
 *	 - optim_driver: penalty function can have additional argument.
 *
 *	Revision 1.12  2002/06/06 16:27:09  eric
 *	 - optim_vmlm: It is now possible to use an active subset of
 *	               parameters, e.g. to apply bound constraints.
 *	 - optim_vmlm: Pairs for which RHO=S'.Y <= 0 are discarded
 *	               (i.e. the inverse Hessian approximation must
 *	               be positive definite).  This avoid potential
 *	               division by zero when numerical rounding
 *	               errors prevent progress.
 *
 *	Revision 1.11  2002/05/31 16:09:04  eric
 *	 - Fix bisection step test in optim_csrch.
 *
 *	Revision 1.10  2002/05/29 15:26:17  eric
 *	 - Fixed setting of STAGE according to TASK in optim_driver.
 *
 *	Revision 1.9  2002/05/29 13:58:16  eric
 *	 - Removed option TASK="stop", i.e. there is no need to
 *	   save ``best solution found so far'' since the routines
 *	   warrant that X, F and G are left unchanged when the
 *	   current step satisfies the sufficient decrease condition
 *	   and is therefore accepted as the new solution (in fact X
 *	   is only changed upon return when TASK is "fg").  This saves
 *	   2N storage.
 *	 - Renamed optim_vmlm_ws and optim_cgmn_ws as optim_vmlm_setup
 *	   and optim_cgmn_setup (in fact this was done in previous
 *	   registered version but I forgot to mention it).
 *
 *	Revision 1.8  2002/05/29 12:00:55  eric
 *	 - Add possibility to use effective step (PRJ option) in VLML.
 *	 - Add option TASK="stop" to get the best solution found
 *	   so far.
 *	 - Fix documentation of optim_csrch.
 *
 *	Revision 1.7  2002/05/24 15:55:02  eric
 *	 - vmlm: some cleanup of the code.
 *	 - vmlm: add a hack to optionally use projected gradient (GP).
 *
 *	Revision 1.6  2002/05/23 22:20:10  eric
 *	 - Complete rewrite of conjugate gradient minimization (CGMN):
 *	     (i) same kind of workspace handling has vmlm;
 *	     (ii) a single routine for preconditioning gradient or not;
 *	     (iii) optim_driver modified to account for these changes.
 *	 - Minor changes in vmlm code and some fixes in the documentation.
 *
 *	Revision 1.5  2002/05/23 08:36:06  eric
 *	 - vmlm: remove the "goto done;" statements and "done:" label.
 *	 - vmlm: store a private copy of FGRAD instead of &FGRAD in case caller
 *	   changes this array in-place later.
 *	 - vmlm: minor (cosmetic) changes to improve readability.
 *
 *	Revision 1.4  2002/03/25 16:33:54  eric
 *	 - Define X0 as local in optim_driver (thanks to Ch. Pichon).
 *
 *	Revision 1.3  2002/03/19 08:06:46  eric
 *	 - fixed optim_driver and add the display of CPU
 *
 *	Revision 1.2  2002/03/15 20:01:39  eric
 *	 - Use array of pointers for Y and S in optim_vmlm.
 *	 - Completely change the interface of optim_vmlm: a single
 *	   workspace structure is used instead of many arrays.
 *
 *	Revision 1.1  2002/03/15 10:09:41  eric
 *	Initial revision
 *
 */

/*---------------------------------------------------------------------------*/
/* DRIVER ROUTINES */

func optim_driver(f, x, &fout, &gout, farg=,
                  method=, maxiter=, maxeval=, fmin=,
                  frtol=, fatol=, verb=, ndirs=, checkderiv=)
/* DOCUMENT xmin = optim_driver(f, x, fmin_out)
     Get the location of a minimum of a multivariate function.
     Arguments are:
       F - user defined function to optimize.  The prototype of F is
           as follow:
             func F(x, &gx)
             {
               fx = ....; // compute function value at X
               gx = ....; // store gradient of F in GX
               return fx; // return F(X)
             }
       X - starting solution (a floating point array).
       FOUT - optional output variable to store the value of F at
           the minimum.
       GOUT - optional output variable to store the value of the
           gradient of F at the minimum.

     The following keywords are available:
       FARG  - supplemental  argument for  F; if  non-nil, F  is  called as
           F(X,GX,FARG) so its prototype must be: func F(x, &gx, arg) ...
       METHOD - scalar integer that define the optimization method to use.
           If METHOD is nil or zero, variable metric limited memory (VMLM)
           method is used (this is the default, also see keyword NDIRS).
           Otherwise conjugate gradient method is used, the bits of METHOD
           are interpreted as follow:
             (METHOD&3) = 1   for Fletcher-Reeves method
                          2   for Polak-Ribiere method
                          3   for Polak-Ribiere method with BETA>0
             (METHOD&4) = 0   to use the previous search direction
                          4   to use effective step
       NDIRS - number of previous directions used in variable metric
           limited memory method (default 5).
       FMIN - expected value at the minimum (default: 0.0).
       MAXITER - maximum number of iterations (default: no limits).
       MAXEVAL - maximum number of function evaluations (default: no limits).
       FRTOL - tolerance for relative function change for convergence
           (default: 1e-8).
       FATOL - tolerance for absolute function change for convergence
           (default: 0.0).
       VERB - verbose mode?  If non-nil and non-zero, print out information
           every VERB iterations.

     Note: in case of early termination, the best solution found so far
     is returned.
           
   SEE ALSO: optim_cgmn, optim_vmlm, optim_vmlm_setup. */
{
  /* Initializes workspace according to optimization method and options. */
  if (is_void(frtol)) frtol = 1e-8;
  if (is_void(fatol)) fatol = 0.0;
  if (method) {
    /* Use conjugate gradients. */
    ws = optim_cgmn_setup(method, fmin=fmin, fatol=fatol, frtol=frtol);
    next_stage = optim_cgmn;
    get_step   = optim_cgmn_step;
  } else {
    /* Use variable metric. */
    if (is_void(ndirs)) ndirs = 5;
    ws = optim_vmlm_setup(ndirs, fmin=fmin, fatol=fatol, frtol=frtol);
    next_stage = optim_vmlm;
    get_step   = optim_vmlm_step;
  }
  task = "start";

  cmplx = (structof(x) == complex);                                                          
  
  msg = string(0);
  if (verb) {
    elapsed = array(double, 3);
    timer, elapsed;
    cpu_start = elapsed(1);
  }
  eval = iter = 0;
  stp = 0.0;
  for (;;) {
    /* Setup STAGE according to TASK. */
    if (task == "fg") {
      stage = 1;
    } else if (task == "start") {
      stage = 0;
    } else if ((ident = strpart(task, 1:4)) == "newx") {
      stage = 2;
    } else if (ident == "conv") {
      stage = 3;
    } else {
      stage = 4;
    }

    if (stage <= 1) {
      /* Compute function value and gradient at X
         (TASK is "fg" or "start"). */
      local gx;
      fx = (is_void(farg) ? f(x, gx) : f(x, gx, farg));
      if (checkderiv) {
        epsilon = 1e-5;
        ax = max(abs(x));
        if (! ax) ax = 1.0;
        for (j=1 ; j<=checkderiv ; ++j) {
          i = 1 + long(numberof(x)*random());
          xsave = x(i);
          if (cmplx) {
            if (! (dx = abs(xsave.re)*epsilon)) dx = ax*epsilon;
            x(i) = xsave + dx;
            _op_checkder_print, "[Re]", i, dx, fx, (is_void(farg) ? f(x) : f(x, , farg)), gx(i).re;
            if (! (dx = abs(xsave.im)*epsilon)) dx = ax*epsilon;
            x(i) = xsave + 1i*dx;
            _op_checkder_print, "[Im]", i, dx, fx, (is_void(farg) ? f(x) : f(x, , farg)), gx(i).im;
          } else {
            if (! (dx = abs(xsave)*epsilon)) dx = ax*epsilon;
            x(i) = xsave + dx;
            _op_checkder_print, "", i, dx, fx, (is_void(farg) ? f(x) : f(x, , farg)), gx(i);
          }
          x(i) = xsave;
        }
      }
      checkderiv = 0;
      ++eval;
    }

    /* Unless TASK is "fg", decide to continue/stop and/or
       display information. */
    if (stage != 1) {
      if (stage == 2 || stage == 3) ++iter;
      if (! (stop = (stage >= 3))) {
        if ((stop = (! is_void(maxiter) && iter > maxiter))) {
          msg = swrite(format="warning: too many iterations (%d)", iter);
        } else if ((stop = (! is_void(maxeval) && eval > maxeval))) {
          msg = swrite(format="warning: too many function evaluations (%d)",
                       eval);
        }
      }
      if (verb) {
        if (eval == 1) {
          write, format="%s  %s\n%s  %s\n",
            " ITER    EVAL     CPU [s]            FUNC             GNORM   ",
            " STEPLEN",
            "------  ------  ----------  -----------------------  ---------",
            "---------";
        }
        if (stop || iter%verb == 0) {
          timer, elapsed;
          cpu = elapsed(1) - cpu_start;
          write, format="%6d  %6d %11.2e  %23.16e  %9.1e  %9.1e\n",
            iter, eval, cpu, (fx), (optim_nrm2(gx)), stp;
        }
      }
      if (stop) {
        if (msg) write, format="%s\n", msg;
        return x;
      }
    }
    
    /* Call optimizer. */
    if (cmplx) {                                                                             
      x = op_z2d(x);                                                                         
      gx = op_z2d(gx);                                                                       
      next_stage, x, fx, gx, task, ws;
      x = op_d2z(x);                                                                         
      gx = op_d2z(gx);                                                                       
    } else {                                                                                 
      next_stage, x, fx, gx, task, ws;
    }                                                                                        
    stp = get_step(ws);
  }
}

func _op_checkder_print(id, i, dx, fx, fxpdx, g)
{
  if (0) {
    gp = (fxpdx - fx)/dx;
    write, format="CHECKDER: %s %8d  comp=%+e  finitdif=%+e  relat=%f\n",
      id, i, g, gp, (g == gp ? 0.0 : 2.0*abs(g - gp)/(abs(g) + abs(gp)));
  } else {
    df1 = fxpdx - fx;
    df2 = g*dx;
    if (df1 && df2) {
      write, format="CHECKDER: %s %8d  df(auto)=%+e  df(code)=%+e  reldif=%f\n",
        id, i, df1, df2, (df1 == df2 ? 0.0 : 2.0*abs(df1 - df2)/(abs(df1) + abs(df2)));
    } else {
      write, format="CHECKDER: %s %8d  df(auto)=%+e  df(code)=%+e  absdif=%f\n",
        id, i, df1, df2, abs(df1 - df2);
    }
  }
}

func op_z2d(z) { return [double(z), z.im]; }                                                 
func op_d2z(d) { return (1i*d(..,2) + d(..,1)); }                                            
/* DOCUMENT d = op_z2d(z);                                                                   
       -or- z = op_d2z(d);                                                                   
     Convert between complex array Z and array of doubles D.  The real part                  
     of Z is D(..,1) and the imaginary part of Z is D(..,2).                                 
                                                                                             
   SEE ALSO: double, complex. */                                                             
                                                                                             
/*---------------------------------------------------------------------------*/
/* VARIABLE METRIC METHOD (LIMITED MEMORY BFGS) */

func optim_vmlm(&x, &f, &g, &task, ws, bound)
/* DOCUMENT optim_vmlmb, x, f, g, task, ws;
       -or- optim_vmlmb, x, f, g, task, ws, bound;
     This  subroutine  computes  a  local  minimizer of  a  function  of  N
     variables by a limited  memory variable metric method; optionally, the
     parameters may  be bounded.  The  user must evaluate the  function and
     the gradient.

     This subroutine  uses reverse communication.  The user  must choose an
     initial approximation  X to the  minimizer, evaluate the  function and
     the gradient at X, and make the initial call with TASK set to "start".
     On exit TASK indicates the required action.

     A typical invocation of optim_vmlm has the following outline:

     |  x = ... ;                                   // choose a starting vector
     |  ws = optim_vmlm_setup(m, fmin=lower_bound); // allocate workspace
     |  task = "start";
     |  for (;;) {
     |    if (task == "fg" || task == "start") {
     |      f = ...;  // evaluate the function at X; store in F
     |      g = ...;  // evaluate the gradient of F at X; store in G
     |    } else if (task == "newx") {
     |       // New successful step: the approximation X, function F, and
     |       // gradient G, are available for inspection.
     |    } else {
     |      break;  // convergence, or error, or warning
     |    }
     |    optim_vmlm, x, f, g, task, ws;
     |  }

     The subroutine arguments are:

       X   is a double precision array.  On entry, X is an approximation to
           the solution.  On exit, X is the current approximation.
       F   is a double precision variable.  On entry, F is the value of the
           function at X.  On final exit, F is the function value at X.
       G   is a double precision array  with same dimension list  as X.  On
           entry,  G is the  value of  the gradient at X.  On final exit, G
           is the value of the gradient at X.
       TASK is  a string variable.  On  initial entry, TASK must  be set to
           "start".  On exit, TASK indicates the required action:
           - If TASK  = "fg" then evaluate  the function and  gradient at X
             and call optim_vmlm again.
           - If strpart(TASK,  1:4) =  "newx" then a  new iterate  has been
             computed.   The  approximation X,  function F, and  gradient G
             are available for examination.
           - If strpart(TASK, 1:4) = "conv" then the search is successful.
           - If strpart(TASK, 1:4) = "warn" then the subroutine is not able
             to satisfy  the convergence conditions.   The exit value  of X
             contains the best approximation found so far.
           - If strpart(TASK, 1:5) = "error"  then there is an error in the
             input arguments.
           On exit  with convergence, a  warning or an error,  the variable
           TASK contains additional information.
       WS  is a workspace array that can be obtained from optim_vmlm_setup.
           The user must not alter  this workspace between calls (except by
           using optim_vmlm_fatol, optim_vmlm_frtol, or optim_vmlm_fmin).
       BOUND  is  an optional  array  with same  dimension  list  as X  and
           provided by the caller if  the values in X has bounds.  Non-zero
           values in BOUND indicate that  the corresponding values in X has
           reached a bound  and should not be changed  during the next step
           because  the gradient  has the  wrong sign  (i.e.   the steepest
           descent direction would violate the bound constraints):
             BOUND(i) = 0 if i-th value has a lower bound XLO(i)
                             and X(i)=XLO(i) and G(i)>0 
                        0 if i-th value has an upper bound XHI(i)
                             and X(i)=XHI(i) and G(i)<0
                        1 (or any non-zero value) otherwise
           If X has  (some) bounds, the caller is  responsible for applying
           the bounds to  X before evaluating the function  value F and the
           gradient G, e.g.:
             if (X(i) < XLO(i)) X(i) = XLO(i)
             if (X(i) > XHI(i)) X(i) = XHI(i)
           If  X has  (some) bounds,  it is  also recommended  to  "use the
           effective step"  in the updating of the  Hessian information (to
           that end, activate option PRJ in optim_vmlm_setup).  BOUND needs
           only  to  be  computed  (and  specified)  when  TASK="start"  or
           TASK="newx" (i.e.  after a  successful step).  BOUND may also be
           specified when strpart(TASK, 1:4)="conv" (i.e. after convergence
           if caller wish to continue with minimization).


   HISTORY:
     MINPACK-2 Project. April 1995.
     Argonne National Laboratory and University of Minnesota.
     Brett M. Averick, Richard G. Carter, and Jorge J. More'.

     Yorick translation an improvements.  October 2001 - June 2002.
     Observatoire de Lyon.
     Eric Thiebaut.

   SEE ALSO: optim_csrch, optim_vmlm_setup, optim_vmlm_iter. */
{
  /*
   * Differences with original FORTRAN version:
   *
   *   (1) I use integer STAGE (saved in ISAVE) instead of string WORK; the
   *       correspondance is:
   *             stage     work
   *               0     "start search"
   *               1     "search"
   *               2     "search direction"
   *
   *   (2) Slots in ISAVE and DSAVE are shifted to avoid extracting parts
   *       of these array before calling DCSRCH and then copying back (this
   *       trick can be avoided in C but not in Yorick).
   *
   *   (3) I use a more consistent number of iterations (or successful steps):
   *       on start, I set ITER=0 (instead of 1); therefore:
   *         mp = min(m, iter);
   *       instead of:
   *         mp = min(m, iter-1);
   *
   *   (4) All work arrays are stored in a single workspace structure.
   *
   *   (5) The effective step may be used instead of the search direction
   *       times the step size (useful, e.g., when parameter constraints
   *       are imposed by projections).
   *
   *   (6) It is possible to use an active subset of parameters, e.g. to apply
   *       bound constraints.
   *
   *   (7) Discard pairs for which RHO = S'.Y <= 0 (i.e. H must be positive
   *       definite).
   */

  /* Make an alias for each workspace arrays and get some "constants". */
  local s, y, rho, isave, dsave, v, y_k;
  eq_nocopy, s,     *ws(1);
  eq_nocopy, y,     *ws(2);
  eq_nocopy, isave, *ws(3);
  eq_nocopy, dsave, *ws(4);
  m = numberof(s);
  prj   = isave( 7);
  sftol = dsave(14);
  sgtol = dsave(15);
  sxtol = dsave(16);
  frtol = dsave(23);
  fatol = dsave(24);
  fmin  = dsave(25);

  if (task == "start") {
    /* Check the input arguments for errors.
       Exit if there are errors on input. */
#if 0
    if (n <= 0)      { task = "error: n <= 0";            return; }
    if (m <= 0)      { task = "error: m <= 0";            return; }
    if (frtol < 0.0) { task = "error: frtol < 0";         return; }
    if (fatol < 0.0) { task = "error: fatol < 0";         return; }
#endif
    if (f <= fmin)   { task = "error: initial f <= fmin"; return; }

    /* Initialize local variables. */
    iter = 0; /* number of successful iterations */
    mark = 1; /* index of current direction */
    mp = 0;   /* number of saved directions */

    /* Initialize step information. */
    v = is_void(bound) ? g : (! bound)*g;
    if (noneof(v)) {
      task = "convergence: local minimum found";
      return;
    }
    scale = optim_nrm2(v);
    s(mark) = &((1.0/scale)*v);
    v = [];
    gd = -scale;

    /* Set stage to start the search. */
    stage = 0;
  } else {
    /* Restore local variables.*/
    stage  = isave( 3);
    iter   = isave( 4);
    mark   = isave( 5);
    mp     = isave( 6);
    sftol  = dsave(14);
    sgtol  = dsave(15);
    sxtol  = dsave(16);
    f0     = dsave(17);
    gd     = dsave(18);
    gd0    = dsave(19);
    stp    = dsave(20);
    stpmin = dsave(21);
    stpmax = dsave(22);

    if (stage == 2) {
      /* Compute H_k.g_k by the two-loop recursion.  H_k is the limited
	 memory BFGS approximation of the inverse Hessian, g_k is the
	 gradient at k-th step.  H_k is approximated by using the M last
         pairs (s,y) where:
           s_k = x_k - x_{k-1}  (in fact -s is stored in this program).
           y_k = g_k - g_{k-1}
         The two-loop recursion algorithm is:
           1- start with: v = g_k   (current gradient)
	   2- for k=m-1, ..., 0
	        rho_k = 1/(y_k'.s_k)
	        alpha_k = rho_k s_k'.v
	        v = v - alpha_k y_k
	   3- r = H0_k q
	      e.g. r = (y_k'.s_k)/(y_k'.y_k) q
	   4- for k=0, ..., m-1
	        beta_k = rho_k y_k'.r
	        r = r + (alpha_k - beta_k) s_k
       */
      if ((bnd = ! is_void(bound))) {
        active = double(! bound);
        v = active*g;
      } else {
        eq_nocopy, v, g;
      }
      if (noneof(v)) {
        task = "convergence: local minimum found";
        return;
      }
      alpha = rho = array(0.0, m);
      k = mark + 1;
      gamma = 0.0;
      for (i=1 ; i<=mp ; ++i) {
        if (--k < 1) k = m;
        if (bnd) y_k = (*y(k))*active;
        else eq_nocopy, y_k, *y(k);
        if ((rho(k) = sum((*s(k))*y_k)) > 0.0) {
          v -= (alpha(k) = sum((*s(k))*v)/rho(k))*y_k;
          if (k == mark) gamma = sum(y_k*y_k);
        }
        y_k = [];
      }
      if (gamma > 0.0) {
        /* We also known that RHO(MARK) > 0 in this case. */
        v *= rho(mark)/gamma;
      }
      for (i=1 ; i<=mp ; ++i) {
        if (rho(k) > 0.0) {
          beta = sum((*y(k))*v)/rho(k);
          if (bnd) v += (alpha(k) - beta)*(*s(k))*active;
          else     v += (alpha(k) - beta)*(*s(k));
        }
        if (++k > m) k = 1;
      }
      if (++mark > m) mark = 1;
      s(mark) = &v;
      active = [];

      /* Set stage to initialize the line search. */
      stage = 0;
    }

    /* Compute derivative with respect to step size STP. */
    gd = -sum(g*(*s(mark)));
  }

  if (stage == 0) {
    /* Initialize the line search subroutine. */
    if (gd >= 0.0) {
      /* BFGS recursion yields a search direction which is not a descent.
         Restart the algorithm with only steepest descent. */
      write, "WARNING: not a descent direction (algorithm restarted)";
      mp = 0;
      mark = 1;
      v = is_void(bound) ? g : (! bound)*g;
      if (noneof(v)) {
        task = "convergence: local minimum found";
        return;
      }
      scale = optim_nrm2(v);
      s(mark) = &((1.0/scale)*v);
      v = [];
      gd = -scale;
    }
    f0 = f;
    gd0 = gd;
    stpmin = 0.0;
    stpmax = (fmin - f0)/(sgtol*gd0);
    stp = min(1.0, stpmax);
    copy = g; y(mark) = &copy; // keep a private copy of G
    copy = x; ws(5)   = &copy; // idem for X0 (X at start of line search)
    task = "start search";
    stage = 1;
  }

  if (stage == 1) {
    /* Determine the line search parameter. */
    if (f < fmin) {
      task = "warning: f < fmin";
    } else {
      optim_csrch,f,gd,stp,sftol,sgtol,sxtol,stpmin,stpmax,task,isave,dsave;

      /* Continue if the line search has converged. */
      if ((ident = strpart(task, 1:4) == "conv") ||
          task == "warning: xtol test satisfied") {

        /* Compute the step and gradient change. */
        ++iter;
        y(mark) = &((*y(mark)) - g);
        s(mark) = &(prj ? *ws(5) - x : (*s(mark))*stp);
        if (mp < m) ++mp;
        if (noneof(*s(mark))) {
          task = "convergence: no parameter change";
        } else if (noneof(*y(mark))) {
          task = "convergence: no gradient change";
        } else {
          /* Test for convergence otherwise set TASK to signal a new iterate.
             Set STAGE to compute a new search direction. */
          change = max(abs(f - f0), abs(stp*gd0));
          if (change <= frtol*abs(f0)) {
            task = "convergence: frtol test satisfied";
          } else if (change <= fatol) {
            task = "convergence: fatol test satisfied";
          } else {
            task = "newx";
          }
        }
        stage = 2;
      } else if (ident == "warn") {
        /* Restore solution at start of line search. */
        x = *ws(5);
        f = f0;
        g = *y(mark);
      } else {
        /* Compute the new iterate. */
        x = (*ws(5)) - stp*(*s(mark));
      }
    }
  }

  /* Save local variables (but the constants). */
  isave( 3) = stage;
  isave( 4) = iter;
  isave( 5) = mark;
  isave( 6) = mp;
  dsave(17) = f0;
  dsave(18) = gd;
  dsave(19) = gd0;
  dsave(20) = stp;
  dsave(21) = stpmin;
  dsave(22) = stpmax;
}

func optim_vmlm_setup(m, frtol=, fatol=, fmin=, sftol=, sgtol=, sxtol=,
                      prj=)
/* DOCUMENT optim_vmlm_setup(m)
     Create workspace  array used by variable metric  limited memory (VMLM)
     method.   M  specifies the  number  of  previous  successful steps  to
     remember.  The higher is M the larger will be the memory footprint but
     the more  efficient should  be the  method.  As a  rule of  thumb, for
     small  size problems,  M  should be  of  the order  of  the number  of
     unknowns; but for large size or higly non-quadratic problems, M=3-5 is
     generally sufficient.   The storage is  about 2*N*(1+M)*sizeof(double)
     bytes for large N = number of unknowns.

     Keyword FRTOL set the relative  error desired in the function (default
     value:  FRTOL=1e-8).   Convergence  occurs  if  the  estimate  of  the
     relative  error  between F(X)  and  F(XSOL),  where  XSOL is  a  local
     minimizer,  is  less  than  FRTOL.   FRTOL must  have  a  non-negative
     floating point value.  During the optimization, the value of FRTOL can
     be get/changed with the function optim_vmlm_frtol.

     Keyword FATOL set the absolute  error desired in the function (default
     value: FATOL=0.0).  Convergence occurs if the estimate of the absolute
     error between  F(X) and F(XSOL), where  XSOL is a  local minimizer, is
     less than FATOL.  FATOL must have a non-negative floating point value.
     During the  optimization, the value  of FATOL can be  get/changed with
     the function optim_vmlm_fatol.

     Keyword  FMIN set  a  lower  bound for  the  function (default  value:
     FMIN=0.0).  The subroutine  exits with a warning if  F < FMIN.  During
     the  optimization, the  value  of  FMIN can  be  get/changed with  the
     function optim_vmlm_fmin.

     If value of keyword PRJ  is true (non-nil and non-zero), the effective
     parameter  change is  registered instead  of the  step size  times the
     search  direction.   This  is  intended  for cases  where  the  caller
     modifies the  parameter value(s) given by optim_vmlm  upon return with
     TASK set to "fg" (e.g. bound constraints imposed by projection).

     Keywords SFTOL,  SGTOL, and SXTOL  set tolerances for the  line search
     subroutine (default values: SFTOL=0.001, SGTOL=0.9 and SXTOL=0.1; also
     see optim_csrch).
     
  
   SEE ALSO: optim_vmlm, optim_vmlm_fatol, optim_vmlm_frtol,
             optim_vmlm_fmin. */
{
  /*
   * The workspace array is:
   *   WS(1) = &S
   *   WS(2) = &Y
   *   WS(3) = &ISAVE
   *   WS(4) = &DSAVE
   *   WS(5) = &X0 or nil (parameters at beginning of line search)
   *
   * ISAVE is an array of 6 long's:
   *   ISAVE(1:2) is used by optim_csrch
   *   ISAVE(3) = STAGE
   *   ISAVE(4) = ITER
   *   ISAVE(5) = MARK
   *   ISAVE(6) = MP
   *   ISAVE(7) = PRJ
   *
   * DSAVE is an array of 26 double's:
   *   DSAVE(1:13) is used by optim_csrch
   *   DSAVE(14) = SFTOL
   *   DSAVE(15) = SGTOL
   *   DSAVE(16) = SXTOL
   *   DSAVE(17) = F0
   *   DSAVE(18) = GD
   *   DSAVE(19) = GD0
   *   DSAVE(20) = STP
   *   DSAVE(21) = STPMIN
   *   DSAVE(22) = STPMAX
   *   DSAVE(23) = FRTOL
   *   DSAVE(24) = FATOL
   *   DSAVE(25) = FMIN
   *
   * S is pointer array with M elements used to store the previous search
   *   directions.
   *
   * Y is pointer array with M elements used to store the previous gradient
   *   changes.
   *
   * RHO is a double precision work array of dimension M.
   *
   * X0  is a double precision work  array with same dimension  list as X
   *   used  to  store X  at the  beginning of  the line  search  (can be
   *   undefined on first entry).
   */
  if (is_void(m) || m <= 0) error, "M must be greater or equal 1";
  if (is_void(frtol)) frtol = 1e-8;
  else if (frtol < 0.0) error, "FRTOL must be non-negative";
  if (is_void(fatol)) fatol = 0.0;
  else if (fatol < 0.0) error, "FATOL must be non-negative";

  if (is_void(sftol)) sftol = 0.001;
  if (is_void(sgtol)) sgtol = 0.9;
  if (is_void(sxtol)) sxtol = 0.1;

  isave = array(long, 7);
  isave(7) = (prj ? 1 : 0);
  
  dsave = array(double, 25);
  dsave(14) = sftol;
  dsave(15) = sgtol;
  dsave(16) = sxtol;
  dsave(23) = frtol;
  dsave(24) = fatol;
  dsave(25) = (is_void(fmin) ? 0.0 : fmin);

  ws = array(pointer, 5);
  ws(1) = &array(pointer, m); // S
  ws(2) = &array(pointer, m); // Y
  ws(3) = &isave;             // ISAVE
  ws(4) = &dsave;             // DSAVE
  return ws;
}

func optim_vmlm_frtol(ws, newvalue)
/* DOCUMENT optim_vmlm_frtol(ws)
       -or- optim_vmlm_frtol(ws, frtol)
     Query/change  value  of parameter  FRTOL  in  variable metric  limited
     memory  (VMLM)  method.   WS   is  the  workspace  array  returned  by
     optim_vmlm_setup and used in  optim_vmlm.  FRTOL is the relative error
     desired in  the function:  convergence occurs if  the estimate  of the
     relative  error  between F(X)  and  F(XSOL),  where  XSOL is  a  local
     minimizer,  is less  than FRTOL.   The return  value is  the parameter
     value prior to the change if any.

   SEE ALSO: optim_vmlm_setup, optim_vmlm. */
{
  oldvalue = (*ws(4))(23);
  if (! is_void(newvalue)) {
    if (newvalue < 0.0) error, "FRTOL must be non-negative";
    (*ws(4))(23) = newvalue;
  }
  return oldvalue;
}

func optim_vmlm_fatol(ws, newvalue)
/* DOCUMENT optim_vmlm_fatol(ws)
       -or- optim_vmlm_fatol(ws, fatol)
     Query/change  value  of parameter  FATOL  in  variable metric  limited
     memory  (VMLM)  method.   WS   is  the  workspace  array  returned  by
     optim_vmlm_setup and used in  optim_vmlm.  FATOL is the absolute error
     desired in  the function:  convergence occurs if  the estimate  of the
     absolute  error  between F(X)  and  F(XSOL),  where  XSOL is  a  local
     minimizer,  is less  than FATOL.   The return  value is  the parameter
     value prior to the change if any.

   SEE ALSO: optim_vmlm_setup, optim_vmlm. */
{
  oldvalue = (*ws(4))(24);
  if (! is_void(newvalue)) {
    if (newvalue < 0.0) error, "FATOL must be non-negative";
    (*ws(4))(24) = newvalue;
  }
  return oldvalue;
}

func optim_vmlm_fmin(ws, newvalue)
/* DOCUMENT optim_vmlm_fmin(ws)
       -or- optim_vmlm_fmin(ws, fmin)
     Query/change value of parameter FMIN in variable metric limited memory
     (VMLM) method.  WS is the workspace array returned by optim_vmlm_setup
     and used  in optim_vmlm.  FMIN  is a lower  bound for the  function to
     minimize.  The return value is the parameter value prior to the change
     if any.

   SEE ALSO: optim_vmlm_setup, optim_vmlm. */
{
  oldvalue = (*ws(4))(25);
  if (! is_void(newvalue)) (*ws(4))(25) = newvalue;
  return oldvalue;
}

func optim_vmlm_step(ws) { return (*ws(4))(20); }
func optim_vmlm_iter(ws) { return (*ws(3))(4); }
/* DOCUMENT optim_vmlm_step(ws) -- current step lenght
       -or- optim_vmlm_iter(ws) -- number of iterations so far
     Query information about the  state/progress of variable metric limited
     memory  (VMLM)   method;  WS  is  the  workspace   array  returned  by
     optim_vmlm_setup and used in optim_vmlm.

   SEE ALSO: optim_vmlm_setup, optim_vmlm. */

/*---------------------------------------------------------------------------*/
/* CONJUGATE GRADIENT METHOD WITH OPTIONAL PRECONDITIONING */

func optim_cgmn(&x, f, g, &task, ws, h)
/* DOCUMENT optim_cgmn, x, f, g, task, ws;
       -or- optim_cgmn, x, f, g, task, ws, h;

     Find a  local minimizer of  a function of  N variables by  a conjugate
     gradient  method.  This  subroutine uses  reverse  communication.  The
     user must choose an initial approximation X to the minimizer, evaluate
     the function  and the gradient  at X, and  make the initial  call with
     TASK set to "start".  On exit TASK indicates the required action.

     A typical invocation of optim_cgmn has the following outline:

     |  x = ... ;                                        // starting vector
     |  ws = optim_cgmn_setup(method, fmin=lower_bound); // allocate workspace
     |  task = "start";
     |  for (;;) {
     |    if (task == "fg" || task == "start") {
     |      fx = ...;  // evaluate the function at X; store in FX
     |      gx = ...;  // evaluate the gradient of F at X; store in GX
     |    } else if (task == "newx") {
     |       // New successful step: the approximation X, function FX, and
     |       // gradient GX, are available for inspection. Optionally compute
     |       // preconditioned gradient HX.
     |    } else {
     |      break;  // convergence, or error, or warning
     |    }
     |    optim_cgmn, x, fx, gx, task, ws, hx;
     |  }

     The subroutine arguments are:
  
       X   is a double precision array.  On entry, X is an approximation to
           the solution.  On exit, X is the current approximation.

       F   is a double precision variable. On  entry, F is the value of the
           function at X.  On final exit, F is the function value at X.

       G   is a double  precision array with same dimension  list as X.  On
           entry, G is  the value of the gradient at X.   On final exit, GX
           is the value of the gradient at X.

       TASK is  a string variable.  On  initial entry, TASK must  be set to
           "start".  On return of  the routine, TASK indicates the required
           action:
           - If TASK  = "fg" then evaluate  the function and  gradient at X
             and call optim_cgmn again.             
           - If  strpart(TASK, 1:4)  = "newx" then  a new iterate  has been
             computed.  The  approximation X,  function F, and  gradient GX
             are  available   for  examination.   Optionally   compute  the
             preconditioned   conjugate   gradient  H   at   X   (i.e.   an
             approximation of the inverse  of the Hessian matrix applied to
             the gradient).             
           - If strpart(TASK, 1:4) = "conv" then the search is successful.
           - If strpart(TASK, 1:4) = "warn" then the subroutine is not able
             to satisfy  the convergence conditions.   The exit value  of X
             contains the best approximation found so far.
           - If strpart(TASK, 1:5) = "error"  then there is an error in the
             input arguments.
           On exit  with convergence, a  warning or an error,  the variable
           TASK  contains additional  information.  Set  TASK to  "stop" in
           order to recover the values of  X, F and G for the best solution
           found so far.
           
       WS  is a workspace array that can be obtained from optim_cgmn_setup.
           The user must not alter  this workspace between calls (except by
           using optim_cgmn_fatol, optim_cgmn_frtol, or optim_cgmn_fmin).

       H   is an  optional double precision array with  same dimension list
           as  X which  is only  required when  TASK="newx".  If  H  is not
           specified,  an ordinary conjugate  gradient search  direction is
           used  for  the next  step;  otherwise, H  is  assumed  to be  an
           approximation   of   the   preconditioned  gradient   (i.e.   an
           approximation of  the inverse of  the Hessian matrix  applied to
           the gradient) and is used to derive a better search direction.
       
   SEE ALSO: optim_cgmn_setup, optim_csrch, optim_vmlm. */
{
  local x0, g0, s, isave, dsave, t;

  /* Restore some constants. */
  eq_nocopy, isave, *ws(4);
  eq_nocopy, dsave, *ws(5);
  sftol  = dsave(14);
  sgtol  = dsave(15);
  sxtol  = dsave(16);
  frtol  = dsave(22);
  fatol  = dsave(23);
  fmin   = dsave(24);
  method = isave(6);

  /* Set/restore local variables. */
  if (task == "start") {
    stage = 0;
    iter  = 0;
    eval  = 1;
    stp = 0.0;
  } else {
    stage  = isave(3);
    if (stage <= 0 || stage > 2) {
      isave(3) = -1;
      task = "error: corrupted workspace (unexpected STAGE)";
      return;
    }
    iter   = isave(4);
    eval   = isave(5);
    f0     = dsave(17);
    sg0    = dsave(18);
    stp    = dsave(19);
    stpmin = dsave(20);
    stpmax = dsave(21);
    eq_nocopy, x0, *ws(1);
    eq_nocopy, g0, *ws(2);
    eq_nocopy, s,  *ws(3);
  }
  
  if (stage == 1) {
    /* Line search in progress. */
    ++eval;
    sg = sum(s*g);
  } else {
    /* Either first entry or previous line search has converged. */

    /* Get various method flags. */
    use_effective_step = (method & 4);
    shanno_pua = ! (method & 8);
    method &= 3;
    
    /* The next search direction is estimated by the conjugate gradient
     * recursion (generalized Polak-Ribiere formula):
     *
     *             (G - G0)'.H
     *   S1 = H + ------------- S0           (use previous search dir.)
     *                 G0'.S0
     *
     *             (G - G0)'.H
     *      = H + ------------- (X - X0)    (use "effective" step)
     *             G0'.(X - X0)
     *
     *   X0 = parameters at previous step
     *   G0 = gradient at X0
     *   S0 = previous search direction
     *
     *   X = parameters at current step
     *   G = gradient at X
     *   S = next search direction
     *   H = approximation of preconditioned anti-gradient at X
     *       (can be -G if no-preconditioning).  Also note that
     *       one should replace H by -H if H is an ascent direction.
     */

    /* Fix preconditioned anti-gradient: use anti-gradient if H not
       provided, otherwise make sure H is a descent direction. */
    if (is_void(h)) {
      h = -g;
      hg = -sum(g*g);
    } else {
      hg = sum(h*g);
      if (hg > 0.0) {
        h = -h;
        hg = -hg;
      }
    }
    
    if (stage == 2) {
      /* Use conjugate gradient recursion to compute new search direction. */
      if (! use_effective_step || ! (tg0 = sum((t = x - x0)*g0))) {
        eq_nocopy, t, s;
        tg0 = sg0;
      }
      if (method == 1) {
        /* Fletcher-Reeves formula. */
        beta = hg / tg0;
      } else {
        /* Polak-Ribiere formula. */
        beta = sum((g - g0)*h) / tg0;
        if (method == 3 && beta < 0.0) {
          /* Reset recursion, same as using BETA = 0 (see below). */
          stage = 0;
        }
      }
      if (stage == 2) {
        s = h + beta*t;
        sg = sum(s*g);
        if (sg >= 0.0) {
          /* Discard non-descent direction produced by conjugate
             gradient recursion (see below). */
          stage = 0;
        }
      }
    }
    if (stage == 0) {
      /* Simply use (preconditioned) anti-gradient as search direction. */
      s = h; /* do _not_ use eq_nocopy here */
      sg = hg;
    }
    h = []; /* no longer needed */
    ws(3) = &s;
    if (sg >= 0.0) {
      if (noneof(g)) {
        task = "convergence: local minimum found";
        isave(3) = 2;
      } else {
        task = "error: search direction perpendicular to gradient";
        isave(3) = -1;
      }
      return;
    }

    /* Initialize for line search. */
    task = "start";
    if (stage == 0) {
      /* FIXME: Use a step proportionnal to 1/|g|^2 (so that the problem is
         invariant by a simple parameter change). */
      scale = optim_nrm2(g);
      // stp = sum(abs(x))/sum(abs(s))*1e-3;
      stp = 1.0/scale;
    } else if (shanno_pua) {
      /* Use Shanno-Phua's formula for length of trial step along the new
         search direction. */
      stp = stp*sg0/sg;
    } else {
      /* Try with previous sucessful step size. */
      stp = stp;
    }
    stpmin = 0.0;
    stpmax = (fmin - f)/(sgtol*sg);
    //write, stp, stpmax, (fmin - f)/sg;
    stp = min(stp, stpmax);

    /* Save variables. */
    x0 = x; ws(1) = &x0;
    g0 = g; ws(2) = &g0;
    f0 = f;
    sg0 = sg;
  }

  /* Call line search routine. */
  optim_csrch,f,sg,stp,sftol,sgtol,sxtol,stpmin,stpmax,task,isave,dsave;
  if (task == "fg") {
    stage = 1;
    x = x0 + stp*s;
  } else {
    ident = strpart(task, 1:4);
    conv = (ident == "conv");
    if (conv || ident == "warn") {
      df = max(f0 - f, stp*abs(sg0));
      if (df <= frtol*abs(f0)) {
        task = "convergence: frtol test satisfied";
      } else if (df <= fatol) {
        task = "convergence: fatol test satisfied";
      } else if (conv) {
        task = "newx";
      }
    
      /* Line search has converged or line search cannot make further
         progress: request preconditioned (anti-)gradient if the caller
         wish to continue. */
      stage = 2; /* expect H on next entry */
      ++iter;
    } else {
      /* Some error occured. */
      stage = -1;
    }
  }

  /* Save local variables. */
  isave(3) = stage;
  isave(4) = iter;
  isave(5) = eval;
  dsave(17) = f0;
  dsave(18) = sg0;
  dsave(19) = stp;
  dsave(20) = stpmin;
  dsave(21) = stpmax;
}

func optim_cgmn_setup(method, frtol=, fatol=, fmin=, sftol=, sgtol=, sxtol=)
/* DOCUMENT optim_cgmn_setup(method)
     Create workspace array used by conjugate gradient optimization method.
     The argument METHOD is a positive integer:
       1   for Fletcher-Reeves method
       2   for Polak-Ribiere method
       3   for Polak-Ribiere method with BETA>0  
     optionally, a value  of 4 can be added to METHOD  to use the effective
     step (instead of the previous search direction).
     
     Keyword FRTOL set the relative  error desired in the function (default
     value:  FRTOL=1e-8).   Convergence  occurs  if  the  estimate  of  the
     relative  error  between F(X)  and  F(XSOL),  where  XSOL is  a  local
     minimizer,  is  less  than  FRTOL.   FRTOL must  have  a  non-negative
     floating point value.  During the optimization, the value of FRTOL can
     be get/changed with the function optim_cgmn_frtol.

     Keyword FATOL set the absolute  error desired in the function (default
     value: FATOL=0.0).  Convergence occurs if the estimate of the absolute
     error between  F(X) and F(XSOL), where  XSOL is a  local minimizer, is
     less than FATOL.  FATOL must have a non-negative floating point value.
     During the  optimization, the value  of FATOL can be  get/changed with
     the function optim_cgmn_fatol.

     Keyword  FMIN set  a  lower  bound for  the  function (default  value:
     FMIN=0.0).  The subroutine  exits with a warning if  F < FMIN.  During
     the  optimization, the  value  of  FMIN can  be  get/changed with  the
     function optim_cgmn_fmin.

     Keywords SFTOL,  SGTOL, and SXTOL  set tolerances for the  line search
     subroutine (default values: SFTOL=0.001, SGTOL=0.9 and SXTOL=0.1; also
     see optim_csrch).
     
  
   SEE ALSO: optim_cgmn, optim_cgmn_fatol, optim_cgmn_frtol,
             optim_cgmn_fmin, optim_csrch. */
{
  /* Workspace arrays:
   *   *WS(1) = X0 parameters at start of line search
   *   *WS(2) = G0 gradient at X0
   *   *WS(3) = S  search direction
   *   *WS(4) = ISAVE array of local integer variables
   *   *WS(5) = DSAVE array of local floating point variables
   *   ISAVE(1:2)  = used by optim_csrch
   *   DSAVE(1:13) = used by optim_csrch
   */

  if (is_void(frtol)) frtol = 1e-8;
  else if (frtol < 0.0) error, "FRTOL must be non-negative";
  if (is_void(fatol)) fatol = 0.0;
  else if (fatol < 0.0) error, "FATOL must be non-negative";
  
  if (is_void(sftol)) sftol = 1e-4;  // 0.001;
  if (is_void(sgtol)) sgtol = 1e-1;  // 0.9;
  if (is_void(sxtol)) sxtol = 1e-17; // 0.1;
  
  isave = array(long, 6);
  isave(6) = (is_void(method) ? 7 : method);
  
  dsave = array(double, 24);
  dsave(14) = sftol;
  dsave(15) = sgtol;
  dsave(16) = sxtol;
  dsave(22) = frtol;
  dsave(23) = fatol;
  dsave(24) = (is_void(fmin) ? 0.0 : fmin);
  
  ws = array(pointer, 5);
  ws(4) = &isave;
  ws(5) = &dsave;
  return ws;
}

func optim_cgmn_frtol(ws, newvalue)
/* DOCUMENT optim_cgmn_frtol(ws)
       -or- optim_cgmn_frtol(ws, frtol)
     Query/change   value  of   parameter  FRTOL   in   conjugate  gradient
     optimization  method.    WS  is   the  workspace  array   returned  by
     optim_cgmn_setup and used in  optim_cgmn.  FRTOL is the relative error
     desired in  the function:  convergence occurs if  the estimate  of the
     relative  error  between F(X)  and  F(XSOL),  where  XSOL is  a  local
     minimizer,  is less  than FRTOL.   The return  value is  the parameter
     value prior to the change if any.

   SEE ALSO: optim_cgmn_setup, optim_cgmn. */
{
  oldvalue = (*ws(5))(22);
  if (! is_void(newvalue)) {
    if (newvalue < 0.0) error, "FRTOL must be non-negative";
    (*ws(5))(22) = newvalue;
  }
  return oldvalue;
}

func optim_cgmn_fatol(ws, newvalue)
/* DOCUMENT optim_cgmn_fatol(ws)
       -or- optim_cgmn_fatol(ws, fatol)
     Query/change   value  of   parameter  FATOL   in   conjugate  gradient
     optimization  method.    WS  is   the  workspace  array   returned  by
     optim_cgmn_setup and used in  optim_cgmn.  FATOL is the absolute error
     desired in  the function:  convergence occurs if  the estimate  of the
     absolute  error  between F(X)  and  F(XSOL),  where  XSOL is  a  local
     minimizer,  is less  than FATOL.   The return  value is  the parameter
     value prior to the change if any.

   SEE ALSO: optim_cgmn_setup, optim_cgmn. */
{
  oldvalue = (*ws(5))(23);
  if (! is_void(newvalue)) {
    if (newvalue < 0.0) error, "FATOL must be non-negative";
    (*ws(5))(23) = newvalue;
  }
  return oldvalue;
}

func optim_cgmn_fmin(ws, newvalue)
/* DOCUMENT optim_cgmn_fmin(ws)
       -or- optim_cgmn_fmin(ws, fmin)
     Query/change   value   of  parameter   FMIN   in  conjugate   gradient
     optimization  method.    WS  is   the  workspace  array   returned  by
     optim_cgmn_setup and  used in optim_cgmn.   FMIN is a lower  bound for
     the function  to minimize.   The return value  is the  parameter value
     prior to the change if any.

   SEE ALSO: optim_cgmn_setup, optim_cgmn. */
{
  oldvalue = (*ws(5))(24);
  if (! is_void(newvalue)) (*ws(5))(24) = newvalue;
  return oldvalue;
}

func optim_cgmn_step(ws) { return (*ws(5))(19); }
func optim_cgmn_iter(ws) { return (*ws(4))(4); }
func optim_cgmn_eval(ws) { return (*ws(4))(5); }
/* DOCUMENT optim_cgmn_step(ws) -- best step lenght so far
       -or- optim_cgmn_iter(ws) -- number of iterations so far
       -or- optim_cgmn_eval(ws) -- number of function evaluations so far
     Query  information  about  the  state/progress of  conjugate  gradient
     optimization   method;  WS   is  the   workspace  array   returned  by
     optim_cgmn_setup and used in optim_cgmn.

   SEE ALSO: optim_cgmn_setup, optim_cgmn. */

/*---------------------------------------------------------------------------*/
/* LINE SEARCH ROUTINES */

func optim_csrch(f, g, &stp, ftol, gtol, xtol, stpmin, stpmax, &task, isave,
                 dsave)
/* DOCUMENT optim_csrch, f, g, stp, ftol, gtol, xtol, stpmin, stpmax, task,
 *                       isave, dsave;
 *
 *   This subroutine  finds a  step that  satisfies  a  sufficient decrease
 *   condition and a curvature condition.
 *
 *   Each call of the subroutine updates an interval with endpoints stx and
 *   sty. The interval is initially  chosen so that it contains a minimizer
 *   of the modified function
 *
 *         psi(stp) = f(stp) - f(0) - ftol*stp*f'(0).
 *
 *   If psi(stp) <= 0 and f'(stp)  >= 0 for some step, then the interval is
 *   chosen  so  that  it contains  a  minimizer of  f.   The algorithm  is
 *   designed  to find  a  step  that  satisfies  the  sufficient  decrease
 *   condition
 *
 *         f(stp) <= f(0) + ftol*stp*f'(0),
 *
 *   and the curvature condition
 *
 *         abs(f'(stp)) <= gtol*abs(f'(0)).
 *
 *   If ftol is less than gtol and if, for example, the function is bounded
 *   below,  then there is always  a step which  satisfies both conditions.
 *   If  no  step can  be found  that satisfies  both conditions,  then the
 *   algorithm stops  with a warning.  In this case  stp only satisfies the
 *   sufficient decrease condition.
 *
 *   A typical invocation of dcsrch has the following outline:
 *
 *   |  task = "start";
 *   |  stp = ...; // Initial STP value
 *   |  for (;;) {
 *   |    if (task == "fg" || task == "start") {
 *   |      // Evaluate the function and the gradient at STP.
 *   |    } else if (strpart(task, 1:4) == "conv") {
 *   |      // Search has converged.
 *   |    } else if (strpart(task, 1:4) == "warn") {
 *   |      // Some problem prevents further progress.
 *   |    } else {
 *   |      // An error occured.
 *   |    }
 *   |    optim_csrch, ...;
 *   |  }
 *
 *   NOTE: The user must not alter work arrays between calls.
 *
 *   The subroutine arguments are:
 *
 *     F is a double precision  variable.  On initial entry, F is the value
 *       of  the function at 0.  On  subsequent entries, F is  the value of
 *       the function at STP.  On exit, F is left unchanged.
 *
 *     G is  a double  precision  variable.  On  initial entry,  G  is  the
 *       derivative of the function at  0.  On subsequent entries, G is the
 *       derivative of the function at STP.  On exit, G is left unchanged.
 *
 *    STP is a double  precision  variable.  On entry, STP  is the  current
 *       estimate  of a  satisfactory step.  On  initial entry,  a positive
 *       initial  estimate must be  provided.   On exit with TASK="fg", STP
 *       is  the  new  estimate  of a  satisfactory  step.   On  exit  with
 *       strpart(TASK,1:3)="conv",   STP is  left  unchanged and  satisfies
 *       the sufficient  decrease and curvature  condition.   On exit  with
 *       TASK not equal to "fg", STP is left unchanged.
 *
 *     FTOL  is a  double precision variable.   On entry, FTOL  specifies a
 *       nonnegative  tolerance for the sufficient  decrease condition.  On
 *       exit, FTOL is unchanged.
 *
 *     GTOL  is a  double precision  variable.  On entry, GTOL  specifies a
 *       nonnegative tolerance for  the curvature condition.  On exit, GTOL
 *       is unchanged.
 *
 *     XTOL  is a  double precision variable.   On entry, XTOL  specifies a
 *       nonnegative  relative  tolerance  for  an  acceptable  step.   The
 *       subroutine exits with a warning if the relative difference between
 *       STY and STX is less than XTOL.  On exit, XTOL is unchanged.
 *
 *     STPMIN  is  a  double precision  variable.   On entry,  STPMIN is  a
 *       nonnegative  lower  bound  for  the step.   On  exit,  STPMIN   is
 *       unchanged.
 *
 *     STPMAX is  a  double precision   variable.   On entry,  STPMAX is  a
 *       nonnegative   upper  bound  for  the step.   On  exit,  STPMAX  is
 *       unchanged.
 *
 *     TASK is a character variable of length at least 60.
 *       On initial entry, task must be set to "start".
 *       On exit, TASK indicates the required action:
 *
 *          If task(1:2) = "FG" then evaluate the function and
 *          derivative at stp and call dcsrch again.
 *
 *          If task(1:4) = "CONV" then the search is successful.
 *
 *          If task(1:4) = "WARN" then the subroutine is not able
 *          to satisfy the convergence conditions. The exit value of
 *          stp contains the best point found during the search.
 *
 *          If task(1:5) = "ERROR" then there is an error in the
 *          input arguments.
 *
 *       On exit with convergence, a warning or an error, the variable TASK
 *       contains additional information.
 *
 *     ISAVE is an integer work array of, at least, 2 elements.
 *
 *     DSAVE is a double precision work array of, at least, 13 elements.
 *
 *
 * HISTORY:
 *   MINPACK-1 Project. June 1983.
 *   Argonne National Laboratory.
 *   Jorge J. More' and David J. Thuente.
 *
 *   MINPACK-2 Project. November 1993.
 *   Argonne National Laboratory and University of Minnesota.
 *   Brett M. Averick, Richard G. Carter, and Jorge J. More'.
 *
 *   Yorick translation an improvements.  October 2001.
 *   Observatoire de Lyon.
 *   Eric Thiebaut.
 *
 *
 * SEE ALSO: optim_cstep, optim_cgmn, optim_vmlm.
 */
{
  /* Initialization block.*/
  zero = 0.0;
  xtrapl = 1.1;
  xtrapu = 4.0;

  if (strpart(task, 1:5) == "start") {

    /* Check the input arguments for errors.
       Exit if there are errors on input. */
    if (stpmax < stpmin) { task = "error: stpmax < stpmin";    return; }
    if (stpmin < zero)   { task = "error: stpmin < zero";      return; }
    if (xtol < zero)     { task = "error: xtol < zero";        return; }
    if (gtol < zero)     { task = "error: gtol < zero";        return; }
    if (ftol < zero)     { task = "error: ftol < zero";        return; }
    if (g >= zero)       { task = "error: initial g >= zero";  return; }
    if (stp > stpmax)    { task = "error: stp > stpmax";       return; }
    if (stp < stpmin)    { task = "error: stp < stpmin";       return; }

    /* Initialize local variables.
       The variables stx, fx, gx contain the values of the step,
       function, and derivative at the best step.
       The variables sty, fy, gy contain the value of the step,
       function, and derivative at sty.
       The variables stp, f, g contain the values of the step,
       function, and derivative at stp. */
    brackt = 0;
    stage = 1;
    finit = f;
    ginit = g;
    gtest = ftol*ginit;
    width = stpmax - stpmin;
    width1 = 2.0*width;
    stx = zero;
    fx = finit;
    gx = ginit;
    sty = zero;
    fy = finit;
    gy = ginit;
    stmin = zero;
    stmax = stp + xtrapu*stp;
    task = "fg";

  } else {

    /* Restore local variables. */
    brackt = isave( 1);
    stage  = isave( 2);
    ginit  = dsave( 1);
    gtest  = dsave( 2);
    gx     = dsave( 3);
    gy     = dsave( 4);
    finit  = dsave( 5);
    fx     = dsave( 6);
    fy     = dsave( 7);
    stx    = dsave( 8);
    sty    = dsave( 9);
    stmin  = dsave(10);
    stmax  = dsave(11);
    width  = dsave(12);
    width1 = dsave(13);

    /* If psi(stp) <= 0 and f'(stp) >= 0 for some step, then the
       algorithm enters the second stage.*/
    ftest = finit + stp*gtest;
    if (stage == 1 && f <= ftest && g >= zero) stage = 2;

    /* Test for termination: convergence or warnings. */
    if (f <= ftest && abs(g) <= gtol*(-ginit)) {
      task = "convergence";
    } else if (stp == stpmin && (f > ftest || g >= gtest)) {
      task = "warning: stp = stpmin";
    } else if (stp == stpmax && f <= ftest && g <= gtest) {
      task = "warning: stp = stpmax";
    } else if (brackt && stmax - stmin <= xtol*stmax) {
      task = "warning: xtol test satisfied";
    } else if (brackt && (stp <= stmin || stp >= stmax)) {
      task = "warning: rounding errors prevent progress";
    } else {
      /* A modified function is used to predict the step during the first
         stage if a lower function value has been obtained but the decrease
         is not sufficient. */
      if (stage == 1 && f <= fx && f > ftest) {
        /* Define the modified function and derivative values. */
        fm  = f  - stp*gtest;
        fxm = fx - stx*gtest;
        fym = fy - sty*gtest;
        gm  = g  - gtest;
        gxm = gx - gtest;
        gym = gy - gtest;

        /* Call dcstep to update stx, sty, and to compute the new step. */
        optim_cstep,stx,fxm,gxm,sty,fym,gym,stp,fm,gm,brackt,stmin,stmax;

        /* Reset the function and derivative values for f. */
        fx = fxm + stx*gtest;
        fy = fym + sty*gtest;
        gx = gxm + gtest;
        gy = gym + gtest;

      } else {
        /* Call dcstep to update stx, sty, and to compute the new step. */
        optim_cstep,stx,fx,gx,sty,fy,gy,stp,f,g,brackt,stmin,stmax;
      }

      if (brackt) {
        /* Decide if a bisection step is needed. */
        if ((wcur = abs(sty - stx)) >= 0.66*width1) stp = stx + 0.5*(sty - stx);
        width1 = width;
        width = wcur;

        /* Set the minimum and maximum steps allowed for stp. */
        stmin = min(stx, sty);
        stmax = max(stx, sty);
      } else {
        /* Set the minimum and maximum steps allowed for stp. */
        stmin = stp + xtrapl*(stp - stx);
        stmax = stp + xtrapu*(stp - stx);
      }

      /* Force the step to be within the bounds stpmax and stpmin. */
      stp = min(max(stp, stpmin), stpmax);

      /* If further progress is not possible, let stp be the best
         point obtained during the search. */
      if (brackt && (stp <= stmin || stp >= stmax)
          || (brackt && stmax-stmin <= xtol*stmax)) stp = stx;

      /* Obtain another function and derivative. */
      task = "fg";
    }
  }

  /* Save local variables. */
  isave( 1) = brackt;
  isave( 2) = stage;
  dsave( 1) = ginit;
  dsave( 2) = gtest;
  dsave( 3) = gx;
  dsave( 4) = gy;
  dsave( 5) = finit;
  dsave( 6) = fx;
  dsave( 7) = fy;
  dsave( 8) = stx;
  dsave( 9) = sty;
  dsave(10) = stmin;
  dsave(11) = stmax;
  dsave(12) = width;
  dsave(13) = width1;
}

func optim_cstep(&stx, &fx, &dx, &sty, &fy, &dy, &stp, fp, dp,
                 &brackt, stpmin, stpmax)
/* DOCUMENT optim_cstep, stx, fx, dx,
 *                       sty, fy, dy,
 *                       stp, fp, dp,
 *                       brackt, stpmin, stpmax;
 *
 *   This subroutine computes a safeguarded step for a search procedure and
 *   updates an  interval that contains a step  that satisfies a sufficient
 *   decrease and a curvature condition.
 *
 *   The parameter STX contains the  step with the least function value. If
 *   brackt  is set  to .true. then  a minimizer  has been bracketed  in an
 *   interval with  endpoints STX and sty.  The  parameter STP contains the
 *   current step.  The subroutine assumes that if BRACKT is true then:
 *
 *         min(STX,STY) < STP < max(STX,STY),
 *
 *   and  that the derivative  at STX is  negative in the direction  of the
 *   step.
 *
 *   The subroutine arguments are:
 *
 *     STX is a double precision  variable.  On entry, STX is the best step
 *       obtained so  far and is an endpoint of  the interval that contains
 *       the minimizer.  On exit, STX is the updated best step.
 *
 *     FX is a double precision variable.  On entry, FX  is the function at
 *       STX.  On exit, FX is the function at STX.
 *
 *     DX is  a double precision variable.  On entry,  DX is the derivative
 *       of  the function at STX.   The derivative must be  negative in the
 *       direction  of  the  step,  that is,  DX  and  STP - STX must  have
 *       opposite signs.  On exit, DX  is the derivative of the function at
 *       STX.
 *
 *     STY  is a double  precision variable.  On  entry, STY is  the second
 *       endpoint  of the interval  that contains the minimizer.   On exit,
 *       STY is the  updated  endpoint  of the  interval that  contains the
 *       minimizer.
 *
 *     FY is a double precision  variable.  On entry, FY is the function at
 *       STY.  On exit, FY is the function at STY.
 *
 *     DY is a double precision variable.
 *       On entry, DY is the derivative of the function at STY.
 *       On exit, DY is the derivative of the function at the exit STY.
 *
 *     STP  is a double precision  variable.  On entry, STP  is the current
 *       step. If  BRACKT is true, then, on input, STP  must be between STX
 *       and STY.  On exit, STP is a new trial step.
 *
 *     FP is a double precision  variable.  On entry, FP is the function at
 *       STP.  On exit, FP is unchanged.
 *
 *     DP is  a double precision variable.  On entry,  DP is the derivative
 *       of the function at STP.  On exit, DP is unchanged.
 *
 *     BRACKT is  a logical  variable.  On  entry, BRACKT  specifies  if  a
 *       minimizer  has been  bracketed.  Initially BRACKT  must be  set to
 *       FALSE.    On  exit,  BRACKT  specifies if  a  minimizer  has  been
 *       bracketed.  When a minimizer is bracketed brackt is set to TRUE.
 *
 *     STPMIN is a double precision  variable.  On entry, STPMIN is a lower
 *       bound for the step.  On exit, STPMIN is unchanged.
 *
 *     STPMAX is a double precision variable.  On entry, STPMAX is an upper
 *       bound for the step.  On exit, STPMAX is unchanged.
 *
 *
 * HISTORY:
 *   MINPACK-1 Project. June 1983
 *   Argonne National Laboratory.
 *   Jorge J. More' and David J. Thuente.
 *
 *   MINPACK-2 Project. November 1993.
 *   Argonne National Laboratory and University of Minnesota.
 *   Brett M. Averick and Jorge J. More'.
 *
 *   Yorick translation an improvements.  October 2001.
 *   Observatoire de Lyon.
 *   Eric Thiebaut.
 *
 *
 * SEE ALSO: optim_csrch.
 */
{
#if 0
  /* Check the input parameters for errors. */
  if ((brackt && (stx < sty ? (stp <= stx || stp >= sty)
		            : (stp >= stx || stp <= sty)))) {
    error, "STP outside bracket (STX,STY)";
  } else if (dx*(stp - stx) >= 0.0) {
    error, "descent condition violated";
  } else if (stpmax < stpmin) {
    error, "STPMAX < STPMIN";
  }
#endif

  /* Determine if the derivatives have opposite sign. */
  sgnd = dp*(dx/abs(dx));

  if (fp > fx) {
    /* First case: A higher function value. The minimum is bracketed.  If
       the cubic step is closer to stx than the quadratic step, the cubic
       step is taken, otherwise the average of the cubic and quadratic
       steps is taken. */
    theta = 3.0*(fx - fp)/(stp - stx) + dx + dp;
    s = max(abs(theta), abs(dx), abs(dp));
    gamma = s*sqrt((theta/s)^2 - (dx/s)*(dp/s));
    if (stp < stx) gamma = -gamma;
    p = (gamma - dx) + theta;
    q = ((gamma - dx) + gamma) + dp;
    r = p/q;
    stpc = stx + r*(stp - stx);
    stpq = stx + ((dx/((fx - fp)/(stp - stx) + dx))/2.0)*(stp - stx);
    if (abs(stpc - stx) < abs(stpq - stx)) {
      stpf = stpc;
    } else {
      /* stpf = (stpq + stpc)/2.0; */
      stpf = stpc + (stpq - stpc)/2.0;
    }
    brackt = 1;
  } else if (sgnd < 0.0) {
    /* Second case: A lower function value and derivatives of opposite
       sign. The minimum is bracketed. If the cubic step is farther from
       stp than the secant step, the cubic step is taken, otherwise the
       secant step is taken. */
    theta = 3.0*(fx - fp)/(stp - stx) + dx + dp;
    s = max(abs(theta),abs(dx),abs(dp));
    gamma = s*sqrt((theta/s)^2 - (dx/s)*(dp/s));
    if (stp > stx) gamma = -gamma;
    p = (gamma - dp) + theta;
    q = ((gamma - dp) + gamma) + dx;
    r = p/q;
    stpc = stp + r*(stx - stp);
    stpq = stp + (dp/(dp - dx))*(stx - stp);
    if (abs(stpc - stp) > abs(stpq - stp)) {
      stpf = stpc;
    } else {
      stpf = stpq;
    }
    brackt = 1;
  } else if (abs(dp) < abs(dx)) {
    /* Third case: A lower function value, derivatives of the same sign,
       and the magnitude of the derivative decreases.  The cubic step is
       computed only if the cubic tends to infinity in the direction of the
       step or if the minimum of the cubic is beyond stp. Otherwise the
       cubic step is defined to be the secant step. */
    theta = 3.0*(fx - fp)/(stp - stx) + dx + dp;
    s = max(abs(theta),abs(dx),abs(dp));
    /* The case gamma = 0 only arises if the cubic does not tend to
       infinity in the direction of the step. */
    gamma = s*sqrt(max(0.0, (theta/s)^2 - (dx/s)*(dp/s)));
    if (stp > stx) gamma = -gamma;
    p = (gamma - dp) + theta;
    q = (gamma + (dx - dp)) + gamma;
    r = p/q;
    if (r < 0.0 && gamma != 0.0) {
      stpc = stp + r*(stx - stp);
    } else if (stp > stx) {
      stpc = stpmax;
    } else {
      stpc = stpmin;
    }
    stpq = stp + (dp/(dp - dx))*(stx - stp);
    if (brackt) {
      /* A minimizer has been bracketed. If the cubic step is closer to stp
	 than the secant step, the cubic step is taken, otherwise the
	 secant step is taken. */
      if (abs(stpc - stp) < abs(stpq - stp)) {
	stpf = stpc;
      } else {
	stpf = stpq;
      }
      if (stp > stx) {
	stpf = min(stp + 0.66*(sty - stp), stpf);
      } else {
	stpf = max(stp + 0.66*(sty - stp), stpf);
      }
    } else {
      /* A minimizer has not been bracketed. If the cubic step is farther
	 from stp than the secant step, the cubic step is taken, otherwise
	 the secant step is taken. */
      if (abs(stpc-stp) > abs(stpq-stp)) {
	stpf = stpc;
      } else {
	stpf = stpq;
      }
      stpf = min(stpmax,stpf);
      stpf = max(stpmin,stpf);
    }
  } else {
    /* Fourth case: A lower function value, derivatives of the same sign,
       and the magnitude of the derivative does not decrease. If the
       minimum is not bracketed, the step is either stpmin or stpmax,
       otherwise the cubic step is taken. */
    if (brackt) {
      theta = 3.0*(fp - fy)/(sty - stp) + dy + dp;
      s = max(abs(theta),abs(dy),abs(dp));
      gamma = s*sqrt((theta/s)^2 - (dy/s)*(dp/s));
      if (stp > sty) gamma = -gamma;
      p = (gamma - dp) + theta;
      q = ((gamma - dp) + gamma) + dy;
      r = p/q;
      stpc = stp + r*(sty - stp);
      stpf = stpc;
    } else if (stp > stx) {
      stpf = stpmax;
    } else {
      stpf = stpmin;
    }
  }

  /* Update the interval which contains a minimizer. */
  if (fp > fx) {
    sty = stp;
    fy = fp;
    dy = dp;
  } else {
    if (sgnd < 0.0) {
      sty = stx;
      fy = fx;
      dy = dx;
    }
    stx = stp;
    fx = fp;
    dx = dp;
  }

  /* Compute the new step. */
  stp = stpf;
}

/*---------------------------------------------------------------------------*/
/* UTILITIES */

func optim_nrm2(x)
/* DOCUMENT optim_nrm2(x)
    Returns the Euclidian norm of X: sqrt(sum(X*X)), taking care of overflows.
*/
{
  if (structof(x) == complex) {
    y = x.im;
    x = double(x);
    s = max(-min(x), max(x), -min(y), max(y));
    if (! s) return 0.0;
    x *= (1.0/s);
    y *= (1.0/s);
    return s*sqrt(sum(x*x) + sum(y*y));
  } else {
    s = max(-min(x), max(x));
    if (! s) return 0.0;
    x *= (1.0/s);
    return s*sqrt(sum(x*x));
  }
}

/*---------------------------------------------------------------------------*/
