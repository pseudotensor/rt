/* Runge-Kutta 2(3), Euler-Cauchy */
/* Author:  G. Jungman */
/* Reference: Abramowitz & Stegun, section 25.5. Runge-Kutta 2nd (25.5.7)
   and 3rd (25.5.8) order methods */
#include "config.h"
#include <stdlib.h>
#include <string.h>
#include "gsl_errno.h"
#include "gsl_odeiv.h"
#include "odeiv_util.h"

typedef struct { doub *k1,*k2,*k3,*ytmp;} rk2_state_t;

static void *
rk2_alloc (size_t dim)
{
  rk2_state_t *state = (rk2_state_t *) malloc (sizeof (rk2_state_t));
  state->k1 = (doub *) malloc (dim * sizeof (doub));
  state->k2 = (doub *) malloc (dim * sizeof (doub));
  state->k3 = (doub *) malloc (dim * sizeof (doub));
  state->ytmp = (doub *) malloc (dim * sizeof (doub));
  return state;}


static int
rk2_apply (void *vstate,
           size_t dim,
           doub t,
           doub h,
           doub y[],
           doub yerr[],
           const doub dydt_in[],
           doub dydt_out[], 
           const gsl_odeiv_system * sys)
{
  rk2_state_t *state = (rk2_state_t *) vstate;

  size_t i;

  doub* const k1 = state->k1;
  doub* const k2 = state->k2;
  doub* const k3 = state->k3;
  doub* const ytmp = state->ytmp;

  /* k1 step */
  /* k1 = f(t,y) */

  if (dydt_in != NULL)
    {DBL_MEMCPY (k1, dydt_in, dim);}
  else GSL_ODEIV_FN_EVAL (sys, t, y, k1);

  /* k2 step */
  /* k2 = f(t + 0.5*h, y + 0.5*k1) */

  for (i = 0; i < dim; i++){ytmp[i] = y[i] + 0.5 * h * k1[i];}
GSL_ODEIV_FN_EVAL (sys, t + 0.5 * h, ytmp, k2);

  /* k3 step */
  /* for 3rd order estimates, is used for error estimation
     k3 = f(t + h, y - k1 + 2*k2) */
 
  for (i = 0; i < dim; i++){ytmp[i] = y[i] + h * (-k1[i] + 2.0 * k2[i]);}
GSL_ODEIV_FN_EVAL (sys, t + h, ytmp, k3);

  /* final sum */
  
  for (i = 0; i < dim; i++)
    {
      /* Save original values if derivative evaluation below fails */
//      ytmp[i] = y[i];

      {
	const doub ksum3 = (k1[i] + 4.0 * k2[i] + k3[i]) / 6.0;
	y[i] += h * ksum3;
      yerr[i] = h * (k2[i] - ksum3);  /* Error estimation */
      }
    }
   return GSL_SUCCESS;
}

static int
rk2_reset (void *vstate, size_t dim)
{
  rk2_state_t *state = (rk2_state_t *) vstate;

  DBL_ZERO_MEMSET (state->k1, dim);
  DBL_ZERO_MEMSET (state->k2, dim);
  DBL_ZERO_MEMSET (state->k3, dim);
  DBL_ZERO_MEMSET (state->ytmp, dim);

  return GSL_SUCCESS;
}

static unsigned int
rk2_order (void *vstate)
{
  rk2_state_t *state = (rk2_state_t *) vstate;
  state = 0; /* prevent warnings about unused parameters */
  return 2;
}

static void
rk2_free (void *vstate)
{
  rk2_state_t *state = (rk2_state_t *) vstate;
  free (state->k1);
  free (state->k2);
  free (state->k3);
  free (state->ytmp);
  free (state);
}

static const gsl_odeiv_step_type rk2_type = { "rk2",    /* name */
  1,                           /* can use dydt_in */
  1,                           /* gives exact dydt_out */
  &rk2_alloc,
  &rk2_apply,
  &rk2_reset,
  &rk2_order,
  &rk2_free
};

const gsl_odeiv_step_type *gsl_odeiv_step_rk2 = &rk2_type;
