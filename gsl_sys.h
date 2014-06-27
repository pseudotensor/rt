/* sys/gsl_sys.h
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Gerard Jungman, Brian Gough
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef __GSL_SYS_H__
#define __GSL_SYS_H__

#undef __BEGIN_DECLS
#undef __END_DECLS
#ifdef __cplusplus
# define __BEGIN_DECLS extern "C" {
# define __END_DECLS }
#else
# define __BEGIN_DECLS /* empty */
# define __END_DECLS /* empty */
#endif

__BEGIN_DECLS

doub gsl_log1p (const doub x);
doub gsl_expm1 (const doub x);
doub gsl_hypot (const doub x, const doub y);
doub gsl_hypot3 (const doub x, const doub y, const doub z);
doub gsl_acosh (const doub x);
doub gsl_asinh (const doub x);
doub gsl_atanh (const doub x);

int gsl_isnan (const doub x);
int gsl_isinf (const doub x);
int gsl_finite (const doub x);

doub gsl_nan (void);
doub gsl_posinf (void);
doub gsl_neginf (void);
doub gsl_fdiv (const doub x, const doub y);

doub gsl_coerce_doub (const doub x);
float gsl_coerce_float (const float x);
long double gsl_coerce_long_double (const long double x);

doub gsl_ldexp(const doub x, const int e);
doub gsl_frexp(const doub x, int * e);

int gsl_fcmp (const doub x1, const doub x2, const doub epsilon);

__END_DECLS

#endif /* __GSL_SYS_H__ */
