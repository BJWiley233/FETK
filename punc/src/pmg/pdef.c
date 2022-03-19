/* ./src_f77/pdef.f -- translated by f2c (version 20200916).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

#include <punc/vf2c.h>

/* * /////////////////////////////////////////////////////////////////////////// */
/* * @file    pded.f */
/* * @author  Michael Holst */
/* * @brief   PDE definition file. */
/* * @version $Id: pdef.f,v 1.2 2010/08/12 05:53:09 fetk Exp $ */
/* * @attention */
/* * @verbatim */
/* * */
/* * PMG -- Parallel algebraic MultiGrid */
/* * Copyright (C) 1994-- Michael Holst. */
/* * */
/* * Michael Holst <mholst@math.ucsd.edu> */
/* * University of California, San Diego */
/* * Department of Mathematics, 5739 AP&M */
/* * 9500 Gilman Drive, Dept. 0112 */
/* * La Jolla, CA 92093-0112 USA */
/* * http://math.ucsd.edu/~mholst */
/* * */
/* * This file is part of PMG. */
/* * */
/* * PMG is free software; you can redistribute it and/or modify */
/* * it under the terms of the GNU General Public License as published by */
/* * the Free Software Foundation; either version 2 of the License, or */
/* * (at your option) any later version. */
/* * */
/* * PMG is distributed in the hope that it will be useful, */
/* * but WITHOUT ANY WARRANTY; without even the implied warranty of */
/* * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the */
/* * GNU General Public License for more details. */
/* * */
/* * You should have received a copy of the GNU General Public License */
/* * along with PMG; if not, write to the Free Software */
/* * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307  USA */
/* * */
/* * Linking PMG statically or dynamically with other modules is making a */
/* * combined work based on PMG. Thus, the terms and conditions of the GNU */
/* * General Public License cover the whole combination. */
/* * */
/* * SPECIAL GPL EXCEPTION */
/* * In addition, as a special exception, the copyright holders of PMG */
/* * give you permission to combine the PMG program with free software */
/* * programs and libraries that are released under the GNU LGPL or with */
/* * code included in releases of ISIM, PMV, PyMOL, SMOL, VMD, and Vision. */
/* * Such combined software may be linked with PMG and redistributed together */
/* * in original or modified form as mere aggregation without requirement that */
/* * the entire work be under the scope of the GNU General Public License. */
/* * This special exception permission is also extended to any software listed */
/* * in the SPECIAL GPL EXCEPTION clauses by the FEtk and APBS libraries. */
/* * */
/* * Note that people who make modified versions of PMG are not obligated */
/* * to grant this special exception for their modified versions; it is */
/* * their choice whether to do so. The GNU General Public License gives */
/* * permission to release a modified version without this exception; this */
/* * exception also makes it possible to release a modified version which */
/* * carries forward this exception. */
/* * */
/* * @endverbatim */
/* * /////////////////////////////////////////////////////////////////////////// */
/* * ZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZ */
/* * */
/* * DIFFERENTIAL EQUATION: Poisson's */
/* * */
/* * BOUNDARY CONDITIONS: */
/* * */
/* *     East   Face (xmin):  Dirichlet, homogeneous */
/* *     West   Face (xmax):  Dirichlet, homogeneous */
/* *     North  Face (ymin):  Dirichlet, homogeneous */
/* *     South  Face (ymax):  Dirichlet, homogeneous */
/* *     Top    Face (zmin):  Dirichlet, homogeneous */
/* *     Bottom Face (zmax):  Dirichlet, homogeneous */
/* * */
/* * MESH:                  hx  = (xmax-xmin) / (nx-1) */
/* *                        hy  = (ymax-ymin) / (ny-1) */
/* *                        hz  = (zmax-zmin) / (nz-1) */
/* *                        xi = xmin + (i-1) * hx,  i=1,...,nx */
/* *                        yi = ymin + (j-1) * hy,  j=1,...,ny */
/* *                        zi = zmin + (k-1) * hk,  k=1,...,nz */
/* * */
/* * CHOSEN TRUE SOLUTION:  u(x,y,z) = Sin(n pi x)*Sin(m pi y)*Sin(l pi z) */
/* * */
/* * author:  michael holst */
/* * ZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZ */
/* * */
doublereal c_scal__(doublereal *coef, doublereal *u, integer *ipkey)
{
    /* System generated locals */
    doublereal ret_val;

/* * ********************************************************************* */
/* * purpose: */
/* * */
/* *    define the nonlinearity (scalar version) */
/* * */
/* * author:  michael holst */
/* * ********************************************************************* */
/* * */
/* *    *** create the nonlinear term *** */
    ret_val = *coef * *u;
/* * */
/* *    *** end it *** */
    return ret_val;
} /* c_scal__ */

doublereal dc_scal__(doublereal *coef, doublereal *u, integer *ipkey)
{
    /* System generated locals */
    doublereal ret_val;

/* * ********************************************************************* */
/* * purpose: */
/* * */
/* *    define the derivative of the nonlinearity (scalar version) */
/* * */
/* * author:  michael holst */
/* * ********************************************************************* */
/* * */
/* *    *** create the nonlinear term *** */
    ret_val = *coef;
/* * */
/* *    *** end it *** */
    return ret_val;
} /* dc_scal__ */

/* Subroutine */ int c_vec__(doublereal *coef, doublereal *uin, doublereal *
	uout, integer *nx, integer *ny, integer *nz, integer *ipkey)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__, n, ii, ipara, ivect;

/* * ********************************************************************* */
/* * purpose: */
/* * */
/* *    define the nonlinearity (vector version) */
/* * */
/* * author:  michael holst */
/* * ********************************************************************* */
/* * */
/* mdir 0 0 */
/* * */
/* *    *** find parallel loops (ipara), remainder (ivect) *** */
    /* Parameter adjustments */
    --uout;
    --uin;
    --coef;

    /* Function Body */
    n = *nx * *ny * *nz;
    ipara = n;
    ivect = n % 1;
/* * */
/* *    *** do parallel loops *** */
/* mdir 2 1 */
    for (ii = 1; ii <= 1; ++ii) {
/* mdir 2 2 */
	i__1 = ipara * ii;
	for (i__ = ipara * (ii - 1) + 1; i__ <= i__1; ++i__) {
	    uout[i__] = coef[i__] * uin[i__];
/* L11: */
	}
/* L10: */
    }
/* * */
/* *    *** do vector loops *** */
/* mdir 1 1 */
    i__1 = n;
    for (i__ = ipara + 1; i__ <= i__1; ++i__) {
	uout[i__] = coef[i__] * uin[i__];
/* L20: */
    }
/* * */
/* *    *** end it *** */
    return 0;
} /* c_vec__ */

/* Subroutine */ int dc_vec__(doublereal *coef, doublereal *uin, doublereal *
	uout, integer *nx, integer *ny, integer *nz, integer *ipkey)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__, n, ii, ipara, ivect;

/* * ********************************************************************* */
/* * purpose: */
/* * */
/* *    define the derivative of the nonlinearity (vector version) */
/* * */
/* * author:  michael holst */
/* * ********************************************************************* */
/* * */
/* mdir 0 0 */
/* * */
/* *    *** find parallel loops (ipara), remainder (ivect) *** */
    /* Parameter adjustments */
    --uout;
    --uin;
    --coef;

    /* Function Body */
    n = *nx * *ny * *nz;
    ipara = n;
    ivect = n % 1;
/* * */
/* *    *** do parallel loops *** */
/* mdir 2 1 */
    for (ii = 1; ii <= 1; ++ii) {
/* mdir 2 2 */
	i__1 = ipara * ii;
	for (i__ = ipara * (ii - 1) + 1; i__ <= i__1; ++i__) {
	    uout[i__] = coef[i__];
/* L11: */
	}
/* L10: */
    }
/* * */
/* *    *** do vector loops *** */
/* mdir 1 1 */
    i__1 = n;
    for (i__ = ipara + 1; i__ <= i__1; ++i__) {
	uout[i__] = coef[i__];
/* L20: */
    }
/* * */
/* *    *** end it *** */
    return 0;
} /* dc_vec__ */

/* Subroutine */ int fillco_(integer *iparm, doublereal *rparm, integer *nx, 
	integer *ny, integer *nz, doublereal *xf, doublereal *yf, doublereal *
	zf, doublereal *gxcf, doublereal *gycf, doublereal *gzcf, doublereal *
	a1cf, doublereal *a2cf, doublereal *a3cf, doublereal *ccf, doublereal 
	*fcf, doublereal *tcf)
{
    /* System generated locals */
    integer a1cf_dim1, a1cf_dim2, a1cf_offset, a2cf_dim1, a2cf_dim2, 
	    a2cf_offset, a3cf_dim1, a3cf_dim2, a3cf_offset, ccf_dim1, 
	    ccf_dim2, ccf_offset, fcf_dim1, fcf_dim2, fcf_offset, tcf_dim1, 
	    tcf_dim2, tcf_offset, gxcf_dim1, gxcf_dim2, gxcf_offset, 
	    gycf_dim1, gycf_dim2, gycf_offset, gzcf_dim1, gzcf_dim2, 
	    gzcf_offset, i__1, i__2, i__3;

    /* Builtin functions */
    double atan(doublereal), sin(doublereal);

    /* Local variables */
    static integer i__, j, k;
    static doublereal pi, hx, hy, hz, hxo2, hyo2, hzo2, xmin, ymin, xmax, 
	    ymax, zmin, zmax;
    static integer iinfo;

/* * ********************************************************************* */
/* * purpose: */
/* * */
/* *   this subroutine defines poisson's equation on the unit cube. */
/* *   boundary conditions are zero dirichlet. */
/* * */
/* * details: */
/* * */
/* *   this routine sets up the coefficients of the following */
/* *   three-dimensional, 2nd order linear elliptic partial */
/* *   differential equation: */
/* * */
/* *      lu = f, u in omega */
/* *       u = g, u on boundary of omega */
/* * */
/* *   where */
/* *      omega = [xmin,xmax]x[ymin,ymax]x[zmin,zmax] */
/* * */
/* *   and l is taken to be in the following form: */
/* * */
/* *      lu = - \nabla \cdot (a \nabla u) + c u */
/* * */
/* *   the coefficients must be supplied as follows: */
/* * */
/* *      a1cf(,,,):  at x-half grid points, hence array is:  (nx-1)*ny*nz */
/* *      a2cf(,,,):  at y-half grid points, hence array is:  nx*(ny-1)*nz */
/* *      a3cf(,,,):  at z-half grid points, hence array is:  nx*ny*(nz-1) */
/* *                  (we index all three as:  nx*ny*nz however) */
/* * */
/* *      ccf:        at grid points, hence array is:  nx*ny*nz */
/* *      fcf:        at grid points, hence array is:  nx*ny*nz */
/* *      u:          must be set to have appropriate boundary values */
/* * */
/* * author:  michael holst */
/* * ********************************************************************* */
/* * */
/* *    *** variable declarations *** */
/* * */
/* *    *** some parameters *** */
    /* Parameter adjustments */
    --iparm;
    --rparm;
    --xf;
    gzcf_dim1 = *nx;
    gzcf_dim2 = *ny;
    gzcf_offset = 1 + gzcf_dim1 * (1 + gzcf_dim2);
    gzcf -= gzcf_offset;
    --yf;
    tcf_dim1 = *nx;
    tcf_dim2 = *ny;
    tcf_offset = 1 + tcf_dim1 * (1 + tcf_dim2);
    tcf -= tcf_offset;
    fcf_dim1 = *nx;
    fcf_dim2 = *ny;
    fcf_offset = 1 + fcf_dim1 * (1 + fcf_dim2);
    fcf -= fcf_offset;
    ccf_dim1 = *nx;
    ccf_dim2 = *ny;
    ccf_offset = 1 + ccf_dim1 * (1 + ccf_dim2);
    ccf -= ccf_offset;
    a3cf_dim1 = *nx;
    a3cf_dim2 = *ny;
    a3cf_offset = 1 + a3cf_dim1 * (1 + a3cf_dim2);
    a3cf -= a3cf_offset;
    a2cf_dim1 = *nx;
    a2cf_dim2 = *ny;
    a2cf_offset = 1 + a2cf_dim1 * (1 + a2cf_dim2);
    a2cf -= a2cf_offset;
    a1cf_dim1 = *nx;
    a1cf_dim2 = *ny;
    a1cf_offset = 1 + a1cf_dim1 * (1 + a1cf_dim2);
    a1cf -= a1cf_offset;
    gycf_dim1 = *nx;
    gycf_dim2 = *nz;
    gycf_offset = 1 + gycf_dim1 * (1 + gycf_dim2);
    gycf -= gycf_offset;
    gxcf_dim1 = *ny;
    gxcf_dim2 = *nz;
    gxcf_offset = 1 + gxcf_dim1 * (1 + gxcf_dim2);
    gxcf -= gxcf_offset;
    --zf;

    /* Function Body */
    pi = atan(1.) * 4.;
/* * */
/* mdir 0 0 */
/* * */
/* *    *** info messages *** */
    iinfo = iparm[12];
/* * */
/* * ********************************************************************* */
/* * *** definition of the domain */
/* * ********************************************************************* */
/* * */
/* *    *** set the correct boundary ranges *** */
    xmin = 0.;
    xmax = 1.;
    ymin = 0.;
    ymax = 1.;
    zmin = 0.;
    zmax = 1.;
/* * */
/* *    *** store in rparm array *** */
    rparm[3] = xmin;
    rparm[4] = xmax;
    rparm[5] = ymin;
    rparm[6] = ymax;
    rparm[7] = zmin;
    rparm[8] = zmax;
/* * */
/* * ********************************************************************* */
/* * *** boundary value definitions */
/* * ********************************************************************* */
/* * */
/* *    *** determine mesh widths *** */
    hx = (xmax - xmin) / (doublereal) (*nx - 1);
    hy = (ymax - ymin) / (doublereal) (*ny - 1);
    hz = (zmax - zmin) / (doublereal) (*nz - 1);
    hxo2 = hx / 2.;
    hyo2 = hy / 2.;
    hzo2 = hz / 2.;
/* * */
/* * ********************************************************************* */
/* * *** pde coefficient value definitions */
/* * ********************************************************************* */
/* * */
/* *    *** fill coefficient arrays *** */
/* mdir 3 1 */
    i__1 = *nz;
    for (k = 1; k <= i__1; ++k) {
	zf[k] = zmin + (doublereal) (k - 1) * hz;
/* mdir 3 2 */
	i__2 = *ny;
	for (j = 1; j <= i__2; ++j) {
	    yf[j] = ymin + (doublereal) (j - 1) * hy;
/* mdir 3 3 */
	    i__3 = *nx;
	    for (i__ = 1; i__ <= i__3; ++i__) {
		xf[i__] = xmin + (doublereal) (i__ - 1) * hx;
/* * */
/* *             *** the diagonal tensor (2nd derivative) entries *** */
		a1cf[i__ + (j + k * a1cf_dim2) * a1cf_dim1] = 1.;
		a2cf[i__ + (j + k * a2cf_dim2) * a2cf_dim1] = 1.;
		a3cf[i__ + (j + k * a3cf_dim2) * a3cf_dim1] = 1.;
/* * */
/* *             *** the scalar (0th derivative) entry *** */
		ccf[i__ + (j + k * ccf_dim2) * ccf_dim1] = 0.;
/* * */
/* *             *** the rhs entry *** */
/* ZZZ           fcf(i,j,k) = ((zn*zn+zm*zm+zl*zl)*pi*pi */
/* ZZZ 2         *dsin(zn*pi*xf(i))*dsin(zm*pi*yf(j))*dsin(zl*pi*zf(k))) */
		fcf[i__ + (j + k * fcf_dim2) * fcf_dim1] = 1.f;
/* * */
/* *             *** the chosen true solution *** */
		tcf[i__ + (j + k * tcf_dim2) * tcf_dim1] = sin(pi * 1. * xf[
			i__]) * sin(pi * 1. * yf[j]) * sin(pi * 1. * zf[k]);
/* L12: */
	    }
/* L11: */
	}
/* L10: */
    }
/* * */
/* *    *** the (i=1) and (i=nx) boundaries (dirichlet) *** */
/* mdir 2 1 */
    i__1 = *nz;
    for (k = 1; k <= i__1; ++k) {
/* mdir 2 2 */
	i__2 = *ny;
	for (j = 1; j <= i__2; ++j) {
	    gxcf[j + (k + gxcf_dim2) * gxcf_dim1] = tcf[(j + k * tcf_dim2) * 
		    tcf_dim1 + 1];
	    gxcf[j + (k + (gxcf_dim2 << 1)) * gxcf_dim1] = tcf[*nx + (j + k * 
		    tcf_dim2) * tcf_dim1];
	    gxcf[j + (k + gxcf_dim2 * 3) * gxcf_dim1] = 0.;
	    gxcf[j + (k + (gxcf_dim2 << 2)) * gxcf_dim1] = 0.;
/* L51: */
	}
/* L50: */
    }
/* * */
/* *    *** the (j=1) and (j=ny) boundaries (dirichlet) *** */
/* mdir 2 1 */
    i__1 = *nz;
    for (k = 1; k <= i__1; ++k) {
/* mdir 2 2 */
	i__2 = *nx;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    gycf[i__ + (k + gycf_dim2) * gycf_dim1] = tcf[i__ + (k * tcf_dim2 
		    + 1) * tcf_dim1];
	    gycf[i__ + (k + (gycf_dim2 << 1)) * gycf_dim1] = tcf[i__ + (*ny + 
		    k * tcf_dim2) * tcf_dim1];
	    gycf[i__ + (k + gycf_dim2 * 3) * gycf_dim1] = 0.;
	    gycf[i__ + (k + (gycf_dim2 << 2)) * gycf_dim1] = 0.;
/* L61: */
	}
/* L60: */
    }
/* * */
/* *    *** the (k=1) and (k=nz) boundaries (dirichlet) *** */
/* mdir 2 1 */
    i__1 = *ny;
    for (j = 1; j <= i__1; ++j) {
/* mdir 2 2 */
	i__2 = *nx;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    gzcf[i__ + (j + gzcf_dim2) * gzcf_dim1] = tcf[i__ + (j + tcf_dim2)
		     * tcf_dim1];
	    gzcf[i__ + (j + (gzcf_dim2 << 1)) * gzcf_dim1] = tcf[i__ + (j + *
		    nz * tcf_dim2) * tcf_dim1];
	    gzcf[i__ + (j + gzcf_dim2 * 3) * gzcf_dim1] = 0.;
	    gzcf[i__ + (j + (gzcf_dim2 << 2)) * gzcf_dim1] = 0.;
/* L71: */
	}
/* L70: */
    }
/* * */
/* *    *** end it *** */
    return 0;
} /* fillco_ */

