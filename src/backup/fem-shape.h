#ifndef _SHAPE_H_
#define _SHAPE_H_

/**
 * @file   shape.h
 * @author Bernardo M. Rocha <bernardo@lncc.br>
 * @date   Wed Mar  3 13:35:27 2010
 * 
 * @brief  FEs shape functions
 *         Change element type names to: 
 *         TRI3N2D, QUA4N2D, TET4N3D, HEX8N3D 
 * 
 */

#include <cstdio>
#include "../util/datatype.h"

void
fem_iso_bilinear_quad(Real ***shl, Real *w);

void
fem_iso_linear_tri(Real ***shl, Real *w);

void
fem_iso_trilinear_hex(Real ***shl, Real *w);

void
fem_quad_basis(Real *x, Real *basis);

void
fem_quad_dbasis(Real *x, Real **dbasis);

#endif
