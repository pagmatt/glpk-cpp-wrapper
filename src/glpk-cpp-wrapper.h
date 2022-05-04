////////////////////////////////////////////////////////////////////////
//
// Copyright (C) 2005-2020 The Octave Project Developers
//
// See the file COPYRIGHT.md in the top-level directory of this
// distribution or <https://octave.org/copyright/>.
//
// This file is part of Octave.
//
// Octave is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Octave is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Octave; see the file COPYING.  If not, see
// <https://www.gnu.org/licenses/>.
//
////////////////////////////////////////////////////////////////////////

#ifndef GLPK_CPP_WRAPPER_H
#define GLPK_CPP_WRAPPER_H

#include <ctime>
#include <limits>

extern "C"
{
  #include <glpk.h>
}

struct optParams
{
  int msglev;
  int dual;
  int price;
  int itlim;
  int outfrq;
  int branch;
  int btrack;
  int presol;
  int rtest;
  int tmlim;
  int outdly;
  double tolbnd;
  double toldj;
  double tolpiv;
  double objll;
  double objul;
  double tolint;
  double tolobj;
};

/**
 * Solves the LP maximization problem, using the GLPK library as solver backend
 * 
 * \param n The number of rows of A
 * \param m The number of columns of A
 * \param c The objective function coefficients vector c
 * \param a The constraint matrix A, encoded as a vector containing its non-zero elements
 * \param ctype A vector of chars, specifying the sense of each constraint in b
 * \param nz The number of non-zero element in the constraint matrix A
 * \param rn The row indexes of the non-zero elements in the constraint matrix A
 * \param cn The column indexes of the non-zero elements in the constraint matrix A
 * \param freeLB A vector which specifies whether the optimization variables are bounded
 *                in the lower bound (1) or not (0)
 * \param lb The vector of lower bounds for the optimization variables
 * \param freeUB A vector which specifies whether the optimization variables are bounded
 *                in the upper bound (1) or not (0)
 * \param ub The vector of upper bounds for the optimization variables
 * \param lpsolver The solver to use 
 * \param scale The scale to use for the solver
 * \param optParams The GLPK solver parameters
 * \param xmin A pointer to where the x solution will be stored
 * \param fOpt A pointer to where the optimal value of the objective function will be stored
 */
int
glpk (int n, int m, double *c, int nz, int *rn, int *cn,
      double *a, double *b, char *ctype, int *freeLB, double *lb,
      int *freeUB, double *ub, double *xOpt, double *fOpt);

#endif
