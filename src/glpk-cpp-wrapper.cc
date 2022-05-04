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

#include "glpk-cpp-wrapper.h"

int
glpk (int n, int m, double *c, int nz, int *rn, int *cn,
      double *a, double *b, char *ctype, int *freeLB, double *lb,
      int *freeUB, double *ub, double *xOpt, double *fOpt)
{
  int typx = 0;
  int errnum = 0;

  // GLPK solver parameters
  optParams par;

  // The type of solver to use
  int lpsolver = 1;
  // The scale to use
  int scale = 16;
  // Level of messages output by the solver
  par.msglev = 1;
  // Dual simplex option
  par.dual = 1;
  // Pricing option
  par.price = 34;
  // Simplex iterations limit
  par.itlim = std::numeric_limits<int>::max ();
  // Output frequency, in iterations
  par.outfrq = 200;
  // Branching heuristic option
  par.branch = 4;
  // Backtracking heuristic option
  par.btrack = 4;
  // Presolver option (1 active, 0 off)
  par.presol = 1; 
  // Ratio test option
  par.rtest = 34;  
  par.tmlim = std::numeric_limits<int>::max ();
  par.outdly = 0;
  // Relative tolerance used to check if the current basic solution
  // is primal feasible
  par.tolbnd = 1e-7;
  // Absolute tolerance used to check if the current basic solution
  // is dual feasible
  par.toldj = 1e-7;
  // Relative tolerance used to choose eligible pivotal elements of
  //  the simplex table in the ratio test
  par.tolpiv = 1e-10;
  par.objll = -std::numeric_limits<double>::max ();
  par.objul = std::numeric_limits<double>::max ();
  par.tolint = 1e-5;
  par.tolobj = 1e-7;


  glp_prob *lp = glp_create_prob ();
  // Set the sense of optimization
  glp_set_obj_dir (lp, GLP_MAX);
  NS_LOG_LOGIC ("Created GLPK LP MAX problem");

  glp_add_cols (lp, n);
  for (int i = 0; i < n; i++)
    {
      // Define type of the structural variables
      if (! freeLB[i] && ! freeUB[i])
        {
          if (lb[i] != ub[i])
          {
            NS_LOG_INFO ("At index: " << std::to_string (i) << " DB var with lower bound: " <<
                          std::to_string(lb[i]) << " and upper bound: " << std::to_string (ub[i]));

            glp_set_col_bnds (lp, i+1, GLP_DB, lb[i], ub[i]);
          }
          else
          {
            NS_LOG_INFO ("At index: " << std::to_string (i) << " FX var with bounds: " <<
                          std::to_string (ub[i]));
            glp_set_col_bnds (lp, i+1, GLP_FX, lb[i], ub[i]);
          }
        }
      else
        {
          if (! freeLB[i] && freeUB[i])
          {
            NS_LOG_INFO ("At index: " << std::to_string (i) << " LB var with lower bound: " <<
                          std::to_string (lb[i]));
            glp_set_col_bnds (lp, i+1, GLP_LO, lb[i], ub[i]);
          }
          else
            {
              if (freeLB[i] && ! freeUB[i])
              {
                NS_LOG_INFO ("At index: " << std::to_string (i) << " UB var with upper bound: " <<
                          std::to_string (ub[i]));
                glp_set_col_bnds (lp, i+1, GLP_UP, lb[i], ub[i]);
              }
              else
              {
                NS_LOG_INFO ("At index: " << std::to_string (i) << " unbounded variable");
                glp_set_col_bnds (lp, i+1, GLP_FR, lb[i], ub[i]);
              }
            }
        }

      // -- Set the objective coefficient of the corresponding
      // -- structural variable.  No constant term is assumed.
      NS_LOG_INFO ("At index: " << std::to_string (i) << 
                   " objective function coefficient: " << std::to_string (c[i]));
      glp_set_obj_coef(lp,i+1,c[i]);
    }

  glp_add_rows (lp, m);

  for (int i = 0; i < m; i++)
    {
      // If the i-th row has no lower bound (types F,U), the
      // corrispondent parameter will be ignored.  If the i-th row has
      // no upper bound (types F,L), the corrispondent parameter will be
      // ignored.  If the i-th row is of S type, the i-th LB is used,
      // but the i-th UB is ignored.

      switch (ctype[i])
        {
        case 'F':
          typx = GLP_FR;
          NS_LOG_INFO ("At index: " << std::to_string (i) << 
                   " constraint of type F");
          break;

        case 'U':
          typx = GLP_UP;
          NS_LOG_INFO ("At index: " << std::to_string (i) << 
                   " constraint of type U");
          break;

        case 'L':
          typx = GLP_LO;
          NS_LOG_INFO ("At index: " << std::to_string (i) << 
                   " constraint of type L");
          break;

        case 'S':
          typx = GLP_FX;
          NS_LOG_INFO ("At index: " << std::to_string (i) << 
                   " constraint of type S");
          break;

        case 'D':
          NS_LOG_INFO ("At index: " << std::to_string (i) << 
                   " constraint of type D");
          typx = GLP_DB;
          break;
        }

      NS_LOG_INFO ("Value of the constraint: " << std::to_string (b[i]));
      glp_set_row_bnds (lp, i+1, typx, b[i], b[i]);
    }

  glp_load_matrix (lp, nz, rn, cn, a);

  // scale the problem data
  if (! par.presol || lpsolver != 1)
    glp_scale_prob (lp, scale);

  // build advanced initial basis (if required)
  if (lpsolver == 1 && ! par.presol)
    glp_adv_basis (lp, 0);

  // For MIP problems without a presolver, a first pass with glp_simplex
  // is required
  if (lpsolver == 1)
    {
      glp_smcp smcp;
      glp_init_smcp (&smcp);
      smcp.msg_lev = par.msglev;
      smcp.meth = par.dual;
      smcp.pricing = par.price;
      smcp.r_test = par.rtest;
      smcp.tol_bnd = par.tolbnd;
      smcp.tol_dj = par.toldj;
      smcp.tol_piv = par.tolpiv;
      smcp.obj_ll = par.objll;
      smcp.obj_ul = par.objul;
      smcp.it_lim = par.itlim;
      smcp.tm_lim = par.tmlim;
      smcp.out_frq = par.outfrq;
      smcp.out_dly = par.outdly;
      smcp.presolve = par.presol;
      errnum = glp_simplex (lp, &smcp);
    }

  if (lpsolver == 2)
    {
      glp_iptcp iptcp;
      glp_init_iptcp (&iptcp);
      iptcp.msg_lev = par.msglev;
      errnum = glp_interior (lp, &iptcp);
    }

 if (errnum == 0)
    {
      if (lpsolver == 1)
        {
          *fOpt = glp_get_obj_val (lp);
        }
      else
        {
          *fOpt = glp_ipt_obj_val (lp);
        }
      // Primal values
      for (int i = 0; i < n; i++)
        {
          if (lpsolver == 1)
            xOpt[i] = glp_get_col_prim (lp, i+1);
          else
            xOpt[i] = glp_ipt_col_prim (lp, i+1);
        }
    }

  glp_delete_prob (lp);
  // Request that GLPK free all memory resources.
  // This prevents reported memory leaks, but isn't strictly necessary.
  // The memory blocks used are allocated once and don't grow with further
  // calls to glpk so they would be reclaimed anyways when Octave exits.
  glp_free_env ();

  return errnum;
}
