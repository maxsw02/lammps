/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS
// clang-format off
FixStyle(sph,FixSDPD);
// clang-format on
#else

#ifndef LMP_FIX_SPH_H
#define LMP_FIX_SPH_H

#include "fix.h"

namespace LAMMPS_NS {

class FixSDPD : public Fix {
 public:
  FixSDPD(class LAMMPS *, int, char **);
  int setmask() override;
  void init() override;
  void setup_pre_force(int) override;
  void initial_integrate(int) override;
  void final_integrate() override;
  void reset_dt() override;

 private:
  class NeighList *list;

 protected:
  //double dtv, dtf;
  double dtv, dt2;
  double *step_respa;
  int mass_require;

  class Pair *pair;
};

}    // namespace LAMMPS_NS

#endif
#endif
