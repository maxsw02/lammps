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

#ifdef COMPUTE_CLASS
// clang-format off
ComputeStyle(sph/entropy/atom,ComputeSPHEntropyAtom);
// clang-format on
#else

#ifndef LMP_COMPUTE_SPH_ENT_ATOM_H
#define LMP_COMPUTE_SPH_ENT_ATOM_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeSPHEntropyAtom : public Compute {
 public:
  ComputeSPHEntropyAtom(class LAMMPS *, int, char **);
  ~ComputeSPHEntropyAtom() override;
  void init() override;
  void compute_peratom() override;
  double memory_usage() override;

 private:
  int nmax;
  double *entropyVector;
};

}    // namespace LAMMPS_NS

#endif
#endif
