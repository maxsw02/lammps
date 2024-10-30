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

#ifdef PAIR_CLASS
// clang-format off
PairStyle(sdpd/ml,PairSDPDML);
// clang-format on
#else

#ifndef LMP_PAIR_SDPD_IDEAL
#define LMP_PAIR_SDPD_IDEAL

#include "pair.h"
//#ifdef USE_ZEST
//#include "zest.hpp"
#include <random>
//#endif

namespace LAMMPS_NS {

class PairSDPDML : public Pair {
 public:

  PairSDPDML(class LAMMPS *);
  ~PairSDPDML() override;
  void compute(int, int) override;
  void settings(int, char **) override;
  void coeff(int, char **) override;
  double init_one(int, int) override;
  void init_style() override;
  void *extract_peratom(const char *, int &) override;

  int pack_reverse_comm(int, int, double *) override;
  void unpack_reverse_comm(int, int *, double *) override;
  int pack_forward_comm(int, int *, double *, int, int *) override;
  void unpack_forward_comm(int, int, double *) override;
  double memory_usage() override;



 protected:
  int nmax;
  double **cut;
  double cutforcesq;
  int first;
  double self_cut;
  double eta, visc, therm;
  double np;

  //per_atom arrays
  double *d;
  virtual void allocate();

  unsigned int seed;
  class RanMars *random;
};

}    // namespace LAMMPS_NS
#endif
#endif