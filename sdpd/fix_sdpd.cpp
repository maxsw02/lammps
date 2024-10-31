// clang-format off
/* ----------------------------------------------------------------------
 LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
 https://www.lammps.org/, Sandia National Laboratories
 LAMMPS development team: developers@lammps.org

 Copyright (2003) Sandia Corporation.  Under the terms of Contract
 DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
 certain rights in this software.  This software is distributed under
 the GNU General Public License.

 See the README file in the top-level LAMMPS directory.
 ------------------------------------------------------------------------- */

#include "fix_sdpd.h"
#include "atom.h"
#include "force.h"
#include "update.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixSDPD::FixSDPD(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg) {

  if ((atom->esph_flag != 1) || (atom->rho_flag != 1))
    error->all(FLERR,
        "Fix sph command requires atom_style with both energy and density");

  if (narg != 3)
    error->all(FLERR,"Illegal number of arguments for fix sph command");

  time_integrate = 1;
}

/* ---------------------------------------------------------------------- */

int FixSDPD::setmask() {
  int mask = 0;
  mask |= INITIAL_INTEGRATE;
  mask |= FINAL_INTEGRATE;
  mask |= PRE_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixSDPD::init() {
  dtv = update->dt;
  dt2 = 0.5 * dtv * force->ftm2v; 
  //dtv = update->dt;
  //dtf = 0.5 * update->dt * force->ftm2v;
}

void FixSDPD::setup_pre_force(int /*vflag*/)
{
  // set vest equal to v
  double **v = atom->v;
  double **vest = atom->vest;
  double *entropy = atom->entropy;
  double *entropy_est = atom->entropyest;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  if (igroup == atom->firstgroup)
    nlocal = atom->nfirst;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      vest[i][0] = v[i][0];
      vest[i][1] = v[i][1];
      vest[i][2] = v[i][2];

      entropy_est[i] = entropy[i];
    }
  }
}

/* ----------------------------------------------------------------------
 allow for both per-type and per-atom mass
 ------------------------------------------------------------------------- */

void FixSDPD::initial_integrate(int /*vflag*/) {
  // update v and x and rho and e of atoms in group

  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  double **vest = atom->vest;
  //double *rho = atom->rho;
  //double *drho = atom->drho;
  //double *esph = atom->esph;
  double *desph = atom->desph;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  double *entropy = atom->entropy;
  double *entropyest = atom->entropyest;
  double *dentropy = atom->dentropy;
  int rmass_flag = atom->rmass_flag;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int i;
  double dtfm;

  if (igroup == atom->firstgroup)
    nlocal = atom->nfirst;

  for (i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      if (rmass_flag) {
        dtfm = dt2 / rmass[i];
      } else {
        dtfm = dt2 / mass[type[i]];
      }

      //esph[i] += dt2 * desph[i]; // half-step update of particle internal energy
            //rho[i] += dt * drho[i]; // ... and density

      // extrapolate velocity for use with velocity-dependent potentials, e.g. SPH
      vest[i][0] = v[i][0] + 2.0 * dtfm * f[i][0];
      vest[i][1] = v[i][1] + 2.0 * dtfm * f[i][1];
      vest[i][2] = v[i][2] + 2.0 * dtfm * f[i][2];

      entropyest[i] = entropy[i] + 2.0 * dtfm * dentropy[i];

      v[i][0] += dtfm * f[i][0];
      v[i][1] += dtfm * f[i][1];
      v[i][2] += dtfm * f[i][2];

      entropy[i] += dtfm * dentropy[i];

      x[i][0] += dtv * v[i][0];
      x[i][1] += dtv * v[i][1];
      x[i][2] += dtv * v[i][2];
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixSDPD::final_integrate() {

  // update v, rho, and e of atoms in group

  double **v = atom->v;
  double **f = atom->f;
  //double *esph = atom->esph;
  //double *desph = atom->desph;
  //double *rho = atom->rho;
  //double *drho = atom->drho;
  double *entropy = atom->entropy;
  double *dentropy = atom->dentropy;
  int *type = atom->type;
  int *mask = atom->mask;
  double *mass = atom->mass;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup)
    nlocal = atom->nfirst;
  double dtfm;
  double *rmass = atom->rmass;
  int rmass_flag = atom->rmass_flag;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {

      if (rmass_flag) {
        dtfm = dt2 / rmass[i];
      } else {
        dtfm = dt2 / mass[type[i]];
      }
      v[i][0] += dtfm * f[i][0];
      v[i][1] += dtfm * f[i][1];
      v[i][2] += dtfm * f[i][2];

      entropy[i] += dtfm * dentropy[i];
      //esph[i] += dt2 * desph[i];
      //rho[i] += dt2 * drho[i];
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixSDPD::reset_dt() {
  dtv = update->dt;
  dt2 = 0.5 * dtv * force->ftm2v;
  //dt2 = std::sqrt(dtv); 
  //dtv = update->dt;
  //dtf = 0.5 * update->dt * force->ftm2v;
}
