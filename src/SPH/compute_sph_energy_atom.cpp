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

#include "compute_sph_energy_atom.h"
#include <cstring>
#include "atom.h"
#include "update.h"
#include "modify.h"
#include "comm.h"
#include "memory.h"
#include "error.h"
#include <iostream>


using namespace LAMMPS_NS;
using namespace std;


/* ---------------------------------------------------------------------- */

ComputeSPHEnergyAtom::ComputeSPHEnergyAtom(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  if (narg != 3) error->all(FLERR,"Illegal compute sph/energy/atom command");
  if (atom->rho_flag != 1)
    error->all(FLERR,"Compute sph/energy/atom command requires atom_style sph");

  peratom_flag = 1;
  size_peratom_cols = 0;

  nmax = 0;
  energyVector = nullptr;
}

/* ---------------------------------------------------------------------- */

ComputeSPHEnergyAtom::~ComputeSPHEnergyAtom()
{
  memory->sfree(energyVector);
}

/* ---------------------------------------------------------------------- */

void ComputeSPHEnergyAtom::init()
{

  int count = 0;
  for (int i = 0; i < modify->ncompute; i++)
    if (strcmp(modify->compute[i]->style,"energyVector/atom") == 0) count++;
  if (count > 1 && comm->me == 0)
    error->warning(FLERR,"More than one compute energyVector/atom");
}

/* ---------------------------------------------------------------------- */

void ComputeSPHEnergyAtom::compute_peratom()
{
  invoked_peratom = update->ntimestep;

  // grow rhoVector array if necessary

  if (atom->nmax > nmax) {
    memory->sfree(energyVector);
    nmax = atom->nmax;
    energyVector = (double *) memory->smalloc(nmax*sizeof(double),"atom:energyVector");
    vector_atom = energyVector;
  }


  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  double *temperature = atom->temperature;
  double **v = atom->v;
  double vdotv;
  double *mass = atom->mass;
  double *esph = atom->esph;
  int *type = atom->type;
  double imass;
  int itype;


    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
              itype = type[i];
              imass = mass[itype];
              vdotv = v[i][0] * v[i][0] + v[i][1] * v[i][1] + v[i][2] * v[i][2];
              //std::cout << "i:" << i << std::endl;
              //std::cout << "vdotv:" << vdotv << std::endl;
             //Ti = 2./3./kb/np * ei;
              energyVector[i] = esph[i] + 0.5 * imass * vdotv;
      }
      else {
              energyVector[i] = 0.0;
      }
    }
    
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double ComputeSPHEnergyAtom::memory_usage()
{
  double bytes = (double)nmax * sizeof(double);
  return bytes;
}
