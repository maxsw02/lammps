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

#include "compute_sph_entropy_atom.h"
#include <cstring>
#include "atom.h"
#include "update.h"
#include "modify.h"
#include "comm.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeSPHEntropyAtom::ComputeSPHEntropyAtom(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  if (narg != 3) error->all(FLERR,"Illegal compute sph/entropy/atom command");
  if (atom->rho_flag != 1)
    error->all(FLERR,"Compute sph/entropy/atom command requires atom_style sph");

  peratom_flag = 1;
  size_peratom_cols = 0;

  nmax = 0;
  entropyVector = nullptr;
}

/* ---------------------------------------------------------------------- */

ComputeSPHEntropyAtom::~ComputeSPHEntropyAtom()
{
  memory->sfree(entropyVector);
}

/* ---------------------------------------------------------------------- */

void ComputeSPHEntropyAtom::init()
{

  int count = 0;
  for (int i = 0; i < modify->ncompute; i++)
    if (strcmp(modify->compute[i]->style,"entropyVector/atom") == 0) count++;
  if (count > 1 && comm->me == 0)
    error->warning(FLERR,"More than one compute entropyVector/atom");
}

/* ---------------------------------------------------------------------- */

void ComputeSPHEntropyAtom::compute_peratom()
{
  invoked_peratom = update->ntimestep;

  // grow rhoVector array if necessary

  if (atom->nmax > nmax) {
    memory->sfree(entropyVector);
    nmax = atom->nmax;
    entropyVector = (double *) memory->smalloc(nmax*sizeof(double),"atom:entropyVector");
    vector_atom = entropyVector;
  }

  // compute entropy for each atom in group

  double *entropy = atom->entropy;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
              entropyVector[i] = entropy[i];
      }
      else {
              entropyVector[i] = 0.0;
      }
    }
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double ComputeSPHEntropyAtom::memory_usage()
{
  double bytes = (double)nmax * sizeof(double);
  return bytes;
}
