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

/* ----------------------------------------------------------------------
   Contributing author:
    Max Win (UPENN)
    references: Espanol and Revenga, Phys Rev E 67, 026705 (2003)
------------------------------------------------------------------------- */

#include "pair_sdpd_ml.h"

#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "memory.h"
#include "neigh_list.h"
#include "neighbor.h"
#include "update.h"
#include "random_park.h"


#include <cmath>
#include <iostream>
//#ifndef USE_ZEST
#include "random_mars.h"

#include <vector>
#include <numeric>
//#endif


using namespace LAMMPS_NS;
using namespace std;
static const double pi=3.141592653589793238462643383279502884197;
//static const double sqrt_2_inv = std::sqrt(0.5);

/* ---------------------------------------------------------------------- */
PairSDPDML::PairSDPDML(LAMMPS *lmp): Pair (lmp)
{
  restartinfo = 0;
  manybody_flag = 1;  
  single_enable = 0;

  //d = F = nullptr;
  nmax = 0;
  d = nullptr;

  random = nullptr;
  first = 1;
  comm_forward = 1;
  //comm_reverse = 0;
}

/* ---------------------------------------------------------------------- */

PairSDPDML::~PairSDPDML() {
  
  if (copymode) return;
  memory->destroy(d);
  //memory->destroy(F);

  if (allocated) {
    memory->destroy (setflag);
    memory->destroy (cutsq);

    memory->destroy (cut);
  

  };
  if (random) delete random;
}

/* ---------------------------------------------------------------------- */

void PairSDPDML::compute (int eflag, int vflag) {
  int i, j, ii, jj, inum, jnum, itype, jtype;
  double xtmp, ytmp, ztmp, delx, dely, delz, ex,ey,ez, h, Pi, Pj;
  double fpgradx, fpgrady,fpgradz, fviscx, fviscy, fviscz, fvisc2x, fvisc2y, fvisc2z;
  double comb_temp, fdd;
  double Aij, Bij, Cij;
  double Ti, Tj;
  double term1, term2, term3, term4, term5, term6, term7 , term9; //Entropy terms
  double ih, ihsq, ihc, velx, vely, velz;
  double rsq, wfd, VijdotEij, W, vdotv, trace, term_rand;
  double wiener[3][3],sym_wiener[3][3], f_random[3], f_rand[3];
  double ei,ej;
  double F,wfdist;


  //ML vars
  int dim;
  double aux_x, aux_y, aux_z, aux2_x, aux2_y, aux2_z;
  double added_temp, gradAij_Ti, gradAij_Tj, gradBij_Ti, gradBij_Tj, gradCij_Ti, gradCij_Tj;
  double prefactor_x, prefactor_y, prefactor_z;
  double D;

  if (atom->nmax > nmax) {
    memory->destroy(d);
    nmax = atom->nmax;
    memory->create(d,nmax,"pair:d");
  }
  
  int *ilist, *jlist, *numneigh, **firstneigh;
  double vxtmp, vytmp, vztmp, imass, jmass, dV;
  
  eflag = 0;
  vflag = 0;
  ev_init(eflag, vflag);

  double **v = atom->v;
  double **x = atom->x;
  double **f = atom->f;
  double *mass = atom->mass;
  double *cv = atom->cv;
  double *esph = atom->esph;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int newton_pair = force->newton_pair;
  double dtinv = 1.0 / update->dt;
  double kb = force->boltz;
  double h_p = 0.1;
  double power_term, exp_term;
  
  int nall = nlocal + atom->nghost;

  double *entropy = atom->entropy;
  double *dentropy = atom->dentropy;
  double *temperature = atom->temperature;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  if (first) {
    for (i = 1; i <= atom->ntypes; i++) {
      for (j = 1; i <= atom->ntypes; i++) {
        if (cutsq[i][j] > 1.e-32) {
          if (!setflag[i][i] || !setflag[j][j]) {
            if (comm->me == 0) {
              printf(
                  "SPH particle types %d and %d interact with cutoff=%g, but not all of their single particle properties are set.\n",
                  i, j, sqrt(cutsq[i][j]));
            }
          }
        }
      }
    }
    first = 0;
  } 


  if (newton_pair) {
    for (i = 0; i < nall; i++) d[i] = 0.0;
  } else for (i = 0; i < nlocal; i++)  d[i] = 0.0;

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    h = self_cut; 
    ih = 1.0 / h;
    ihsq = ih * ih;
    ihc = ih * ih * ih;
  
    wfd = h;
    if (domain->dimension == 3) {
          W = 2.08890862808 * ihc ;  
          d[i] += W; 
        } else {
          W = 2.08890862808 * ihc * W * wfd * wfd * ihc * ihc ; //THIS IS A PLACEHOLDER WRONG
          d[i] += W; 
          }
    
    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];     
      jtype = type[j];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];

      rsq = delx*delx + dely*dely + delz*delz;
      

      if (rsq < cutsq[itype][jtype]) {
        
        h = cut[itype][jtype]; 
        ih = 1.0 / h;
        ihsq = ih * ih;
        ihc = ih * ih * ih;

        wfd = h - sqrt(rsq);
        W = h + 3*sqrt(rsq);
        if (domain->dimension == 3) {
          W = 2.08890862808 * W * wfd * wfd * wfd * ihc * ihc * ih;  
        } else {
          W = 2.08890862808 * ihc * W * wfd * wfd * ihc * ihc ;
        }
        d[i] += W;
        if (newton_pair || j < nlocal) {
          d[j] += W;
        }
        }
      }
    }
  if (newton_pair) comm->reverse_comm(this);
  comm->forward_comm(this);

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    vxtmp = v[i][0];
    vytmp = v[i][1];
    vztmp = v[i][2];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];
    imass = mass[itype];
    
    power_term = pow(np* d[i],2./3.);
    exp_term = exp((2.*entropy[i]/kb/np-5.)/3.);
    ei = 3.*h_p*h_p*np*np/4./pi/imass * power_term * exp_term;
    
    Ti = 2./3./kb/np * ei;
    Pi = d[i] * Ti * kb * np;
    temperature[i] = Ti;
    esph[i] = ei;

    for (jj = 0; jj < jnum; jj++) {

      j = jlist[jj];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];

      rsq = delx * delx + dely * dely + delz * delz;
      jtype = type[j];
      jmass = mass[jtype];

      if (rsq < cutsq[itype][jtype]) {
        h = cut[itype][jtype];
        ih = 1.0 / h;
        ihsq = ih * ih;
        
        double r = sqrt (rsq);
        if (delx == 0.) {

          ex = 0;
        } else {
          ex = delx / r;
        }
        if (dely == 0.) {
          ey = 0;
        } else {
          ey = dely / r;
        }     
        if (delz == 0.) {
          ez = 0;
        } else {
          ez = delz / r;
        } 


        wfdist = h - r;
        dim = domain -> dimension;
        if (dim == 3) {
          D = 3.0;
          // Lucy Kernel, 3d
          // This F function in the espanol paper
          F = 25.066903536973515383e0 * wfdist * wfdist * ihsq * ihsq * ihsq * ih;
          // Regular lucy function not the derivative
        } else {
          // Lucy Kernel, 2d
          wfd = -19.098593171027440292e0 * wfd * wfd * ihsq * ihsq * ihsq;
        }

        power_term = pow(np* d[j],2./3.);
        exp_term = exp((2.*entropy[j]/kb/np-5.)/3.);
        
        ej = 3.*h_p*h_p*np*np/4./pi/jmass * power_term * exp_term;
        Tj = 2./3. / kb / np * ej;
        Pj = d[j] * Tj * kb*np;
        comb_temp = (Ti * Tj) /(Ti + Tj);
        added_temp = Ti + Tj;
        
        velx= vxtmp - v[j][0] ;
        vely= vytmp - v[j][1];
        velz= vztmp - v[j][2];

        VijdotEij = ex * velx + ey * vely + ez * velz;
        fdd = F / (d[i]* d[j]);
        
        //double D = D;
        // Stochastic constants
        Aij = sqrt(8. * comb_temp * kb * (((5./3.0) * visc - eta))*fdd);
        gradAij_Ti = 8*kb*fdd*((5./3.0) * visc - eta)* Tj * Tj / added_temp / added_temp;
        gradAij_Tj = 8*kb*fdd*((5./3.0) * visc - eta)* Ti * Ti / added_temp / added_temp;
        //aij = (Aij * Aij / (8. * kb)) *(1./Ti + 1./Tj);
        
        Bij = sqrt(8. * comb_temp * kb * (((5./3.) * visc + 8. * eta))*fdd);
        gradBij_Ti = 8*kb*fdd*((5./3.) * visc + 8. * eta)* Tj * Tj / added_temp / added_temp;
        gradBij_Tj = 8*kb*fdd*((5./3.) * visc + 8. * eta)* Ti * Ti / added_temp / added_temp;
        //bij = (Bij * Bij / (12. * kb)) *(1./Ti  + 1./Tj);

        Cij = sqrt(4. * therm * kb * Ti * Tj * fdd);
        gradCij_Ti = 4*kb*therm*Tj*fdd;
        gradCij_Tj = 4*kb*therm*Ti*fdd;
        //aux = A_ij**2/2 * v_ij + (A_ij**2/2 + (B_ij**2 - A_ij**2) / self.D) * dot(e_ij, v_ij) * e_ij

        aux_x = Aij * Aij / 2 * velx + (Aij * Aij / 2 + (Bij * Bij - Aij * Aij) / D) * VijdotEij * ex;
        aux_y = Aij * Aij / 2 * vely + (Aij * Aij / 2 + (Bij * Bij - Aij * Aij) / D) * VijdotEij * ey;
        aux_z = Aij * Aij / 2 * velz + (Aij * Aij / 2 + (Bij * Bij - Aij * Aij) / D) * VijdotEij * ez;
        
        prefactor_x = -1.0 * aux_x * (1/cv[i]/Ti + 1/cv[j]/Tj); //term1
        prefactor_y = -1.0 * aux_y * (1/cv[i]/Ti + 1/cv[j]/Tj);
        prefactor_z = -1.0 * aux_z * (1/cv[i]/Ti + 1/cv[j]/Tj);

        aux2_x = -aux_x*(1/Ti +1/Tj);
        aux2_y = -aux_y*(1/Ti +1/Tj);
        aux2_z = -aux_z*(1/Ti +1/Tj);

        //MgradS_v 
        fviscx = 0.25/kb * aux2_x;
        fviscy = 0.25/kb * aux2_y;
        fviscz = 0.25/kb * aux2_z;
        //MgradS_v
        
        fvisc2x = fvisc2y = fvisc2z = 0.0;
        fvisc2x += (gradAij_Ti/2. * velx + (gradAij_Ti/2. + ((gradBij_Ti - gradAij_Ti)/ D)) * VijdotEij * ex) / cv[i];
        fvisc2y += (gradAij_Ti/2. * vely + (gradAij_Ti/2. + ((gradBij_Ti - gradAij_Ti)/ D)) * VijdotEij * ey) / cv[i];
        fvisc2z += (gradAij_Ti/2. * velz + (gradAij_Ti/2. + ((gradBij_Ti - gradAij_Ti)/ D)) * VijdotEij * ez) / cv[i];

        fvisc2x += (gradAij_Tj/2. * velx + (gradAij_Tj/2. + (gradBij_Tj - gradAij_Tj)/ D) * VijdotEij * ex) / cv[j];
        fvisc2y += (gradAij_Tj/2. * vely + (gradAij_Tj/2. + (gradBij_Tj - gradAij_Tj)/ D) * VijdotEij * ey) / cv[j];
        fvisc2z += (gradAij_Tj/2. * velz + (gradAij_Tj/2. + (gradBij_Tj - gradAij_Tj)/ D) * VijdotEij * ez) / cv[j];

        //term1 =
        fvisc2x = -0.25/kb * (fvisc2x + prefactor_x);
        fvisc2y = -0.25/kb * (fvisc2y + prefactor_y);
        fvisc2z = -0.25/kb * (fvisc2z + prefactor_z);

        fpgradx = delx * F * (Pi / (d[i]*d[i]) + Pj / (d[j]*d[j]));
        fpgrady = dely * F * (Pi / (d[i]*d[i]) + Pj / (d[j]*d[j]));
        fpgradz = delz * F * (Pi / (d[i]*d[i]) + Pj / (d[j]*d[j]));

        vdotv = velx * velx + vely * vely + velz * velz;

        // random force and entropy calculation
        // independent increments of a Wiener process matrix
        wiener[0][0] = random->gaussian(0.0,1.0);
        wiener[1][1] = random->gaussian(0.0,1.0);
        wiener[2][2] = random->gaussian(0.0,1.0);
        
        wiener[0][1] = random->gaussian(0.0,1.0);
        wiener[0][2] = random->gaussian(0.0,1.0);
        wiener[1][2] = random->gaussian(0.0,1.0);
        
        wiener[2][0] = random->gaussian(0.0,1.0);
        wiener[2][1] = random->gaussian(0.0,1.0);
        wiener[1][0] = random->gaussian(0.0,1.0);

        dV = random->gaussian(0.0,1.0);
        
        trace = wiener[0][0] + wiener[1][1] + wiener[2][2];
        sym_wiener[0][0] = 0.5 * (wiener[0][0] + wiener[0][0]) - 1./3. * trace;  
        sym_wiener[1][1] = 0.5 * (wiener[1][1] + wiener[1][1]) - 1./3. * trace;
        sym_wiener[2][2] = 0.5 * (wiener[2][2] + wiener[2][2]) - 1./3. * trace;
        
        sym_wiener[1][0] = sym_wiener[0][1] = 0.5 * (wiener[0][1] + wiener[1][0]);
        sym_wiener[2][1] = sym_wiener[1][2] = 0.5 * (wiener[1][2] + wiener[2][1]);
        sym_wiener[2][0] = sym_wiener[0][2] = 0.5 * (wiener[0][2] + wiener[2][0]);

        f_rand[0] = (Aij * sym_wiener[0][0] + Bij * (1./3.) * trace) * ex + (Aij * sym_wiener[0][1]) * ey + (Aij * sym_wiener[0][2]) * ez;
        f_rand[1] = (Aij * sym_wiener[1][0] ) * ex + (Aij * sym_wiener[1][1] + Bij * (1./3.) * trace) * ey + (Aij * sym_wiener[1][2]) * ez;
        f_rand[2] = (Aij * sym_wiener[2][0]) * ex + (Aij * sym_wiener[2][1] ) * ey + (Aij * sym_wiener[2][2] + Bij * (1./3.) * trace) * ez;

        f_random[0] =  f_rand[0];
        f_random[1] =  f_rand[1];
        f_random[2] =  f_rand[2];
        //f_random[0] = f_random[1] = f_random[2] = 0.0;

        term_rand = -0.5 * (f_rand[0] * velx + f_rand[1] * vely + f_rand[2] * velz) + Cij * dV;
        //term_rand = 0.0;

        f[i][0] += fpgradx + fviscx + fvisc2x + sqrt(dtinv) * jmass * f_random[0];
        f[i][1] += fpgrady + fviscy + fvisc2y + sqrt(dtinv) * jmass * f_random[1];
        f[i][2] += fpgradz + fviscz + fvisc2z + sqrt(dtinv) * jmass * f_random[2];

        //f[i][0] += fviscx + fvisc2x; //-fviscx - fvisc2x; //- fviscx - fvisc2x + sqrt(dtinv) * jmass * f_rand[0];
        //f[i][1] += fviscy + fvisc2y; //-fviscy - fvisc2y; //- fviscy - fvisc2y + sqrt(dtinv) * jmass * f_rand[1];
        //f[i][2] += fviscz + fvisc2z; //-fviscz - fvisc2z; //- fviscz - fvisc2z + sqrt(dtinv) * jmass * f_rand[2];
        
        //MgradS_S 
        term1 = (Aij*Aij/2 * vdotv + (Aij*Aij/2 + (Bij*Bij - Aij*Aij)/D) * VijdotEij * VijdotEij) /4.0;
        term2 = (1/Ti + 1/Tj)*term1 + (1/Ti - 1/Tj)*Cij*Cij;
        //MgradS_S 

        //divM_S
        term3 = - 1.0 * (2./ cv[i]/Ti + 1./ cv[j]/Tj) * term1;
        term4 = (gradAij_Ti/2. * vdotv + (gradAij_Ti/2. + (gradBij_Ti - gradAij_Ti)/ D) * VijdotEij * VijdotEij) / cv[i] /4.0;

        term5 = (gradAij_Tj/2. * vdotv + (gradAij_Tj/2. + (gradBij_Tj - gradAij_Tj)/ D) * VijdotEij * VijdotEij) / cv[j]/4.0;

        term6 = - 1.0 * (2./ cv[i]/Ti - 1./ cv[j]/Tj) * Cij * Cij;

        term7 = gradCij_Ti / cv[i] - gradCij_Tj / cv[j];

        term9 = - 1.0 * ((D + 1.) * Aij * Aij / 2. + (Bij * Bij - Aij * Aij) / D);   

        //dentropy[i] += term9/imass;
        dentropy[i] +=  0.5/kb/Tj * term2 + 0.5 * (term3 + term4 +term5 +term6+ term7 +term9/imass)/Ti + sqrt(dtinv) * term_rand;


        if (newton_pair || j < nlocal) {
        
          fvisc2x = fvisc2y = fvisc2z = 0.0;
          fvisc2x += (gradAij_Tj/2. * velx + (gradAij_Tj/2. + (gradBij_Tj - gradAij_Tj)/ D) * VijdotEij * ex) / cv[j];
          fvisc2y += (gradAij_Tj/2. * vely + (gradAij_Tj/2. + (gradBij_Tj - gradAij_Tj)/ D) * VijdotEij * ey) / cv[j];
          fvisc2z += (gradAij_Tj/2. * velz + (gradAij_Tj/2. + (gradBij_Tj - gradAij_Tj)/ D) * VijdotEij * ez) / cv[j];

          fvisc2x += (gradAij_Ti/2. * velx + (gradAij_Ti/2. + (gradBij_Ti - gradAij_Ti)/ D) * VijdotEij * ex) / cv[i];
          fvisc2y += (gradAij_Ti/2. * vely + (gradAij_Ti/2. + (gradBij_Ti - gradAij_Ti)/ D) * VijdotEij * ey) / cv[i];
          fvisc2z += (gradAij_Ti/2. * velz + (gradAij_Ti/2. + (gradBij_Ti - gradAij_Ti)/ D) * VijdotEij * ez) / cv[i];

        //term1 =
          fvisc2x = -0.25/kb * (fvisc2x + prefactor_x);
          fvisc2y = -0.25/kb * (fvisc2y + prefactor_y);
          fvisc2z = -0.25/kb * (fvisc2z + prefactor_z);
          
          f[j][0] -= fpgradx + fviscx + fvisc2x + sqrt(dtinv) * imass * f_random[0];
          f[j][1] -= fpgrady + fviscy + fvisc2y + sqrt(dtinv) * imass * f_random[1];
          f[j][2] -= fpgradz + fviscz + fvisc2z + sqrt(dtinv) * imass * f_random[2];

          //f[j][0] -= fviscx + fvisc2x; //- fviscx - fvisc2x + sqrt(dtinv) * imass * f_rand[0];
          //f[j][1] -= fviscy + fvisc2y; //- fviscy - fvisc2y + sqrt(dtinv) * imass * f_rand[1];
          //f[j][2] -= fviscz + fvisc2z; //- fviscz - fvisc2z + sqrt(dtinv) * imass * f_rand[2];

        //term1 = (Aij*Aij/2 * vdotv + (Aij*Aij/2 + (Bij*Bij - Aij*Aij)/D) * VijdotEij * VijdotEij) /4.0;
        term1 = (Aij*Aij/2 * vdotv + (Aij*Aij/2 + (Bij*Bij - Aij*Aij)/D) * VijdotEij * VijdotEij) /4.0;
        term2 = (1/Tj + 1/Ti)*term1 + (1/Tj - 1/Ti)*Cij*Cij;
        //MgradS_S 

        //divM_S
        term3 = - 1.0 * (2./ cv[j]/Tj + 1./ cv[i]/Ti) * term1;
        term4 = (gradAij_Tj/2. * vdotv + (gradAij_Tj/2. + (gradBij_Tj - gradAij_Tj)/ D) * VijdotEij * VijdotEij) / cv[j]/4.0;

        term5 = (gradAij_Ti/2. * vdotv + (gradAij_Tj/2. + (gradBij_Ti - gradAij_Ti)/ D) * VijdotEij * VijdotEij) / cv[i]/4.0;

        term6 = - 1.0 * (2./ cv[j]/Tj - 1./ cv[i]/Ti) * Cij * Cij;

        term7 = gradCij_Tj / cv[j] - gradCij_Ti / cv[i];

        dentropy[j] +=  0.5/kb/Tj *term2 + 0.5 * (term3 +term4+term5+term6+term7+ term9/jmass)/Tj + sqrt(dtinv) * term_rand;
        }
      }

    }
  }
}

/* ----------------------------------------------------------------------
 allocate all arrays
 ------------------------------------------------------------------------- */

void PairSDPDML::allocate () {
  allocated = 1;
  int n = atom->ntypes;

  memory->create (setflag, n + 1, n + 1, "pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create (cutsq, n + 1, n + 1, "pair:cutsq");
  memory->create (cut, n + 1, n + 1, "pair:cut");
}

/* ----------------------------------------------------------------------
 global settings
 ------------------------------------------------------------------------- */

void PairSDPDML::settings (int narg, char **arg) {
  if (narg != 2)
    error->all (FLERR, "Illegal number of arguments for "
                "pair_style sdpd/ideal");
  
  // seed is immune to underflow/overflow because it is unsigned
  seed = comm->nprocs + comm->me + atom->nlocal;
  
  int seed = utils::inumeric(FLERR, arg[0], false, lmp);
  double num_per_part = utils::inumeric(FLERR, arg[1], false, lmp);
  np = num_per_part;

  if (seed <= 0) error->all(FLERR,"Invalid random number seed");

  delete random;

  random = new RanMars(lmp,(seed + comm->me) % 900000000);
}

/* ----------------------------------------------------------------------
 set coeffs for one or more type pairs
 ------------------------------------------------------------------------- */

void PairSDPDML::coeff (int narg, char **arg) {
  if (narg != 6)
    error->all (FLERR, "Incorrect args for pair_style "
                "sdpd/ideal coefficients");
  if (!allocated) allocate();


  int ilo, ihi, jlo, jhi;
  utils::bounds(FLERR, arg[0], 1, atom->ntypes, ilo, ihi, error);
  utils::bounds(FLERR, arg[1], 1, atom->ntypes, jlo, jhi, error);


  double cut_one = utils::numeric(FLERR,arg[2], false, lmp);
  double visc_one = utils::numeric(FLERR,arg[3], false, lmp);
  double therm_one = utils::numeric(FLERR,arg[4], false, lmp);
  double eta_one = utils::numeric(FLERR,arg[5], false, lmp);


  therm = therm_one;
  eta = eta_one;
  visc = visc_one;
  self_cut = cut_one;
  if (cut_one <= 0) error->all (FLERR, "Cutoff must be positive");

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      printf("setting cut[%d][%d] = %f\n", i, j, cut_one);
      cut[i][j] = cut_one;
      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0)
    error->all(FLERR,"Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
 init specific to this pair style
------------------------------------------------------------------------- */

void PairSDPDML::init_style()
{
  if ((!atom->dentropy_flag) || (atom->dentropy == nullptr))
    error->all(FLERR,"Pair style sdpd/ideal requires atom attributes entropy and dentropy");

  neighbor->add_request(this);
}

/* ----------------------------------------------------------------------
 init for one type pair i,j and corresponding j,i
 ------------------------------------------------------------------------- */

double PairSDPDML::init_one (int i, int j) {
  if (setflag[i][j] == 0)
    error->all(FLERR,"Not all pair sdpd/ideal coeffs are set");

  cut[j][i] = cut[i][j];
  cutforcesq = cut[i][j] * cut[i][j];
  return cut[i][j];
}


void *PairSDPDML::extract_peratom(const char *str, int &ncol)
{
  if (strcmp(str,"d") == 0) {
    ncol = 0;
    return (void *) d;
  }

  return nullptr;
}

int PairSDPDML::pack_reverse_comm(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) buf[m++] = d[i];
  return m;
}

void PairSDPDML::unpack_reverse_comm(int n, int *list, double *buf)
{
  int i,j,m;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    d[j] += buf[m++];
  }
}

int PairSDPDML::pack_forward_comm(int n, int *list, double *buf, int /*pbc_flag*/, int * /*pbc*/)
{
  int i,j,m;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    buf[m++] = d[j];
  }
  return m;
}

void PairSDPDML::unpack_forward_comm(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) d[i] = buf[m++];
}

double PairSDPDML::memory_usage()
{
  double bytes = 2* (double)nmax * sizeof(double);
  return bytes;
}