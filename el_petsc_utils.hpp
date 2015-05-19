#include "El.hpp"
#include <petscvec.h>
#include <vector>

using namespace El;

/*
 * Takes an elemental distributed vector and places it into a 
 * petsc distributed vector. For now, the el vector MUST be  
 * VC,STAR distributed, and the petsc vector is assumed to have 
 * memory stored in contiguous blocks
 */
void El2Petsc_vec(DistMatrix<double,VC,STAR>& el_vec, Vec& pt_vec);

/*
 * Takes a petsc distributed vector and places it into an 
 * el distributed vector. For now, the el vector MUST be  
 * VC,STAR distributed, and the petscvector is assumed to have 
 * memory stored in contiguous blocks
 */
void Petsc2El_vec(Vec& pt_vec,DistMatrix<double,VC,STAR>& el_vec);

/*
 * Simple exclusive scan routine for displacements
 */
std::vector<int> exscan(std::vector<int> x);


