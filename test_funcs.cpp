#include "El.hpp"
#include "el_petsc_utils.hpp"
#include <petscvec.h>

int main(int argc, char* argv []){
	PetscInt rankp;
	//PetscInt N = 30; //Size of vectors to use
	//bool print = false;
	static char help[] = "Converts between Elemental and Petsc vectors.\n\n";

	// Initialize mpi
	El::Initialize(argc,argv);
	PetscInitialize(&argc,&argv,(char*)0,help);
	
	// Input args
	const int N      = Input("--N","size of vector",30);
	const bool print = Input("--p","print vectors?",false);
	
	// Finish up with inputs
	El::ProcessInput();
	El::PrintInputReport();
	
	// MPI_Comm
	MPI_Comm comm = MPI_COMM_WORLD;
	El::mpi::Comm el_comm = El::mpi::Comm(comm);
	MPI_Comm_rank(comm,&rankp);
	int ranke = El::mpi::WorldRank();

	int mpicomms = El::mpi::Size(el_comm);
	const El::Grid grid(el_comm);
	
	// Check that ranks are allocated the same way
	int rank;
	if(rankp != ranke){
		if(ranke == 0){
			std::cout << "Ranks are not distributed in the same way! Be aware" << std::endl;
		}
	}else{
		rank = ranke;
	}

	
	// Print grid and ownership
	//DistMatrix<double,VR,STAR> test(grid);
	//test.Resize(N,1);
	//for(int p = 0; p<mpicomms;p++){
	//	if(rank == p){std::cout << "Process: " << p << " (" << grid.Row() << ", " <<grid.Col() << ")" <<std::endl;}
	//	MPI_Barrier(comm);
	//}
	//for(int p = 0; p<mpicomms;p++){
	//	if(rank == p){std::cout << "VR: " << grid.VRRank() << " VC: " <<grid.VCRank()<<std::endl;}
	//	MPI_Barrier(comm);
	//}
	//if(rank == 2){
	//	for(int i=0;i<N;i++){
	//		std::cout << test.Owner(i,0) <<std::endl;
	//	}
	//}
	

	//-------- Test Petsc --> El ---------//
	// Create Petsc vec of 0,1,2 ..
	PetscScalar *array;
	PetscInt nlocal;
	Vec x_pt;

	VecCreate(PETSC_COMM_WORLD,&x_pt);
	VecSetSizes(x_pt,PETSC_DECIDE,N);
	VecSetFromOptions(x_pt);
	VecSet(x_pt,0.0);
	VecAssemblyBegin(x_pt);
	VecAssemblyEnd(x_pt);

	PetscInt offset;
	VecGetOwnershipRange(x_pt,&offset,NULL);
	VecGetLocalSize(x_pt,&nlocal);
	VecGetArray(x_pt,&array);
	for(int i=0;i<nlocal;i++){array[i] = i + offset;}
	VecRestoreArray(x_pt,&array);

	// View to make sure
	if(print){VecView(x_pt,PETSC_VIEWER_STDOUT_WORLD);}
	
	// Create el vec
	El::DistMatrix<double,VC,STAR> x_el(grid);
	x_el.Resize(N,1);
	El::Fill(x_el,1.5);

	// Call function
	double start = El::mpi::Time(); 
	Petsc2El_vec(x_pt,x_el);
	double p2e_time = El::mpi::Time() - start;
	
	// Get times
	double max_p2e_time;
	El::mpi::Reduce(&p2e_time,&max_p2e_time,1,El::mpi::MAX,0,El::mpi::COMM_WORLD);
	PetscPrintf(PETSC_COMM_WORLD,"Petsc -> El time: %g \n", max_p2e_time);

	// Print vec
	if(print){El::Print(x_el);}


	// ---------- Test El --> Petsc ----------- //
	// Pt vec
	Vec x_pt2;
	VecCreate(PETSC_COMM_WORLD,&x_pt2);
	VecSetSizes(x_pt2,PETSC_DECIDE,N);
	VecSetFromOptions(x_pt2);
	VecSet(x_pt2,0.0);
	VecAssemblyBegin(x_pt2);
	VecAssemblyEnd(x_pt2);

	// Call function
	start = El::mpi::Time();
	El2Petsc_vec(x_el,x_pt2);	
	double e2p_time = El::mpi::Time() - start;
	
	// Get times
	double max_e2p_time;
	El::mpi::Reduce(&e2p_time,&max_e2p_time,1,El::mpi::MAX,0,El::mpi::COMM_WORLD);
	PetscPrintf(PETSC_COMM_WORLD,"El -> Petsc time: %g \n", max_e2p_time);
	
	// Print vec
	if(print){VecView(x_pt2,PETSC_VIEWER_STDOUT_WORLD);}



	// ------------ Report errors  ----------------//
	VecAXPY(x_pt,-1.0,x_pt2);
	PetscReal norm;
	VecNorm(x_pt,NORM_2,&norm);
	if(norm >-PETSC_SMALL && norm < PETSC_SMALL){norm = 0.0;}
	PetscPrintf(PETSC_COMM_WORLD, "Petsc -> El -> Petsc norm: %g \n",(double)norm);
	
	
	El::Finalize();

	PetscFinalize();
}
