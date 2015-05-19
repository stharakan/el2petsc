#include "el_petsc_utils.hpp"

template<typename T>
std::ostream &operator<<(std::ostream &stream, std::vector<T> ob)
{
	for(int i=0;i<ob.size();i++){
	  stream << ob[i] << '\n';
	}
	return stream;
}


std::vector<int> exscan(std::vector<int> x){
	int ll = x.size();
	std::vector<int> y(ll);
	y[0] = 0;

	for(int i=1;i<ll;i++){
		y[i] = y[i-1] + x[i-1];
	}
	
	return y;
}

void El2Petsc_vec(El::DistMatrix<double,VC,STAR>& el_vec,Vec& pt_vec){
	PetscInt nlocal, nstart; // petsc vec info
	PetscScalar *pt_array,*pt_perm_array;
	int r,q,ll,rq; // el vec info
	int nbigs; //Number of large recv (i.e. recv 1 extra data point)
	int pstart; // p_id of nstart
	int p = El::mpi::WorldRank(); //p_id
	int recv_size; // base recv size
	bool print = p == -1; 

	// Get el vec info
	ll = el_vec.Height();
	const El::Grid* g = &(el_vec.Grid());
	r = g->Height();
	q = g->Width();
	MPI_Comm comm = (g->Comm()).comm;
	
	// Get petsc vec params
	VecGetLocalSize(pt_vec,&nlocal);
	VecGetArray(pt_vec,&pt_array);
	VecGetOwnershipRange(pt_vec,&nstart,NULL);

	// Determine who owns the first element we want
	rq = r * q;
	pstart = nstart % rq;
	nbigs = nlocal % rq;
	recv_size = nlocal / rq;
	
	if(print){
		std::cout << "r: " << r << " q: " << q <<std::endl;
		std::cout << "nstart: " << nstart << std::endl;
		std::cout << "ps: " << pstart << std::endl;
		std::cout << "nbigs: " << nbigs << std::endl;
		std::cout << "recv_size: " << recv_size << std::endl;
	}

	// Make recv sizes
	std::vector<int> recv_lengths(rq);
	std::fill(recv_lengths.begin(),recv_lengths.end(),recv_size);
	if(nbigs >0){
		for(int i=0;i<nbigs;i++){
			recv_lengths[(pstart + i) % rq] += 1;
		}
	}

	// Make recv disps
	std::vector<int> recv_disps = exscan(recv_lengths);

	// All2all to get send sizes
	std::vector<int> send_lengths(rq);
	MPI_Alltoall(&recv_lengths[0], 1, MPI_INT, &send_lengths[0], 1, MPI_INT,comm);

	// Scan to get send_disps
	std::vector<int> send_disps = exscan(send_lengths);

	// Do all2allv to get data on correct processor
	std::vector<double> recv_data(nlocal);
	MPI_Alltoallv(el_vec.Buffer(),&send_lengths[0],&send_disps[0],MPI_DOUBLE, \
			&recv_data[0],&recv_lengths[0],&recv_disps[0],MPI_DOUBLE,comm);
	
	if(print){
		//std::cout << "Send data: " <<std::endl << *el_vec.Buffer() <<std::endl;
		std::cout << "Send lengths: " <<std::endl << send_lengths <<std::endl;
		std::cout << "Send disps: " <<std::endl << send_disps <<std::endl;
		std::cout << "Recv data: " <<std::endl << recv_data <<std::endl;
		std::cout << "Recv lengths: " <<std::endl << recv_lengths <<std::endl;
		std::cout << "Recv disps: " <<std::endl << recv_disps <<std::endl;
	}
	
	// Allocate the memory? just so we only do it once
	for(int p=0;p<rq;p++){
		int base_idx = (p - pstart + rq) % rq;
		int offset = recv_disps[p];
		for(int i=0;i<recv_lengths[p];i++){
			pt_array[base_idx + rq*i] = recv_data[offset + i];
		}
	}
	if(print){std::cout <<"here?"<<std::endl;}

	// Copy into array
	VecRestoreArray(pt_vec,&pt_array);
}

void Petsc2El_vec(Vec& pt_vec,El::DistMatrix<double,VC,STAR>& el_vec){
	PetscInt nlocal,nstart,gsize; //local elements, start p_id, global size
	PetscScalar *pt_array; // will hold local array
	int r,q,rq; //Grid sizes
	int nbigs; //Number of large sends (i.e. send 1 extra data point)
	int pstart; // p_id of nstart
	int p = El::mpi::WorldRank(); //p_id
	int send_size; // base send size
	bool print = p == -1; 


	// Get Grid and associated params
	const El::Grid* g = &(el_vec.Grid());
	r = g->Height();
	q = g->Width();
	MPI_Comm comm = (g->Comm()).comm;

	// Get sizes, array in petsc 
	VecGetSize(pt_vec,&gsize);
	VecGetLocalSize(pt_vec,&nlocal);
	VecGetArray(pt_vec,&pt_array);
	VecGetOwnershipRange(pt_vec,&nstart,NULL);

	//Find processor that nstart belongs to, number of larger sends
	rq = r * q;
	pstart = nstart % rq; //int div
	nbigs = nlocal % rq;
	send_size = nlocal/rq;
	
	if(print){
		std::cout << "r: " << r << " q: " << q <<std::endl;
		std::cout << "nstart: " << nstart << std::endl;
		std::cout << "ps: " << pstart << std::endl;
		std::cout << "nbigs: " << nbigs << std::endl;
		std::cout << "send_size: " << send_size << std::endl;
	}

	// Make send_lengths
	std::vector<int> send_lengths(rq);
	std::fill(send_lengths.begin(),send_lengths.end(),send_size);
	if(nbigs >0){
		for(int j=0;j<nbigs;j++){
			send_lengths[(pstart + j) % rq] += 1;
		}
	}

	// Make send_disps
	std::vector<int> send_disps = exscan(send_lengths);

	// Make send_data
	std::vector<double> send_data(nlocal);
	for(int proc=0;proc<rq;proc++){
		int offset = send_disps[proc];
		int base_idx = (proc - pstart + rq) % rq; 
		for(int j=0; j<send_lengths[proc]; j++){
			int idx = base_idx + (j * rq);
			send_data[offset + j] = pt_array[idx];
		}
	}

	// Do all2all to get recv_lengths
	std::vector<int> recv_lengths(rq);
	MPI_Alltoall(&send_lengths[0], 1, MPI_INT, &recv_lengths[0], 1, MPI_INT,comm);

	// Scan to get recv_disps
	std::vector<int> recv_disps = exscan(recv_lengths);

	// Do all2allv to get data on correct processor
	double * recv_data = el_vec.Buffer();
	MPI_Alltoallv(&send_data[0],&send_lengths[0],&send_disps[0],MPI_DOUBLE, \
			&recv_data[0],&recv_lengths[0],&recv_disps[0],MPI_DOUBLE,comm);

	if(print){
		std::cout << "Send data: " <<std::endl << send_data <<std::endl;
		std::cout << "Send lengths: " <<std::endl << send_lengths <<std::endl;
		std::cout << "Send disps: " <<std::endl << send_disps <<std::endl;
		std::cout << "Recv data: " <<std::endl << recv_data <<std::endl;
		std::cout << "Recv lengths: " <<std::endl << recv_lengths <<std::endl;
		std::cout << "Recv disps: " <<std::endl << recv_disps <<std::endl;
	}
}


