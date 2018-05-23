#include <cstdint> /* int*_t types */
#include <cstdlib> /* atoi */
#include <iostream> /* for console output */
#include <fstream> /* for file I/O */

#include <NTL/GF2.h>
#include <NTL/GF2E.h>
#include <NTL/GF2X.h>
#include <NTL/GF2XFactoring.h>
#include <NTL/mat_GF2.h>
#include <NTL/vec_GF2E.h>
#include <NTL/vec_vec_GF2E.h>
#include <NTL/BasicThreadPool.h>

//This version of the code is a non-MPIed version.
#ifndef NO_MPI
#define NO_MPI
#endif

//In the future may provide 2 versions depending
//on the NO_MPI flag
#ifndef NO_MPI
#include <mpi.h> /* For extra-node parallelism */
#endif

#include <pthread.h> /* For intra-node parallelism */

//Deprecated: used in the initial pthreads version.
typedef struct thread_data {
	int n;
	int MPI_rank;
	int MPI_size;
	int thread_id;
	int NUM_THREADS;
} thread_data_t;

//Return data for main outer loop
typedef struct return_data {
	int status;
	long min_cplex;
	NTL::GF2XModulus modulus;
	NTL::GF2E min_element;
} return_data_t;

//Key function definitions
uint64_t integer_representation(NTL::GF2E* a);
NTL::vec_vec_GF2E precompute(int* traces, const int n);
int is_lex_first_int_rep(NTL::GF2E* current, const int n);
int calculate_complexity(NTL::GF2E *current, const int n);
int graycode_unrank(NTL::GF2E* current, unsigned long long index, const int n);
return_data_t exhaustive_search(const int n, unsigned long long first, unsigned long long last);


//Default trailing zeroes function. Uses the builting __builtin_ctzll
//inline int get_graycode_index(uint64_t c);
//Placeholder: to be defined if compiler builtin does not exist.
inline int __deprecated_get_graycode_index2(long c);
