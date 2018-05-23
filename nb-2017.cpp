#include "nb-2017.h" //Function definitions and structures
#include <time.h>
#include <chrono>

thread_local 	NTL::Vec< NTL::Vec<NTL::GF2> > vec_P;
thread_local	NTL::Vec< NTL::Vec<NTL::GF2> > vec_B;
thread_local	NTL::mat_GF2 P;
thread_local	NTL::mat_GF2 B;
thread_local	NTL::GF2 det;
thread_local	NTL::GF2E bi;
thread_local	NTL::GF2X rep;

thread_local	int BYTES_PER_INT = NTL_BITS_PER_INT/8;

// This function accepts a GF2E element expressed in a basis,
// and the second argument can be dereferenced to find an
// integer representation of that element.

uint64_t integer_representation(NTL::GF2E *a) {
	// This is a direct copy of sage's integer representation
	// calculator for GF2E elements

	uint64_t i = 0;
	uint64_t ret = 0;
	int shift = 0;
	rep = a->_GF2E__rep; //This is a (necesary) copy

	if (NTL::IsZero(rep)) {
		return 0;
	}

	while(NTL::deg(rep) >= NTL_BITS_PER_INT) {
		NTL::BytesFromGF2X( (unsigned char*)&i, rep, BYTES_PER_INT);
		ret += i<<shift;
		shift += NTL_BITS_PER_INT;
		NTL::RightShift(rep, rep, NTL_BITS_PER_INT);
	}
	NTL::BytesFromGF2X((unsigned char*)&i, rep, BYTES_PER_INT);
	ret += i<<shift;
	return ret;
}

// This function returns the matrix e_i^{2^j} for 0 <= i,j <= n-1
// and places Tr(e_j) in the j-th element of the first argument.
NTL::vec_vec_GF2E precompute(int *traces, const int n) {
	int i,j;
	NTL::GF2X rep;
	NTL::vec_vec_GF2E stdbasis;
	stdbasis.SetLength(n);

	for(i = 0; i < n; i++) {
		rep = NTL::GF2X(NTL::INIT_MONO, i); // Initialize rep = x^i
		stdbasis[i].SetLength(n);
		NTL::conv(stdbasis[i][0], rep);
		for (j = 1; j < n; j++) {
			stdbasis[i][j] = NTL::sqr(stdbasis[i][j-1]);
		}
		conv(traces[i], NTL::trace(stdbasis[i][0]));
	}

	return stdbasis;
}

//Returns 1 if the current element ordering is canonical Otherwise returns 0;
int is_lex_first_int_rep(NTL::GF2E *current, const int n) {
	int i;
	uint64_t test, next;
int check = 0;
	//read integer representation of a polynomial
	test = integer_representation(current);
	for(i = 1; i < n; i++) {
		next = integer_representation(current+i);
		if (next >= test) {
			return 0;
		}
	}
	return 1;
}

//Calculates the complexity of the putative basis defined by
//the elements of current. Returns the (positive) complexity of the
//basis if the putative basis is indeed a basis, otherwise returns 0;
int calculate_complexity(NTL::GF2E *current, const int n) {
	int cplex = 0;
	int i,j;

	//std::cout << "Constructing the matrix P" << std::endl; //VERBOSE
	//Construct the matrix P = (current_ij)
	for (i = 0; i < n; i++) {
		//convert from GF2X to vec_GF2 -- a necessary copy
		NTL::VectorCopy(*(vec_P.elts()+i), (current + i)->_GF2E__rep, n);
		//for (j = 0; j < n; j++) { //SLOW
		//	vec_P[i][j] = current[i]._GF2E__rep[j];
		//}
	}

	//NTL::MakeMatrix(P, vec_P); //SLOW: the NTL `safe' way - makes a copy
	P._mat__rep.swap(vec_P); //Hack, cheaper (less safe) than MakeMatrix

	//Stores P^{-1} = P if P is invertible.
	//std::cout << "Done. Computing P_inv" << std::endl; //VERBOSE
	NTL::inv(det, P, P);
	//std::cout << "Done. Checking determinant" << std::endl; //VERBOSE

	//If P has determinant 0, break out of the loop
	if (NTL::IsZero(det)) {
		//std::cout << "Determinant is " << det << std::endl; //VERBOSE
		return 0;
	}

	//Construct the matrix B = [ current[0]*current[i]_j ]

	//Save a multiplication because current[0]current[0] = current[1]
	//std::cout << "Done. Constructing the matrix B" << std::endl; //VERBOSE
	NTL::VectorCopy(*vec_B.elts(), (current + 1)->_GF2E__rep, n);
	//for (j = 0; j < n; j++) vec_B[0][j] = current[1]._GF2E__rep[j]; //SLOW

	for (i = 1; i < n; i++) {
		NTL::mul(bi, *current, *(current+i));
		NTL::VectorCopy(*(vec_B.elts()+i), bi._GF2E__rep, n); //this is a copy.
		//for (j = 0; j < n; j++) vec_B[i][j] = bi._GF2E__rep[j]; //SLOW
	}

	//NTL::MakeMatrix(B,vec_B); //SLOW the NTL `safe' way - makes a copy
	B._mat__rep.swap(vec_B); //Hack -- cheaper than MakeMatrix

	//Set B = BP^{-1}
	//std::cout << "Done, Computing B in the normal basis." << std::endl; //VERBOSE
	NTL::mul(B,B,P);
	//std::cout << "Done. Counting the complexity." << std::endl; //VERBOSE
	for (i = 0; i < n; i++) {
		//Point rather than use accessor
		//cplex += NTL::weight(*(B._mat__rep.elts()+i)); //SLOWish
		//cplex += NTL::weight(B[i]); //SLOW
		//	This is a total hack
		//	B._mat__rep.elts()           is a pointer to row 0
		//	                  +i         is a pointer to row i
		//	                    .rep     is the underlying WordVector
		//	                        .rep is the _ntl_ulong
		//	and dereference everything to
		//
		cplex += __builtin_popcountll( *(*(B._mat__rep.elts()+i)).rep.rep );
	}
	//std::cout << "Done. Complexity is: " << cplex << std::endl;

//Should probably only be used with 1 thread for file concurrency!!
#if defined(ALL_TABLES) || defined(ALL_BASES)
	char str[256];
	sprintf(str, "all-normals-%d", n);
	std::ofstream outfile;
	outfile.open(str, std::ios::app);
	outfile << current[0] << " " << cplex << std::endl;
#if defined(ALL_TABLES)
  	outfile << B << std::endl;
#endif
	outfile.close();
#endif
	return cplex;
}

 /*
 * TODO: Check for builtin compiler intrinsics and use the correct
 *       function only if found. Not sure of a way to do this other
 *       than to check for compiler versions and write a proper
 *       configuration script.
 */

//Counts the number of trailing zero bits of a 64 bit number
inline int get_graycode_index(unsigned long c) { return __builtin_ctzl(c); }

inline int __deprecated_get_graycode_index2(long c) {
	if (c == 0) return 0;
	int zero_bits = 64;
	//c &= -signed long (c);
	if (c) zero_bits--;
	if (c & 0x00000000FFFFFFFF) zero_bits -= 32;
	if (c & 0x0000FFFF0000FFFF) zero_bits -= 16;
	if (c & 0x00FF00FF00FF00FF) zero_bits -= 8;
	if (c & 0x0F0F0F0F0F0F0F0F) zero_bits -= 4;
	if (c & 0x3333333333333333) zero_bits -= 2;
	if (c & 0x5555555555555555) zero_bits -= 1;

	return zero_bits;
}

int graycode_unrank(NTL::GF2E *current, long index, const int n) {
	//Use the code from Kreher-Stinson GrayCodeUnrank
	//to place current[0] = graycode_unrank(index)
	//and current[j] = sqr(current[j-1]) for j > 0
	//
	//Write index-1 = (b_n-1, b_n-2, ..., b_1, b_0)
	//unranked: c_n-1 = b_n-1, c_j = b_j XOR b_j+1
	long i;
	uint64_t b, bb = 0;
	NTL::GF2X T;
	T.SetLength(n);
	//T = 0;

	if (index > 1ULL<<n) std::cout << "WTF" << std::endl;

	for (i = n-1; i >= 0; i--) {
		b = index>>i; // get 2^i-the coefficient of b
		if (b != bb) {
			NTL::SetCoeff(T, i);
		} else {
			NTL::SetCoeff(T,i,0); //XOR
		}
		bb = b; //bb becomes b_j+1
		index = index - b * (1ULL<<i); //Set 2^i-th coeff to 0, iterate
	}

	if (T == 0) {
		for (i = 0; i < n; i++) *(current+i) = 0;
		return 0;
	}


	(*current)._GF2E__rep = T;

	for (i = 1; i < n; i++) {
		*(current+i) = NTL::sqr(*(current + i-1));
	}

	return 0;
}

return_data_t exhaustive_search(const int n, long first, long last) {
	//required local variables
	int i;
	long c;
	int gc_idx;
	return_data_t ret;

 	//Persistent intermediate values
	int *traces; // traces[i] = Tr(e_i)
	traces = (int*)calloc(n, sizeof(int));
	uint64_t *complexity_histogram;
	complexity_histogram = (uint64_t*)calloc(n*n, sizeof(uint64_t));

	NTL::GF2E min_element;
	int min_complexity = n*n;
	NTL::vec_vec_GF2E stdbasis; // stdbasis[i][j] = e_i^{q^j}
	NTL::GF2X modulus;

	//Transient values
	int current_complexity = n*n;
	int current_trace;

	//NTL has a table of sparse irreducibles, so this is a constant
	//for any n
	NTL::vec_GF2E current;

	//Build GF(2^n)
	NTL::BuildSparseIrred(modulus, n);
	NTL::GF2E::init(modulus);

	//Initialize sizes of local structures
	current.SetLength(n);
	//Initialize sizes of thread_local global structures
	vec_P.SetLength(n);
	vec_B.SetLength(n);
	P.SetDims(n,n);
	B.SetDims(n,n);
	rep.SetLength(n);
	bi._GF2E__rep.SetLength(n);
	//Don't forget substructures
	for (i = 0; i < n; i++) {
		current[i]._GF2E__rep.SetLength(n);
		vec_P[i].SetLength(n);
		vec_B[i].SetLength(n);
	}

	//std::cout << "Precomputing" << std::endl; //VERBOSE
	stdbasis = precompute(traces,n);

	//Initialize the search and checkpoint data
	graycode_unrank(current.elts(), first, n);
	conv(current_trace, NTL::trace(current[0]));

	if (!first) first++; //Element starts at 0; index starts at 1

	//Run the main loop (for this rank)
	for (c = first; c < last; c++) {
		//Update the current element considered and all its conjugates (and its trace)
		//std::cout << "Done. Checking trace." << std::endl; //VERBOSE
		//Normal elements must have trace 1
		if(current_trace) { //Cutdown by a factor of 2
			if (is_lex_first_int_rep(current.elts(), n)) { //Time cutdown by a factor of n
				current_complexity = calculate_complexity(current.elts(), n);
				if (current_complexity > 0) {
					//Keep a histogram of complexities seen.
					//Currently unused, but may be useful in
					//the future and essentially free
					complexity_histogram[current_complexity]++;
					if (current_complexity < min_complexity) {
						min_complexity = current_complexity;
						min_element = current[0];
					} // End minimality check
				} //End complexity check
			} // End canonical element check
		} // End trace check
		gc_idx = get_graycode_index(c);
		NTL::add(current, current, stdbasis[gc_idx]);
		current_trace ^= traces[gc_idx];

	} //End main outer loop

	//Cleanup
	free(complexity_histogram);
	ret.min_cplex = min_complexity;
	ret.min_element = min_element;
	ret.modulus = modulus;
	ret.status = 0;
	return ret;
}

int main(int argc, char **argv) {
	const int n = std::atoi(argv[1]);
	const int num_pthreads = std::atoi(argv[2]);
	const int num_tasks = std::atoi(argv[3]);
	const int task_index = std::atoi(argv[4]);

	if (argc != 5) {
		std::cout << "Invalid execution: <n> <num_threads> <num_tasks> <task_index>" << std::endl << std::flush;
		return 1;
	}

	if (task_index >= num_tasks) {
		std::cout << "Task index must be smaller than number of tasks! Exiting." << std::endl << std::flush;
		return 1;
	}

	//Instruct NTL to use this many threads
	NTL::SetNumThreads(num_pthreads);
	//Set up return data
	return_data_t ret;
	//Indices - Save the name "rank" from legacy MPI code
	long rank_first, rank_last;
	int i;

	NTL::PartitionInfo global_pinfo(1ULL<<n, num_tasks);
	//Since the task_index is user-input, can already set up rank_first and last
	global_pinfo.interval(rank_first, rank_last, task_index);

	std::cout << "n: " << n << std::endl << std::flush;
	std::cout << "Task index: " << task_index << " of " << num_tasks << " tasks." << std::endl << std::flush;
	std::cout << "Number of threads: " << num_pthreads << std::endl << std::flush;

	//Use NTL's system clock
	double my_time = NTL::GetTime();
	//Use std::chrono's wall-time
	auto start = std::chrono::system_clock::now();

	//Set up thread partition info
	NTL::PartitionInfo thread_pinfo(rank_last - rank_first, num_pthreads);

	NTL_EXEC_INDEX(num_pthreads, thread_index)
		//Execute each thread synchronously
		long thread_first, thread_last;
		thread_pinfo.interval(thread_first, thread_last, thread_index);
		thread_first += rank_first; thread_last += rank_first;

		ret = exhaustive_search(n, thread_first, thread_last);
		//Not the most graceful outputting, but effective: Print to stdout
		//Always remember to flush: each thread is printing
		std::cout << "n: " << n << ", Modulus: " << ret.modulus << std::endl << std::flush;
		std::cout << "Min element: " << ret.min_element << ", Min complexity: " <<
			ret.min_cplex << std::endl << std::flush;
		std::cout << "Complexity: " << ret.min_cplex << " Task index: " << task_index <<
			" Task First: " << rank_first << " Task Last: " << rank_last <<
			" Thread index: " << thread_index << " Thread First: " << thread_first <<
			" Thread Last: " << thread_last << std::endl << std::flush;
		//End thread execution
	NTL_EXEC_INDEX_END

	//Collect runtime information

	my_time = NTL::GetTime() - my_time;
	auto end = std::chrono::system_clock::now();
	std::chrono::duration<double> elapsed_seconds = end-start;
	std::cout << "Wall time on task: " << task_index << ": " << elapsed_seconds.count() << std::endl << std::flush;
	std::cout << "CPU time on task: " << task_index << ": " << my_time << std::endl << std::flush;
	//Success!
	return 0;
}
