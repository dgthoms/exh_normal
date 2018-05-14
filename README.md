# exh_normal
Version 0.0.0

Exhaustive search of GF(2^n) over GF(2) for normal bases of low complexity

This code is a single-node, multi-threaded version of our exhaustive search. Each run is performed on a Task, which is a contiguous region of outer loop indices. Each Task is split into \<number of threads\> subtasks and run synchronously one-per-thread. At the conclusion of its subtask, the thread returns the smallest complexity element found to master, which prints (to stdout) a summary of the information.

It is not too hard to alter the existing code into a primitive MPI program. Another option is to write a simple controlling function (e.g., using MPI4Py) to perform adaptive Tasking and job management. We may eventually provide some of these functionalities if there is enough interest.

Note that this program is highly compute-bound with negligble memory usage. The synchronous execution of threads in subtasks can create overhead on the Tasks with the largest amount of compute if the task size is too large. There are various work-arounds, but the simplest is to ensure that number of tasks is large enough that this overhead is minimal.

## Getting Started

To get started, we assume you have a built binary (see [Building](#building)) as ./foo. To run the search algorithm, use
```
  ./foo <n> <number of threads> <number of tasks> <task index>
```

To run the entire search on a single thread, use
```
  ./foo <n> 1 1 0
```

To run the search for n=34 on 2 cores/threads using 100 Tasks on a single node and store the output of each task to the file 34-\<task index\>-of-100, use:
```
  for i in {0..99}; do ./foo 34 2 100 ${i} >>34-${i}-of-100; done
```

### Prerequisites
This package has the following pre-requisites:

* [NTL](http://www.shoup.net/ntl)
    - Compile with thread-safety (NTL_THREADS=on) in configuration file
    - Highly recommended to check and tune your build of NTL

### Building
The easiest way to build this package is to compile directly on the command line. This can be done using a derivative of the following command:
```
  g++ -o foo ./nb-2017.cpp -lm -lntl -lpthread
```
If NTL is installed in a non-standard location (e.g., /home/dgthoms/NTL), then specifying the -L and -I flags may be necessary; e.g.,:
```
  g++ -o foo ./nb-2017.cpp -lm -lntl -lpthread -L/home/dgthoms/NTL/lib -I/home/dgthoms/NTL/include
```
If NTL is installed with GMP back-end, then ensure the -lgmp flag is passed (and library and include directories are indicated, if necessary). Note that compiling NTL with GMP is not required for this program, though it is highly recommended for performance for general NTL programs requiring big integers and multiprecision.

### Profiling
We recommend using the [gperftools](https://github.com/gperftools) (originally Google Performance Tools) for profiling. Once gperftools is installed, simply add the
```
  -g
```
flag and the
```
  -lprofiler
```
linker flag to your compile line. See the documentation [here](https://github.com/gperftools/gperftools) for more information. Note that linking the profiler does not incur a performance hit (aside from passing the -g flag, which downgrades compiler optimizations) unless the environment variable $CPUPROFILE is set.

## Contributing
All contributions, comments and suggestions welcome!

## Authors
* Lucia Moura (University of Ottawa)
* Daniel Panario (Carleton University)
* David Thomson (Carleton University) - **Main developer**

## License and Disclaimer
Coming soon.

The original authors do not hold any responsibilty for the use of this code or any of its derivatives.

## Acknowledgements
* We would like to thank Ariane Masuda for her contributions to the original paper which kicked off this work in 2008, and to David Marquis for some initial discussions when we considered rebooting this search.
* We would also like to thank Leland McInnes and Chris North for some helpful discussions on parallel programming.

## References
* L. Moura, D. Panario and D. Thomson, Normal basis exhaustive search: 10 years later, Accepted to the [2018 WAIFI conference](http://www.waifi.org).
* A. Masuda, L. Moura, D. Panario and D. Thomson, Low complexity normal elements over finite fields of characteristic two, IEEE Transactions on Computers, volume 57 (2008), 990-1001.
