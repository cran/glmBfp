/*
 *
 * Definition of types used in the code.
 *
 */

#ifndef TYPES_H_
#define TYPES_H_

#include <set>
#include <vector>

#define ARMA_BLAS_CAPITALS
#define ARMA_LAPACK_CAPITALS
#include <RcppArmadillo.h>


// Integer types ***************************************

// general signed integer
typedef int Int;
// general unsigned integer
typedef unsigned int PosInt;
// general large unsigned integer
typedef unsigned long int PosLargeInt;



// STL Set types ***************************************

// a set of ints
typedef std::set<Int> IntSet;
// a set of unsigned ints
typedef std::set<PosInt> PosIntSet;


// STL Vector types ************************************

// (1) Double vectors:

// here is a double vector
typedef std::vector<double> MyDoubleVector;
// and a long double vector
typedef std::vector<long double> LongDoubleVector;
// and a string vector
typedef std::vector<std::string> StrVector;

// (2) Integer vectors:

// here is an integer vector
typedef std::vector<Int> IntVector;
// and an unsigned int vector
typedef std::vector<PosInt> PosIntVector;
// the machine precision
static const double EPS = sqrt(DOUBLE_EPS);


// Armadillo types ************************************

// a vector
typedef arma::colvec AVector;
// a matrix
typedef arma::mat AMatrix;


// Special types ***************************************

// array of double vectors
typedef std::vector<std::vector<AVector> > AVectorArray;


// one power set
typedef std::multiset<Int> Powers;
// vector of powers
typedef std::vector<Powers> PowersVector;


// the "index" type
typedef LongDoubleVector::size_type Ind;
// one index set
typedef std::set<Ind> IndSet;






#endif /* TYPES_H_ */
