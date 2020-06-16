/*
 * linalgInterface.cpp
 *
 *  Created on: 10.12.2010
 *      Author: daniel
 */


#define USE_FC_LEN_T
#include <Rconfig.h>

#include <cassert>
#include <linalgInterface.h>

#include <R_ext/Lapack.h>
#include <R_ext/BLAS.h>

#ifndef FCONE
# define FCONE
#endif

// constants:
static const double doubleOne = 1.0;
static const int intOne = 1;



// triangular solve of L * x = R
// where L can be provided lower or upper-triangular and be transposed.
// the solution is directly written into R.
void
  trs(const bool upper,
      const bool transpose,
      const AMatrix& L,
      AMatrix& R)
  {
    const char* side = "L";
    const char* uplo = upper ? "U" : "L";
    const char* transa = transpose ? "T" : "N";
    const char* diag = "N";
    
    const int m = R.n_rows;
    const int n = R.n_cols;
    assert(static_cast<PosInt>(m) == L.n_rows);
    
    F77_CALL(dtrsm)(side,
             uplo,
             transa,
             diag,
             & m,
             & n,
             & doubleOne,
             L.memptr(),
             & m,
             R.memptr(),
             & m FCONE FCONE FCONE FCONE);
  }

// Cholesky decomposition of (symmetric) A,
// which can either be provided in lower or upper-triangular storage.
// Result will be in the same layout.
// Returns an error code which should be zero if all went well.
int
  potrf(const bool upper,
        AMatrix& A)
  {
    const char* uplo = upper ? "U" : "L";
    const int n = A.n_rows;
    assert(A.is_square());
    int info = 0;
    
    F77_CALL(dpotrf)(uplo,
             & n,
             A.memptr(),
             & n,
             & info FCONE);
             
             return info;
  }

// Solve the system LL' * x = R
// Instead of lower-triangular L also the upper-triangular form can be provided.
// The solution is directly written into R.
// Returns an error code which should be zero if all went well.
int
  potrs(const bool upper,
        const AMatrix& L,
        AMatrix& R)
  {
    const char* uplo = upper ? "U" : "L";
    const int n = L.n_rows;
    assert(L.is_square());
    assert(static_cast<PosInt>(n) == R.n_rows);
    const int nrhs = R.n_cols;
    int info = 0;
    
    F77_CALL(dpotrs)(uplo,
             & n,
             & nrhs,
             L.memptr(),
             & n,
             R.memptr(),
             & n,
             & info FCONE);
             
             return info;
  }

// Do a symmetric rank k update of C:
// C := A * A' + beta * C
// where C must be symmetric and can be provided in lower or upper-triangular storage.
// Result will be in the same layout.
// If transpose == true, then the crossproduct will be added, i.e., C := A' * A + beta*C.
void
  syrk(const bool upper,
       const bool transpose,
       const AMatrix& A,
       const double beta,
       AMatrix& C)
  {
    const char* uplo = upper ? "U" : "L";
    const char* trans = transpose ? "T" : "N";
    const int n = C.n_rows;
    assert(C.is_square());
    const int k = transpose ? A.n_rows : A.n_cols;
    const int lda = transpose ? k : n;
    
    F77_CALL(dsyrk)(uplo,
             trans,
             & n,
             & k,
             & doubleOne,
             A.memptr(),
             & lda,
             & beta,
             C.memptr(),
             & n FCONE FCONE);
  }

// Do a triangular matrix-vector multiplication:
// x := A*x,   or   x := A'*x,
//
// where x is an n element vector and  A is an n by n upper or lower triangular matrix.
void
  trmv(const bool upper,
       const bool transpose,
       const AMatrix& A,
       AVector& x)
  {
    const char* uplo = upper ? "U" : "L";
    const char* trans = transpose ? "T" : "N";
    const char* diag = "N";
    const int n = A.n_rows;
    
    F77_CALL(dtrmv)(uplo,
             trans,
             diag,
             & n,
             A.memptr(),
             & n,
             x.memptr(),
             & intOne FCONE FCONE FCONE);
  }
