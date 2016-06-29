/*
 * linalgInterface.h
 *
 *  Created on: 10.12.2010
 *      Author: daniel
 */

#ifndef LINALGINTERFACE_H_
#define LINALGINTERFACE_H_

#include <types.h>


// triangular solve of L * x = R
// where L can be provided lower or upper-triangular and be transposed.
// the solution is directly written into R.
void
trs(const bool upper,
    const bool transpose,
    const AMatrix& L,
    AMatrix& R);

// Cholesky decomposition of (symmetric) A,
// which can either be provided in lower or upper-triangular storage.
// Result will be in the same layout.
// Returns an error code which should be zero if all went well.
int
potrf(const bool upper,
      AMatrix& A);

// Solve the system LL' * x = R
// Instead of lower-triangular L also the upper-triangular form can be provided.
// The solution is directly written into R.
// Returns an error code which should be zero if all went well.
int
potrs(const bool upper,
      const AMatrix& L,
      AMatrix& R);


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
     AMatrix& C);

// Do a triangular matrix-vector multiplication:
// x := A*x,   or   x := A'*x,
//
// where x is an n element vector and  A is an n by n upper or lower triangular matrix.
void
trmv(const bool upper,
     const bool transpose,
     const AMatrix& A,
     AVector& x);


#endif /* LINALGINTERFACE_H_ */
