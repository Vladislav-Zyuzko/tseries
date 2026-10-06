// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, Vladislav Zyuzko
//
// Replacement for ctsa's `minfit()` (upstream `src/lls.c`), whose code was
// tagged MIT upstream without any named copyright holder and has been removed
// from the vendored tree. See third_party/ctsa/PROVENANCE.md.
//
// Contract (as documented at the call site, lls_svd2() in lls.c):
//   AB  max(M, N) x (N + P), row-major. On input AB = [A : B] with A M x N and
//       B M x P (rows M..max(M,N)-1, if any, are ignored and treated as 0).
//       On output the first N rows hold [V : C], where A = U diag(q) V' is the
//       singular value decomposition and C = U' B.
//   q   N singular values.
// Returns 0 on success, or the non-zero code of the SVD.
//
// Implemented on top of ctsa's own BSD-3-Clause `svd()` (matrix.c, Golub &
// Reinsch, Num. Math. 14, 403-420 (1970)) rather than as a second SVD: A is
// padded with zero rows to at least N rows (which leaves its singular values
// and the least-squares solution unchanged) because svd() needs M >= N.
// Singular values come back sorted by svd(); lls_svd2() does not depend on
// their order.

#include <stdlib.h>
#include <string.h>

#include "lls.h"
#include "matrix.h"

int minfit(double *AB, int M, int N, int P, double *q) {
  const int rows = M > N ? M : N;
  const int np = N + P;
  double *A, *U, *V, *C;
  int i, j, k, r, ret;

  if (M <= 0 || N <= 0 || P < 0) return -1;
  A = (double *)calloc((size_t)rows * N, sizeof(double));
  U = (double *)malloc(sizeof(double) * (size_t)rows * N);
  V = (double *)malloc(sizeof(double) * (size_t)N * N);
  C = (double *)malloc(sizeof(double) * ((size_t)N * P + 1));
  if (A == NULL || U == NULL || V == NULL || C == NULL) {
    free(A);
    free(U);
    free(V);
    free(C);
    return -1;
  }

  for (i = 0; i < M; ++i) {
    memcpy(A + (size_t)i * N, AB + (size_t)i * np, sizeof(double) * N);
  }

  ret = svd(A, rows, N, U, V, q);
  if (ret == 0) {
    // C = U' B, using only the M real rows of B.
    for (i = 0; i < N; ++i) {
      for (k = 0; k < P; ++k) {
        double s = 0.0;
        for (r = 0; r < M; ++r) {
          s += U[(size_t)r * N + i] * AB[(size_t)r * np + N + k];
        }
        C[(size_t)i * P + k] = s;
      }
    }
    for (i = 0; i < N; ++i) {
      for (j = 0; j < N; ++j) AB[(size_t)i * np + j] = V[(size_t)i * N + j];
      for (k = 0; k < P; ++k) {
        AB[(size_t)i * np + N + k] = C[(size_t)i * P + k];
      }
    }
  }

  free(A);
  free(U);
  free(V);
  free(C);
  return ret;
}
