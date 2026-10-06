// SPDX-License-Identifier: BSD-3-Clause
/*
 * lls.c
 *
 *  Created on: Apr 14, 2014
 *      Author: Rafat Hussain
 */

#include "lls.h"

int lls_normal(double *A,double *b,int M,int N,double *x) {
	int retcode,sc,i,j,c1,l;
	double sum;
	double *AT,*d,*C,*y,*CT;
	// M - data points
	// N - Number of parameters
	// A - MXN; b size - M vector; AT - NXM

	AT = (double*) malloc(sizeof(double) * M * N);
	d = (double*) malloc(sizeof(double) * N);
	C = (double*) malloc(sizeof(double) * N * N);
	y = (double*) malloc(sizeof(double) * N);
	CT = (double*) malloc(sizeof(double) * N * N);

	retcode = 0;

	mtranspose(A,M,N,AT);
	mmult(AT,b,d,N,M,1);
	mmult(AT,A,C,N,M,N);

	sc = chol(C,N);
	if (sc == -1) {
		return -1;
	}

	mtranspose(C,N,N,CT);
	//Forward Substitution

	y[0] = d[0]/CT[0];
	for(i = 1; i < N; ++i) {
		sum = 0.;
		c1 = i*N;
		for(j = 0; j < i; ++j) {
			sum += y[j] * CT[c1 + j];
		}
		y[i] = (d[i] - sum)/CT[c1+i];
	}

	//Back Substitution

	x[N - 1] = y[N - 1]/C[N * N - 1];

	for (i = N - 2; i >= 0; i--) {
		sum = 0.;
		c1 = i*(N+1);
		l=0;
		for(j = i+1; j < N;j++) {
			l++;
			sum += C[c1 + l] * x[j];
		}
		x[i] = (y[i] - sum) / C[c1];
	}

	free(AT);
	free(d);
	free(C);
	free(y);
	return retcode;
}

static int lls_qr_determined(double *Ai,double *bi,int M,int N,double *xo) {
	int j,i,k,u,t,retcode,c1,l;
	double *x,*v,*AT,*w,*bvec,*b,*A,*R;
	double beta,sum;

	retcode = 0;

	if (M < N) {
			printf("M should be greater than or equal to N");
			exit(-1);
	}
	x = (double*) malloc(sizeof(double) * M);
	b = (double*) malloc(sizeof(double) * M);
	bvec = (double*) malloc(sizeof(double) * N);
	v = (double*) malloc(sizeof(double) * M);
	AT = (double*) malloc(sizeof(double) * M * N);
	A = (double*) malloc(sizeof(double) * M * N);
	w = (double*) malloc(sizeof(double) * N);
	R = (double*) malloc(sizeof(double) * N * N);

	for(j = 0; j < M;++j) {
		b[j] = bi[j];
	}
	for(j = 0; j < M*N;++j) {
		A[j] = Ai[j];
	}

	for(j = 0; j < N;++j) {
		for(i=j;i < M;++i) {
			x[i-j] = A[i*N+j];

		}

		beta = house(x,M-j,v);
		bvec[j] = beta;

		for (i=j; i < M; i++) {
			t = i * N;
			u = 0;
			for (k=j; k < N; k++) {
				AT[u+i-j] = A[k+t];
				u+=(M-j);

			}

		}


		mmult(AT,v,w,N-j,M-j,1);
		scale(w,N-j,1,beta);
		mmult(v,w,AT,M-j,1,N-j);
		for (i=j; i < M; i++) {
			t = i *N;
			for (k=j; k < N; k++) {
				A[t+k] -= AT[(i-j)*(N-j) + k - j];
			}
		}
		if (j < M) {

			for(i=j+1;i < M;++i) {
				A[i*N+j] = v[i-j];
			}
		}

	}

	for(i = 0; i < N;++i) {
		t = i *N;
		for(j = 0; j < N;++j) {
			if (i > j) {
				R[t+j] = 0.;
			} else {
				R[t+j] = A[t+j];
			}
		}
	}

	//Calculate Rank of R

	retcode = rank(R,N,N);


	for(j = 0; j < N;++j) {
		v[j] = 1;
		for(i = j+1; i < M;++i) {
			v[i] = A[i * N + j];//edit
		}
		mmult(b+j,v+j,w,1,M-j,1);
		*w = *w * bvec[j];
		for(i = j; i < M;++i) {
			v[i] = *w * v[i];
		}
		for(i = j; i < M;++i) {
			b[i] = b[i] - v[i];
		}


	}

	//mdisplay(b,1,M);

	//back substitution

	xo[N - 1] = b[N - 1]/R[N * N - 1];

	for (i = N - 2; i >= 0; i--) {
		sum = 0.;
		c1 = i*(N+1);
		l=0;
		for(j = i+1; j < N;j++) {
			l++;
			sum += R[c1 + l] * xo[j];
		}
		xo[i] = (b[i] - sum) / R[c1];
	}

	free(x);
	free(v);
	free(AT);
	free(w);
	free(bvec);
	free(R);
	free(b);
	free(A);

	return retcode;
}

static int lls_qr_undetermined(double *Ai,double *bi,int M,int N,double *xo) {
	int retval, i,c1,l,j;
	double *A, *bvec, *Q, *R, *z, *b;
	double sum;

	retval = 0;
	A = (double*) malloc(sizeof(double) * N * M);
	bvec = (double*) malloc(sizeof(double) * M);
	Q = (double*) malloc(sizeof(double) * N * M);
	R = (double*) malloc(sizeof(double) * M * M);
	z = (double*) malloc(sizeof(double) * M);
	b = (double*) malloc(sizeof(double) * M);

	for(j = 0; j < M;++j) {
		b[j] = bi[j];
	}

	mtranspose(Ai,M,N,A);
	// A NXM
	qrdecomp(A,N,M,bvec);

	getQR(A,N,M,bvec,Q,R);

	retval = rank(R,M,M);


	itranspose(R,M,M);

	// Forward Substitution

	z[0] = b[0] / R[0];

	for(i = 1; i < M;i++) {
		sum = 0.;
		c1 = i*M;
		l = i*M+i;

		for(j = 0; j < i;j++) {
			sum += R[c1 + j] * z[j];
		}
		z[i] = (b[i] - sum) / R[l];
	}


	mmult(Q,z,xo,N,M,1);

	free(A);
	free(bvec);
	free(Q);
	free(R);
	free(z);
	free(b);
	return retval;
}

int lls_qr(double *Ai,double *bi,int M,int N,double *xo) {
	int ret;

	if (M >= N) {
		ret = lls_qr_determined(Ai,bi,M,N,xo);
	} else {
		ret = lls_qr_undetermined(Ai,bi,M,N,xo);
	}

	return ret;
}

void bidiag(double *A, int M, int N) {
	int j,i,k,u,t;
	double *x,*v,*AT,*w;
	double beta;

	if (M < N) {
			printf("M should be greater than or equal to N");
			exit(1);
	}
	x = (double*) malloc(sizeof(double) * M);
	v = (double*) malloc(sizeof(double) * M);
	AT = (double*) malloc(sizeof(double) * M * N);
	w = (double*) malloc(sizeof(double) * M * M);


	for(j = 0; j < N;++j) {
		for(i=j;i < M;++i) {
			x[i-j] = A[i*N+j];

		}

		beta = house(x,M-j,v);

		for (i=j; i < M; i++) {
			t = i * N;
			u = 0;
			for (k=j; k < N; k++) {
				AT[u+i-j] = A[k+t];
				u+=(M-j);

			}

		}


		mmult(AT,v,w,N-j,M-j,1);
		scale(w,N-j,1,beta);
		mmult(v,w,AT,M-j,1,N-j);
		for (i=j; i < M; i++) {
			t = i *N;
			for (k=j; k < N; k++) {
				A[t+k] -= AT[(i-j)*(N-j) + k - j];
			}
		}

		for(i=j+1;i < M;++i) {
			A[i*N+j] = v[i-j];
		}

		if (j < N - 2) {
			for(i=j+1;i < N;++i) {
				x[i-j-1] = A[j*N+i];
			}
			beta = house(x,N-j-1,v);

			for (i=j; i < M; i++) {
				t = i * N;
				u = (i-j) *(N-j-1);
				for (k=j+1; k < N; k++) {
					AT[u+k-j-1] = A[k+t];
				}
			}

			mmult(AT,v,w,M-j,N-j-1,1);
			scale(w,M-j,1,beta);
			mmult(w,v,AT,M-j,1,N-j-1);

			for (i=j; i < M; i++) {
				t = i * N;
				u = (i-j) *(N-j-1);
				for (k=j+1; k < N; k++) {
					A[k+t] -= AT[u+k-j-1];
				}
			}
			u = 1;
			t = j*N;
			for(i = j + 2; i < N;++i) {
				A[t+i] = v[u];
				u++;
			}

		}

	}


	free(x);
	free(v);
	free(AT);

}

void bidiag_orth(double *A, int M, int N,double *U,double *V) {
	int j,i,k,u,t;
	double *x,*v,*AT,*w;
	double beta;

	if (M < N) {
			printf("M should be greater than or equal to N");
			exit(1);
	}
	x = (double*) malloc(sizeof(double) * M);
	v = (double*) malloc(sizeof(double) * M);
	AT = (double*) malloc(sizeof(double) * M * M);
	w = (double*) malloc(sizeof(double) * M * M);


	eye(U,M);
	eye(V,N);

	for(j = 0; j < N;++j) {
		for(i=j;i < M;++i) {
			x[i-j] = A[i*N+j];

		}

		beta = house(x,M-j,v);
		//mdisplay(v,M-j,1);

		for (i=j; i < M; i++) {
			t = i * N;
			u = 0;
			for (k=j; k < N; k++) {
				AT[u+i-j] = A[k+t];
				u+=(M-j);

			}

		}


		mmult(AT,v,w,N-j,M-j,1);
		scale(w,N-j,1,beta);
		mmult(v,w,AT,M-j,1,N-j);
		for (i=j; i < M; i++) {
			t = i *N;
			for (k=j; k < N; k++) {
				A[t+k] -= AT[(i-j)*(N-j) + k - j];
			}
		}

		for(i=j+1;i < M;++i) {
			A[i*N+j] = v[i-j];
		}


		for (i=j; i < M; i++) {
			t = i * M;
			u = (i-j) *(M-j);
			for (k=j; k < M; k++) {
				AT[u+k-j] = U[k+t];
			}
		}

		mmult(AT,v,w,M-j,M-j,1);
		scale(w,M-j,1,beta);
		mmult(w,v,AT,M-j,1,M-j);

		for (i=j; i < M; i++) {
			t = i * M;
			u = (i-j) *(M-j);
			for (k=j; k < M; k++) {
				U[k+t] -= AT[u+k-j];
			}
		}

		//mdisplay(U,M,M);

		if (j < N - 2) {
			for(i=j+1;i < N;++i) {
				x[i-j-1] = A[j*N+i];
			}
			beta = house(x,N-j-1,v);

			for (i=j; i < M; i++) {
				t = i * N;
				u = (i-j) *(N-j-1);
				for (k=j+1; k < N; k++) {
					AT[u+k-j-1] = A[k+t];
				}
			}

			mmult(AT,v,w,M-j,N-j-1,1);
			scale(w,M-j,1,beta);
			mmult(w,v,AT,M-j,1,N-j-1);

			for (i=j; i < M; i++) {
				t = i * N;
				u = (i-j) *(N-j-1);
				for (k=j+1; k < N; k++) {
					A[k+t] -= AT[u+k-j-1];
				}
			}
			u = 1;
			t = j*N;
			for(i = j + 2; i < N;++i) {
				A[t+i] = v[u];
				u++;
			}

		}

	}


	free(x);
	free(v);
	free(AT);

}


/* tseries patch (license): svd_gr(), svd_gr2() and minfit() stood here
 * upstream, tagged "SPDX-License-Identifier: MIT" with no copyright holder
 * named anywhere. svd_gr/svd_gr2 had no callers and were dropped; minfit (only
 * reached by lls_svd2() below, i.e. regress() with llsmethod "svd", which
 * nothing in ctsa or tseries selects) is replaced by a BSD-3-Clause wrapper
 * over ctsa's own svd() (matrix.c) in
 * third_party/ctsa_replacements/minfit.c. Declarations in lls.h are kept.
 * See third_party/ctsa/PROVENANCE.md, "Local modifications". */

int lls_svd2(double *Ai,double *bi,int M,int N,double *xo) {
	int retcode,P,mnmax,np,i,j,t,t2;
	double *AB,*q,*V,*C;
	double eps;

	if (M > N) {
		mnmax = M;
	} else {
		mnmax = N;
	}
	retcode = 0;
	P = 1;
	np = N + P;
	eps = macheps();

	AB = (double*) malloc(sizeof(double) * mnmax * np);
	q = (double*) malloc(sizeof(double) * N);
	C = (double*) malloc(sizeof(double) * N);
	V = (double*) malloc(sizeof(double) * N * N);

	for(i = 0; i < M;++i) {
		t = i * N;
		t2 = i * np;
		for(j = 0; j < N;++j) {
			AB[t2 + j] = Ai[t+j];
		}
	}

	for(i = 0; i < M;++i) {
		t2 = i * np;
		AB[t2+N] = bi[i];
	}

	minfit(AB,M,N,P,q);

	for(i = 0; i < N;++i) {
		t = i * N;
		t2 = i * np;
		for(j = 0; j < N;++j) {
			V[t+j] = AB[t2 + j];
		}
	}

	for(i = 0; i < N;++i) {
		t = i *N;
		for(j = 0; j < N;++j) {
			if (fabs(q[j]) > eps) {
				V[t+j] /= q[j];
			} else {
				V[t+j] = 0.0;
			}
		}
	}

	for(i = 0; i < N;++i) {
		t2 = i * np;
		C[i] = AB[t2+N];
	}
	mmult(V,C,xo,N,N,1);

	free(AB);
	free(q);
	free(V);
	free(C);
	return retcode;
}

int lls_svd(double *Ai,double *bi,int M,int N,double *xo) {
	int retcode,i,j,t;
	double *U,*V,*q,*C;
	double eps;

	U = (double*) malloc(sizeof(double) * M*N);
	V = (double*) malloc(sizeof(double) * N*N);
	q = (double*) malloc(sizeof(double) * N);
	C = (double*) malloc(sizeof(double) * N);

	retcode = 0;
	eps = macheps();

	retcode = svd(Ai,M,N,U,V,q);;
	if (retcode != 0) {
		return retcode;
	}

	mmult(bi,U,C,1,M,N);

	for(i = 0; i < N;++i) {
		t = i *N;
		for(j = 0; j < N;++j) {
			if (fabs(q[j]) > eps) {
				V[t+j] /= q[j];
			} else {
				V[t+j] = 0.0;
			}
		}
	}

	mmult(V,C,xo,N,N,1);

	free(U);
	free(V);
	free(q);

	return retcode;
}

