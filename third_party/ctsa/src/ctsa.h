// SPDX-License-Identifier: BSD-3-Clause
/*
* arma.h
*
*  Created on: Sep 11, 2014
*      Author: Rafat Hussain
*/

#ifndef CTSA_H_
#define CTSA_H_

#include "emle.h"
#include "errors.h"  // tseries patch (license): was "autoutils.h" (removed GPL cluster); errors.h is what it pulled in transitively

#define MAX_STR_LEN 10

#ifdef __cplusplus
extern "C" {
#endif
	/*
	ARIMA vs SARIMA
	Typically you should use SARIMA if you are calculating seasonal models and ARIMA for non-seasonal models. However, ARIMA wil not work if the number of parameters p and q exceed 100.
	For these extreme-case non-seasonal models, use SARIMA and set seasonal paramters to zero. SARIMA and ARIMA work identically for non-seasonal models otherwise.
	*/

typedef struct arima_set* arima_object;

arima_object arima_init(int p, int d, int q, int N);

struct arima_set{
	int N;// length of time series
	int Nused;//length of time series after differencing, Nused = N - d
	int method;
	int optmethod;
	int p;// size of phi
	int d;// Number of times the series is to be differenced
	int q;//size of theta
	int M; // M = 0 if mean is 0.0 else M = 1
	int ncoeff;// Total Number of Coefficients to be estimated
	int cssml;// Uses CSS before MLE if 1 else uses MLE only if cssml is 0
	double *phi;
	double *theta;
	double *vcov;// Variance-Covariance Matrix Of length lvcov
	int lvcov; //length of VCOV
	double *res;
	double mean;
	double var;
	double loglik;
	double aic;
	int retval;
	double params[];
};

typedef struct sarima_set* sarima_object;

sarima_object sarima_init(int p, int d, int q,int s,int P,int D,int Q, int N);

struct sarima_set{
	int N;// length of time series
	int Nused;//length of time series after differencing, Nused = N - d - s*D
	int method;
	int optmethod;
	int p;// size of phi
	int d;// Number of times the series is to be differenced
	int q;//size of theta
	int s;// Seasonality/Period
	int P;//Size of seasonal phi
	int D;// The number of times the seasonal series is to be differenced
	int Q;//size of Seasonal Theta
	int M; // M = 0 if mean is 0.0 else M = 1
	int ncoeff;// Total Number of Coefficients to be estimated
	int cssml;// Uses CSS before MLE if 1 else uses MLE only if cssml is 0
	double *phi;
	double *theta;
	double *PHI;
	double *THETA;
	double *vcov;// Variance-Covariance Matrix Of length lvcov
	int lvcov; //length of VCOV
	double *res;
	double mean;
	double var;
	double loglik;
	double aic;
	int retval;
	double params[];
};

typedef struct ar_set* ar_object;

ar_object ar_init(int method, int N);

struct ar_set{
	int N;// length of time series
	int method;
	int optmethod;// Valid only for MLE estimation
	int p;// size of phi
	int order; // order = p
	int ordermax; // Set Maximum order to be fit
	double *phi;
	double *res;
	double mean;
	double var;
	double aic;
	int retval;
	double params[];
};

typedef struct sarimax_set* sarimax_object;

sarimax_object sarimax_init(int p, int d, int q,int P, int D, int Q,int s, int r,int imean, int N);

struct sarimax_set{
	int N;// length of time series
	int Nused;//length of time series after differencing, Nused = N - d
	int method;
	int optmethod;
	int p;// size of phi
	int d;// Number of times the series is to be differenced
	int q;//size of theta
	int s;// Seasonality/Period
	int P;//Size of seasonal phi
	int D;// The number of times the seasonal series is to be differenced
	int Q;//size of Seasonal Theta
	int r;// Number of exogenous variables
	int M; // M = 0 if mean is 0.0 else M = 1
	int ncoeff;// Total Number of Coefficients to be estimated
	double *phi;
	double *theta;
	double *PHI;
	double *THETA;
	double *exog;
	double *vcov;// Variance-Covariance Matrix Of length lvcov
	int lvcov; //length of VCOV
	double *res;
	double mean;
	double var;
	double loglik;
	double aic;
	int retval;
	int start;
	int imean;
	double params[];
};

// tseries patch (license): types and prototypes of the auto-ARIMA cluster removed
// (auto_arima_*, myarima*, aa_ret*, auto_arima1, search_arima) together with their
// implementations in ctsa.c -- a C port of the GPL-3 R package forecast. See PROVENANCE.md.


// tseries patch (license): the sarimax_wrapper* group (type sarimax_wrapper_object,
// struct sarimax_wrapper_set and the sarimax_wrapper, sarimax_wrapper_predict,
// sarimax_wrapper_summary and sarimax_wrapper_free prototypes) removed together with
// their implementations in ctsa.c -- unaudited and unreachable from tseries. See
// PROVENANCE.md, Local modifications 11.



void arima_exec(arima_object obj, double *x);

void sarimax_exec(sarimax_object obj, double *inp,double *xreg) ;

void sarima_exec(sarima_object obj, double *x);

void arima_predict(arima_object obj, double *inp, int L, double *xpred, double *amse);

void sarima_predict(sarima_object obj, double *inp, int L, double *xpred, double *amse);

void sarimax_predict(sarimax_object obj, double *inp, double *xreg, int L,double *newxreg, double *xpred, double *amse);

void ar_predict(ar_object obj, double *inp, int L, double *xpred, double *amse);

void arima_setMethod(arima_object obj, int value);

void sarima_setMethod(sarima_object obj, int value);

void sarimax_setMethod(sarimax_object obj, int value);

void arima_setCSSML(arima_object obj, int cssml);

void sarima_setCSSML(sarima_object obj, int cssml);

void arima_setOptMethod(arima_object obj, int value);

void sarima_setOptMethod(sarima_object obj, int value);

void sarimax_setOptMethod(sarimax_object obj, int value);

void arima_vcov(arima_object obj, double *vcov);

void sarima_vcov(sarima_object obj, double *vcov);

void sarimax_vcov(sarimax_object obj, double *vcov);

void arima_summary(arima_object obj);

void sarimax_summary(sarimax_object obj);

void sarima_summary(sarima_object obj);

void ar_summary(ar_object obj);

void acvf(double* vec, int N, double* par, int M);

void acvf_opt(double* vec, int N, int method, double* par, int M);

void acvf2acf(double *acf, int M);

void arima_free(arima_object object);

void sarimax_free(sarimax_object object);

void sarima_free(sarima_object object);

void ar_free(ar_object object);

// tseries patch (license): yw(), burg(), hr(), ar(), ar_exec(), ar_estimate(), model_estimate(), pacf(), pacf_opt() removed with initest.c (original licence not found). See PROVENANCE.md, Local modifications 17.


#ifdef __cplusplus
}
#endif


#endif /* CTSA_H_ */
