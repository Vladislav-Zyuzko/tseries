// SPDX-License-Identifier: BSD-3-Clause
/*
 * pred.h
 *
 *  Created on: Jul 14, 2014
 *      Author: Rafat Hussaint
 */

#ifndef PRED_H_
#define PRED_H_

#include "boxjenkins.h" /* tseries patch (license): was "initest.h", deleted with initest.c (PROVENANCE.md, Local modifications 17); this is what it included */

#ifdef __cplusplus
extern "C" {
#endif

void predictarima(double *zt,int lenzt,int p,int d, int q,double *phi, double *theta,
		double constant,int forlength,double *oup);

#ifdef __cplusplus
}
#endif

#endif /* PRED_H_ */
