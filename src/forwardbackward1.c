#include <R.h>

void calAlpha (double *alpha, double *c0, double *A, double *f0x, double *f1x, int *NUM)
{
	int k;
	double c0denom;
	
	for(k = 0; k < (*NUM - 1); k++){
		alpha[k+1] = (alpha[k]*A[k*4] + alpha[*NUM+k]*A[k*4+1])*f0x[k+1];
		alpha[*NUM+k+1] = (alpha[k]*A[k*4+2] + alpha[*NUM+k]*A[k*4+3])*f1x[k+1];
	/* rescaling alpha */
		c0denom = alpha[k+1]+alpha[*NUM+k+1];
		c0[k+1] = 1/c0denom;
		alpha[k+1] = c0[k+1]*alpha[k+1];
		alpha[*NUM+k+1] = c0[k+1]*alpha[*NUM+k+1];
	}

}

void calBeta (double *beta, double *c0, double *A, double *f0x, double *f1x, int *NUM)
{
	int k;

	for(k = (*NUM - 2); k >= 0; k--){
		beta[k] = A[k*4]*f0x[k+1]*beta[k+1] + A[k*4+2]*f1x[k+1]*beta[k+1+*NUM];
		beta[k+*NUM] = A[k*4+1]*f0x[k+1]*beta[k+1] + A[k*4+3]*f1x[k+1]*beta[k+1+*NUM];
		beta[k] = c0[k]*beta[k];
		beta[k+*NUM] = c0[k]*beta[k+*NUM];
	}
}


void calLfdr (double *alpha, double *beta, double *lfdr, int *NUM)
{
	int k;
	double q1, q2;

	for(k = 0; k < *NUM; k++){
		q1 = alpha[k]*beta[k];
		q2 = alpha[k+*NUM]*beta[k+*NUM];
		lfdr[k] = q1/(q1+q2);
	}
}


void calGamma (double *alpha, double *beta, double *A, double *f0x, double *f1x, double *gamma, double *dgamma, int *NUM )
{
	int i, j, k;
	double denom, fx;
	
	for(k = 0; k < (*NUM - 1); k++){
		denom = 0;	
		for(i = 0; i < 2; i++){
			for(j = 0; j < 2; j++){
				fx = (1-j)*f0x[k+1]+j*f1x[k+1];
				denom = denom + alpha[k+*NUM*i]*A[k*4+i+j*2]*fx*beta[k+1+*NUM*j];
			}
		}

		for(i = 0; i < 2; i++){
			gamma[k+*NUM*i] = 0;
			for(j = 0; j < 2; j++){
				fx = (1-j)*f0x[k+1]+j*f1x[k+1];
				dgamma[k*4+i+j*2] = alpha[k+*NUM*i]*A[k*4+i+j*2]*fx*beta[k+1+*NUM*j]/denom;
				gamma[k+*NUM*i] = gamma[k+*NUM*i] + dgamma[k*4+i+2*j];
			}	
		}	
	}
}

