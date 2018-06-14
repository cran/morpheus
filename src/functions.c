#include <stdlib.h>

//index matrix (by columns)
int mi(int i, int j, int d1, int d2)
{
	return j*d1+i;
}

//index 3-tensor (by columns, matrices ordered by last dim)
int ti(int i, int j, int k, int d1, int d2, int d3)
{
	return k*d1*d2 + j*d1 + i;
}

// Emprical cross-moment of order 2 between X size nxd and Y size n
void Moments_M2(double* X, double* Y, int* pn, int* pd, double* M2)
{
	int n=*pn, d=*pd;
	//double* M2 = (double*)calloc(d*d,sizeof(double));

	// M2 = E[Y*X^*2] - E[Y*e^*2] = E[Y (X^*2 - I)]
	for (int j=0; j<d; j++)
	{
		for (int i=0; i<n; i++)
		{
			M2[mi(j,j,d,d)] -= Y[i] / n;
			for (int k=0; k<d; k++)
				M2[mi(j,k,d,d)] += Y[i] * X[mi(i,j,n,d)]*X[mi(i,k,n,d)] / n;
		}
	}
}

// Emprical cross-moment of order 3 between X size nxd and Y size n
void Moments_M3(double* X, double* Y, int* pn, int* pd, double* M3)
{
	int n=*pn, d=*pd;
	//double* M3 = (double*)calloc(d*d*d,sizeof(double));

	// M3 = E[Y*X^*3] - E[Y*e*X*e] - E[Y*e*e*X] - E[Y*X*e*e]
	for (int j=0; j<d; j++)
	{
		for (int k=0; k<d; k++)
		{
			for (int i=0; i<n; i++)
			{
				double tensor_elt = Y[i]*X[mi(i,k,n,d)] / n;
				M3[ti(j,k,j,d,d,d)] -= tensor_elt;
				M3[ti(j,j,k,d,d,d)] -= tensor_elt;
				M3[ti(k,j,j,d,d,d)] -= tensor_elt;
				for (int o=0; o<d; o++)
					M3[ti(j,k,o,d,d,d)] += Y[i] * X[mi(i,j,n,d)]*X[mi(i,k,n,d)]*X[mi(i,o,n,d)] / n;
			}
		}
	}
}
