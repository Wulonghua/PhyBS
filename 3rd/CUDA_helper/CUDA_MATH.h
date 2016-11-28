///////////////////////////////////////////////////////////////////////////////////////////
//  Copyright (C) 2016 - 2016, Huamin Wang
//  All rights reserved.
//  Redistribution and use in source and binary forms, with or without
//  modification, are permitted provided that the following conditions
//  are met:
//     1. Redistributions of source code must retain the above copyright
//        notice, this list of conditions and the following disclaimer.
//     2. Redistributions in binary form must reproduce the above copyright
//        notice, this list of conditions and the following disclaimer in the
//        documentation and/or other materials provided with the distribution.
//     3. The names of its contributors may not be used to endorse or promote
//        products derived from this software without specific prior written
//        permission.
//  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
//  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
//  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
//  A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
//  CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
//  EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
//  PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
//  PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
//  LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
//	NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
//	SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
///////////////////////////////////////////////////////////////////////////////////////////
//  CUDA_MATH library.
///////////////////////////////////////////////////////////////////////////////////////////


#pragma once
#include <cuda_runtime.h>
#include "svd3_cuda.h"

#define	MIN(a,b)		((a)<(b)?(a):(b))
#define	MAX(a,b)		((a)>(b)?(a):(b))
#define SIGN(a)			((a)<0?-1:1)

namespace cuMath{




///////////////////////////////////////////////////////////////////////////////////////////
//  R = A - B
///////////////////////////////////////////////////////////////////////////////////////////
__host__ __device__ __forceinline__
 void Matrix_Substract_3(float *A, float *B, float *R)								
{
	for (int i = 0; i < 9; i++)
		R[i]=A[i]-B[i];
}

///////////////////////////////////////////////////////////////////////////////////////////
//  R = A * B
///////////////////////////////////////////////////////////////////////////////////////////
__host__ __device__ __forceinline__
 void Matrix_Product(float *A, float *B, float *R, int nx, int ny, int nz)			
{
	memset(R, 0, sizeof(float)*nx*nz);
	for (int i = 0; i < nx; i++)
	for (int j = 0; j < nz; j++)
	for (int k = 0; k < ny; k++)
		R[i*nz + j] += A[i*ny + k] * B[k*nz + j];
}

__host__ __device__ __forceinline__
 void Matrix_Product_3(const float *A, const float *B, float *R)					
{
	R[0]=A[0]*B[0]+A[1]*B[3]+A[2]*B[6];
	R[1]=A[0]*B[1]+A[1]*B[4]+A[2]*B[7];
	R[2]=A[0]*B[2]+A[1]*B[5]+A[2]*B[8];
	R[3]=A[3]*B[0]+A[4]*B[3]+A[5]*B[6];
	R[4]=A[3]*B[1]+A[4]*B[4]+A[5]*B[7];
	R[5]=A[3]*B[2]+A[4]*B[5]+A[5]*B[8];
	R[6]=A[6]*B[0]+A[7]*B[3]+A[8]*B[6];
	R[7]=A[6]*B[1]+A[7]*B[4]+A[8]*B[7];
	R[8]=A[6]*B[2]+A[7]*B[5]+A[8]*B[8];
}

///////////////////////////////////////////////////////////////////////////////////////////
//  R = A^T * B
///////////////////////////////////////////////////////////////////////////////////////////
__host__ __device__ __forceinline__
 void Matrix_T_Product_3(const float *A, const float *B, float *R)
{
	R[0]=A[0]*B[0]+A[3]*B[3]+A[6]*B[6];
	R[1]=A[0]*B[1]+A[3]*B[4]+A[6]*B[7];
	R[2]=A[0]*B[2]+A[3]*B[5]+A[6]*B[8];
	R[3]=A[1]*B[0]+A[4]*B[3]+A[7]*B[6];
	R[4]=A[1]*B[1]+A[4]*B[4]+A[7]*B[7];
	R[5]=A[1]*B[2]+A[4]*B[5]+A[7]*B[8];
	R[6]=A[2]*B[0]+A[5]*B[3]+A[8]*B[6];
	R[7]=A[2]*B[1]+A[5]*B[4]+A[8]*B[7];
	R[8]=A[2]*B[2]+A[5]*B[5]+A[8]*B[8];
}


///////////////////////////////////////////////////////////////////////////////////////////
//  R = A * B^T
///////////////////////////////////////////////////////////////////////////////////////////
__host__ __device__ __forceinline__
 void Matrix_Product_T_3(const float *A, const float *B, float *R)
{	
	R[0]=A[0]*B[0]+A[1]*B[1]+A[2]*B[2];
	R[1]=A[0]*B[3]+A[1]*B[4]+A[2]*B[5];
	R[2]=A[0]*B[6]+A[1]*B[7]+A[2]*B[8];
	R[3]=A[3]*B[0]+A[4]*B[1]+A[5]*B[2];
	R[4]=A[3]*B[3]+A[4]*B[4]+A[5]*B[5];
	R[5]=A[3]*B[6]+A[4]*B[7]+A[5]*B[8];
	R[6]=A[6]*B[0]+A[7]*B[1]+A[8]*B[2];
	R[7]=A[6]*B[3]+A[7]*B[4]+A[8]*B[5];
	R[8]=A[6]*B[6]+A[7]*B[7]+A[8]*B[8];
}

///////////////////////////////////////////////////////////////////////////////////////////
// r = A * x
///////////////////////////////////////////////////////////////////////////////////////////
__host__ __device__ __forceinline__
void Matrix_Vector_Product_3(float *A, float *x, float *r)
{
	r[0]=A[0]*x[0]+A[1]*x[1]+A[2]*x[2];
	r[1]=A[3]*x[0]+A[4]*x[1]+A[5]*x[2];
	r[2]=A[6]*x[0]+A[7]*x[1]+A[8]*x[2];
}

///////////////////////////////////////////////////////////////////////////////////////////
//  return (A * B')[i, j]
///////////////////////////////////////////////////////////////////////////////////////////
__host__ __device__ __forceinline__
float Matrix_Product_T(float *A, float *B, int nx, int ny, int nz, int i, int j)
{
    float ret=0;
    for(int k=0; k<ny; k++)
        ret+=A[i*ny+k]*B[j*ny+k];
    return ret;
}

__host__ __device__ __forceinline__
void Matrix_Inverse_3(float *A, float *R)
{
	R[0] = A[4] * A[8] - A[7] * A[5];
	R[1] = A[7] * A[2] - A[1] * A[8];
	R[2] = A[1] * A[5] - A[4] * A[2];
	R[3] = A[5] * A[6] - A[3] * A[8];
	R[4] = A[0] * A[8] - A[2] * A[6];
	R[5] = A[2] * A[3] - A[0] * A[5];
	R[6] = A[3] * A[7] - A[4] * A[6];
	R[7] = A[1] * A[6] - A[0] * A[7];
	R[8] = A[0] * A[4] - A[1] * A[3];
	float det = A[0] * R[0] + A[3] * R[1] + A[6] * R[2];
	float inv_det = 1 / det;
	for (int i = 0; i<9; i++)	
		R[i] *= inv_det;
}

__host__ __device__ __forceinline__
float Matrix_Determinant_3(float *x)
{
	return x[0] * (x[4] * x[8] - x[7] * x[5]) + x[3] * (x[7] * x[2] - x[1] * x[8]) + x[6] * (x[1] * x[5] - x[4] * x[2]);
}

__host__ __device__ __forceinline__
void Matrix_Transose_3(float *A, float *R)
{
	R[0] = A[0];
	R[1] = A[3];
	R[2] = A[6];
	R[3] = A[1];
	R[4] = A[4];
	R[5] = A[7];
	R[6] = A[2];
	R[7] = A[5];
	R[8] = A[8];
}

/**************More helper functions************************************************/
// added by Longhua Wu
//
/***********************************************************************************/

__host__ __device__ __forceinline__
void Matrix_Diagonal_Matrix_T_Product_3(float *U, float *S, float *V, float *F)
{
	float tmp[9];
	for (int i = 0; i < 3; ++i)
	{
		tmp[i] = U[i] * S[i];
		tmp[i + 3] = U[i + 3] * S[i];
		tmp[i + 6] = U[i + 6] * S[i];
	}

	Matrix_Product_T_3(tmp, V, F);
}

__host__ __device__ __forceinline__
void Matrix_Diagonal_Matrix_Product_3(float *U, float *S, float *V, float *F)
{
	float tmp[9];
	for (int i = 0; i < 3; ++i)
	{
		tmp[i] = U[i] * S[i];
		tmp[i + 3] = U[i + 3] * S[i];
		tmp[i + 6] = U[i + 6] * S[i];
	}

	Matrix_Product_3(tmp, V, F);
}

// solve Ax = b in which A is 2*2 matrix and x, b are 2 vetors
__host__ __device__ __forceinline__
void solveLinearSystem22(float *A, float *b, float *x)
{
	x[0] = (b[0] * A[3] - b[1] * A[1]) / (A[0] * A[3] - A[1] * A[2]);
	x[1] = (b[1] - A[2] * x[0]) / A[3];
}

__host__ __device__ __forceinline__
float pythag(float a, float b)
{
	float at = fabs(a), bt = fabs(b), ct, result;
	if (at > bt)       { ct = bt / at; result = at * sqrt(1.0 + ct * ct); }
	else if (bt > 0.0) { ct = at / bt; result = bt * sqrt(1.0 + ct * ct); }
	else result = 0.0;
	return(result);
}

__host__ __device__ __forceinline__
void SVD3(float *u, float *w, float *v)
{

	float	anorm, c, f, g, h, s, scale;
	float	x, y, z;
	float	rv1[3];
	g = scale = anorm = 0.0; //Householder reduction to bidiagonal form.

	for (int i = 0; i<3; i++)
	{
		int l = i + 1;
		rv1[i] = scale*g;
		g = s = scale = 0.0;
		if (i<3)
		{
			for (int k = i; k<3; k++) scale += fabsf(u[k * 3 + i]);
			if (scale != 0)
			{
				for (int k = i; k<3; k++)
				{
					u[k * 3 + i] /= scale;
					s += u[k * 3 + i] * u[k * 3 + i];
				}
				f = u[i * 3 + i];
				g = -sqrtf(s)*SIGN(f);
				h = f*g - s;
				u[i * 3 + i] = f - g;
				for (int j = l; j<3; j++)
				{
					s = 0;
					for (int k = i; k<3; k++)	s += u[k * 3 + i] * u[k * 3 + j];
					f = s / h;
					for (int k = i; k<3; k++)	u[k * 3 + j] += f*u[k * 3 + i];
				}
				for (int k = i; k<3; k++)		u[k * 3 + i] *= scale;
			}
		}
		w[i] = scale*g;

		g = s = scale = 0.0;
		if (i <= 2 && i != 2)
		{
			for (int k = l; k<3; k++)	scale += fabsf(u[i * 3 + k]);
			if (scale != 0)
			{
				for (int k = l; k<3; k++)
				{
					u[i * 3 + k] /= scale;
					s += u[i * 3 + k] * u[i * 3 + k];
				}
				f = u[i * 3 + l];
				g = -sqrtf(s)*SIGN(f);
				h = f*g - s;
				u[i * 3 + l] = f - g;
				for (int k = l; k<3; k++) rv1[k] = u[i * 3 + k] / h;
				for (int j = l; j<3; j++)
				{
					s = 0;
					for (int k = l; k<3; k++)	s += u[j * 3 + k] * u[i * 3 + k];
					for (int k = l; k<3; k++)	u[j * 3 + k] += s*rv1[k];
				}
				for (int k = l; k<3; k++) u[i * 3 + k] *= scale;
			}
		}
		anorm = MAX(anorm, (fabs(w[i]) + fabs(rv1[i])));
	}

	for (int i = 2, l; i >= 0; i--) //Accumulation of right-hand transformations.
	{
		if (i<2)
		{
			if (g != 0)
			{
				for (int j = l; j<3; j++) //Double division to avoid possible under
					v[j * 3 + i] = (u[i * 3 + j] / u[i * 3 + l]) / g;
				for (int j = l; j<3; j++)
				{
					s = 0;
					for (int k = l; k<3; k++)	s += u[i * 3 + k] * v[k * 3 + j];
					for (int k = l; k<3; k++)	v[k * 3 + j] += s*v[k * 3 + i];
				}
			}
			for (int j = l; j<3; j++)	v[i * 3 + j] = v[j * 3 + i] = 0.0;
		}
		v[i * 3 + i] = 1.0;
		g = rv1[i];
		l = i;
	}

	for (int i = 2; i >= 0; i--) //Accumulation of left-hand transformations.
	{
		int l = i + 1;
		g = w[i];
		for (int j = l; j<3; j++) u[i * 3 + j] = 0;
		if (g != 0)
		{
			g = 1 / g;
			for (int j = l; j<3; j++)
			{
				s = 0;
				for (int k = l; k<3; k++)	s += u[k * 3 + i] * u[k * 3 + j];
				f = (s / u[i * 3 + i])*g;
				for (int k = i; k<3; k++)	u[k * 3 + j] += f*u[k * 3 + i];
			}
			for (int j = i; j<3; j++)		u[j * 3 + i] *= g;
		}
		else for (int j = i; j<3; j++)		u[j * 3 + i] = 0.0;
		u[i * 3 + i]++;
	}

	for (int k = 2; k >= 0; k--)				//Diagonalization of the bidiagonal form: Loop over
	{
		for (int its = 0; its<30; its++)	//singular values, and over allowed iterations.
		{
			bool flag = true;
			int  l;
			int	 nm;
			for (l = k; l >= 0; l--)			//Test for splitting.
			{
				nm = l - 1;
				if ((float)(fabs(rv1[l]) + anorm) == anorm)
				{
					flag = false;
					break;
				}
				if ((float)(fabs(w[nm]) + anorm) == anorm)	break;
			}
			if (flag)
			{
				c = 0.0; //Cancellation of rv1[l], if l > 0.
				s = 1.0;
				for (int i = l; i<k + 1; i++)
				{
					f = s*rv1[i];
					rv1[i] = c*rv1[i];
					if ((float)(fabs(f) + anorm) == anorm) break;
					g = w[i];
					h = pythag(f, g);
					w[i] = h;
					h = 1.0 / h;
					c = g*h;
					s = -f*h;
					for (int j = 0; j<3; j++)
					{
						y = u[j * 3 + nm];
						z = u[j * 3 + i];
						u[j * 3 + nm] = y*c + z*s;
						u[j * 3 + i] = z*c - y*s;
					}
				}
			}
			z = w[k];
			if (l == k)		// Convergence.
			{
				if (z<0.0)	// Singular value is made nonnegative.
				{
					w[k] = -z;
					for (int j = 0; j<3; j++) v[j * 3 + k] = -v[j * 3 + k];
				}
				break;
			}
			if (its == 29) { printf("Error: no convergence in 30 svdcmp iterations"); getchar(); }
			x = w[l]; //Shift from bottom 2-by-2 minor.
			nm = k - 1;
			y = w[nm];
			g = rv1[nm];
			h = rv1[k];
			f = ((y - z)*(y + z) + (g - h)*(g + h)) / (2.0*h*y);
			g = pythag(f, (float)1.0);
			f = ((x - z)*(x + z) + h*((y / (f + fabs(g)*SIGN(f))) - h)) / x;
			c = s = 1.0; //Next QR transformation:

			for (int j = l; j <= nm; j++)
			{
				int i = j + 1;
				g = rv1[i];
				y = w[i];
				h = s*g;
				g = c*g;
				z = pythag(f, h);
				rv1[j] = z;
				c = f / z;
				s = h / z;
				f = x*c + g*s;
				g = g*c - x*s;
				h = y*s;
				y *= c;
				for (int jj = 0; jj<3; jj++)
				{
					x = v[jj * 3 + j];
					z = v[jj * 3 + i];
					v[jj * 3 + j] = x*c + z*s;
					v[jj * 3 + i] = z*c - x*s;
				}
				z = pythag(f, h);
				w[j] = z; //Rotation can be arbitrary if z D 0.
				if (z)
				{
					z = 1.0 / z;
					c = f*z;
					s = h*z;
				}
				f = c*g + s*y;
				x = c*y - s*g;
				for (int jj = 0; jj<3; jj++)
				{
					y = u[jj * 3 + j];
					z = u[jj * 3 + i];
					u[jj * 3 + j] = y*c + z*s;
					u[jj * 3 + i] = z*c - y*s;
				}
			}
			rv1[l] = 0.0;
			rv1[k] = f;
			w[k] = x;
		}
	}
}

__host__ __device__ __forceinline__
void Matrix_SVD_3(float *A, float *U, float *S, float *V)
{
	float T[6]; // parameters holder
	svd(A[0], A[1], A[2], A[3], A[4], A[5], A[6], A[7], A[8],
		U[0], U[1], U[2], U[3], U[4], U[5], U[6], U[7], U[8],
		S[0], T[0], T[1], T[2], S[1], T[3], T[4], T[5], S[2],
		V[0], V[1], V[2], V[3], V[4], V[5], V[6], V[7], V[8]);

	//memcpy(U, A, sizeof(float) * 9);
	//cuMath::SVD3(U, S, V);
}


}