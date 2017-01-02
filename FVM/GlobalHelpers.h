#pragma once

#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include <device_functions.h>
#include <cstdio>

struct CSRmatrix
{
	int nnz;   // number of non-zero elements
	int m;     // rows' number
	float *d_Valptr;
	int *d_Rowptr;
	int *d_Colptr;
	//   helper indices
	int *d_diagonalIdx;
};

__global__ void setDeviceArray(float *devicePtr, float v, int num);
__global__ void setDeviceArray(float *devicePtr, int * mask, float v, int num);

// to prevent precision losing, split scale to s1 * s2
__global__ void scaleDeviceArray(float *devicePtr, float s1, float s2, int num);
__global__ void scaleDeviceArray(float *devicePtr, float * deviceScaledPtr, float s, int n);

// this is helper function that is used for NSIGHT debug.
__global__ void checkData(float *devicePtr);

// helper function to print array in m*n format.
__global__ void printData(float *devicePtr, int m, int n);
__global__ void printData(int *devicePtr, int m, int n);