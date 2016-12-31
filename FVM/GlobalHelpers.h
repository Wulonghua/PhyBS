#pragma once

#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include <device_functions.h>
#include <cstdio>

__global__ void setDeviceArray(float *devicePtr, float v, int num)
{
	int i = blockDim.x * blockIdx.x + threadIdx.x;
	if (i >= num)	return;

	devicePtr[i] = v;
}

__global__ void setDeviceArray(float *devicePtr, int * mask, float v, int num)
{
	int i = blockDim.x * blockIdx.x + threadIdx.x;
	if (i >= num)	return;

	if (mask[i] > 0)
		devicePtr[i] = v;
}

// to prevent precision losing, split scale to s1 * s2
__global__ void scaleDeviceArray(float *devicePtr, float s1, float s2, int num)
{
	int i = blockDim.x * blockIdx.x + threadIdx.x;
	if (i >= num)	return;

	devicePtr[i] = devicePtr[i] * s1*s2;
}

__global__ void scaleDeviceArray(float *devicePtr, float * deviceScaledPtr, float s, int n)
{
	int i = blockDim.x * blockIdx.x + threadIdx.x;
	if (i >= n)	return;

	deviceScaledPtr[i] = devicePtr[i] * s;
}



// this is helper function that is used for NSIGHT debug.
__global__ void checkData(float *devicePtr)
{
	int i = blockDim.x * blockIdx.x + threadIdx.x;
	//printf("v: %f\n",devicePtr[7199]);
}

// helper function to print array in m*n format.
__global__ void printData(float *devicePtr, int m, int n)
{
	int i = blockDim.x * blockIdx.x + threadIdx.x;
	if (i >= n)	return;
	for (int i = 0; i < m; ++i)
	{
		for (int j = 0; j < n; ++j)
		{
			printf(" %f", devicePtr[i * n + j]);
		}
		printf("\n");
	}
}

__global__ void printData(int *devicePtr, int m, int n)
{
	int i = blockDim.x * blockIdx.x + threadIdx.x;
	if (i >= n)	return;
	for (int i = 0; i < m; ++i)
	{
		for (int j = 0; j < n; ++j)
		{
			printf(" %d", devicePtr[i * n + j]);
		}
		printf("\n");
	}
}