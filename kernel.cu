#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <malloc.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

#define E 0.001
#define N  (512*512)
#define BLOCKS 256

__global__ void Kernel(double* present, double* past, double* grid, double* f, int SIZE)
{
	double h = 1.0 / (SIZE - 1);
	int id = blockIdx.x * blockDim.x + threadIdx.x;

	if (id >= SIZE && id <= (SIZE * SIZE - SIZE - 1) && id % SIZE != 0 && (id + 1) % SIZE != 0)
		present[id] = 0.25 * (past[id - SIZE] + past[id - 1] + past[id + SIZE] + past[id + 1] - pow(h, 2) * f[id]);
	else
		present[id] = past[id];
}

__global__ void Errors(double* present, double* grid, double* error)
{
	int id = blockIdx.x * blockDim.x + threadIdx.x;
	error[id] = fabs(present[id] - grid[id]);
}

double fr(double a, double b) {
	return 4 + 2 * pow(a, 2) - 2 * a + 2 * pow(b, 2) - 2 * b;
}

double u(double a, double b) {
	double c = 0;
	c = (pow(a, 2) - a + 1) * (pow(b, 2) - b + 1);
	return c;
}

void gran(double* arr, double* Z, int X0, int Xn, int Y0, int Yn, int SIZE)
{
	int k = 0;
	for (int i = X0; i < Xn; ++i)
	{
		for (int j = Y0; j < Yn; ++j)
		{
			k = i - X0 + j - Y0;
			arr[i * SIZE + j] = Z[k] * Z[k] - Z[k] + 1;
		}
	}

}

bool Writetofile(char* filename, double* X, double* Y, double* Z, int SIZE)
{

	FILE* fp = fopen(filename, "w");
	if (fp == NULL)
	{
		printf("Can not open file %s\n", filename);
		return false;
	}
	for (int i = 0; i < SIZE; i++)
		for (int j = 0; j < SIZE; j++)
			fprintf(fp, "%7.8f %7.8f %7.8f\n", X[i], Y[j], Z[i * SIZE + j]);
	fclose(fp);
	return true;
}

int main()
{
	int SIZE = 512;
	int steperror = 1000;
	double* present = (double*)malloc((SIZE * SIZE) * sizeof(double));
	double* grid = (double*)malloc((SIZE * SIZE) * sizeof(double));
	double* past = (double*)malloc((SIZE * SIZE) * sizeof(double));


	double* X = (double*)malloc(SIZE * sizeof(double));
	double* Y = (double*)malloc(SIZE * sizeof(double));
	double* f = (double*)malloc(SIZE * SIZE * sizeof(double));
	double* errors = (double*)malloc(SIZE * SIZE * sizeof(double));

	for (int i = 0; i < SIZE; i++) {
		X[i] = double(i) / (SIZE - 1);
		Y[i] = double(i) / (SIZE - 1);
	}

	for (int i = 0; i < SIZE; i++) {
		for (int j = 0; j < SIZE; ++j) {
			past[i * SIZE + j] = 0;
			grid[i * SIZE + j] = u(X[i], Y[j]);
			f[i * SIZE + j] = fr(X[i], Y[j]);
		}
	}

	clock_t start, end;
	double time;

	start = clock();

	gran(past, X, 0, SIZE, 0, 1, SIZE);
	gran(past, X, 0, SIZE, SIZE - 1, SIZE, SIZE);
	gran(past, Y, 0, 1, 0, SIZE, SIZE);
	gran(past, Y, SIZE - 1, SIZE, 0, SIZE, SIZE);

	double localMax = 1;
	int iter = 0;

	double* dev_arrpresent = 0;
	double* dev_arrpast = 0;
	double* dev_arrgrid = 0;
	double* dev_f = 0;
	double* dev_error = 0;

	cudaSetDevice(0);
	cudaMalloc((void**)&dev_arrpresent, N * sizeof(double));
	cudaMalloc((void**)&dev_arrpast, N * sizeof(double));
	cudaMalloc((void**)&dev_arrgrid, N * sizeof(double));
	cudaMalloc((void**)&dev_f, N * sizeof(double));
	cudaMalloc((void**)&dev_error, N * sizeof(double));

	cudaMemcpy(dev_arrpast, past, N * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(dev_arrgrid, grid, N * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(dev_f, f, N * sizeof(double), cudaMemcpyHostToDevice);

	while (localMax > E)
	{
		Kernel << < BLOCKS, N / BLOCKS >> > (dev_arrpresent, dev_arrpast, dev_arrgrid, dev_f, SIZE);
		if (iter % steperror == 0) {
			localMax = 0;
			Errors << < BLOCKS, N / BLOCKS >> > (dev_arrpresent, dev_arrgrid, dev_error);
			cudaMemcpy(errors, dev_error, N * sizeof(double), cudaMemcpyDeviceToHost);
			for (int i = 0; i < SIZE * SIZE; i++)
				if (localMax < errors[i])
					localMax = errors[i];
			printf("%d, %lf\n", iter, localMax);
		}
		double* temp = dev_arrpast;
		dev_arrpast = dev_arrpresent;
		dev_arrpresent = temp;
		iter++;
	}
	cudaMemcpy(past, dev_arrpresent, N * sizeof(double), cudaMemcpyDeviceToHost);
	Writetofile("Jakobi.dat", X, Y, past, SIZE);

	end = clock();
	time = ((double)(end - start)) / CLOCKS_PER_SEC;

	printf("Number of iterations:\n");
	printf("%d\n", iter - 1);
	printf("Time, s:\n");
	printf("%4.2lf\n", time);

	cudaFree(dev_arrpresent);
	cudaFree(dev_arrpast);
	cudaFree(dev_arrgrid);
	cudaFree(dev_f);
	cudaFree(dev_error);

	free(X);
	free(Y);
	free(present);
	free(grid);
	free(past);

	return 0;
}
