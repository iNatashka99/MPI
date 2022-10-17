#include "Rects.h"
#include "Somefunc.h"
#include <math.h>
#include <mpi.h>


// �������������� ��������������
int in_rect(int x, int y, int x1, int x2, int y1, int y2)
{
	if ((x < x2) & (x > x1) & (y < y2) & (y > y1))
		return 1;
	if ((x == x1) & (y <= y2) & (y >= y1))
		return 2;
	if ((x == x2) & (y <= y2) & (y >= y1))
		return 2;
	if ((y == y1) & (x <= x2) & (x >= x1))
		return 2;
	if ((y == y2) & (x <= x2) & (x >= x1))
		return 2;
	return 0;
}

// �������� ��������������
void send_rect(double **data, int x1, int x2, int y1, int y2, int proc, int shape)
{
	for (int i = x1; i <= x2; i++)
	{
		for (int j = y1; j <= y2; j++)
		{
			MPI_Send(&data[j][i], 1, MPI_DOUBLE, proc, shape, MPI_COMM_WORLD);
		}
	}
}

// ����� ��������������
void recive_rect(double **data, int x1, int x2, int y1, int y2, int proc, int shape)
{
	MPI_Status stat;
	double el;
	for (int i = x1; i <= x2; i++)
	{
		for (int j = y1; j <= y2; j++)
		{
			MPI_Recv(&el, 1, MPI_DOUBLE, proc, shape, MPI_COMM_WORLD, &stat);
			data[j][i] = el;
		}
	}
}

// ������ ��� �������� ����������
void set_rect(double **data, double **layer, int x1, int x2, int y1, int y2)
{
	for (int i = x1; i <= x2; i++)
	{
		for (int j = y1; j <= y2; j++)
		{
			layer[j][i] = data[j][i];
		}
	}
}

// ����� � ��������������
double norma_rect(double **layer1, double **layer2, int x1, int x2, int y1, int y2)
{
	double nrm = 0;
	for (int i = x1; i <= x2; i++)
	{
		for (int j = y1; j <= y2; j++)
		{
			nrm = nrm + fabs(layer1[j][i] - layer2[j][i]);
		}
	}
	return nrm;
}