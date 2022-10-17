#include "MPI_funcsr.h"
#include "SSOR.h"
#include "ForMatr.h"
#include "Somefunc.h"
#include <mpi.h>
#include <iostream>
#include <iomanip>
using namespace std;

void SendShapes(double **matr, double **local_layer, int size, double R, double r,
	int circ1_x0, int circ1_y0,
	int circ2_x0, int circ2_y0,
	int circ3_x0, int circ3_y0,
	int rect_circ1_x1, int rect_circ1_x2, int rect_circ1_y1, int rect_circ1_y2,
	int rect_circ2_x1, int rect_circ2_x2, int rect_circ2_y1, int rect_circ2_y2,
	int rect_circ3_x1, int rect_circ3_x2, int rect_circ3_y1, int rect_circ3_y2,
	int rect1_x1, int rect1_x2, int rect1_y1, int rect1_y2,
	int rect2_x1, int rect2_x2, int rect2_y1, int rect2_y2,
	int numprocs, int **proc_shapes)
{
	for (int i = 1; i < numprocs; i++)
	{
		for (int j = 0; j < 8 / numprocs; j++)
		{
			if (proc_shapes[i][j] == 0)
			{
				send_circle(matr, circ1_x0, circ1_y0, R, r, i, 0);
			}
			if (proc_shapes[i][j] == 1)
			{
				send_circle(matr, circ2_x0, circ2_y0, R, r, i, 1);
			}
			if (proc_shapes[i][j] == 2)
			{
				send_circle(matr, circ3_x0, circ3_y0, R, r, i, 2);
			}
			if (proc_shapes[i][j] == 3)
			{
				send_rect(matr, rect_circ1_x1, rect_circ1_x2, rect_circ1_y1, rect_circ1_y2, i, 3);
			}
			if (proc_shapes[i][j] == 4)
			{
				send_rect(matr, rect_circ2_x1, rect_circ2_x2, rect_circ2_y1, rect_circ2_y2, i, 4);
			}
			if (proc_shapes[i][j] == 5)
			{
				send_rect(matr, rect_circ3_x1, rect_circ3_x2, rect_circ3_y1, rect_circ3_y2, i, 5);
			}
			if (proc_shapes[i][j] == 6)
			{
				send_rect(matr, rect1_x1, rect1_x2, rect1_y1, rect1_y2, i, 6);
			}
			if (proc_shapes[i][j] == 7)
			{
				send_rect(matr, rect2_x1, rect2_x2, rect2_y1, rect2_y2, i, 7);
			}
		}
	}
	for (int j = 0; j < 8 / numprocs; j++)
	{
		if (proc_shapes[0][j] == 0)
		{
			set_circle(matr, local_layer, circ1_x0, circ1_y0, R, r);
		}
		if (proc_shapes[0][j] == 1)
		{
			set_circle(matr, local_layer, circ2_x0, circ2_y0, R, r);
		}
		if (proc_shapes[0][j] == 2)
		{
			set_circle(matr, local_layer, circ3_x0, circ3_y0, R, r);
		}
		if (proc_shapes[0][j] == 3)
		{
			set_rect(matr, local_layer, rect_circ1_x1, rect_circ1_x2, rect_circ1_y1, rect_circ1_y2);
		}
		if (proc_shapes[0][j] == 4)
		{
			set_rect(matr, local_layer, rect_circ2_x1, rect_circ2_x2, rect_circ2_y1, rect_circ2_y2);
		}
		if (proc_shapes[0][j] == 5)
		{
			set_rect(matr, local_layer, rect_circ3_x1, rect_circ3_x2, rect_circ3_y1, rect_circ3_y2);
		}
		if (proc_shapes[0][j] == 6)
		{
			set_rect(matr, local_layer, rect1_x1, rect1_x2, rect1_y1, rect1_y2);
		}
		if (proc_shapes[0][j] == 7)
		{
			set_rect(matr, local_layer, rect2_x1, rect2_x2, rect2_y1, rect2_y2);
		}
	}
}

void SendBorders(double **matr, double **local_layer, int size, double R, double r,
	int circ1_x0, int circ1_y0,
	int circ2_x0, int circ2_y0,
	int circ3_x0, int circ3_y0,
	int rect_circ1_x1, int rect_circ1_x2, int rect_circ1_y1, int rect_circ1_y2,
	int rect_circ2_x1, int rect_circ2_x2, int rect_circ2_y1, int rect_circ2_y2,
	int rect_circ3_x1, int rect_circ3_x2, int rect_circ3_y1, int rect_circ3_y2,
	int rect1_x1, int rect1_x2, int rect1_y1, int rect1_y2,
	int rect2_x1, int rect2_x2, int rect2_y1, int rect2_y2,
	int numprocs, int **proc_shapes,
	int rank)
{
	int proc;
	if (numprocs > 1)
	{
		int i = rank;
		int x1, x2, y1, y2;
		for (int j = 0; j < 8 / numprocs; j++)
		{
			if (proc_shapes[i][j] == 0)
			{
				x1 = circ1_x0 - R;
				x2 = circ1_x0 + R;
				y1 = circ1_y0 - R;
				y2 = circ1_y0 + R;
				proc = in_proc(3, numprocs, proc_shapes);
				if (rank != proc)
				{
					for (int i = x1; i <= x2; i++)
					{
						for (int j = y1; j <= y2; j++)
						{
							if ((in_circle_dec(i, j, R, r, circ1_x0, circ1_y0) == 1) &
								(in_rect(i, j, rect_circ1_x1, rect_circ1_x2, rect_circ1_y1, rect_circ1_y2) == 2))
							{
								MPI_Send(&local_layer[j][i], 1, MPI_DOUBLE, proc, 0, MPI_COMM_WORLD);
							}
						}
					}
				}
				proc = in_proc(6, numprocs, proc_shapes);
				if (rank != proc)
				{
					for (int i = x1; i <= x2; i++)
					{
						for (int j = y1; j <= y2; j++)
						{
							if ((in_circle_dec(i, j, R, r, circ1_x0, circ1_y0) == 1) &
								(in_rect(i, j, rect1_x1, rect1_x2, rect1_y1, rect1_y2) == 2))
							{
								MPI_Send(&local_layer[j][i], 1, MPI_DOUBLE, proc, 0, MPI_COMM_WORLD);
							}
						}
					}
				}
			}
			if (proc_shapes[i][j] == 1)
			{
				x1 = circ2_x0 - R;
				x2 = circ2_x0 + R;
				y1 = circ2_y0 - R;
				y2 = circ2_y0 + R;
				proc = in_proc(4, numprocs, proc_shapes);
				if (rank != proc)
				{
					for (int i = x1; i <= x2; i++)
					{
						for (int j = y1; j <= y2; j++)
						{
							if ((in_circle_dec(i, j, R,r, circ2_x0, circ2_y0) == 1) &
								(in_rect(i, j, rect_circ2_x1, rect_circ2_x2, rect_circ2_y1, rect_circ2_y2) == 2))
							{
								MPI_Send(&local_layer[j][i], 1, MPI_DOUBLE, proc, 0, MPI_COMM_WORLD);
							}
						}
					}
				}
				proc = in_proc(6, numprocs, proc_shapes);
				if (rank != proc)
				{
					for (int i = x1; i <= x2; i++)
					{
						for (int j = y1; j <= y2; j++)
						{
							if ((in_circle_dec(i, j, R,r, circ2_x0, circ2_y0) == 1) &
								(in_rect(i, j, rect1_x1, rect1_x2, rect1_y1, rect1_y2) == 2))
							{
								MPI_Send(&local_layer[j][i], 1, MPI_DOUBLE, proc, 0, MPI_COMM_WORLD);
							}
						}
					}
				}
				proc = in_proc(7, numprocs, proc_shapes);
				if (rank != proc)
				{
					for (int i = x1; i <= x2; i++)
					{
						for (int j = y1; j <= y2; j++)
						{
							if ((in_circle_dec(i, j, R, r, circ2_x0, circ2_y0) == 1) &
								(in_rect(i, j, rect2_x1, rect2_x2, rect2_y1, rect2_y2) == 2))
							{
								MPI_Send(&local_layer[j][i], 1, MPI_DOUBLE, proc, 0, MPI_COMM_WORLD);
							}
						}
					}
				}
			}
			if (proc_shapes[i][j] == 2)
			{
				x1 = circ3_x0 - R;
				x2 = circ3_x0 + R;
				y1 = circ3_y0 - R;
				y2 = circ3_y0 + R;
				proc = in_proc(5, numprocs, proc_shapes);
				if (rank != proc)
				{
					for (int i = x1; i <= x2; i++)
					{
						for (int j = y1; j <= y2; j++)
						{
							if ((in_circle_dec(i, j, R, r, circ3_x0, circ3_y0) == 1) &
								(in_rect(i, j, rect_circ3_x1, rect_circ3_x2, rect_circ3_y1, rect_circ3_y2) == 2))
							{
								MPI_Send(&local_layer[j][i], 1, MPI_DOUBLE, proc, 0, MPI_COMM_WORLD);
							}
						}
					}
				}
				proc = in_proc(7, numprocs, proc_shapes);
				if (rank != proc)
				{
					for (int i = x1; i <= x2; i++)
					{
						for (int j = y1; j <= y2; j++)
						{
							if ((in_circle_dec(i, j, R, r, circ3_x0, circ3_y0) == 1) &
								(in_rect(i, j, rect2_x1, rect2_x2, rect2_y1, rect2_y2) == 2))
							{
								MPI_Send(&local_layer[j][i], 1, MPI_DOUBLE, proc, 0, MPI_COMM_WORLD);
							}
						}
					}
				}
			}
			if (proc_shapes[i][j] == 3)
			{
				proc = in_proc(0, numprocs, proc_shapes);
				if (rank != proc)
				{
					for (int i = rect_circ1_x1; i <= rect_circ1_x2; i++)
					{
						for (int j = rect_circ1_y1; j <= rect_circ1_y2; j++)
						{
							if ((in_circle_dec(i, j, R, r, circ1_x0, circ1_y0) == 2) &
								(in_rect(i, j, rect_circ1_x1, rect_circ1_x2, rect_circ1_y1, rect_circ1_y2) == 1))
							{
								MPI_Send(&local_layer[j][i], 1, MPI_DOUBLE, proc, 0, MPI_COMM_WORLD);
							}
						}
					}
				}
			}
			if (proc_shapes[i][j] == 4)
			{
				proc = in_proc(1, numprocs, proc_shapes);
				if (rank != proc)
				{
					for (int i = rect_circ2_x1; i <= rect_circ2_x2; i++)
					{
						for (int j = rect_circ2_y1; j <= rect_circ2_y2; j++)
						{
							if ((in_circle_dec(i, j, R, r, circ2_x0, circ2_y0) == 2) &
								(in_rect(i, j, rect_circ2_x1, rect_circ2_x2, rect_circ2_y1, rect_circ2_y2) == 1))
							{
								MPI_Send(&local_layer[j][i], 1, MPI_DOUBLE, proc, 0, MPI_COMM_WORLD);
							}
						}
					}
				}
			}
			if (proc_shapes[i][j] == 5)
			{
				proc = in_proc(2, numprocs, proc_shapes);
				if (rank != proc)
				{
					for (int i = rect_circ3_x1; i <= rect_circ3_x2; i++)
					{
						for (int j = rect_circ3_y1; j <= rect_circ3_y2; j++)
						{
							if ((in_circle_dec(i, j, R, r, circ3_x0, circ3_y0) == 2) &
								(in_rect(i, j, rect_circ3_x1, rect_circ3_x2, rect_circ3_y1, rect_circ3_y2) == 1))
							{
								MPI_Send(&local_layer[j][i], 1, MPI_DOUBLE, proc, 0, MPI_COMM_WORLD);
							}
						}
					}
				}
			}
			if (proc_shapes[i][j] == 6)
			{
				proc = in_proc(0, numprocs, proc_shapes);
				if (rank != proc)
				{
					for (int i = rect1_x1; i <= rect1_x2; i++)
					{
						for (int j = rect1_y1; j <= rect1_y2; j++)
						{
							if ((in_circle_dec(i, j, R, r, circ1_x0, circ1_y0) == 2) &
								(in_rect(i, j, rect1_x1, rect1_x2, rect1_y1, rect1_y2) == 1))
							{
								MPI_Send(&local_layer[j][i], 1, MPI_DOUBLE, proc, 0, MPI_COMM_WORLD);
							}
						}
					}
				}
				proc = in_proc(1, numprocs, proc_shapes);
				if (rank != proc)
				{
					for (int i = rect1_x1; i <= rect1_x2; i++)
					{
						for (int j = rect1_y1; j <= rect1_y2; j++)
						{
							if ((in_circle_dec(i, j, R, r, circ2_x0, circ2_y0) == 2) &
								(in_rect(i, j, rect1_x1, rect1_x2, rect1_y1, rect1_y2) == 1))
							{
								MPI_Send(&local_layer[j][i], 1, MPI_DOUBLE, proc, 0, MPI_COMM_WORLD);
							}
						}
					}
				}
			}
			if (proc_shapes[i][j] == 7)
			{
				proc = in_proc(1, numprocs, proc_shapes);
				if (rank != proc)
				{
					for (int i = rect2_x1; i <= rect2_x2; i++)
					{
						for (int j = rect2_y1; j <= rect2_y2; j++)
						{
							if ((in_circle_dec(i, j, R, r, circ2_x0, circ2_y0) == 2) &
								(in_rect(i, j, rect2_x1, rect2_x2, rect2_y1, rect2_y2) == 1))
							{
								MPI_Send(&local_layer[j][i], 1, MPI_DOUBLE, proc, 0, MPI_COMM_WORLD);
							}
						}
					}
				}
				proc = in_proc(2, numprocs, proc_shapes);
				if (rank != proc)
				{
					for (int i = rect2_x1; i <= rect2_x2; i++)
					{
						for (int j = rect2_y1; j <= rect2_y2; j++)
						{
							if ((in_circle_dec(i, j, R, r, circ3_x0, circ3_y0) == 2) &
								(in_rect(i, j, rect2_x1, rect2_x2, rect2_y1, rect2_y2) == 1))
							{
								MPI_Send(&local_layer[j][i], 1, MPI_DOUBLE, proc, 0, MPI_COMM_WORLD);
							}
						}
					}
				}
			}
		}
	}
}

void ReciveShapes(double **local_layer, int size, double R, double r,
	int circ1_x0, int circ1_y0,
	int circ2_x0, int circ2_y0,
	int circ3_x0, int circ3_y0,
	int rect_circ1_x1, int rect_circ1_x2, int rect_circ1_y1, int rect_circ1_y2,
	int rect_circ2_x1, int rect_circ2_x2, int rect_circ2_y1, int rect_circ2_y2,
	int rect_circ3_x1, int rect_circ3_x2, int rect_circ3_y1, int rect_circ3_y2,
	int rect1_x1, int rect1_x2, int rect1_y1, int rect1_y2,
	int rect2_x1, int rect2_x2, int rect2_y1, int rect2_y2,
	int numprocs, int **proc_shapes,
	int rank)
{
	int i = rank;
	for (int j = 0; j < 8 / numprocs; j++)
	{
		if (proc_shapes[i][j] == 0)
		{
			recive_circle(local_layer, circ1_x0, circ1_y0, R, r, 0, 0);
		}
		if (proc_shapes[i][j] == 1)
		{
			recive_circle(local_layer, circ2_x0, circ2_y0, R, r, 0, 1);
		}
		if (proc_shapes[i][j] == 2)
		{
			recive_circle(local_layer, circ3_x0, circ3_y0, R, r, 0, 2);
		}
		if (proc_shapes[i][j] == 3)
		{
			recive_rect(local_layer, rect_circ1_x1, rect_circ1_x2, rect_circ1_y1, rect_circ1_y2, 0, 3);
		}
		if (proc_shapes[i][j] == 4)
		{
			recive_rect(local_layer, rect_circ2_x1, rect_circ2_x2, rect_circ2_y1, rect_circ2_y2, 0, 4);
		}
		if (proc_shapes[i][j] == 5)
		{
			recive_rect(local_layer, rect_circ3_x1, rect_circ3_x2, rect_circ3_y1, rect_circ3_y2, 0, 5);
		}
		if (proc_shapes[i][j] == 6)
		{
			recive_rect(local_layer, rect1_x1, rect1_x2, rect1_y1, rect1_y2, 0, 6);
		}
		if (proc_shapes[i][j] == 7)
		{
			recive_rect(local_layer, rect2_x1, rect2_x2, rect2_y1, rect2_y2, 0, 7);
		}
	}
}


// SSOR 1 без внутр распараллеливания
double ** SSOR_1(double **local_layer, int size, double R, double r, double h, double w,
	int circ1_x0, int circ1_y0,
	int circ2_x0, int circ2_y0,
	int circ3_x0, int circ3_y0,
	int rect_circ1_x1, int rect_circ1_x2, int rect_circ1_y1, int rect_circ1_y2,
	int rect_circ2_x1, int rect_circ2_x2, int rect_circ2_y1, int rect_circ2_y2,
	int rect_circ3_x1, int rect_circ3_x2, int rect_circ3_y1, int rect_circ3_y2,
	int rect1_x1, int rect1_x2, int rect1_y1, int rect1_y2,
	int rect2_x1, int rect2_x2, int rect2_y1, int rect2_y2,
	int rank, int numprocs, int **proc_shapes)
{
	double **prom_layer = local_layer;
	int i = rank;
	for (int j = 0; j < 8 / numprocs; j++)
	{
		if (proc_shapes[i][j] == 0)
		{
			prom_layer = SSOR_in_circle_dec(prom_layer, size, circ1_x0, circ1_y0, R,r, h, w);
		}
		if (proc_shapes[i][j] == 1)
		{
			prom_layer = SSOR_in_circle_dec(prom_layer, size, circ2_x0, circ2_y0, R, r, h, w);
		}
		if (proc_shapes[i][j] == 2)
		{
			prom_layer = SSOR_in_circle_dec(prom_layer, size, circ3_x0, circ3_y0, R, r, h, w);
		}
		if (proc_shapes[i][j] == 3)
		{
			prom_layer = SSOR_in_rect(prom_layer, size, rect_circ1_x1, rect_circ1_x2, rect_circ1_y1, rect_circ1_y2, h, w);
		}
		if (proc_shapes[i][j] == 4)
		{
			prom_layer = SSOR_in_rect(prom_layer, size, rect_circ2_x1, rect_circ2_x2, rect_circ2_y1, rect_circ2_y2, h, w);
		}
		if (proc_shapes[i][j] == 5)
		{
			prom_layer = SSOR_in_rect(prom_layer, size, rect_circ3_x1, rect_circ3_x2, rect_circ3_y1, rect_circ3_y2, h, w);
		}
		if (proc_shapes[i][j] == 6)
		{
			prom_layer = SSOR_in_rect(prom_layer, size, rect1_x1, rect1_x2, rect1_y1, rect1_y2, h, w);
		}
		if (proc_shapes[i][j] == 7)
		{
			prom_layer = SSOR_in_rect(prom_layer, size, rect2_x1, rect2_x2, rect2_y1, rect2_y2, h, w);
		}
	}
	return prom_layer;
}

void RecieveBorders(double **matr, double **local_layer, int size, double R, double r,
	int circ1_x0, int circ1_y0,
	int circ2_x0, int circ2_y0,
	int circ3_x0, int circ3_y0,
	int rect_circ1_x1, int rect_circ1_x2, int rect_circ1_y1, int rect_circ1_y2,
	int rect_circ2_x1, int rect_circ2_x2, int rect_circ2_y1, int rect_circ2_y2,
	int rect_circ3_x1, int rect_circ3_x2, int rect_circ3_y1, int rect_circ3_y2,
	int rect1_x1, int rect1_x2, int rect1_y1, int rect1_y2,
	int rect2_x1, int rect2_x2, int rect2_y1, int rect2_y2,
	int numprocs, int **proc_shapes,
	int rank)
{
	MPI_Status stat;
	int i = rank;
	double el;
	int proc;
	if (numprocs > 1)
	{
		int i = rank;
		int x1, x2, y1, y2;
		for (int j = 0; j < 8 / numprocs; j++)
		{
			if (proc_shapes[i][j] == 0)
			{
				x1 = circ1_x0 - R;
				x2 = circ1_x0 + R;
				y1 = circ1_y0 - R;
				y2 = circ1_y0 + R;
				proc = in_proc(3, numprocs, proc_shapes);
				if (rank != proc)
				{
					for (int i = x1; i <= x2; i++)
					{
						for (int j = y1; j <= y2; j++)
						{
							if ((in_circle_dec(i, j, R,r, circ1_x0, circ1_y0) == 2) &
								(in_rect(i, j, rect_circ1_x1, rect_circ1_x2, rect_circ1_y1, rect_circ1_y2) == 1))
							{
								MPI_Recv(&el, 1, MPI_DOUBLE, proc, 0, MPI_COMM_WORLD, &stat);
								local_layer[j][i] = el;
							}
						}
					}
				}
				proc = in_proc(6, numprocs, proc_shapes);
				if (rank != proc)
				{
					for (int i = x1; i <= x2; i++)
					{
						for (int j = y1; j <= y2; j++)
						{
							if ((in_circle_dec(i, j, R, r,circ1_x0, circ1_y0) == 2) &
								(in_rect(i, j, rect1_x1, rect1_x2, rect1_y1, rect1_y2) == 1))
							{
								MPI_Recv(&el, 1, MPI_DOUBLE, proc, 0, MPI_COMM_WORLD, &stat);
								local_layer[j][i] = el;
							}
						}
					}
				}
			}
			if (proc_shapes[i][j] == 1)
			{
				x1 = circ2_x0 - R;
				x2 = circ2_x0 + R;
				y1 = circ2_y0 - R;
				y2 = circ2_y0 + R;
				proc = in_proc(4, numprocs, proc_shapes);
				if (rank != proc)
				{
					for (int i = x1; i <= x2; i++)
					{
						for (int j = y1; j <= y2; j++)
						{
							if ((in_circle_dec(i, j, R, r, circ2_x0, circ2_y0) == 2) &
								(in_rect(i, j, rect_circ2_x1, rect_circ2_x2, rect_circ2_y1, rect_circ2_y2) == 1))
							{
								MPI_Recv(&el, 1, MPI_DOUBLE, proc, 0, MPI_COMM_WORLD, &stat);
								local_layer[j][i] = el;
							}
						}
					}
				}
				proc = in_proc(6, numprocs, proc_shapes);
				if (rank != proc)
				{
					for (int i = x1; i <= x2; i++)
					{
						for (int j = y1; j <= y2; j++)
						{
							if ((in_circle_dec(i, j, R, r, circ2_x0, circ2_y0) == 2) &
								(in_rect(i, j, rect1_x1, rect1_x2, rect1_y1, rect1_y2) == 1))
							{
								MPI_Recv(&el, 1, MPI_DOUBLE, proc, 0, MPI_COMM_WORLD, &stat);
								local_layer[j][i] = el;
							}
						}
					}
				}
				proc = in_proc(7, numprocs, proc_shapes);
				if (rank != proc)
				{
					for (int i = x1; i <= x2; i++)
					{
						for (int j = y1; j <= y2; j++)
						{
							if ((in_circle_dec(i, j, R, r, circ2_x0, circ2_y0) == 2) &
								(in_rect(i, j, rect2_x1, rect2_x2, rect2_y1, rect2_y2) == 1))
							{
								MPI_Recv(&el, 1, MPI_DOUBLE, proc, 0, MPI_COMM_WORLD, &stat);
								local_layer[j][i] = el;
							}
						}
					}
				}
			}
			if (proc_shapes[i][j] == 2)
			{
				x1 = circ3_x0 - R;
				x2 = circ3_x0 + R;
				y1 = circ3_y0 - R;
				y2 = circ3_y0 + R;
				proc = in_proc(5, numprocs, proc_shapes);
				if (rank != proc)
				{
					for (int i = x1; i <= x2; i++)
					{
						for (int j = y1; j <= y2; j++)
						{
							if ((in_circle_dec(i, j, R, r, circ3_x0, circ3_y0) == 2) &
								(in_rect(i, j, rect_circ3_x1, rect_circ3_x2, rect_circ3_y1, rect_circ3_y2) == 1))
							{
								MPI_Recv(&el, 1, MPI_DOUBLE, proc, 0, MPI_COMM_WORLD, &stat);
								local_layer[j][i] = el;
							}
						}
					}
				}
				proc = in_proc(7, numprocs, proc_shapes);
				if (rank != proc)
				{
					for (int i = x1; i <= x2; i++)
					{
						for (int j = y1; j <= y2; j++)
						{
							if ((in_circle_dec(i, j, R, r, circ3_x0, circ3_y0) == 2) &
								(in_rect(i, j, rect2_x1, rect2_x2, rect2_y1, rect2_y2) == 1))
							{
								MPI_Recv(&el, 1, MPI_DOUBLE, proc, 0, MPI_COMM_WORLD, &stat);
								local_layer[j][i] = el;
							}
						}
					}
				}
			}
			if (proc_shapes[i][j] == 3)
			{
				proc = in_proc(0, numprocs, proc_shapes);
				if (rank != proc)
				{
					for (int i = rect_circ1_x1; i <= rect_circ1_x2; i++)
					{
						for (int j = rect_circ1_y1; j <= rect_circ1_y2; j++)
						{
							if ((in_circle_dec(i, j, R, r, circ1_x0, circ1_y0) == 1) &
								(in_rect(i, j, rect_circ1_x1, rect_circ1_x2, rect_circ1_y1, rect_circ1_y2) == 2))
							{
								MPI_Recv(&el, 1, MPI_DOUBLE, proc, 0, MPI_COMM_WORLD, &stat);
								local_layer[j][i] = el;
							}
						}
					}
				}
			}
			if (proc_shapes[i][j] == 4)
			{
				proc = in_proc(1, numprocs, proc_shapes);
				if (rank != proc)
				{
					for (int i = rect_circ2_x1; i <= rect_circ2_x2; i++)
					{
						for (int j = rect_circ2_y1; j <= rect_circ2_y2; j++)
						{
							if ((in_circle_dec(i, j, R, r, circ2_x0, circ2_y0) == 1) &
								(in_rect(i, j, rect_circ2_x1, rect_circ2_x2, rect_circ2_y1, rect_circ2_y2) == 2))
							{
								MPI_Recv(&el, 1, MPI_DOUBLE, proc, 0, MPI_COMM_WORLD, &stat);
								local_layer[j][i] = el;
							}
						}
					}
				}
			}
			if (proc_shapes[i][j] == 5)
			{
				proc = in_proc(2, numprocs, proc_shapes);
				if (rank != proc)
				{
					for (int i = rect_circ3_x1; i <= rect_circ3_x2; i++)
					{
						for (int j = rect_circ3_y1; j <= rect_circ3_y2; j++)
						{
							if ((in_circle_dec(i, j, R, r, circ3_x0, circ3_y0) == 1) &
								(in_rect(i, j, rect_circ3_x1, rect_circ3_x2, rect_circ3_y1, rect_circ3_y2) == 2))
							{
								MPI_Recv(&el, 1, MPI_DOUBLE, proc, 0, MPI_COMM_WORLD, &stat);
								local_layer[j][i] = el;
							}
						}
					}
				}
			}
			if (proc_shapes[i][j] == 6)
			{
				proc = in_proc(0, numprocs, proc_shapes);
				if (rank != proc)
				{
					for (int i = rect1_x1; i <= rect1_x2; i++)
					{
						for (int j = rect1_y1; j <= rect1_y2; j++)
						{
							if ((in_circle_dec(i, j, R, r, circ1_x0, circ1_y0) == 1) &
								(in_rect(i, j, rect1_x1, rect1_x2, rect1_y1, rect1_y2) == 2))
							{
								MPI_Recv(&el, 1, MPI_DOUBLE, proc, 0, MPI_COMM_WORLD, &stat);
								local_layer[j][i] = el;
							}
						}
					}
				}
				proc = in_proc(1, numprocs, proc_shapes);
				if (rank != proc)
				{
					for (int i = rect1_x1; i <= rect1_x2; i++)
					{
						for (int j = rect1_y1; j <= rect1_y2; j++)
						{
							if ((in_circle_dec(i, j, R, r, circ2_x0, circ2_y0) == 1) &
								(in_rect(i, j, rect1_x1, rect1_x2, rect1_y1, rect1_y2) == 2))
							{
								MPI_Recv(&el, 1, MPI_DOUBLE, proc, 0, MPI_COMM_WORLD, &stat);
								local_layer[j][i] = el;
							}
						}
					}
				}
			}
			if (proc_shapes[i][j] == 7)
			{
				proc = in_proc(1, numprocs, proc_shapes);
				if (rank != proc)
				{
					for (int i = rect2_x1; i <= rect2_x2; i++)
					{
						for (int j = rect2_y1; j <= rect2_y2; j++)
						{
							if ((in_circle_dec(i, j, R, r, circ2_x0, circ2_y0) == 1) &
								(in_rect(i, j, rect2_x1, rect2_x2, rect2_y1, rect2_y2) == 2))
							{
								MPI_Recv(&el, 1, MPI_DOUBLE, proc, 0, MPI_COMM_WORLD, &stat);
								local_layer[j][i] = el;
							}
						}
					}
				}
				proc = in_proc(2, numprocs, proc_shapes);
				if (rank != proc)
				{
					for (int i = rect2_x1; i <= rect2_x2; i++)
					{
						for (int j = rect2_y1; j <= rect2_y2; j++)
						{
							if ((in_circle_dec(i, j, R, r, circ3_x0, circ3_y0) == 1) &
								(in_rect(i, j, rect2_x1, rect2_x2, rect2_y1, rect2_y2) == 2))
							{
								MPI_Recv(&el, 1, MPI_DOUBLE, proc, 0, MPI_COMM_WORLD, &stat);
								local_layer[j][i] = el;
							}
						}
					}
				}
			}
		}
	}
}

double norma(double **prom_layer, double **local_layer, int size, double R, double r,
	int circ1_x0, int circ1_y0,
	int circ2_x0, int circ2_y0,
	int circ3_x0, int circ3_y0,
	int rect_circ1_x1, int rect_circ1_x2, int rect_circ1_y1, int rect_circ1_y2,
	int rect_circ2_x1, int rect_circ2_x2, int rect_circ2_y1, int rect_circ2_y2,
	int rect_circ3_x1, int rect_circ3_x2, int rect_circ3_y1, int rect_circ3_y2,
	int rect1_x1, int rect1_x2, int rect1_y1, int rect1_y2,
	int rect2_x1, int rect2_x2, int rect2_y1, int rect2_y2,
	int numprocs, int **proc_shapes,
	int rank)
{
	MPI_Status stat;
	double norm = 0;
	int i = rank;
	// считаем свою
	for (int j = 0; j < 8 / numprocs; j++)
	{
		if (proc_shapes[i][j] == 0)
		{
			norm = norm + norma_circle(local_layer, prom_layer, circ1_x0, circ1_y0, R, r);
		}
		if (proc_shapes[i][j] == 1)
		{
			norm = norm + norma_circle(local_layer, prom_layer, circ2_x0, circ2_y0, R, r);
		}
		if (proc_shapes[i][j] == 2)
		{
			norm = norm + norma_circle(local_layer, prom_layer, circ3_x0, circ3_y0, R, r);
		}
		if (proc_shapes[i][j] == 3)
		{
			norm = norm + norma_rect(local_layer, prom_layer, rect_circ1_x1, rect_circ1_x2, rect_circ1_y1, rect_circ1_y2);
		}
		if (proc_shapes[i][j] == 4)
		{
			norm = norm + norma_rect(local_layer, prom_layer, rect_circ2_x1, rect_circ2_x2, rect_circ2_y1, rect_circ2_y2);
		}
		if (proc_shapes[i][j] == 5)
		{
			norm = norm + norma_rect(local_layer, prom_layer, rect_circ3_x1, rect_circ3_x2, rect_circ3_y1, rect_circ3_y2);
		}
		if (proc_shapes[i][j] == 6)
		{
			norm = norm + norma_rect(local_layer, prom_layer, rect1_x1, rect1_x2, rect1_y1, rect1_y2);
		}
		if (proc_shapes[i][j] == 7)
		{
			norm = norm + norma_rect(local_layer, prom_layer, rect2_x1, rect2_x2, rect2_y1, rect2_y2);
		}
	}
	double el;
	//отправляем и добираем
	
	if (numprocs > 1)
	{
		for (int i = 0; i < numprocs; i++)
		{
			if (i != rank)
			{
				MPI_Send(&norm, 1, MPI_DOUBLE, i, 1, MPI_COMM_WORLD);
			}
		}
		for (int i = 0; i < numprocs; i++)
		{
			if (i != rank)
			{
				MPI_Recv(&el, 1, MPI_DOUBLE, i, 1, MPI_COMM_WORLD, &stat);
				norm = norm + el;
			}
		}
	}
	return norm;
}

void SendShapes_back(double **matr, int size, double R, double r,
	int circ1_x0, int circ1_y0,
	int circ2_x0, int circ2_y0,
	int circ3_x0, int circ3_y0,
	int rect_circ1_x1, int rect_circ1_x2, int rect_circ1_y1, int rect_circ1_y2,
	int rect_circ2_x1, int rect_circ2_x2, int rect_circ2_y1, int rect_circ2_y2,
	int rect_circ3_x1, int rect_circ3_x2, int rect_circ3_y1, int rect_circ3_y2,
	int rect1_x1, int rect1_x2, int rect1_y1, int rect1_y2,
	int rect2_x1, int rect2_x2, int rect2_y1, int rect2_y2,
	int numprocs, int **proc_shapes,
	int rank)
{
	int i = rank;
	for (int j = 0; j < 8 / numprocs; j++)
	{
		if (proc_shapes[i][j] == 0)
		{
			send_circle(matr, circ1_x0, circ1_y0, R, r, 0, 0);
		}
		if (proc_shapes[i][j] == 1)
		{
			send_circle(matr, circ2_x0, circ2_y0, R, r, 0, 1);
		}
		if (proc_shapes[i][j] == 2)
		{
			send_circle(matr, circ3_x0, circ3_y0, R, r, 0, 2);
		}
		if (proc_shapes[i][j] == 3)
		{
			send_rect(matr, rect_circ1_x1, rect_circ1_x2, rect_circ1_y1, rect_circ1_y2, 0, 3);
		}
		if (proc_shapes[i][j] == 4)
		{
			send_rect(matr, rect_circ2_x1, rect_circ2_x2, rect_circ2_y1, rect_circ2_y2, 0, 4);
		}
		if (proc_shapes[i][j] == 5)
		{
			send_rect(matr, rect_circ3_x1, rect_circ3_x2, rect_circ3_y1, rect_circ3_y2, 0, 5);
		}
		if (proc_shapes[i][j] == 6)
		{
			send_rect(matr, rect1_x1, rect1_x2, rect1_y1, rect1_y2, 0, 6);
		}
		if (proc_shapes[i][j] == 7)
		{
			send_rect(matr, rect2_x1, rect2_x2, rect2_y1, rect2_y2, 0, 7);
		}
	}
}

void ReciveShapes_back(double **local_layer, int size, double R, double r,
	int circ1_x0, int circ1_y0,
	int circ2_x0, int circ2_y0,
	int circ3_x0, int circ3_y0,
	int rect_circ1_x1, int rect_circ1_x2, int rect_circ1_y1, int rect_circ1_y2,
	int rect_circ2_x1, int rect_circ2_x2, int rect_circ2_y1, int rect_circ2_y2,
	int rect_circ3_x1, int rect_circ3_x2, int rect_circ3_y1, int rect_circ3_y2,
	int rect1_x1, int rect1_x2, int rect1_y1, int rect1_y2,
	int rect2_x1, int rect2_x2, int rect2_y1, int rect2_y2,
	int numprocs, int **proc_shapes)
{
	for (int i = 1; i < numprocs; i++)
	{
		for (int j = 0; j < 8 / numprocs; j++)
		{
			if (proc_shapes[i][j] == 0)
			{
				recive_circle(local_layer, circ1_x0, circ1_y0, R, r, i, 0);
			}
			if (proc_shapes[i][j] == 1)
			{
				recive_circle(local_layer, circ2_x0, circ2_y0, R, r, i, 1);
			}
			if (proc_shapes[i][j] == 2)
			{
				recive_circle(local_layer, circ3_x0, circ3_y0, R, r, i, 2);
			}
			if (proc_shapes[i][j] == 3)
			{
				recive_rect(local_layer, rect_circ1_x1, rect_circ1_x2, rect_circ1_y1, rect_circ1_y2, i, 3);
			}
			if (proc_shapes[i][j] == 4)
			{
				recive_rect(local_layer, rect_circ2_x1, rect_circ2_x2, rect_circ2_y1, rect_circ2_y2, i, 4);
			}
			if (proc_shapes[i][j] == 5)
			{
				recive_rect(local_layer, rect_circ3_x1, rect_circ3_x2, rect_circ3_y1, rect_circ3_y2, i, 5);
			}
			if (proc_shapes[i][j] == 6)
			{
				recive_rect(local_layer, rect1_x1, rect1_x2, rect1_y1, rect1_y2, i, 6);
			}
			if (proc_shapes[i][j] == 7)
			{
				recive_rect(local_layer, rect2_x1, rect2_x2, rect2_y1, rect2_y2, i, 7);
			}
		}
	}
}

// Параллельный метод
void SSOR_paral(double **matr, int size, double R, double r, double h, double w, double eps,
	int circ1_x0, int circ1_y0,
	int circ2_x0, int circ2_y0,
	int circ3_x0, int circ3_y0,
	int circ4_x0, int circ4_y0,
	int rect1_x1, int rect1_x2, int rect1_y1, int rect1_y2,
	int rect2_x1, int rect2_x2, int rect2_y1, int rect2_y2,
	int rect3_x1, int rect3_x2, int rect3_y1, int rect3_y2,
	int rect4_x1, int rect4_x2, int rect4_y1, int rect4_y2)
{
	int rank, numprocs;
	MPI_Status stat;
	MPI_Init(NULL, NULL);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	int shapes_count = 8 / numprocs;
	double **shapes[8];
	double ***local_layer = new double **[shapes_count];
	int sizes_x[8], sizes_y[8];
	//double **prom_layer = new double *[size];
	int ** proc_shapes = MakeShapes(numprocs);

	double step = dec_to_pol_fi(R+1, 0.5);
	int N_fi = 360 / step;
	
	if (rank == 0)
	{
		start_fill_matr(matr, size, R, r,
			circ1_x0, circ1_y0,
			circ2_x0, circ2_y0,
			circ3_x0, circ3_y0,
			circ4_x0, circ4_y0,
			rect1_x1, rect1_x2, rect1_y1, rect1_y2,
			rect2_x1, rect2_x2, rect2_y1, rect2_y2,
			rect3_x1, rect3_x2, rect3_y1, rect3_y2,
			rect4_x1, rect4_x2, rect4_y1, rect4_y2);
		shapes[0] = new double *[R + 1];
		sizes_x[0] = R;
		sizes_y[0] = N_fi;
		for (int i = 0; i < R + 1; i++)
		{
			shapes[0][i] = new double[N_fi];
			for (int j = 0; j < N_fi; j++)
			{
				shapes[0][i][j] = 0;
			}
		}
		for (int i = 0; i < R + 1; i++)
		{
			for (int j = 0; j < N_fi; j++)
			{
				shapes[0][i][j] = 0;
				int x = circ1_x0 - pol_to_dec_x(i, j* step);
				int y = circ1_y0 - pol_to_dec_y(i, j* step);
				int jj = j * step;
				if (in_circle_pol(i, j*step, R, r, circ1_x0, circ1_y0) == 1)
				{
					shapes[0][i][j] = matr[y][x];
				}
				if ((in_circle_pol(i, j*step, R, r, circ1_x0, circ1_y0) == 2) & !(in_rect(i, j, rect1_x1, rect1_x2, rect1_y1, rect1_y2) == 1))
				{
					shapes[0][i][j] = matr[y][x];
				}
			}
			//printShape(shapes[0], sizes_x[0], sizes_y[0]);
		}
		cout << endl << endl;

	}

	MPI_Barrier(MPI_COMM_WORLD);

	// ----------------------------   Рассылка кусков по процессорам    --------------------
	if (numprocs > 1)
	{
		/*if (rank == 0)
		{
			SendShapes(matr, local_layer, size, R, r,
				circ1_x0, circ1_y0,
				circ2_x0, circ2_y0,
				circ3_x0, circ3_y0,
				rect_circ1_x1, rect_circ1_x2, rect_circ1_y1, rect_circ1_y2,
				rect_circ2_x1, rect_circ2_x2, rect_circ2_y1, rect_circ2_y2,
				rect_circ3_x1, rect_circ3_x2, rect_circ3_y1, rect_circ3_y2,
				rect1_x1, rect1_x2, rect1_y1, rect1_y2,
				rect2_x1, rect2_x2, rect2_y1, rect2_y2,
				numprocs, proc_shapes);
		}
		else
		{
			ReciveShapes(local_layer, size, R, r,
				circ1_x0, circ1_y0,
				circ2_x0, circ2_y0,
				circ3_x0, circ3_y0,
				rect_circ1_x1, rect_circ1_x2, rect_circ1_y1, rect_circ1_y2,
				rect_circ2_x1, rect_circ2_x2, rect_circ2_y1, rect_circ2_y2,
				rect_circ3_x1, rect_circ3_x2, rect_circ3_y1, rect_circ3_y2,
				rect1_x1, rect1_x2, rect1_y1, rect1_y2,
				rect2_x1, rect2_x2, rect2_y1, rect2_y2,
				numprocs, proc_shapes,
				rank);
		}*/
	}
	else
	{
		
	}
	MPI_Barrier(MPI_COMM_WORLD);

	// ----------------------------   Конец рассылки кусков по процессорам    --------------------

/*
	// ----------------------------   SSOR    --------------------
	// Первая итерация
	prom_layer = SSOR_1(local_layer, size, R,r, h, w,
		circ1_x0, circ1_y0,
		circ2_x0, circ2_y0,
		circ3_x0, circ3_y0,
		rect_circ1_x1, rect_circ1_x2, rect_circ1_y1, rect_circ1_y2,
		rect_circ2_x1, rect_circ2_x2, rect_circ2_y1, rect_circ2_y2,
		rect_circ3_x1, rect_circ3_x2, rect_circ3_y1, rect_circ3_y2,
		rect1_x1, rect1_x2, rect1_y1, rect1_y2,
		rect2_x1, rect2_x2, rect2_y1, rect2_y2,
		rank, numprocs, proc_shapes);


	// Норма
	double norm = norma(prom_layer, local_layer, size, R, r,
		circ1_x0, circ1_y0,
		circ2_x0, circ2_y0,
		circ3_x0, circ3_y0,
		rect_circ1_x1, rect_circ1_x2, rect_circ1_y1, rect_circ1_y2,
		rect_circ2_x1, rect_circ2_x2, rect_circ2_y1, rect_circ2_y2,
		rect_circ3_x1, rect_circ3_x2, rect_circ3_y1, rect_circ3_y2,
		rect1_x1, rect1_x2, rect1_y1, rect1_y2,
		rect2_x1, rect2_x2, rect2_y1, rect2_y2,
		numprocs, proc_shapes,
		rank);
	local_layer = prom_layer;

	// Отправка вычисленных границ
	if (eps < norm)
	{
		SendBorders(matr, local_layer, size, R, r,
			circ1_x0, circ1_y0,
			circ2_x0, circ2_y0,
			circ3_x0, circ3_y0,
			rect_circ1_x1, rect_circ1_x2, rect_circ1_y1, rect_circ1_y2,
			rect_circ2_x1, rect_circ2_x2, rect_circ2_y1, rect_circ2_y2,
			rect_circ3_x1, rect_circ3_x2, rect_circ3_y1, rect_circ3_y2,
			rect1_x1, rect1_x2, rect1_y1, rect1_y2,
			rect2_x1, rect2_x2, rect2_y1, rect2_y2,
			numprocs, proc_shapes,
			rank);
	}
	
	while (eps < norm)
	{
		// Получаем границы
		RecieveBorders(matr, local_layer, size, R, r,
			circ1_x0, circ1_y0,
			circ2_x0, circ2_y0,
			circ3_x0, circ3_y0,
			rect_circ1_x1, rect_circ1_x2, rect_circ1_y1, rect_circ1_y2,
			rect_circ2_x1, rect_circ2_x2, rect_circ2_y1, rect_circ2_y2,
			rect_circ3_x1, rect_circ3_x2, rect_circ3_y1, rect_circ3_y2,
			rect1_x1, rect1_x2, rect1_y1, rect1_y2,
			rect2_x1, rect2_x2, rect2_y1, rect2_y2,
			numprocs, proc_shapes,
			rank);
		
		// Подсчитываем слой
		prom_layer = SSOR_1(local_layer, size, R, r, h, w,
			circ1_x0, circ1_y0,
			circ2_x0, circ2_y0,
			circ3_x0, circ3_y0,
			rect_circ1_x1, rect_circ1_x2, rect_circ1_y1, rect_circ1_y2,
			rect_circ2_x1, rect_circ2_x2, rect_circ2_y1, rect_circ2_y2,
			rect_circ3_x1, rect_circ3_x2, rect_circ3_y1, rect_circ3_y2,
			rect1_x1, rect1_x2, rect1_y1, rect1_y2,
			rect2_x1, rect2_x2, rect2_y1, rect2_y2,
			rank, numprocs, proc_shapes);

		norm = norma(prom_layer, local_layer, size, R, r,
			circ1_x0, circ1_y0,
			circ2_x0, circ2_y0,
			circ3_x0, circ3_y0,
			rect_circ1_x1, rect_circ1_x2, rect_circ1_y1, rect_circ1_y2,
			rect_circ2_x1, rect_circ2_x2, rect_circ2_y1, rect_circ2_y2,
			rect_circ3_x1, rect_circ3_x2, rect_circ3_y1, rect_circ3_y2,
			rect1_x1, rect1_x2, rect1_y1, rect1_y2,
			rect2_x1, rect2_x2, rect2_y1, rect2_y2,
			numprocs, proc_shapes,
			rank);

		local_layer = prom_layer;

		if (eps < norm)
		{
			SendBorders(matr, local_layer, size, R, r,
				circ1_x0, circ1_y0,
				circ2_x0, circ2_y0,
				circ3_x0, circ3_y0,
				rect_circ1_x1, rect_circ1_x2, rect_circ1_y1, rect_circ1_y2,
				rect_circ2_x1, rect_circ2_x2, rect_circ2_y1, rect_circ2_y2,
				rect_circ3_x1, rect_circ3_x2, rect_circ3_y1, rect_circ3_y2,
				rect1_x1, rect1_x2, rect1_y1, rect1_y2,
				rect2_x1, rect2_x2, rect2_y1, rect2_y2,
				numprocs, proc_shapes,
				rank);
		}
	}

	MPI_Barrier(MPI_COMM_WORLD);
	// ----------------------------   Конец SSOR    --------------------
	*/
	// ----------------------------   Сборка кусок в 1 процессор   --------------------
	
	if (numprocs > 1)
	{
		if (rank != 0)
		{
			/*SendShapes_back(local_layer, size, R, r,
				circ1_x0, circ1_y0,
				circ2_x0, circ2_y0,
				circ3_x0, circ3_y0,
				rect_circ1_x1, rect_circ1_x2, rect_circ1_y1, rect_circ1_y2,
				rect_circ2_x1, rect_circ2_x2, rect_circ2_y1, rect_circ2_y2,
				rect_circ3_x1, rect_circ3_x2, rect_circ3_y1, rect_circ3_y2,
				rect1_x1, rect1_x2, rect1_y1, rect1_y2,
				rect2_x1, rect2_x2, rect2_y1, rect2_y2,
				numprocs, proc_shapes, rank);*/
		}
		else
		{
			/*
			ReciveShapes_back(local_layer, size, R, r,
				circ1_x0, circ1_y0,
				circ2_x0, circ2_y0,
				circ3_x0, circ3_y0,
				rect_circ1_x1, rect_circ1_x2, rect_circ1_y1, rect_circ1_y2,
				rect_circ2_x1, rect_circ2_x2, rect_circ2_y1, rect_circ2_y2,
				rect_circ3_x1, rect_circ3_x2, rect_circ3_y1, rect_circ3_y2,
				rect1_x1, rect1_x2, rect1_y1, rect1_y2,
				rect2_x1, rect2_x2, rect2_y1, rect2_y2,
				numprocs, proc_shapes);
				*/
		}
	}

	if (rank == 0)
	{
		for (int i = 0; i < R + 4; i++)
		{
			for (int j = 0; j < N_fi; j++)
			{
				int x = circ1_x0 - pol_to_dec_x(i, j* step);
				int y = circ1_y0 - pol_to_dec_y(i, j* step);
				if (in_circle_pol(i, j*step, R, r, circ1_x0, circ1_y0) == 1)
				{
					matr[y][x] = 0;
				}
				int jj = j * step;
				if (in_circle_pol(i, j*step, R, r, circ1_x0, circ1_y0) == 1)
				{
					//matr[y][x]=shapes[0][i][j];
				}
				if ((in_circle_pol(i, j*step, R, r, circ1_x0, circ1_y0) == 2) & !(in_rect(x, y, rect1_x1, rect1_x2, rect1_y1, rect1_y2) == 1))
				{
					//matr[y][x] = shapes[0][i][j];
				}

			}
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);
	
	

	// ----------------------------   Конец    --------------------

	if (rank == 0)
	{
		printMatr(matr, size);
		printShape(shapes[0], sizes_x[0], sizes_y[0]);
	}

	
	delete local_layer;
	delete proc_shapes;



	MPI_Finalize();
}