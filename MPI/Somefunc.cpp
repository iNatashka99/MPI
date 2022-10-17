#include "Somefunc.h"
#include <math.h>

// функция u
double U(double x, double y)
{
	return pow(x,2) + pow(y,2);
}

// функция f
double F(double x, double y)
{
	return 4;
}

// Граница
double border(double x, double y)
{
	return 1;
}

// Начальное приближение
double startfill(double x, double y)
{
	if (x > 1.0 / 2)
		return int((x + y) * 10000000) % 9 + 1;
	else
		return 2;
}

// Распределение кусочков по процессорам
int ** MakeShapes(int num_procs)
{
	int ** procs_shapes = new int *[num_procs];
	for (int i = 0; i < num_procs; i++)
	{
		procs_shapes[i] = new int[8 / num_procs];
	}
	if (num_procs == 1)
	{
		for (int i = 0; i < 8; i++)
		{
			procs_shapes[0][i] = i + 1;
		}
	}
	if (num_procs == 2)
	{
		procs_shapes[0][0] = 1;
		procs_shapes[0][1] = 4;
		procs_shapes[0][2] = 2;
		procs_shapes[0][3] = 3;
		procs_shapes[1][0] = 5;
		procs_shapes[1][1] = 8;
		procs_shapes[1][2] = 6;
		procs_shapes[1][3] = 7;
	}
	if (num_procs == 4)
	{
		procs_shapes[0][0] = 1;
		procs_shapes[0][1] = 2;
		procs_shapes[1][0] = 4;
		procs_shapes[1][1] = 3;
		procs_shapes[2][0] = 5;
		procs_shapes[2][1] = 6;
		procs_shapes[3][0] = 8;
		procs_shapes[3][1] = 7;
	}
	return procs_shapes;
}

// Есть ли пересечение
bool have_intersection(int shape1, int shape2)
{
	if (((shape1 == 0) & (shape2 == 3)) || ((shape1 == 3) & (shape2 == 0)))
		return 1;
	if (((shape1 == 0) & (shape2 == 6)) || ((shape1 == 6) & (shape2 == 0)))
		return 1;
	if (((shape1 == 1) & (shape2 == 6)) || ((shape1 == 6) & (shape2 == 1)))
		return 1;
	if (((shape1 == 1) & (shape2 == 4)) || ((shape1 == 4) & (shape2 == 1)))
		return 1;
	if (((shape1 == 1) & (shape2 == 7)) || ((shape1 == 7) & (shape2 == 1)))
		return 1;
	if (((shape1 == 2) & (shape2 == 5)) || ((shape1 == 5) & (shape2 == 2)))
		return 1;
	if (((shape1 == 2) & (shape2 == 7)) || ((shape1 == 7) & (shape2 == 2)))
		return 1;
	return 0;
}

// В каком процесоре кусок
int in_proc(int shape, int num_procs, int ** procs_shapes)
{
	for (int i = 0; i < num_procs; i++)
	{
		for (int j = 0; j < 8 / num_procs; j++)
		{
			if (procs_shapes[i][j] == shape)
			{
				return i;
			}
		}
	}
}