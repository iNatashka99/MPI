#include "Circles.h"
#include "Somefunc.h"
#include <math.h>
#include <mpi.h>

# define M_PI           3.14159265358979323846 

// Уравнение окружности
int circle_x(int x0, double R, double fi)
{
	return int(x0 + R * cos(fi));
}
int circle_y(int y0, double R, double fi)
{
	return int(y0 + R * sin(fi));
}

// Принадлежность кольцу в декартовой
int in_circle_dec(int x, int y, double R, double r, int x0, int y0)
{
	if ((pow(x - x0, 2) + pow(y - y0, 2) >= R * R) & (pow(x - x0, 2) + pow(y - y0, 2) <= pow(R + 1, 2)))
		return 2;
	if ((pow(x - x0, 2) + pow(y - y0, 2) >= r * r) & (pow(x - x0, 2) + pow(y - y0, 2) <= pow(r + 1, 2)))
		return 2;
	if ((pow(x - x0, 2) + pow(y - y0, 2) < R * R) & (pow(x - x0, 2) + pow(y - y0, 2) > R * R / 4))
		return 1;
	return 0;
}


//Перевод из декартовых в полярные
double dec_to_pol_r(double x, double y)
{
	double r = pow(pow(x, 2) + pow(y, 2), 0.5);
	return r;
}
double dec_to_pol_fi(double x, double y)
{
	double fi = 0;
	if ((x > 0) && (y >= 0))
	{
		fi = atan(y / x)*(180 / M_PI);
	}
	if ((x > 0) && (y < 0))
	{
		fi = atan(y / x)*(180 / M_PI) + 360;
	}
	if (x < 0)
	{
		fi = atan(y / x)*(180 / M_PI) + 180;
	}
	if ((x == 0) && (y > 0))
	{
		fi = 90;
	}
	if ((x == 0) && (y < 0))
	{
		fi = 270;
	}
	return fi;
}

//Перевод из полярных в декартовые
double pol_to_dec_x(double r, double fi)
{
	double x = r * cos(fi / 180 * M_PI);
	return x;
}
double pol_to_dec_y(double r, double fi)
{
	double y = r * sin(fi / 180 * M_PI);
	return y;
}

// Принадлежность кольцу в полярной
int in_circle_pol(int rr, int fi, int R, int r, int x0, int y0)
{
	if (rr == R)
		return 2;
	if (rr == r)
		return 2;
	if ((rr < R) && (rr > r))
		return 1;
	return 0;
}



// Отправка кольца
void send_circle(double **data, int x0, int y0, double R, double r, int proc, int shape)
{
	int x1 = x0 - R;
	int x2 = x0 + R;
	int y1 = y0 - R;
	int y2 = y0 + R;
	for (int i = x1; i <= x2; i++)
	{
		for (int j = y1; j <= y2; j++)
		{
			if (in_circle_dec(i, j, R,r, x0, y0))
			{
				MPI_Send(&data[j][i], 1, MPI_DOUBLE, proc, shape, MPI_COMM_WORLD);
			}
		}
	}
}

// Прием кольца
void recive_circle(double **data, int x0, int y0, double R, double r, int proc, int shape)
{
	MPI_Status stat;
	int x1 = x0 - R;
	int x2 = x0 + R;
	int y1 = y0 - R;
	int y2 = y0 + R;
	double el;
	for (int i = x1; i <= x2; i++)
	{
		for (int j = y1; j <= y2; j++)
		{
			if (in_circle_dec(i, j, R, r, x0, y0))
			{
				MPI_Recv(&el, 1, MPI_DOUBLE, proc, shape, MPI_COMM_WORLD, &stat);
				data[j][i] = el;
			}
		}
	}
}

// Запись для нулевого процессора
void set_circle(double **data, double **layer, int x0, int y0, double R, double r)
{
	int x1 = x0 - R;
	int x2 = x0 + R;
	int y1 = y0 - R;
	int y2 = y0 + R;
	for (int i = x1; i <= x2; i++)
	{
		for (int j = y1; j <= y2; j++)
		{
			if (in_circle_dec(i, j, R, r, x0, y0))
			{
				layer[j][i] = data[j][i];
			}
		}
	}
}

// Норма в кольце
double norma_circle(double **layer1, double **layer2, int x0, int y0, double R, double r )
{
	int x1 = x0 - R;
	int x2 = x0 + R;
	int y1 = y0 - R;
	int y2 = y0 + R;
	double nrm = 0;
	for (int i = x1; i <= x2; i++)
	{
		for (int j = y1; j <= y2; j++)
		{
			if (in_circle_dec(i, j, R, r, x0, y0))
			{
				nrm = nrm + fabs(layer1[j][i] - layer2[j][i]);
			}
		}
	}
	return nrm;
}