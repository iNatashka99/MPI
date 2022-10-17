#include "SSOR.h"
#include "Circles.h"
#include "Rects.h"
#include <math.h>
#include "Somefunc.h"

// Гаус в узле
double gaus(double **old_layer, double **new_layer, int i, int j, double h)
{
	return 1.0 / 4 * (old_layer[j][i + 1] + old_layer[j + 1][i] + new_layer[j][i - 1] + new_layer[j - 1][i] - pow(h, 2) * F(j * h, i * h));
}

//Коррекция в узле
double corr(double **old_layer, double **new_layer, int i, int j, double w)
{
	return (1 - w) * old_layer[j][i] + w * new_layer[j][i];
}

// Гаус в кольце 
double ** GAUS_in_circle_dec1(double **old_layer, int size, int x0, int y0, double R, double r, double h)
{
	double **new_layer;
	new_layer = new double *[size / 3];
	for (int i = 0; i < size / 3; i++)
	{
		new_layer[i] = new double[size];
		for (int j = 0; j < size; j++)
		{
			new_layer[i][j] = old_layer[i][j];
		}
	}
	int x1 = x0 - R;
	int x2 = x0 + R;
	int y1 = y0 - R;
	int y2 = y0 + R;
	for (int i = x1; i <= x2; i++)
	{
		for (int j = y1; j <= y2; j++)
		{
			int incirc = in_circle_dec(i, j, R, r, x0, y0);
			if (incirc == 1)
			{
				new_layer[j][i] = gaus(old_layer, new_layer, i, j, h);
			}
			else if (incirc == 2)
			{
				new_layer[j][i] = old_layer[j][i];
			}
		}
	}
	return new_layer;
}
double ** GAUS_in_circle_dec2(double **old_layer, int size, int x0, int y0, double R, double r, double h)
{
	double **new_layer;
	new_layer = new double *[size / 3];
	for (int i = 0; i < size / 3; i++)
	{
		new_layer[i] = new double[size];
		for (int j = 0; j < size; j++)
		{
			new_layer[i][j] = old_layer[i][j];
		}
	}
	int x1 = x0 - R;
	int x2 = x0 + R;
	int y1 = y0 - R;
	int y2 = y0 + R;
	for (int i = x2; i >= x1; i--)
	{
		for (int j = y2; j >= y1; j--)
		{
			int incirc = in_circle_dec(i, j, R, r, x0, y0);
			if (incirc == 1)
			{
				new_layer[j][i] = gaus(new_layer, old_layer, i, j, h);
			}
			else if (incirc == 2)
			{
				new_layer[j][i] = old_layer[j][i];
			}
		}
	}
	return new_layer;
}

// Коррекция в кольце
double ** corr_in_circle_dec(double **old_layer, double **new_layer, int size, int x0, int y0, double R, double r, double w)
{
	double **corr_layer;
	corr_layer = new double *[size / 3];
	for (int i = 0; i < size / 3; i++)
	{
		corr_layer[i] = new double[size];
		for (int j = 0; j < size; j++)
		{
			corr_layer[i][j] = old_layer[i][j];
		}
	}
	int x1 = x0 - R;
	int x2 = x0 + R;
	int y1 = y0 - R;
	int y2 = y0 + R;
	for (int i = x1; i <= x2; i++)
	{
		for (int j = y1; j <= y2; j++)
		{
			int incirc = in_circle_dec(i, j, R, r, x0, y0);
			if (incirc == 1)
			{
				corr_layer[j][i] = corr(old_layer, new_layer, i, j, w);
			}
			else if (incirc == 2)
			{
				corr_layer[j][i] = new_layer[j][i];
			}
		}
	}
	return corr_layer;
}

// SSOR в кольце 
double ** SSOR_in_circle_dec(double **old_layer, int size, int x0, int y0, double R, double r, double h, double w)
{
	double **new_layer;
	new_layer = new double *[size / 3];
	for (int i = 0; i < size / 3; i++)
	{
		new_layer[i] = new double[size];
		for (int j = 0; j < size; j++)
		{
			new_layer[i][j] = old_layer[i][j];
		}
	}
	new_layer = GAUS_in_circle_dec1(old_layer, size, x0, y0, R, r, h);
	new_layer = corr_in_circle_dec(old_layer, new_layer, size, x0, y0, R, r, w);
	old_layer = new_layer;
	new_layer = GAUS_in_circle_dec2(old_layer, size, x0, y0, R,r, h);
	new_layer = corr_in_circle_dec(old_layer, new_layer, size, x0, y0, R, r, w);
	return new_layer;
}

// Гаус в прямоугольнике 
double ** GAUS_in_rect1(double **old_layer, int size, int x1, int x2, int y1, int y2, double h)
{
	double **new_layer;
	new_layer = new double *[size / 3];
	for (int i = 0; i < size / 3; i++)
	{
		new_layer[i] = new double[size];
		for (int j = 0; j < size; j++)
		{
			new_layer[i][j] = old_layer[i][j];
		}
	}
	for (int i = x1; i <= x2; i++)
	{
		for (int j = y1; j <= y2; j++)
		{
			int inrect = in_rect(i, j, x1, x2, y1, y2);
			if (inrect == 1)
			{
				new_layer[j][i] = gaus(old_layer, new_layer, i, j, h);
			}
			else if (inrect == 2)
			{
				new_layer[j][i] = old_layer[j][i];
			}
		}
	}
	return new_layer;
}
double ** GAUS_in_rect2(double **old_layer, int size, int x1, int x2, int y1, int y2, double h)
{
	double **new_layer;
	new_layer = new double *[size / 3];
	for (int i = 0; i < size / 3; i++)
	{
		new_layer[i] = new double[size];
		for (int j = 0; j < size; j++)
		{
			new_layer[i][j] = old_layer[i][j];
		}
	}
	for (int i = x2; i >= x1; i--)
	{
		for (int j = y2; j >= y1; j--)
		{
			int inrect = in_rect(i, j, x1, x2, y1, y2);
			if (inrect == 1)
			{
				new_layer[j][i] = gaus(old_layer, new_layer, i, j, h);
			}
			else if (inrect == 2)
			{
				new_layer[j][i] = old_layer[j][i];
			}
		}
	}
	return new_layer;
}

// Коррекция в прямоугольнике
double ** corr_in_rect(double **old_layer, double **new_layer, int size, int x1, int x2, int y1, int y2, double w)
{
	double **corr_layer;
	corr_layer = new double *[size / 3];
	for (int i = 0; i < size / 3; i++)
	{
		corr_layer[i] = new double[size];
		for (int j = 0; j < size; j++)
		{
			corr_layer[i][j] = old_layer[i][j];
		}
	}
	for (int i = x1; i <= x2; i++)
	{
		for (int j = y1; j <= y2; j++)
		{
			int inrect = in_rect(i, j, x1, x2, y1, y2);
			if (inrect == 1)
			{
				corr_layer[j][i] = corr(old_layer, new_layer, i, j, w);
			}
			else if (inrect == 2)
			{
				corr_layer[j][i] = new_layer[j][i];
			}
		}
	}
	return corr_layer;
}

// SSOR в прямоугольнике 
double ** SSOR_in_rect(double **old_layer, int size, int x1, int x2, int y1, int y2, double h, double w)
{
	double **new_layer;
	new_layer = new double *[size / 3];
	for (int i = 0; i < size / 3; i++)
	{
		new_layer[i] = new double[size];
		for (int j = 0; j < size; j++)
		{
			new_layer[i][j] = old_layer[i][j];
		}
	}
	new_layer = GAUS_in_rect1(old_layer, size, x1, x2, y1, y2, h);
	new_layer = corr_in_rect(old_layer, new_layer, size, x1, x2, y1, y2, w);
	old_layer = new_layer;
	new_layer = GAUS_in_rect1(old_layer, size, x1, x2, y1, y2, h);
	new_layer = corr_in_rect(old_layer, new_layer, size, x1, x2, y1, y2, w);
	return new_layer;
}