#include "ForMatr.h"
#include "Circles.h"
#include "Rects.h"
#include "Somefunc.h"
#include <iostream>
#include <iomanip>

using namespace std;

//Вывод матрицы A
void printMatr(double **data, int size)
{
	//cout << setprecision(1);
	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < size; j++)
		{
			cout << setw(3);
			//cout << setw(10);
			if (data[i][j] == 0)
				cout << " ";
			else
				cout << data[i][j];
		}
		cout << "|" << endl;
	}
	cout << endl;
}

//Вывод куска
void printShape(double **data, int size1, int size2)
{
	//cout << setprecision(1);
	for (int i = 0; i < size1; i++)
	{
		for (int j = 0; j < size2; j++)
		{
			cout << setw(3);
			//cout << setw(10);
			if (data[i][j] == 0)
				cout << " ";
			else
				cout << data[i][j];
		}
		cout << "|" << endl;
	}
	cout << endl;
}

//Очистка матрицы A
void clearMatr(double **data, int size)
{
	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < size; j++)
		{
			data[i][j] = 0;
		}
	}
}

// копирование матрицы
double ** copymatr(double **A1, int size)
{
	double **A2;
	A2 = new double *[size];
	for (int i = 0; i < size; i++)
	{
		A2[i] = new double[size];
		for (int j = 0; j < size; j++)
		{
			A2[i][j] = A1[i][j];
		}
	}
	return A2;
}

//Начальное приближение
void start_fill_matr(double **data, int size, double R, double r,
	int circ1_x0, int circ1_y0,
	int circ2_x0, int circ2_y0,
	int circ3_x0, int circ3_y0,
	int circ4_x0, int circ4_y0,
	int rect1_x1, int rect1_x2, int rect1_y1, int rect1_y2,
	int rect2_x1, int rect2_x2, int rect2_y1, int rect2_y2,
	int rect3_x1, int rect3_x2, int rect3_y1, int rect3_y2,
	int rect4_x1, int rect4_x2, int rect4_y1, int rect4_y2)
{
	int rr, fi;
	for (int i = 0; i < size ; i++)
	{
		for (int j = 0; j < size; j++)
		{
			rr = dec_to_pol_r(j - circ1_y0, i - circ1_x0);
			fi = dec_to_pol_fi(j - circ1_y0, i - circ1_x0);
			data[j][i] = 0;
			if (in_circle_pol(rr, fi, R, r, circ1_x0, circ1_y0) == 1)
			{
				data[j][i] = startfill(j / double(size), i / double(size));
			}
			if ((in_circle_pol(rr, fi, R, r, circ1_x0, circ1_y0) == 2) & !(in_rect(i, j, rect1_x1, rect1_x2, rect1_y1, rect1_y2) == 1))
			{
				data[j][i] = border(j / double(size), i / double(size));
			}
			if ((in_rect(i, j, rect1_x1, rect1_x2, rect1_y1, rect1_y2) == 2) & !(in_circle_pol(rr, fi, R, r, circ1_x0, circ1_y0) == 1) & !(in_rect(i, j, rect2_x1, rect2_x2, rect2_y1, rect2_y2) == 1))
			{
				data[j][i] = border(j / double(size), i / double(size));
			}
			if (in_rect(i, j, rect1_x1, rect1_x2, rect1_y1, rect1_y2) == 1)
			{
				data[j][i] = startfill(j / double(size), i / double(size));
			}

			rr = dec_to_pol_r(j - circ2_y0, i - circ2_x0);
			fi = dec_to_pol_fi(j - circ2_y0, i - circ2_x0);
			if (in_circle_pol(rr, fi, R, r, circ2_x0, circ2_y0) == 1)
			{
				data[j][i] = startfill(j / double(size), i / double(size));
			}
			if ((in_circle_pol(rr, fi, R, r, circ2_x0, circ2_y0) == 2) & !(in_rect(i, j, rect2_x1, rect2_x2, rect2_y1, rect2_y2) == 1))
			{
				data[j][i] = border(j / double(size), i / double(size));
			}
			if ((in_rect(i, j, rect2_x1, rect2_x2, rect2_y1, rect2_y2) == 2) & !(in_circle_pol(rr, fi, R, r, circ2_x0, circ2_y0) == 1) & !(in_rect(i, j, rect1_x1, rect1_x2, rect1_y1, rect1_y2) == 1) & !(in_rect(i, j, rect3_x1, rect3_x2, rect3_y1, rect3_y2) == 1))
			{
				data[j][i] = border(j / double(size), i / double(size));
			}
			if (in_rect(i, j, rect2_x1, rect2_x2, rect2_y1, rect2_y2) == 1)
			{
				data[j][i] = startfill(j / double(size), i / double(size));
			}

			rr = dec_to_pol_r(j - circ3_y0, i - circ3_x0);
			fi = dec_to_pol_fi(j - circ3_y0, i - circ3_x0);
			if (in_circle_pol(rr, fi, R, r, circ3_x0, circ3_y0) == 1)
			{
				data[j][i] = startfill(j / double(size), i / double(size));
			}
			if ((in_circle_pol(rr, fi, R, r, circ3_x0, circ3_y0) == 2) & !(in_rect(i, j, rect3_x1, rect3_x2, rect3_y1, rect3_y2) == 1))
			{
				data[j][i] = border(j / double(size), i / double(size));
			}
			if ((in_rect(i, j, rect3_x1, rect3_x2, rect3_y1, rect3_y2) == 2) & !(in_circle_pol(rr, fi, R, r, circ3_x0, circ3_y0) == 1) & !(in_rect(i, j, rect2_x1, rect2_x2, rect2_y1, rect2_y2) == 1) & !(in_rect(i, j, rect4_x1, rect4_x2, rect4_y1, rect4_y2) == 1))
			{
				data[j][i] = border(j / double(size), i / double(size));
			}
			if (in_rect(i, j, rect3_x1, rect3_x2, rect3_y1, rect3_y2) == 1)
			{
				data[j][i] = startfill(j / double(size), i / double(size));
			}

			rr = dec_to_pol_r(j - circ4_y0, i - circ4_x0);
			fi = dec_to_pol_fi(j - circ4_y0, i - circ4_x0);
			if (in_circle_pol(rr, fi, R, r, circ4_x0, circ4_y0) == 1)
			{
				data[j][i] = startfill(j / double(size), i / double(size));
			}
			if ((in_circle_pol(rr, fi, R, r, circ4_x0, circ4_y0) == 2) & !(in_rect(i, j, rect4_x1, rect4_x2, rect4_y1, rect4_y2) == 1))
			{
				data[j][i] = border(j / double(size), i / double(size));
			}
			if ((in_rect(i, j, rect4_x1, rect4_x2, rect4_y1, rect4_y2) == 2) & !(in_circle_pol(rr, fi, R, r, circ4_x0, circ4_y0) == 1) & !(in_rect(i, j, rect3_x1, rect3_x2, rect3_y1, rect3_y2) == 1))
			{
				data[j][i] = border(j / double(size), i / double(size));
			}
			if (in_rect(i, j, rect4_x1, rect4_x2, rect4_y1, rect4_y2) == 1)
			{
				data[j][i] = startfill(j / double(size), i / double(size));
			}

		}
	}
}



//Начальное разбиение на 8 частей
void start_shape_matr(double **data, int size, double R, double r,
	int circ1_x0, int circ1_y0,
	int circ2_x0, int circ2_y0,
	int circ3_x0, int circ3_y0,
	int circ4_x0, int circ4_y0,
	int rect1_x1, int rect1_x2, int rect1_y1, int rect1_y2,
	int rect2_x1, int rect2_x2, int rect2_y1, int rect2_y2,
	int rect3_x1, int rect3_x2, int rect3_y1, int rect3_y2, 
	int rect4_x1, int rect4_x2, int rect4_y1, int rect4_y2)
{
	int rr, fi;
	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < size; j++)
		{
			rr = dec_to_pol_r(j - circ1_y0, i - circ1_x0);
			fi = dec_to_pol_fi(j - circ1_y0, i - circ1_x0);
			data[j][i] = 0;
			if (in_circle_pol(rr, fi, R, r, circ1_x0, circ1_y0) == 2)
			{
				data[j][i] = border(j / double(size), i / double(size));
			}
			rr = dec_to_pol_r(j - circ2_y0, i - circ2_x0);
			fi = dec_to_pol_fi(j - circ2_y0, i - circ2_x0);
			if (in_circle_pol(rr, fi, R, r, circ2_x0, circ2_y0) == 2)
			{
				data[j][i] = border(j / double(size), i / double(size));
			}
			rr = dec_to_pol_r(j - circ3_y0, i - circ3_x0);
			fi = dec_to_pol_fi(j - circ3_y0, i - circ3_x0);
			if (in_circle_pol(rr, fi, R, r, circ3_x0, circ3_y0) == 2)
			{
				data[j][i] = border(j / double(size), i / double(size));
			}
			rr = dec_to_pol_r(j - circ4_y0, i - circ4_x0);
			fi = dec_to_pol_fi(j - circ4_y0, i - circ4_x0);
			if (in_circle_pol(rr, fi, R, r, circ4_x0, circ4_y0) == 2)
			{
				data[j][i] = border(j / double(size), i / double(size));
			}
			if (in_rect(i, j, rect1_x1, rect1_x2, rect1_y1, rect1_y2) == 2)
			{
				data[j][i] = border(j / double(size), i / double(size));
			}
			if (in_rect(i, j, rect2_x1, rect2_x2, rect2_y1, rect2_y2) == 2)
			{
				data[j][i] = border(j / double(size), i / double(size));
			}
			if (in_rect(i, j, rect3_x1, rect3_x2, rect3_y1, rect3_y2) == 2)
			{
				data[j][i] = border(j / double(size), i / double(size));
			}
			if (in_rect(i, j, rect4_x1, rect4_x2, rect4_y1, rect4_y2) == 2)
			{
				data[j][i] = border(j / double(size), i / double(size));
			}
		}
	}
}

