#define _CRT_SECURE_NO_WARNINGS


#include <iostream>
#include <iomanip>
#include <mpi.h>
#include <string>
#include "math.h"
#include <stdio.h>

#include "ForMatr.h"
#include "Circles.h"
#include "Rects.h"
#include "SSOR.h"
#include "MPI_funcsr.h"

using namespace std;

const int Tag = 0;
const int root = 0;


int main()
{
	int size = 120;

	double w = 0.5;
	double eps = 0.001*pow(size, 2);
	double h = 1.0 / size;

	size = size + 2;
	double R = (size - 2) / 12.0;
	double r = R / 3.0;

	// circle 1
	int circ1_x0 = 1 + R;
	int circ1_y0 = size / 2.0;

	// circle 2
	int circ2_x0 = 1 + 4 * R;
	int circ2_y0 = size / 2.0;

	// circle 3
	int circ3_x0 = 1 + 7 * R;
	int circ3_y0 = size / 2.0;

	// circle 4
	int circ4_x0 = 1 + 10 * R;
	int circ4_y0 = size / 2.0;

	int delt = 2;

	// rect 1
	int rect1_x1 = circ1_x0 - (r + (R-r)/3);
	int rect1_x2 = circ2_x0 - (R + delt);
	int rect1_y1 = circ2_y0 - (r + (R - r) / 3);
	int rect1_y2 = circ2_y0 + (r + (R - r) / 3);

	// rect 2
	int rect2_x1 = rect1_x2 - delt;
	int rect2_x2 = circ2_x0 + (circ3_x0 - circ2_x0)/2 + delt;
	int rect2_y1 = circ2_y0 - (r + (R - r) / 3);
	int rect2_y2 = circ2_y0 + (r + (R - r) / 3);

	// rect 3
	int rect3_x1 = rect2_x2  -  2* delt;
	int rect3_x2 = circ3_x0 + (R + 2* delt);
	int rect3_y1 = circ3_y0 - (r + (R - r) / 3);
	int rect3_y2 = circ3_y0 + (r + (R - r) / 3);

	// rect 4
	int rect4_x1 = rect3_x2 - delt;
	int rect4_x2 = circ4_x0 + (r + (R - r) / 3);
	int rect4_y1 = circ4_y0 - (r + (R - r) / 3);
	int rect4_y2 = circ4_y0 + (r + (R - r) / 3);



	double **matr;
	matr = new double *[size];
	for (int i = 0; i < size; i++)
	{
		matr[i] = new double[size];
	}

	
	SSOR_paral(matr, size, R, r, h, w, eps,
		circ1_x0, circ1_y0,
		circ2_x0, circ2_y0,
		circ3_x0, circ3_y0,
		circ4_x0, circ4_y0,
		rect1_x1, rect1_x2, rect1_y1, rect1_y2,
		rect2_x1, rect2_x2, rect2_y1, rect2_y2,
		rect3_x1, rect3_x2, rect3_y1, rect3_y2,
		rect4_x1, rect4_x2, rect4_y1, rect4_y2);
	
	/*
	double **Proba, **Proba2;
	Proba = new double * [11];

	double step = dec_to_pol_fi(4.0, 1.0);
	int N_fi = 360 / step;
	Proba2 = new double *[N_fi];

	for (int i = 0; i < 11; i++)
	{
		Proba[i] = new double [11];
		for (int j = 0; j < 11; j++)
		{
			Proba[i][j] = 0;
			int rr = dec_to_pol_r(i - 5, j - 5);
			int fi = dec_to_pol_fi(i - 5, j - 5);
			fi = fi / 360 * step;
			if (in_circle_pol(rr, fi, 4, 1, 5, 5)==1)
			{
				Proba[i][j] = 1;
			}
			if (in_circle_pol(rr, fi, 4, 1, 5, 5) == 2)
			{
				Proba[i][j] = 2;
			}
		}
	}
	printMatr(Proba, 11);

	for (int i = 0; i < 5; i++)
	{
		Proba2[i] = new double[N_fi];
		for (int j = 0; j < N_fi; j++)
		{
			int x = 5 - pol_to_dec_x(i, j);
			int y = 5 - pol_to_dec_y(i, j);
			Proba2[i][j] = 0;
			Proba2[i][j] = Proba[x][y];
		}
	}
	for (int i = 5; i < N_fi; i++)
	{
		Proba2[i] = new double[N_fi];
		for (int j = 0; j < N_fi; j++)
		{
			Proba2[i][j] = 0;
		}
	}
	//printMatr(Proba2, N_fi);
	clearMatr(Proba, 11);

	for (int i = 0; i < 11; i++)
	{
		for (int j = 0; j < 11; j++)
		{
			int rr = dec_to_pol_r(i - 5, j - 5);
			int fi = dec_to_pol_fi(i - 5, j - 5);
			fi = fi / 360 * step;
			if (in_circle_pol(rr, fi, 4, 1, 5, 5))
			{
				Proba[i][j] = Proba2[rr][fi];
			}
		}
	}
	printMatr(Proba, 11);
	*/

	return 0;
}