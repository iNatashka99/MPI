#pragma once
#include "SSOR.h"
#include "Circles.h"
#include <math.h>

#include "SSOR.h"
#include "Circles.h"
#include "Rects.h"
#include <math.h>


// ���� � ����
double gaus(double **old_layer, double **new_layer, int i, int j, double h);
// ���� � ������ 
double ** GAUS_in_circle_dec1(double **old_layer, int size, int x0, int y0, double R, double r, double h);
double ** GAUS_in_circle_dec2(double **old_layer, int size, int x0, int y0, double R, double r, double h);

// ��������� � ������
double ** corr_in_circle_dec(double **old_layer, double **new_layer, int size, int x0, int y0, double R, double r, double w);

// SSOR � ������ 
double ** SSOR_in_circle_dec(double **old_layer, int size, int x0, int y0, double R, double r, double h, double w);

// ���� � �������������� 
double ** GAUS_in_rect1(double **old_layer, int size, int x1, int x2, int y1, int y2, double h);
double ** GAUS_in_rect2(double **old_layer, int size, int x1, int x2, int y1, int y2, double h);

// ��������� � ��������������
double ** corr_in_rect(double **old_layer, double **new_layer, int size, int x1, int x2, int y1, int y2, double w);

// SSOR � �������������� 
double ** SSOR_in_rect(double **old_layer, int size, int x1, int x2, int y1, int y2, double h, double w);