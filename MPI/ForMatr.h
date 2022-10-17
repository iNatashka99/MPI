#pragma once

//����� ������� A
void printMatr(double **data, int size);

//����� �����
void printShape(double **data, int size1, int size2);

//������� ������� A
void clearMatr(double **data, int size);

//����������� �������
double ** copymatr(double **A1, int size);

//��������� �����������
void start_fill_matr(double **data, int size, double R, double r,
	int circ1_x0, int circ1_y0,
	int circ2_x0, int circ2_y0,
	int circ3_x0, int circ3_y0,
	int circ4_x0, int circ4_y0,
	int rect1_x1, int rect1_x2, int rect1_y1, int rect1_y2,
	int rect2_x1, int rect2_x2, int rect2_y1, int rect2_y2,
	int rect3_x1, int rect3_x2, int rect3_y1, int rect3_y2,
	int rect4_x1, int rect4_x2, int rect4_y1, int rect4_y2);

//��������� ��������� �� 8 ������
void start_shape_matr(double **data, int size, double R, double r,
	int circ1_x0, int circ1_y0,
	int circ2_x0, int circ2_y0,
	int circ3_x0, int circ3_y0,
	int circ4_x0, int circ4_y0,
	int rect1_x1, int rect1_x2, int rect1_y1, int rect1_y2,
	int rect2_x1, int rect2_x2, int rect2_y1, int rect2_y2,
	int rect3_x1, int rect3_x2, int rect3_y1, int rect3_y2,
	int rect4_x1, int rect4_x2, int rect4_y1, int rect4_y2);