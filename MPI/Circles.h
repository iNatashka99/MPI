#pragma once

// ��������� ����������
int circle_x(int x0, double R, double fi);
int circle_y(int y0, double R, double fi);

// �������������� ������
int in_circle_dec(int x, int y, double R, double r, int x0, int y0);

// �������������� ������ � ��������
int in_circle_pol(int rr, int fi, int R, int r, int x0, int y0);


//������� �� ���������� � ��������
double dec_to_pol_r(double x, double y);
double dec_to_pol_fi(double x, double y);

//������� �� �������� � ����������
double pol_to_dec_x(double r, double fi);
double pol_to_dec_y(double r, double fi);

// �������� ������
void send_circle(double **data, int x0, int y0, double R, double r, int proc, int shape);

// ����� ������
void recive_circle(double **data, int x0, int y0, double R, double r, int proc, int shape);

// ������ ��� �������� ����������
void set_circle(double **data, double **layer, int x0, int y0, double R, double r);

// ����� � ������
double norma_circle(double **layer1, double **layer2, int x0, int y0, double R, double r);