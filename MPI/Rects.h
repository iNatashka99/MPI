#pragma once

// �������������� ��������������
int in_rect(int x, int y, int x1, int x2, int y1, int y2);

// �������� ��������������
void send_rect(double **data, int x1, int x2, int y1, int y2, int proc, int shape);

// ����� ��������������
void recive_rect(double **data, int x1, int x2, int y1, int y2, int proc, int shape);

// ������ ��� �������� ����������
void set_rect(double **data, double **layer, int x1, int x2, int y1, int y2);

// ����� � ��������������
double norma_rect(double **layer1, double **layer2, int x1, int x2, int y1, int y2);

