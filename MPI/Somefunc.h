#pragma once

// ������� u
double U(double x, double y);

// ������� f
double F(double x, double y);

// �������
double border(double x, double y);

// ��������� �����������
double startfill(double x, double y);

// ������������� �������� �� �����������
int ** MakeShapes(int num_procs);

// ���� �� �����������
bool have_intersection(int shape1, int shape2);

// � ����� ��������� �����
int in_proc(int shape, int num_procs, int ** procs_shapes);