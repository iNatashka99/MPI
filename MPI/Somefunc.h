#pragma once

// функция u
double U(double x, double y);

// функция f
double F(double x, double y);

// Граница
double border(double x, double y);

// Начальное приближение
double startfill(double x, double y);

// Распределение кусочков по процессорам
int ** MakeShapes(int num_procs);

// Есть ли пересечение
bool have_intersection(int shape1, int shape2);

// В каком процесоре кусок
int in_proc(int shape, int num_procs, int ** procs_shapes);