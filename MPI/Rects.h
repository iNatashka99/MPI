#pragma once

// Принадлежность прямоугольнику
int in_rect(int x, int y, int x1, int x2, int y1, int y2);

// Отправка прямоугольника
void send_rect(double **data, int x1, int x2, int y1, int y2, int proc, int shape);

// Прием прямоугольника
void recive_rect(double **data, int x1, int x2, int y1, int y2, int proc, int shape);

// Запись для нулевого процессора
void set_rect(double **data, double **layer, int x1, int x2, int y1, int y2);

// Норма в прямоугольнике
double norma_rect(double **layer1, double **layer2, int x1, int x2, int y1, int y2);

