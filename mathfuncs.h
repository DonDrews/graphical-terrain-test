#ifndef MATHFUNCS_H
#define MATHFUNCS_H

int coord(int x, int y, int size);
void adjustArray(float* graph, int size, float* newArray);
void printArray(float* matrix, int size);
float cubicInterpolate(float x, float y0, float y1, float y2, float y3);
void bicubicInterpolate(float* original, int originalSize, float* smoothed, int size);

#endif
