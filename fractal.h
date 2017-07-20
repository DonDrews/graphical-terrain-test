#ifndef FRACTAL_H
#define FRACTAL_H

bool isValid(int x, int y, int size);
float randomRange(float start, float end);
void makeFractalArray(float* starting, int startSize, float* &finished, int finishSize, int iterations);

#endif
