#include "fractal.h"
#include "mathfuncs.h"
#include "Erosion.h"
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <iostream>
#include <vector>
#include <fstream>
#include <string>

using namespace std;

void writeImage(string name, float* data, int size)
{
  ofstream image;
  image.open(name);
  //write header
  image << "P3\n";
  image << size << " " << size << "\n";
  image << "255\n";

  for(int i = 0; i < size; i++)
  {
    for(int j = 0; j < size; j++)
    {
      /*
      unsigned short color = data[coord(i, j, size)] * 65535;
      if(color < 0)
        color = 0;
      if(color > 65535)
        color = 65535;

      unsigned char msb = color & 0xFF;
      unsigned char lsb = (color >> 8) & 0xFF;

      image << to_string(lsb) + " " + to_string(msb) + " " + to_string(0) + "  ";
      */
      unsigned short color = data[coord(i, j, size)] * 255;
      if(color < 0)
        color = 0;
      if(color > 255)
        color = 255;

      image << to_string(color) << " " << to_string(color) << " " << to_string(color) << " ";
    }
    image << "\n";
  }

  image.close();
}

int main()
{
  const int SIZE = 513;

/*

  float* startPerlin = new float[6 * 6];

  //seed random gen
  srand(time(NULL));

  for(int i = 0; i < 6 * 6; i++)
  {
    startPerlin[i] = randomRange(0.0, 1.0);
  }
  float* finishedPerlin = new float[SIZE * SIZE];
  bicubicInterpolate(&startPerlin[0], 6, &finishedPerlin[0], SIZE);

  float* startPerlinSmall = new float[24 * 24];

  for(int i = 0; i < 24 * 24; i++)
  {
    startPerlin[i] = randomRange(-0.15, 0.15);
  }
  float* finishedPerlinSmall = new float[SIZE * SIZE];
  bicubicInterpolate(&startPerlinSmall[0], 24, &finishedPerlinSmall[0], SIZE);

  //combine
  float* finished = new float[SIZE * SIZE];
  for(int x = 0; x < SIZE; x++)
  {
    for(int y = 0; y < SIZE; y++)
    {
      float coefficient = (finishedPerlin[y * SIZE + x] - 0.1) * 1.5;
      if(coefficient < 0)
        coefficient = 0;

      finished[y * SIZE + x] = finishedPerlin[y * SIZE + x] + finishedPerlinSmall[y * SIZE + x] + finishedFractal[coord(x, y, SIZE)] * coefficient;
    }
  }
  */
  float* startFractal = new float[4];

  for(int i = 0; i < 2 * 2; i++)
  {
    startFractal[i] = randomRange(0.5, 0.7);
  }

  float temp = 0;
  float* finishedFractal = &temp;

  makeFractalArray(&startFractal[0], 2, finishedFractal, SIZE, 9);

  writeImage("preerode.ppm", &finishedFractal[0], SIZE);

  float* eroded = erodeField(&finishedFractal[0], SIZE);

  writeImage("bigfinal.ppm", &eroded[0], SIZE);
}
