#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include "mathfuncs.h"

using namespace std;

const float START_HARMONIC = 0.5;

bool isValid(int x, int y, int size)
{
  return (x > -1 && x < size && y > -1 && y < size);
}

float randomRange(float start, float end)
{
  float range = end - start;
  return start + ((float)rand() / (float)RAND_MAX) * range;
}

void makeFractalArray(float* starting, int startSize, float* &finished, int finishSize, int iterations)
{
  //seed random gen
  srand(time(NULL));

  float* currentArray = starting;
  float harmonic = START_HARMONIC;

  for(int i = 0; i < iterations; i++)
  {
    int currentSize = int(pow(2, log2(startSize - 1) + i) + 1);
    int newSize = currentSize + currentSize - 1;

    float* newArray = (float*)malloc(newSize * newSize * sizeof(float));

    //diamond and copy step
    for(int y = 0; y < newSize; y++)
    {
      for(int x = 0; x < newSize; x++)
      {
        if(x & 1 && y & 1) //if x and y are both odd, diamond
        {
          int xDown = floor(x / 2);
          int yDown = floor(y / 2);

          float rands[4];
          rands[0] = currentArray[coord(xDown, yDown, currentSize)];
          rands[1] = currentArray[coord(xDown + 1, yDown, currentSize)];
          rands[2] = currentArray[coord(xDown, yDown + 1, currentSize)];
          rands[3] = currentArray[coord(xDown + 1, yDown + 1, currentSize)];

          float average = (rands[0] + rands[1] + rands[2] + rands[3]) / 4;

          newArray[coord(x, y, newSize)] = average + randomRange(-harmonic, harmonic);
          //cout << randomRange(-harmonic, harmonic);

          //cout << "Diamond: " << x << " " << y << endl;
        }
        else if(!(x & 1) && !(y & 1)) //if x and y are both even, copy
        {
          newArray[coord(x, y, newSize)] = currentArray[coord(x / 2, y / 2, currentSize)];
          //cout << "Copy: " << x << " " << y << endl;
        }
      }
    }

    //printArray(newArray, newSize);

    //square step
    for(int y = 0; y < newSize; y++)
    {
      for(int x = 0; x < newSize; x++)
      {
        if(!(x & 1) != !(y & 1)) //x or y is odd
        {

          float summation = 0;
          int amount = 0;

          if(isValid(x - 1, y, newSize))
          {
            summation += newArray[coord(x - 1, y, newSize)];
            ++amount;
          }
          if(isValid(x + 1, y, newSize))
          {
            summation += newArray[coord(x + 1, y, newSize)];
            ++amount;
          }
          if(isValid(x, y - 1, newSize))
          {
            summation += newArray[coord(x, y - 1, newSize)];
            ++amount;
          }
          if(isValid(x, y + 1, newSize))
          {
            summation += newArray[coord(x, y + 1, newSize)];
            ++amount;
          }

          newArray[coord(x, y, newSize)] = (summation / amount) + randomRange(-harmonic, harmonic);

          //cout << "Square: " << x << " " << y << endl;
        }
      }
    }

    //printArray(newArray, newSize);

    if(currentArray != starting)
    {
      delete[] currentArray;
    }
    currentArray = newArray;
    harmonic *= 0.5;
  }

  finished = currentArray;
}
