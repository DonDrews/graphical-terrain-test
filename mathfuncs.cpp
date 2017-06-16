#include <iostream>
#include <math.h>

//helper function
int coord(int x, int y, int size)
{
  return x + (y * size);
}

//helper function to bicubicInterpolate()
//makes a new array with a 1 unit border of approximated outside values
void adjustArray(float* graph, int size, float* newArray)
{
  int nSize = size + 2;
  //fill corners
  newArray[coord(0, 0, nSize)] = 2 * graph[coord(0, 0, size)] - graph[coord(1, 1, size)];
  newArray[coord(0, nSize - 1, nSize)] = 2 * graph[coord(0, size - 1, size)] - graph[coord(1, size - 2, size)];
  newArray[coord(nSize - 1, 0, nSize)] = 2 * graph[coord(size - 1, 0, size)] - graph[coord(size - 2, 1, size)];
  newArray[coord(nSize - 1, nSize - 1, nSize)] = 2 * graph[coord(size - 1, size - 1, size)] - graph[coord(size - 2, size - 2, size)];

  //fill sides
  for(int i = 1; i < size + 1; i++)
  {
    //LEFT
    newArray[coord(0, i, nSize)] = 2 * graph[coord(0, i - 1, size)] - graph[coord(1, i - 1, size)];

    //RIGHT
    newArray[coord(size + 1, i, nSize)] = 2 * graph[coord(size - 1, i - 1, size)] - graph[coord(size - 2, i - 1, size)];

    //BOTTOM
    newArray[coord(i, 0, nSize)] = 2 * graph[coord(i - 1, 0, size)] - graph[coord(i - 1, 1, size)];

    //TOP
    newArray[coord(i, size + 1, nSize)] = 2 * graph[coord(i - 1, size - 1, size)] - graph[coord(i - 1, size - 2, size)];
  }

  //fill center
  for(int x = 1; x < size + 1; x++)
  {
    for(int y = 1; y < size + 1; y++)
    {
      newArray[coord(x, y, nSize)] = graph[coord(x - 1, y - 1, size)];
    }
  }
}

void printArray(float* matrix, int size)
{
  for(int i = 0; i < size; i++)
  {
    for(int j = 0; j < size; j++)
    {
      std::cout << matrix[i * size + j] << " ";
    }
    std::cout << std::endl;
  }
}

//returns p(x) where x is the normalized position between y1 and y2
//thanks to Paul Bourke (the following code is taken from his website)
float cubicInterpolate(float x, float y0, float y1, float y2, float y3)
{
  float a0,a1,a2,a3,mu2;

  mu2 = x*x;
  a0 = y3 - y2 - y0 + y1;
  a1 = y0 - y1 - a0;
  a2 = y2 - y0;
  a3 = y1;

  return (a0*x*mu2+a1*mu2+a2*x+a3);
}

void bicubicInterpolate(float* original, int originalSize, float* smoothed, int size)
{
  //this is the change factor
  //when a coordinate for the new grid is multiplied by this coefficient, it is converted to an old coordinate
  float cF = (float(originalSize) - 1) / (float(size) - 1);

  //make array with approximations
  float adjusted[originalSize + 2][originalSize + 2];
  adjustArray(&original[0], originalSize, &adjusted[0][0]);
  printArray(&adjusted[0][0], originalSize + 2);

  //iterate through each index on the scaled array
  for(int y = 0; y < size; y++)
  {
    for(int x = 0; x < size; x++)
    {
      //find nearest old coordinates that correspond to these (rounded down)
      int xDown = floor(x * cF + 1);
      if(xDown == originalSize)
        xDown = originalSize - 1;

      int yDown = floor(y * cF + 1);
      if(yDown == originalSize)
        yDown = originalSize - 1;

      float interps[4];
      int index = 0;

      //for each new index, four parallel cubic interpolations are done
      //coordinate should lie between second and third interp
      for(int i = yDown - 1; i <= yDown + 2; i++)
      {
        float pVals[4];
        pVals[0] = adjusted[i][xDown - 1];
        pVals[1] = adjusted[i][xDown];
        pVals[2] = adjusted[i][xDown + 1];
        pVals[3] = adjusted[i][xDown + 2];

        interps[index] = cubicInterpolate((x * cF + 1) - xDown, pVals[0], pVals[1], pVals[2], pVals[3]);
        ++index;
      }

      //interpolate between previous interpolations
      smoothed[coord(x, y, size)] = cubicInterpolate((y * cF + 1) - yDown, interps[0], interps[1], interps[2], interps[3]);
    }
  }
}

//temporary
int main()
{
  float org[3][3] = {2, 0, 3, 0, 1, 0, 0, 0, 1};
  printArray(&org[0][0], 3);
  std::cout << std::endl;
  float parsed[9][9];
  bicubicInterpolate(&org[0][0], 3, &parsed[0][0], 9);
  printArray(&parsed[0][0], 9);
  return 1;
}
