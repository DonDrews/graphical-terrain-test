#include <algorithm>
#include <iostream>
#include <vector>
#include <cmath>

#include "mathfuncs.h"
#include "Erosion.h"

const int ITERATIONS = 1000;
const float TIME_STEP = 0.0002;
const float RAINDROP_SIZE = 0.1;
const float RAIN_PROB = 0.05;
const float PIPE_CROSS_SECTION = 0.05;
const float GRAVITY = 0.05;
const float PIPE_LENGTH = 0.001;
const float TILT_MIN = 0.0000;
const float SEDIMENT_CAP = 0.07;
const float DISSOLVE_COEFF = 0.0002;
const float DEP_COEFF = 0.0008;
const float EVAP_COEFF = 0.001;

struct Crd
{
  Crd(){x = 0; y = 0;}
  Crd(int i, int j){x = i; y = j;}

  int x;
  int y;
};

struct Cell
{
  float b;
  float b1;
  float prevD;
  float d;
  float d1;
  float d2;
  float s;
  float s1;
  float f[4];
  float u;
  float v;
};

const int LEFT = 0;
const int RIGHT = 1;
const int TOP = 2;
const int BOTTOM = 3;

void checkStability(Cell& c)
{
  bool unDef = false;

  unDef = (unDef || std::isnan(c.b) || std::isinf(c.b));
  unDef = (unDef || std::isnan(c.prevD) || std::isinf(c.prevD));
  unDef = (unDef || std::isnan(c.d) || std::isinf(c.d));
  unDef = (unDef || std::isnan(c.s) || std::isinf(c.s));
  unDef = (unDef || std::isnan(c.f[0]) || std::isinf(c.f[0]));
  unDef = (unDef || std::isnan(c.f[1]) || std::isinf(c.f[1]));
  unDef = (unDef || std::isnan(c.f[2]) || std::isinf(c.f[2]));
  unDef = (unDef || std::isnan(c.f[3]) || std::isinf(c.f[3]));
  unDef = (unDef || std::isnan(c.u) || std::isinf(c.u));
  unDef = (unDef || std::isnan(c.v) || std::isinf(c.v));
  unDef = (unDef || abs(c.f[0]) > 10000);
  unDef = (unDef || abs(c.f[1]) > 10000);
  unDef = (unDef || abs(c.f[2]) > 10000);
  unDef = (unDef || abs(c.f[3]) > 10000);
  unDef = (unDef || c.d < 0.0);
  unDef = (unDef || c.b < -1);

  if(unDef)
  {
    std::cout << "Listing Diagnostic:" << std::endl;
    std::cout << c.b << std::endl;
    std::cout << c.d1 << std::endl;
    std::cout << c.d2 << std::endl;
    std::cout << c.d << std::endl;
    std::cout << c.s << std::endl;
    std::cout << c.f[0] << std::endl;
    std::cout << c.f[1] << std::endl;
    std::cout << c.f[2] << std::endl;
    std::cout << c.f[3] << std::endl;
    std::cout << c.u << std::endl;
    std::cout << c.v << std::endl;

    exit(-1);
  }
}

Crd getCoordAtDir(Crd c, int dir, int size, bool& isNull)
{
  Crd val;
  isNull = false;

  if(dir == LEFT)
    val = Crd(c.x - 1, c.y);
  else if(dir == RIGHT)
    val = Crd(c.x + 1, c.y);
  else if(dir == TOP)
    val = Crd(c.x, c.y + 1);
  else if(dir == BOTTOM)
    val = Crd(c.x, c.y - 1);
  else
    isNull = true;

  //check for out of bounds
  if(val.x < 0 || val.x >= size || val.y < 0 || val.y >= size)
	 isNull = true;
  else
    return val;
}

int getOppositeDirection(int dir)
{
  if(dir == LEFT)
    return RIGHT;
  else if(dir == RIGHT)
    return LEFT;
  else if(dir == TOP)
    return BOTTOM;
  else if(dir == BOTTOM)
    return TOP;
  else
    return -1;
}

float getInterpValue(float ll, float lr, float ul, float ur, float x, float y)
{
  float lLerp = x * (lr - ll) + ll;
  float uLerp = x * (ur - ul) + ul;
  return y * (uLerp - lLerp) + lLerp;
}

float* erodeField(float* field, float*& water, int size)
{
  Cell defaultCell;
  defaultCell.b = 0.0f;
  defaultCell.b1 = 0.0f;
  defaultCell.d = 0.0f;
  defaultCell.d1 = 0.0f;
  defaultCell.d2 = 0.0f;
  defaultCell.s = 0.0f;
  defaultCell.s1 = 0.0f;
  defaultCell.f[0] = {0.0f};
  defaultCell.f[1] = { 0.0f };
  defaultCell.f[2] = { 0.0f };
  defaultCell.f[3] = { 0.0f };
  defaultCell.u = 0.0f;
  defaultCell.v = 0.0f;

  //create 2D map of cells
  std::vector<std::vector<Cell>> sim(size, std::vector<Cell>(size, defaultCell));

  //Set terrain height to values stored in field
  for(int x = 0; x < size; x++)
  {
    for(int y = 0; y < size; y++)
    {
      sim[x][y].b = field[coord(x, y, size)];
    }
  }

  //main loop
  for(int i = 0; i < ITERATIONS; i++)
  {
    //std::cout << "Iteration " << i << std::endl;

    //Step 1: Add water through rainfall
    for(int x = 0; x < size; x++)
    {
      for(int y = 0; y < size; y++)
      {
        if(rand() % int(size * size * RAIN_PROB) == 0)
        {
          sim[x][y].d1 = sim[x][y].d + RAINDROP_SIZE;
        }
        else
        {
          sim[x][y].d1 = sim[x][y].d;
        }
      }
    }

    //Step 2: Calculate movement of water
    for(int x = 0; x < size; x++)
    {
      for(int y = 0; y < size; y++)
      {
        float fluxSum = 0.0;
        //for each direction
        for(int j = 0; j < 4; j++)
        {
          //get adjacent cell at direction
    			bool n;
    			Crd side = getCoordAtDir(Crd(x, y), j, size, n);

          if(!n)
          {
            //calculate height difference (including water)
            float deltaHeight = (sim[x][y].b + sim[x][y].d1) - (sim[side.x][side.y].b + sim[side.x][side.y].d1);

            //find new flux value for direction
            sim[x][y].f[j] = std::max(0.0f, sim[x][y].f[j] + (TIME_STEP * PIPE_CROSS_SECTION * GRAVITY * deltaHeight) / PIPE_LENGTH);
          }
          else
          {
            sim[x][y].f[j] = 0.0;
          }

          fluxSum += sim[x][y].f[j];
        }

        float scalingFactor = 1.0f;
        if(fluxSum > 0.000001)
          scalingFactor = (PIPE_LENGTH * PIPE_LENGTH) / (fluxSum * TIME_STEP);

        scalingFactor = std::max(std::min(1.0f, scalingFactor * sim[x][y].d1), 0.0f);

        //for each direction
        for(int j = 0; j < 4; j++)
        {
          //adjust based on scaling factor
          sim[x][y].f[j] *= scalingFactor;
        }

        checkStability(sim[x][y]);
      }
    }

    //Step 3: Apply calculated flux amounts
    for(int x = 0; x < size; x++)
    {
      for(int y = 0; y < size; y++)
      {
        float totalDelta = 0.0;

        for(int j = 0; j < 4; j++)
        {
    			bool n;
    			Crd side = getCoordAtDir(Crd(x, y), j, size, n);
          if(!n)
          {
            //inflow
            totalDelta += sim[side.x][side.y].f[getOppositeDirection(j)];
            //outflow
            totalDelta -= sim[x][y].f[j];
          }
        }

        sim[x][y].d2 = sim[x][y].d1 + ((TIME_STEP * totalDelta) / (PIPE_LENGTH * PIPE_LENGTH));

        if(sim[x][y].d2 < 0.000001)
        {
          if(sim[x][y].d2 < 0.0)
            sim[x][y].d2 = 0.0;
        }

        checkStability(sim[x][y]);
      }
    }

    //Step 4: Adjust velocity field
    for(int x = 0; x < size; x++)
    {
      for(int y = 0; y < size; y++)
      {
        float avgWater = (sim[x][y].d2 + sim[x][y].d1) / 2.0f;

        //horizontal side
        float lContrib = 0.0;
    		bool l;
    		Crd lSide = getCoordAtDir(Crd(x, y), LEFT, size, l);
        if(!l)
        {
          lContrib = sim[lSide.x][lSide.y].f[RIGHT] - sim[x][y].f[LEFT];
        }

        float rContrib = 0.0;
    		bool r;
    		Crd rSide = getCoordAtDir(Crd(x, y), RIGHT, size, r);
        if(!r)
        {
          rContrib = sim[x][y].f[RIGHT] - sim[rSide.x][rSide.y].f[LEFT];
        }

        if(avgWater > 0.0000001)
        {
          sim[x][y].u = ((lContrib + rContrib) / 2.0) / (avgWater * PIPE_LENGTH);
        }
        else
        {
          sim[x][y].u = 0.0;
        }

        //vertical side
        float bContrib = 0.0;
    		bool b;
    		Crd bSide = getCoordAtDir(Crd(x, y), BOTTOM, size, b);
        if(!b)
        {
          bContrib = sim[bSide.x][bSide.y].f[TOP] - sim[x][y].f[BOTTOM];
        }

        float tContrib = 0.0;
    		bool t;
    		Crd tSide = getCoordAtDir(Crd(x, y), TOP, size, t);
        if(!t)
        {
          tContrib = sim[x][y].f[TOP] - sim[tSide.x][tSide.y].f[BOTTOM];
        }

        if(avgWater > 0.0000001)
          sim[x][y].v = ((bContrib + tContrib) / 2.0) / (avgWater * PIPE_LENGTH);
        else
          sim[x][y].v = 0.0;

        checkStability(sim[x][y]);
      }
    }

    //Step 5: Erode and Deposit
    for(int x = 0; x < size; x++)
    {
      for(int y = 0; y < size; y++)
      {
        //find tilt angle

        //horizontal side
        float hHeightDelta = 0.0;
        int index = 0;
    		bool l;
    		Crd lSide = getCoordAtDir(Crd(x, y), LEFT, size, l);
        if(!l)
        {
          hHeightDelta += sim[x][y].b - sim[lSide.x][lSide.y].b;
          ++index;
        }

    		bool r;
    		Crd rSide = getCoordAtDir(Crd(x, y), RIGHT, size, r);
        if(!r)
        {
          hHeightDelta += sim[rSide.x][rSide.y].b - sim[x][y].b;
          ++index;
        }

        //adjust for boundary condition
        if(index != 2)
          hHeightDelta *= 2;

        //vertical side
        float vHeightDelta = 0.0;
        index = 0;
    		bool b;
    		Crd bSide = getCoordAtDir(Crd(x, y), BOTTOM, size, b);
        if(!b)
        {
          vHeightDelta += sim[x][y].b - sim[bSide.x][bSide.y].b;
          ++index;
        }

    		bool t;
    		Crd tSide = getCoordAtDir(Crd(x, y), TOP, size, t);
        if(!t)
        {
          vHeightDelta += sim[tSide.x][tSide.y].b - sim[x][y].b;
          ++index;
        }

        //adjust for boundary condition
        if(index != 2)
          vHeightDelta *= 2;

        //find normal (and normalize)
        float normal[3] = {hHeightDelta, PIPE_LENGTH, vHeightDelta};
        float magnitude = sqrt(pow(normal[0], 2) + pow(normal[1], 2) + pow(normal[2], 2));
        normal[0] /= magnitude;
        normal[1] /= magnitude;
        normal[2] /= magnitude;

        float sinOfAngle = std::max(TILT_MIN, float(sqrt(1.0 - pow(normal[1], 2))));

        float velMagnitude = sqrt(pow(sim[x][y].u, 2) + pow(sim[x][y].v, 2));

        float transCapacity = SEDIMENT_CAP * sinOfAngle * velMagnitude;

        if(transCapacity > sim[x][y].s)
        {
          //erode
			    float sedChange = DISSOLVE_COEFF * (transCapacity - sim[x][y].s);

          sim[x][y].b1 = std::max(0.0f, sim[x][y].b - sedChange);
          sim[x][y].s1 = sim[x][y].s + sedChange;
        }
        else
        {
          //deposit
          float sedChange = DEP_COEFF * (sim[x][y].s - transCapacity);

          sim[x][y].b1 = sim[x][y].b + sedChange;
          sim[x][y].s1 = std::max(0.0f, sim[x][y].s - sedChange);
        }

        checkStability(sim[x][y]);
      }
    }

    //Step 6: Transport Sediment
    for(int x = 0; x < size; x++)
    {
      for(int y = 0; y < size; y++)
      {
        float xSed = x - (sim[x][y].u * TIME_STEP);
        float ySed = y - (sim[x][y].v * TIME_STEP);
        int xDown = floor(xSed);
        int yDown = floor(ySed);

        if(xDown >= size - 1 || xDown < 0 || yDown >= size - 1 || yDown < 0)
        {
          //do not move sediment
        }
        else
        {
          sim[x][y].s = getInterpValue(sim[xDown][yDown].s1, sim[xDown + 1][yDown].s1, sim[xDown][yDown + 1].s1, sim[xDown + 1][yDown + 1].s1, xSed - xDown, ySed - yDown);
        }
      }
    }

    //Step 7: Evaporate Water
    for(int x = 0; x < size; x++)
    {
      for(int y = 0; y < size; y++)
      {
        sim[x][y].d *= 1 - (EVAP_COEFF * TIME_STEP);
      }
    }

    //Step 8: Move all changes back to center
    for(int x = 0; x < size; x++)
    {
      for(int y = 0; y < size; y++)
      {
        sim[x][y].b = sim[x][y].b1;
        sim[x][y].d = sim[x][y].d2;
        /*
        if(x == 100 && y == 100)
        {
          std::cout << "WAT- " << sim[x][y].d << std::endl;
          std::cout << "Left- " << sim[x][y].f[LEFT] << std::endl;
          std::cout << "leftIn- " << sim[x - 1][y].f[RIGHT] << std::endl;
          std::cout << "Right- " << sim[x][y].f[1] << std::endl;
          std::cout << "Top- " << sim[x][y].f[2] << std::endl;
          std::cout << "Bot- " << sim[x][y].f[3] << std::endl;
        }
        */
      }
    }
  }

  //convert water
  water = new float[size * size];
  for(int x = 0; x < size; x++)
  {
    for(int y = 0; y < size; y++)
    {
      water[coord(x, y, size)] = sim[x][y].d;
    }
  }

  //Convert back to float array
  float* eroded = new float[size * size];
  for(int x = 0; x < size; x++)
  {
    for(int y = 0; y < size; y++)
    {
      eroded[coord(x, y, size)] = sim[x][y].b;
    }
  }

  return eroded;
}
