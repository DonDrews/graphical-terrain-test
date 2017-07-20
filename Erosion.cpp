#include <algorithm>
#include <iostream>
#include <vector>
#include <cmath>

#include "mathfuncs.h"
#include "Erosion.h"

const int ITERATIONS = 50;
const float TIME_STEP = 0.0001;
const float RAINDROP_SIZE = 0.05;
const float PIPE_CROSS_SECTION = 0.1;
const float GRAVITY = 0.1;
const float PIPE_LENGTH = 0.2;
const float TILT_MIN = 0.001;
const float SEDIMENT_CAP = 0.001;
const float DISSOLVE_COEFF = 0.1;
const float DEP_COEFF = 0.1;
const float EVAP_COEFF = 0.1;

struct Crd
{
  Crd(){x = 0; y = 0;}
  Crd(int i, int j){x = i; y = j;}

  int x;
  int y;
};

const int LEFT = 0;
const int RIGHT = 1;
const int TOP = 2;
const int BOTTOM = 3;

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

struct Cell
{
  float b;
  float prevD;
  float d;
  float s;
  float f[4];
  float u;
  float v;
};

float* erodeField(float* field, int size)
{
  Cell defaultCell;
  defaultCell.b = 0.0f;
  defaultCell.d = 0.0f;
  defaultCell.s = 0.0f;
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
    std::cout << "Iteration " << i << std::endl;

    //Step 1: Add water through rainfall
    sim[rand() % size][rand() % size].d += RAINDROP_SIZE;

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
            float deltaHeight = (sim[x][y].b + sim[x][y].d) - (sim[side.x][side.y].b + sim[side.x][side.y].d);

            //find new flux value for direction
            sim[x][y].f[j] = std::max(0.0f, sim[x][y].f[j] + (TIME_STEP * PIPE_CROSS_SECTION * GRAVITY) / PIPE_LENGTH);
          }
          else
          {
            sim[x][y].f[j] = 0.0;
          }

          fluxSum += sim[x][y].f[j];
        }

        float scalingFactor = std::min(1.0f, (sim[x][y].d * PIPE_LENGTH * PIPE_LENGTH) / (fluxSum * TIME_STEP));

        //for each direction
        for(int j = 0; j < 4; j++)
        {
          //adjust based on scaling factor
          sim[x][y].f[j] *= scalingFactor;
        }
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
        //save previous water height for velocity calcs
        sim[x][y].prevD = sim[x][y].d;

        sim[x][y].d += totalDelta / (PIPE_LENGTH * PIPE_LENGTH);
      }
    }

    //Step 4: Adjust velocity field
    for(int x = 0; x < size; x++)
    {
      for(int y = 0; y < size; y++)
      {
        float avgWater = (sim[x][y].d + sim[x][y].prevD) / 2.0f;

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

        if(avgWater != 0.0)
          sim[x][y].u = ((lContrib + rContrib) / 2.0) / (avgWater * PIPE_LENGTH);
        else
          sim[x][y].u = 0.0;

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

        if(avgWater != 0.0)
          sim[x][y].v = ((bContrib + tContrib) / 2.0) / (avgWater * PIPE_LENGTH);
        else
          sim[x][y].v = 0.0;
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
        float normal[3] = {hHeightDelta, 2.0, vHeightDelta};
        float magnitude = pow(normal[0], 2) + pow(normal[1], 2) + pow(normal[2], 2);
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

          sim[x][y].b -= sedChange;
          sim[x][y].s += sedChange;
        }
        else
        {
          //deposit
          float sedChange = DEP_COEFF * (sim[x][y].s - transCapacity);

          sim[x][y].b += sedChange;
          sim[x][y].s -= sedChange;
        }
      }
    }

    //Step 6: Transport Sediment
    for(int x = 0; x < size; x++)
    {
      for(int y = 0; y < size; y++)
      {
        float xSed = x - (sim[x][y].u * TIME_STEP);
        float ySed = y - (sim[x][y].v * TIME_STEP);

        if(xSed >= size - 1 || xSed < 0 || ySed >= size - 1 || ySed < 0)
        {
          //do not move sediment
        }
        else
        {
          sim[x][y].s = getInterpValue(sim[xSed][ySed].s, sim[xSed + 1][ySed].s, sim[xSed][ySed + 1].s, sim[xSed + 1][ySed + 1].s, xSed - floor(xSed), ySed - floor(ySed));
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
