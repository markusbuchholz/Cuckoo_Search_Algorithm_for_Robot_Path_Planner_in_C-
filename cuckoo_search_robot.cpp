// Markus Buchholz
// g++ cuckoo_search_robot.cpp -o t -I/usr/include/python3.8 -lpython3.8

#include <iostream>
#include <vector>
#include <tuple>
#include <algorithm>
#include <math.h>
#include <random>

#include "matplotlibcpp.h"

namespace plt = matplotlibcpp;

//--------Path Planner--------------------------------------------------------------

float xmin = 0.0;
float xmax = 50.0;
float ymin = 0.0;
float ymax = 50.0;

float obsX = 25.0;
float obsY = 25.0;
float obsR = 3.0;

float goalX = 45.0;
float goalY = 45.0;

float startX = 2.0;
float startY = 2.0;

float K1 = 0.01;    // / obsR; // fitting parameter table 1
float K2 = 0.0001; // fitting parameter table 2

//--------------------------------------------------------------------------------
int EVOLUTIONS = 20;
int NESTS = 50;
float Pa = 0.3;

//--------------------------------------------------------------------------------

struct Pos
{

    float x;
    float y;
};
//--------------------------------------------------------------------------------

float generateNormalRandom()
{

    std::random_device engine;
    std::mt19937 gen(engine());
    std::normal_distribution<float> distrib(0, 1.0);
    return distrib(gen);
}
//--------------------------------------------------------------------------------

float generateGammaRand()
{

    std::random_device engine;
    std::mt19937 gen(engine());
    std::gamma_distribution<double> distribution(1, 2);

    return distribution(gen);
}
//--------------------------------------------------------------------------------

float generateSigmaU()
{

    float B = 1.5;
    float G1 = generateGammaRand();
    float G2 = generateGammaRand();
    float sigmaU_n = G1 * (1 + B) * std::sin(M_PI * B / 2);
    float sigmaU_d = G2 * ((1 + B) / 2) * B * std::pow(2, ((B - 1) / 2));

    return std::pow(sigmaU_n / sigmaU_d, 1 / B);
}
//--------------------------------------------------------------------------------

float generateS()
{

    float B = 1.5;

    return generateNormalRandom() * generateSigmaU() / std::pow(std::abs(generateNormalRandom()), 1 / B);
}
//--------------------------------------------------------------------------------

float euclid(Pos a, Pos b)
{

    return std::sqrt(std::pow(a.x - b.x, 2) + std::pow(a.y - b.y, 2));
}
//--------------------------------------------------------------------------------

float generateRandom()
{

    std::random_device engine;
    std::uniform_real_distribution<float> distrib(0, 1.0);
    return distrib(engine);
}

//--------------------------------------------------------------------------------

float valueGenerator(float low, float high)
{

    return low + generateRandom() * (high - low);
}

//--------------------------------------------------------------------------------

std::vector<float> function(std::vector<Pos> pos)
{

    std::vector<float> funcValue;
    Pos Obs{obsX, obsY};
    Pos Goal{goalX, goalY};

    for (auto &ii : pos)
    {

        funcValue.push_back(K1 * (1 / euclid(Obs, ii)) + K2 * euclid(Goal, ii));
    }

    return funcValue;
}

//--------------------------------------------------------------------------------

float func(Pos pos)
{
    Pos Obs{obsX, obsY};
    Pos Goal{goalX, goalY};

    return K1 * (1 / euclid(Obs, pos)) + K2 * euclid(Goal, pos);
}

//--------------------------------------------------------------------------------

Pos positionUpdateCheck(Pos actualPos)
{

    Pos Pnew = actualPos;

    if (Pnew.x < xmin)
    {
        Pnew.x = xmin;
    }

    if (Pnew.x > xmax)
    {
        Pnew.x = xmax;
    }

    if (Pnew.y < ymin)
    {
        Pnew.y = ymin;
    }

    if (Pnew.y > ymax)
    {
        Pnew.y = ymax;
    }

    return Pnew;
}

//--------------------------------------------------------------------------------

Pos posUpdateA(Pos act, Pos best)
{

    Pos Xnew;

    Xnew.x = act.x * generateRandom() * 0.01 * generateS() * (act.x - best.x);
    Xnew.y = act.y * generateRandom() * 0.01 * generateS() * (act.y - best.y);

    return positionUpdateCheck(Xnew);
}

//--------------------------------------------------------------------------------

Pos posUpdateB(Pos act, Pos d1, Pos d2)
{

    Pos newPos;

    float r1 = generateRandom();
    float r2 = generateRandom();

    if (r1 < Pa)
    {

        newPos.x = act.x + generateRandom() * (d1.x - d2.x);
    }

    else
    {

        newPos.x = act.x;
    }

    if (r2 < Pa)
    {

        newPos.y = act.y + generateRandom() * (d1.y - d2.y);
    }
    else
    {

        newPos.y = act.y;
    }

    return positionUpdateCheck(newPos);
}

//--------------------------------------------------------------------------------

std::vector<Pos> initPosXY()
{

    std::vector<Pos> pos;

    for (int ii = 0; ii < NESTS; ii++)
    {

        pos.push_back({valueGenerator(xmin, xmax), valueGenerator(ymin, ymax)});
    }

    return pos;
}

//-------------------------------------------------------------------------------
bool compareMin(std::pair<Pos, float> a, std::pair<Pos, float> b)
{

    return a.second < b.second;
}

//-------------------------------------------------------------------------------

// min
std::tuple<Pos, float> findBestPosFuncValue(std::vector<Pos> positions, std::vector<float> func)
{

    std::vector<std::pair<Pos, float>> best;

    for (int ii = 0; ii < func.size(); ii++)
    {

        best.push_back(std::pair<Pos, float>(positions[ii], func[ii]));
    }

    std::sort(best.begin(), best.end(), compareMin);

    return best[0];
}

//-------------------------------------------------------------------------------

int chooseNest()
{

    std::random_device engine;

    std::uniform_int_distribution<int> distribution(0, NESTS);

    return distribution(engine);
}

//-------------------------------------------------------------------------------

std::vector<Pos> runCuckoo()
{

    std::vector<Pos> currentPositions = initPosXY();
    std::vector<float> currentValueFunction = function(currentPositions);

    for (int ii = 0; ii < EVOLUTIONS; ii++)
    {
    std::tuple<Pos, float> bestPosFuncValueIx = findBestPosFuncValue(currentPositions, currentValueFunction);
    Pos bestPosIx = std::get<0>(bestPosFuncValueIx);
    float bestFuncValueIx = std::get<1>(bestPosFuncValueIx);

        for (int jj = 0; jj < NESTS; jj++)
        {

            Pos newPosA = posUpdateA(currentPositions[jj], bestPosIx);
            float newFunctionValueA = func(newPosA);

            if (newFunctionValueA < currentValueFunction[jj])
            {

                currentPositions[jj] = newPosA;
                currentValueFunction[jj] = newFunctionValueA;
            }

            Pos newPosB = posUpdateB(currentPositions[jj], currentPositions[chooseNest()], currentPositions[chooseNest()]);
            float newFunctionValueB = func(newPosB);

            if (newFunctionValueB < currentValueFunction[jj])
            {

                currentPositions[jj] = newPosB;
                currentValueFunction[jj] = newFunctionValueB;
            }
        }
    }


    return currentPositions;
}

//-------------------------------------------------------------------------------
std::tuple<std::vector<float>, std::vector<float>> gen_circle(float a, float b, float r)
{

    std::vector<float> xX;
    std::vector<float> yY;

    for (float dt = -M_PI; dt < M_PI; dt += 0.01)
    {

        xX.push_back(a + r * std::cos(dt));
        yY.push_back(b + r * std::sin(dt));
    }
    return std::make_tuple(xX, yY);
}

//-----------------------------------------------------------------------------------------

void plot2D(std::vector<float> xX, std::vector<float> yY)
{
    std::sort(xX.begin(), xX.end());
    std::sort(yY.begin(), yY.end());

    std::tuple<std::vector<float>, std::vector<float>> circle = gen_circle(obsX, obsY, obsR);

    std::vector<float> xObs = std::get<0>(circle);
    std::vector<float> yObs = std::get<1>(circle);

    plt::plot(xX, yY);
    plt::plot(xObs, yObs);
    plt::xlabel("X");
    plt::ylabel("Y");
    plt::show();
}

//------------------------------------------------------------------------------------------

int main()
{

    std::vector<Pos> path = runCuckoo();

    std::vector<float> xX;
    std::vector<float> yY;

    for (auto &ii : path)
    {
        xX.push_back(ii.x);
        yY.push_back(ii.y);

        std::cout << ii.x << " ," << ii.y << "\n";
    }

    plot2D(xX,yY);
}
