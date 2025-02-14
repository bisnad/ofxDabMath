/** \file dab_math_field_algorithm.cpp
*/

#include "dab_math_field_algorithm.h"

using namespace dab;
using namespace dab::math;

float AbstractFieldAlgorithm::sTime = 0.0;
float AbstractFieldAlgorithm::sTimeStep = 0.001;
float AbstractFieldAlgorithm::sTimeLimit = 1.0;

AbstractFieldAlgorithm::AbstractFieldAlgorithm()
: mTime(sTime)
, mTimeStep(sTimeStep)
, mTimeLimit(sTimeLimit)
{}

AbstractFieldAlgorithm::~AbstractFieldAlgorithm()
{}

void
AbstractFieldAlgorithm::resetTime(float pTime)
{
    mTime = pTime;
}

void
AbstractFieldAlgorithm::setTimeStep(float pTimeStep)
{
    mTimeStep = pTimeStep;
}

void
AbstractFieldAlgorithm::setTimeLimit(float pTimeLimit)
{
    mTimeLimit = pTimeLimit;
}