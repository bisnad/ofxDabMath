/** \file dab_math.h
*/

#pragma once

#include <iostream>
#include <random>
#include <vector>
#include "ofUtils.h"
#include "dab_singleton.h"
#include "dab_exception.h"

namespace dab
{

namespace math
{

template<typename IntType = int, typename RealType = float>
class Math : public Singleton<Math<IntType, RealType> >
{
public:
    Math();
    ~Math();
    
    /**
     \brief integer power function
     \param pBase base
     \param pExponent exponent
     \return pBase ^ pExponent
     */
	inline IntType powI(IntType pBase, IntType pExponent) const;
    
    /**
     \brief uniform random function
     \param pMinValue minimum
     \param pMaxValue maximum
     \return uniform random value
     */
    RealType random( RealType pMin, RealType pMax );
    
    /**
     \brief uniform random function
     \param pMinValue minimum
     \param pMaxValue maximum
     \return uniform random value
     */
    RealType uniformRandom( RealType pMin, RealType pMax );
    
    /**
     \brief normal random function
     \param pMeanValue mean
     \param pStdDev standard deviation
     \return uniform random value
     */
    RealType normalRandom( RealType pMean, RealType pStdDev );
    
    /**
     \brief returns factorial of an positive integer
     \param pValue positive integer
     \return factorial
     */
	IntType
	factorial(IntType pValue)
	{
		IntType value = pValue;
		IntType returnValue = 1;
		
		for(; value > 0; --value) returnValue *= value;
		
		return returnValue;
	}
	
	/**
     \brief checks whether number is power of two
     \return true if number is power of two, false otherwise
     */
	bool checkPowerOfTwo(IntType pValue);
	
	/**
     \brief transforms second input value into a (larger) number that's related to the first input value by a power of two
     \param pFixValue fixed input value
     \param pVarValue variable input value
     */
	void powerOfTwoRatio(RealType pFixValue, RealType& pVarValue);
    
	/**
     \brief converts degrees to radians
     */
	RealType radian(RealType pValue);
    
	/**
     \brief converts radians to degrees
     */
	RealType degree(RealType pValue);
	
	/**
     \brief converts degrees to radians
     */
	RealType degreeToRadian(RealType pValue);
	
	/**
     \brief converts radians to degrees
     */
	RealType radianToDegree(RealType pValue);
	
	/**
     \brief transforms complex number into polar angle
     \param aValue aValue
     \param bValue bValue
     \return angle in radians
     */
	RealType complexToPhase(RealType aValue, RealType bValue);
	
	/**
     \brief gaussian function
     \param pXValue x value
     \param pMean mean
     \param pStd standard deviation
     /return y value
     */
	RealType gaussian(RealType pXValue, RealType pMean, RealType pStd);
    
	/**
     \brief hamming function
     \param pXValue x value
     /return y value
     
     runs from 0.0 to 1.0
     */
	RealType hamming(RealType pXValue);
    
	/**
     \brief hann function
     \param pXValue x value
     /return y value
     
     runs from 0.0 to 1.0
     */
	RealType hann(RealType pXValue);
    
	/**
     \brief hann function
     \param pXValue x value
     /return y value
     
     runs from 0.0 to 1.0
     */
	RealType kaiser(RealType pXValue);
	
	/**
     \brief sigmoid function
     \param pXValue x value
     \param pA a value
     \return y value
     */
	RealType sigmoid(RealType pXValue, RealType pA);
	
	/**
     \brief normalized sigmoid function
     \param pXValue x value
     \param pA a value
     \return y value
     
     the normal function input range from -6 to +6 has been rescaled and shifted to the range of 0 to 1
     */
	RealType nSigmoid(RealType pXValue, RealType pA);
	
	/**
     \brief calculate index into n-dimensional array
     \param pSubdivisions array subdivisions
     \param pPosition position into array
     \return index
     \exception Exception dimensions of subdivisons and position don't match or if positition is out of bounds
     */
	IntType index( const std::vector<IntType>& pSubdivisions, const std::vector<IntType>& pPosition ) throw (Exception);
    
protected:
    std::ranlux48_base mRandomEngine;
    std::uniform_real_distribution<RealType> mUniformRealDistribution;
    std::normal_distribution<RealType> mNormalRealDistribution;
};

template<typename IntType, typename RealType>
IntType
Math<IntType, RealType>::powI(IntType pBase, IntType pExponent) const
{
    IntType result = 1;
    
	while(--pExponent >= 0) result *= pBase;
	
	return result;
}
    
template<typename IntType, typename RealType>
Math<IntType, RealType>::Math()
{
    uint64_t seed = ofGetSystemTimeMicros();
    mRandomEngine.seed(seed);
    mUniformRealDistribution = std::uniform_real_distribution<RealType>(0.0, 1.0);
    mNormalRealDistribution = std::normal_distribution<RealType>(0.0, 1.0);
};
    
template<typename IntType, typename RealType>
RealType
Math<IntType, RealType>::random( RealType pMin, RealType pMax )
{
    return uniformRandom(pMin, pMax);
}

template<typename IntType, typename RealType>
RealType
Math<IntType, RealType>::uniformRandom( RealType pMin, RealType pMax )
{
    typename std::uniform_real_distribution<RealType>::param_type new_param(pMin, pMax);
    mUniformRealDistribution.param(new_param);
    return mUniformRealDistribution(mRandomEngine);
}
    
template<typename IntType, typename RealType>
RealType
Math<IntType, RealType>::normalRandom( RealType pMean, RealType pStdDev )
{
    typename std::normal_distribution<RealType>::param_type new_param(pMean, pStdDev);
    mNormalRealDistribution.param(new_param);
    return mNormalRealDistribution(mRandomEngine);
}
    
template<typename IntType, typename RealType>
bool
Math<IntType, RealType>::checkPowerOfTwo(IntType pValue)
{
    return !((pValue - 1) & pValue);
}

template<typename IntType, typename RealType>
void
Math<IntType, RealType>::powerOfTwoRatio(RealType pFixValue, RealType& pVarValue)
{
    if(pVarValue <= pFixValue)
    {
        RealType currentRatio = pFixValue / pVarValue;
        unsigned int divisionCount = 0;
        
        while(currentRatio >= 2.0)
        {
            divisionCount++;
            currentRatio /= 2.0;
        }
        
        IntType power2Ratio = powI(2, divisionCount);
        pVarValue = pFixValue / static_cast<RealType>(power2Ratio);
    }
    else
    {
        RealType currentRatio = pVarValue / pFixValue;
        IntType divisionCount = 1;
        
        while(currentRatio >= 2.0)
        {
            divisionCount++;
            currentRatio /= 2.0;
        }
        
        IntType power2Ratio = powI(2, divisionCount);
        pVarValue = pFixValue * static_cast<RealType>(power2Ratio);
    }
}

template<typename IntType, typename RealType>
RealType
Math<IntType, RealType>::radian(RealType pValue)
{
    return pValue * PI / 180.0;
}

template<typename IntType, typename RealType>
RealType
Math<IntType, RealType>::degree(RealType pValue)
{
    return pValue * 180.0 / PI;
}

template<typename IntType, typename RealType>
RealType
Math<IntType, RealType>::degreeToRadian(RealType pValue)
{
    return pValue * PI / 180.0;
}

template<typename IntType, typename RealType>
RealType
Math<IntType, RealType>::radianToDegree(RealType pValue)
{
    return pValue * 180.0 / PI;
}

template<typename IntType, typename RealType>
RealType
Math<IntType, RealType>::complexToPhase(RealType aValue, RealType bValue)
{
    if(aValue == 0.0) aValue = 0.00001;
    
    if(aValue >= 0.0 && bValue >= 0.0) return atan(bValue / aValue);
    if(aValue < 0.0) return atan(bValue / aValue) + PI;
    if(bValue < 0.0) return atan(bValue / aValue) + 2.0 * PI;
}

template<typename IntType, typename RealType>
RealType
Math<IntType, RealType>::gaussian(RealType pXValue, RealType pMean, RealType pStd)
{
    return 1.0 / (pStd * sqrt(2.0 * PI)) * exp(- ( pow ( ( pXValue - pMean), static_cast<RealType>( 2.0 ) )  ) / ( 2.0 * pow (pStd, static_cast<RealType>( 2.0 ) )  ) );
}

template<typename IntType, typename RealType>
RealType
Math<IntType, RealType>::hamming(RealType pXValue)
{
    return 0.53836 - 0.46164 * cos( 2.0 * PI * pXValue );
}

template<typename IntType, typename RealType>
RealType
Math<IntType, RealType>::hann(RealType pXValue)
{
    return 0.5 * ( 1.0 - cos( 2.0 * PI * pXValue ) );
}

template<typename IntType, typename RealType>
RealType
Math<IntType, RealType>::sigmoid(RealType pXValue, RealType pA)
{
    return 1.0 / (1.0 + pow(pA, -pXValue));
}

template<typename IntType, typename RealType>
RealType
Math<IntType, RealType>::nSigmoid(RealType pXValue, RealType pA)
{
    return 1.0 / (1.0 + pow(pA, static_cast<RealType>( -(pXValue * 12.0 - 6.0) ) ) );
}

template<typename IntType, typename RealType>
IntType
Math<IntType, RealType>::index( const std::vector<IntType>& pSubdivisions, const std::vector<IntType>& pPosition ) throw (Exception)
{
    int dim = pSubdivisions.size();
    
    if(dim != pPosition.size() ) throw Exception( "MATH ERROR: subdivision dim " + std::to_string(dim) + " does not match position dim " + std::to_string(pPosition.size()), __FILE__, __FUNCTION__, __LINE__ );
    for(unsigned int d=0; d<dim; ++d)
    {
        if(pPosition.c[d] >= pSubdivisions.c[d]) throw MathException( "MATH ERROR: for dimension " + std::to_string(d) + " position " + std::to_string(pPosition[d]) + " exceeds subdivision count " + std::to_string(pSubdivisions[d]), __FILE__, __FUNCTION__, __LINE__ );
    }
    
    unsigned int resultIndex = pPosition.c[0];
    unsigned int offset = 1;
    
    for(int d=1; d<dim; ++d)
    {
        //std::cout << "d " << d << " sub " << pSubdivisions.c[d] << " pos " << pPosition.c[d] << " offset " << offset << " res " << resultIndex << "\n";
        
        offset *= pSubdivisions[d-1];
        resultIndex += pPosition[d] * offset;
        
        //std::cout << "d " << d << " sub " << pSubdivisions.c[d] << " pos " << pPosition.c[d] << " offset " << offset << " res " << resultIndex << "\n"; 
    }
    
    return resultIndex;
}

};
    
};