/** \file dab_math_field_algorithm.h
*/

#pragma once

#include <iostream>
#include <Eigen/Dense>
#include "dab_math_vector_field.h"

namespace dab
{
    
namespace math
{
    
#pragma mark AbstractFieldAlgorithm definition

class AbstractFieldAlgorithm
{
public:
    AbstractFieldAlgorithm();
    virtual ~AbstractFieldAlgorithm();
    
    virtual void update() throw (Exception) = 0;
    
    void resetTime(float pTime);
    void setTimeStep(float pTimeStep);
    void setTimeLimit(float pTimeLimit);
    
protected:
    static float sTime;
    static float sTimeStep;
    static float sTimeLimit;
    
    float mTime;
    float mTimeStep;
    float mTimeLimit;
    
    unsigned int mFieldDim;
    unsigned int mParameterDim;
};
    
#pragma mark FieldAlgorithm definition

// TODO: at the moment: ParameterDim and FieldDim must be identical, figure out how to relax this condition
    
template<typename ValueType>
class FieldAlgorithm : public AbstractFieldAlgorithm
{
friend class VectorField<ValueType>;
    
public:
    FieldAlgorithm();
	FieldAlgorithm(const Eigen::Matrix<ValueType, Eigen::Dynamic, 1>& pMinParValue, const Eigen::Matrix<ValueType, Eigen::Dynamic, 1>& pMaxParValue);
	~FieldAlgorithm();
    
    void setVectorField(VectorField<ValueType>* pVectorField) throw (Exception);
    
protected:
    Eigen::Matrix<ValueType, Eigen::Dynamic, 1> mMinParValue;
    Eigen::Matrix<ValueType, Eigen::Dynamic, 1> mMaxParValue;
    Eigen::Matrix<ValueType, Eigen::Dynamic, 1> mParValueStep;
    VectorField<ValueType>* mVectorField;
};
    
#pragma mark FieldAlgorithm implementation
    
template<typename ValueType>
FieldAlgorithm<ValueType>::FieldAlgorithm()
    : mVectorField(nullptr)
{}
    
template<typename ValueType>
FieldAlgorithm<ValueType>::FieldAlgorithm(const Eigen::Matrix<ValueType, Eigen::Dynamic, 1>& pMinParValue, const Eigen::Matrix<ValueType, Eigen::Dynamic, 1>& pMaxParValue)
    : mVectorField(nullptr)
    , mMinParValue(pMinParValue)
    , mMaxParValue(pMaxParValue)
    , mParValueStep(pMinParValue.rows())
{
    mParameterDim = pMinParValue.rows();
    
    for(int d=0; d<mParameterDim; ++d)
    {
        mParValueStep[d] = 0.0;
    }
}
    
template<typename ValueType>
FieldAlgorithm<ValueType>::~FieldAlgorithm()
{}
    
template<typename ValueType>
void
FieldAlgorithm<ValueType>::setVectorField(VectorField<ValueType>* pVectorField) throw (Exception)
{
    if(pVectorField->VectorDim != mParameterDim) throw Exception("MATH ERROR:  field vector dim " + std::to_string(pVectorField->VectorDim) + " doesn't match algorithm parameter dim " + std::to_string(mParameterDim), __FILE__, __FUNCTION__, __LINE__);
    
    mVectorField = pVectorField;
}
    
};

};