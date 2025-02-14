/** \file dab_math_roesseler_field_algorithm.cpp
*/

#include "dab_math_roesseler_field_algorithm.h"

using namespace dab;
using namespace dab::math;

RoesselerFieldAlgorithm::RoesselerFieldAlgorithm()
: FieldAlgorithm<float>( Eigen::Vector3f(-10.0, -10.0, -5.0), Eigen::Vector3f( 10.0, 10.0, 20.0 ) )
{}

RoesselerFieldAlgorithm::RoesselerFieldAlgorithm(const Eigen::Vector3f& pMinParValue, const Eigen::Vector3f& pMaxParValue)
: FieldAlgorithm<float>(pMinParValue,pMaxParValue)
{}

RoesselerFieldAlgorithm::~RoesselerFieldAlgorithm()
{}

void
RoesselerFieldAlgorithm::setVectorField(VectorField<float>* pVectorField) throw (Exception)
{
    if(pVectorField->fieldDim() != 3) throw Exception("MATH ERROR:  field field dim " + std::to_string(pVectorField->fieldDim()) + " doesn't match required dim " + std::to_string(3), __FILE__, __FUNCTION__, __LINE__);
    if(pVectorField->vectorDim() != 3) throw Exception("MATH ERROR:  field vector dim " + std::to_string(pVectorField->vectorDim()) + " doesn't match required dim " + std::to_string(3), __FILE__, __FUNCTION__, __LINE__);
    
    mVectorField = pVectorField;
    
    const dab::Array<unsigned int>& fieldSize = mVectorField->size();
    for(int d=0; d<3; ++d)
    {
		mParValueStep[d] = (mMaxParValue[d] - mMinParValue[d]) / (static_cast<float>(fieldSize.operator[](d)) - 1.0);
    }
}

void
RoesselerFieldAlgorithm::update() throw(dab::Exception)
{
    VectorField<float>* vectorField = FieldAlgorithm<float>::mVectorField;
    
    if(vectorField == nullptr) return;
    
    const dab::Array<unsigned int>& vectorFieldSize = vectorField->size();
    std::vector< Eigen::Matrix<float, Eigen::Dynamic, 1> >& fieldVectors = vectorField->vectors();
    Eigen::Vector3f parValue;
    
    for(unsigned int z=0, fieldIndex=0; z<vectorFieldSize[2]; ++z)
    {
        parValue[2] = mMinParValue[2] + mParValueStep[2] * static_cast<float>(z);
        
        for(unsigned int y=0; y<vectorFieldSize[1]; ++y)
        {
            parValue[1] = mMinParValue[1] + mParValueStep[1] * static_cast<float>(y);
            
            for(unsigned int x=0; x<vectorFieldSize[0]; ++x, ++fieldIndex)
            {
                //std::cout << "x " << x << " y " << y << " z " << z << " fieldIndex " << fieldIndex << "\n";
                
                parValue[0] = mMinParValue[0] + mParValueStep[0] * static_cast<float>(x);
                
                Eigen::Matrix<float, Eigen::Dynamic, 1>& fieldVector = fieldVectors[fieldIndex];
                
                fieldVector[0] = -parValue[1] - parValue[2];
                fieldVector[1] = parValue[0] + 0.2 * parValue[1];
                fieldVector[2] = 0.2 + parValue[0] * parValue[2] - 4.7 * parValue[2];
                
                //std::cout << "i " << fieldIndex << " vec " << fieldVector[0] << " " << fieldVector[1] << " " << fieldVector[2] << "\n";
            }
        }
    }
}

