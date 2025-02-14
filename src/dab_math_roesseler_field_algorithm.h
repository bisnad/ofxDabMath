/** \file dab_math_roesseler_field_algorithm.h
*/

#pragma once

#include "dab_math_field_algorithm.h"
#include "dab_math_vector_field.h"

namespace dab
{

namespace math
{

class RoesselerFieldAlgorithm : public FieldAlgorithm<float>
{
public:
    RoesselerFieldAlgorithm();
    RoesselerFieldAlgorithm(const Eigen::Vector3f& pMinParValue, const Eigen::Vector3f& pMaxParValue);
    ~RoesselerFieldAlgorithm();
    
    void setVectorField(VectorField<float>* pVectorField) throw (Exception);
    
    void update() throw (Exception);
    
protected:
};
        
};
    
};