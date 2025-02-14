/** \file dab_math_vector_field.h
 */

#pragma once

#include <iostream>
#include <Eigen/Dense>
#include "dab_exception.h"
#include "dab_array.h"

namespace dab
{

namespace math
{

template<typename ValueType>
class VectorField
{
public:
    /**
     \brief create vector field
     \param pSize size of vector field (i.e. number of vectors)
     */
    VectorField( const dab::Array<unsigned int>& pSize, const Eigen::Matrix<ValueType, Eigen::Dynamic, 1>& pRefVector);
    
    /**
     \brief create vector field
     \param pVectorField vector field
     */
    VectorField(const VectorField<ValueType>& pVectorField);

    /**
     \brief destructor
     */
    ~VectorField();

    /**
     \brief copy operator
     \param pVectorField vectorfield to copy
     \return vectorfield
     */
    const VectorField<ValueType>& operator= (const VectorField<ValueType>& pVectorField);

    /**
     \brief return field dimension
     \return field dimension
     */
    unsigned int fieldDim() const;
    
    /**
     \brief return vector dimension
     \return vector dimension
     */
    unsigned int vectorDim() const;
    
    /**
     \brief return size of vector field
     \return size of vector field
     */
    const dab::Array<unsigned int>& size() const;

    /**
     \brief return number of vectors
     \return number of vectors
     */
    unsigned int vectorCount() const;

    /**
     \brief return index offset
     \return index offset
     */
    const dab::Array<unsigned int>& indexOffset() const;

    /**
     \brief return vectors
     \return vectors
     */
    std::vector< Eigen::Matrix<ValueType, Eigen::Dynamic, 1> >& vectors();

    /**
     \brief return vectors
     \return vectors
     */
    const std::vector< Eigen::Matrix<ValueType, Eigen::Dynamic, 1> >& vectors() const;

    /**
     \brief return sum of all vectors
     \return sum of all vectors
     */
    Eigen::Matrix<ValueType, Eigen::Dynamic, 1> sum() const;
    
    /**
     \brief retrieve subfield
     \param pIndex position of subfield (top left)
     \param pSubField subfield
     \exception Exception out of bounds error
     */
    void subField( const dab::Array<unsigned int>& pIndex, VectorField<ValueType>& pSubField ) throw (Exception);
    
    /**
     \brief return vector
     \param pIndex scalar index
     \return vector
     \exception Exception if index is out of bounds
     */
    inline const Eigen::Matrix<ValueType, Eigen::Dynamic, 1>& operator[] (unsigned int pIndex) const throw (Exception);

    /**
     \brief return vector
     \param pIndex scalar index
     \return vector
     \exception MathException index is out of bounds
     */
    inline const Eigen::Matrix<ValueType, Eigen::Dynamic, 1>& get(unsigned int pIndex) const throw (Exception);
    
    /**
     \brief return vector
     \param pIndex nD index
     \return vector
     \exception Exception index is out of bounds
     */
    const Eigen::Matrix<ValueType, Eigen::Dynamic, 1>& get(const dab::Array<unsigned int>& pIndex) const throw (Exception);

    /**
     \brief return interpolated vector
     \param pIndex nD index
     \return vector
     \exception Exception index is out of bounds
     
     perform linear interpolation
     */
    Eigen::Matrix<ValueType, Eigen::Dynamic, 1> get(const dab::Array<float>& pIndex) const throw (Exception);

    /**
     \brief return interpolated vector
     \param pIndex nD index
     \param pResult interpolated vector
     \exception Exception index is out of bounds
     
     perform linear interpolation
     */
    void get(const dab::Array<float>& pIndex, Eigen::Matrix<ValueType, Eigen::Dynamic, 1>& pResult) const throw (Exception);
    
    /**
     \brief set vector
     \param pIndex scalar index
     \param pValue vector value
     \exception Exception out of bounds error
     */
    void set(unsigned int pIndex, const Eigen::Matrix<ValueType, Eigen::Dynamic, 1>& pValue) throw (Exception);
    
    /**
     \brief set vector
     \param pIndex vector index
     \param pValue vector value
     \exception Exception index out of bounds
     */
    void set(const dab::Array<unsigned int>& pIndex, const Eigen::Matrix<ValueType, Eigen::Dynamic, 1>& pValue) throw (Exception);

    /**
     \brief set vector (interpolated mode)
     \param pIndex vector index
     \param pValue vector value
     \exception Exception index out of bounds
     */
    void set(const dab::Array<float>& pIndex, const Eigen::Matrix<ValueType, Eigen::Dynamic, 1>& pValue) throw (Exception);
    
    /**
     \brief set vectors
     \param pStartIndex start vector index
     \param pEndIndex end vector index
     \param pValue vector value
     \exception Exception index out of bounds or start index larger than end index
     */
    void set(const dab::Array<unsigned int>& pStartIndex, const dab::Array<unsigned int>& pEndIndex, const Eigen::Matrix<ValueType, Eigen::Dynamic, 1>& pValue) throw (Exception);

    /**
     \brief set all vector values
     \param pValue vector value
     \exception Exception incompatible value dimension
     */
    void set(const Eigen::Matrix<ValueType, Eigen::Dynamic, 1>& pValue) throw (Exception);
    
    /**
     \brief add vector
     \param pIndex scalar index
     \param pValue vector value
     \exception Exception out of bounds error or incompatible value dimension
     */
    void add(unsigned int pIndex, const Eigen::Matrix<ValueType, Eigen::Dynamic, 1>& pValue) throw (Exception);

    /**
     \brief add vector
     \param pIndex vector index
     \param pValue vector value
     \exception Exception index out of bounds or incompatible value dimension
     */
    void add(const dab::Array<unsigned int>& pIndex, const Eigen::Matrix<ValueType, Eigen::Dynamic, 1>& pValue) throw (Exception);

    /**
     \brief add vector (interpolated mode)
     \param pIndex vector index
     \param pValue vector value
     \exception Exception index out of bounds or incompatible value dimension
     */
    void add(const dab::Array<float>& pIndex, const Eigen::Matrix<ValueType, Eigen::Dynamic, 1>& pValue) throw (Exception);

    /**
     \brief add vector
     \param pStartIndex start vector index
     \param pEndIndex end vector index
     \param pValue vector value
     \exception Exception index out of bounds or start index larger than end index or incompatible value dimension
     */
    void add(const dab::Array<unsigned int>& pStartIndex, const dab::Array<unsigned int>& pEndIndex, const Eigen::Matrix<ValueType, Eigen::Dynamic, 1>& pValue) throw (Exception);

    /**
     \brief add all vector values
     \param pValue vector value
     \exception Exception incompatible value dimension
     */
    void add(const Eigen::Matrix<ValueType, Eigen::Dynamic, 1>& pValue) throw (Exception);
    
    /**
     \brief vectorfield += vector
     \param pValue value
     \return vectorfield
     \exception Exception incompatible value dimension
     */
    inline const VectorField<ValueType>& operator+=( const Eigen::Matrix<ValueType, Eigen::Dynamic, 1>& pValue ) throw (Exception);
    
    /**
     \brief vectorfield -= vector
     \param pValue value
     \return vectorfield
     \exception Exception incompatible value dimension
     */
    inline const VectorField<ValueType>& operator-=( const Eigen::Matrix<ValueType, Eigen::Dynamic, 1>& pValue ) throw (Exception);
    
    /**
     \brief vectorfield *= vector
     \param pValue value
     \return vectorfield
     \exception Exception incompatible value dimension
     */
    inline const VectorField<ValueType>& operator*=( const Eigen::Matrix<ValueType, Eigen::Dynamic, 1>& pValue ) throw (Exception);
    
    /**
     \brief vectorfield /= vector
     \param pValue value
     \return vectorfield
     \exception Exception incompatible value dimension
     */
    inline const VectorField<ValueType>& operator/=( const Eigen::Matrix<ValueType, Eigen::Dynamic, 1>& pValue ) throw (Exception);

    /**
     \brief vectorfield + vector
     \param pValue value
     \return vectorfield
     \exception Exception incompatible value dimension
     */
    inline VectorField<ValueType> operator+( const Eigen::Matrix<ValueType, Eigen::Dynamic, 1>& pValue ) throw (Exception);

    /**
     \brief vectorfield - vector
     \param pValue value
     \return vectorfield
     \exception Exception incompatible value dimension
     */
    inline VectorField<ValueType> operator-( const Eigen::Matrix<ValueType, Eigen::Dynamic, 1>& pValue ) throw (Exception);

    /**
     \brief vectorfield * vector
     \param pValue value
     \return vectorfield
     \exception Exception incompatible value dimension
     */
    inline VectorField<ValueType> operator*( const Eigen::Matrix<ValueType, Eigen::Dynamic, 1>& pValue ) throw (Exception);

    /**
     \brief vectorfield / vector
     \param pValue value
     \return vectorfield
     \exception Exception incompatible value dimension
     */
    inline VectorField<ValueType> operator/( const Eigen::Matrix<ValueType, Eigen::Dynamic, 1>& pValue ) throw (Exception);

    /**
     \brief vectorfield += vectorfield
     \param pVectorField vectorfield
     \return vectorfield
     \exception Exception field size mismatch
     */
    inline const VectorField<ValueType>& operator+=( const VectorField<ValueType>& pVectorField ) throw (Exception);

    /**
     \brief vectorfield -= vectorfield
     \param pVectorField vectorfield
     \return vectorfield
     \exception Exception ield size mismatch
     */
    inline const VectorField<ValueType>& operator-=( const VectorField<ValueType>& pVectorField ) throw (Exception);
    
    /**
     \brief vectorfield *= vectorfield
     \param pVectorField vectorfield
     \return vectorfield
     \exception Exception field size mismatch
     */
    inline const VectorField<ValueType>& operator*=( const VectorField<ValueType>& pVectorField ) throw (Exception);
    

    /**
     \brief vectorfield /= vectorfield
     \param pVectorField vectorfield
     \return vectorfield
     \exception Exception field size mismatch
     */
    inline const VectorField<ValueType>& operator/=( const VectorField<ValueType>& pVectorField ) throw (Exception);

    /**
     \brief vectorfield + vectorfield
     \param pVectorField vectorfield
     \return vectorfield
     \exception Exception field size mismatch
     */
    inline VectorField<ValueType> operator+( const VectorField<ValueType>& pVectorField ) throw (Exception);

    /**
     \brief vectorfield - vectorfield
     \param pVectorField vectorfield
     \return vectorfield
     \exception Exception field size mismatch
     */
    inline VectorField<ValueType> operator-( const VectorField<ValueType>& pVectorField ) throw (Exception);
    
    /**
     \brief vectorfield * vectorfield
     \param pVectorField vectorfield
     \return vectorfield
     \exception Exception field size mismatch
     */
    inline VectorField<ValueType> operator*( const VectorField<ValueType>& pVectorField ) throw (Exception);
    
    /**
     \brief vectorfield / vectorfield
     \param pVectorField vectorfield
     \return vectorfield
     \exception Exception field size mismatch
     */
    inline VectorField<ValueType> operator/( const VectorField<ValueType>& pVectorField ) throw (Exception);
    
    /**
     \brief calculate index
     \param pIndex nD index
     \return 1D index
     \exception Exception dimension mismatch or index out of bounds
     
     calculates scalar index from an n-dimensional index
     */
    inline unsigned int calcIndex(const dab::Array<unsigned int>& pIndex) const throw (Exception);

    /**
     \brief calculate index
     \param pIndex index
     \return nD index
     \exception Exception index out of bounds
     
     calculates n-dimensional index from scalar index
     */
    inline dab::Array<unsigned int> calcIndex(unsigned int pIndex) const throw (Exception);
    
    /**
     \brief calculate index
     \param pIndex scalar index
     \param pVecIndex n-dimensional index
     \exception Exception dimension mismatch or index out of bounds
     
     calculates n-dimensional index from scalar index
     */
    inline void calcIndex(unsigned int pIndex, dab::Array<unsigned int>& pVecIndex) const throw (Exception);

    /**
     \brief convolve vector field with supplied kernel
     \param pKernel convolution kernel
     \exception Exception kernel size error (size not odd)
     */
    void convolve( const VectorField<ValueType>& pKernel ) throw (Exception);

    /**
     \brief obtain textual vector field information
     */
    operator std::string() const;

    /**
     \brief output textual vector field information
     \param os output stream
     \param pVectorField vector field
     \return output stream
     */
    friend std::ostream& operator << ( std::ostream& pOstream, const VectorField<ValueType>& pVectorField )
    {
        pOstream << std::string(pVectorField);
        return pOstream;
    };
    
protected:
    /**
     \brief dimension of vector field
     */
    unsigned int mFieldDim;
    
    /**
     \brief dimension of vector values
     */
    unsigned int mVectorDim;
    
    /**
     \brief total number of vectors
     */
    unsigned int mVectorCount;

    /**
     \brief size of vectors field
     */
    dab::Array<unsigned int> mSize;
    
    /**
     \brief offsets for converting a relative index into an absolute index
     */
    dab::Array<unsigned int> mIndexOffset;
    dab::Array<float> mInvIndexOffset;
    
    /**
     \brief vectors
     */
    std::vector< Eigen::Matrix<ValueType, Eigen::Dynamic, 1> > mVectors;
    
    unsigned int mInterpolVectorCount; ///\brief number of vectors to test for linear interpolation
};

template<typename ValueType>
    VectorField<ValueType>::VectorField(const dab::Array<unsigned int>& pSize, const Eigen::Matrix<ValueType, Eigen::Dynamic, 1>& pRefVector)
    : mSize(pSize)
    , mFieldDim(pSize.size())
    , mVectorDim(pRefVector.rows())
    , mIndexOffset(pSize.size())
    , mInvIndexOffset(pSize.size())
{
    mVectorCount = 1;
    for(unsigned int d=0; d<mFieldDim; ++d)
    {
        mVectorCount *= mSize[d];
        mIndexOffset[d] = 1;
        for(unsigned int j=0; j<d; ++j) mIndexOffset[d] *= mSize[j];
        mInvIndexOffset[d] = 1.0 / static_cast<float>( mIndexOffset[d] );
        
        //std::cout << "d " << d << " size " << mSize[d] << " count " << mVectorCount << " offset " << mIndexOffset[d] << " invOffset " << mInvIndexOffset[d] << "\n";
    }
    
    mVectors.resize(mVectorCount, pRefVector);
    mInterpolVectorCount = static_cast<unsigned int>( pow( 2.0f, static_cast<float>( mFieldDim ) ) );
}
    
template<typename ValueType>
VectorField<ValueType>::VectorField(const VectorField<ValueType>& pVectorField)
    : mFieldDim(pVectorField.mFieldDim)
    , mVectorDim(pVectorField.mVectorDim)
    , mVectorCount(pVectorField.mVectorCount)
    , mSize(pVectorField.mSize)
    , mIndexOffset(pVectorField.mIndexOffset)
    , mInvIndexOffset(pVectorField.mInvIndexOffset)
    , mVectors(pVectorField.mVectors)
    , mInterpolVectorCount(pVectorField.mInterpolVectorCount)
{}

template<typename ValueType>
VectorField<ValueType>::~VectorField()
{}
   
template<typename ValueType>
const VectorField<ValueType>&
VectorField<ValueType>::operator= (const VectorField<ValueType>& pVectorField)
{
    mFieldDim = pVectorField.mFieldDim;
    mVectorDim = pVectorField.mVectorDim;
    mVectorCount = pVectorField.mVectorCount;
    mSize = pVectorField.mSize;
    mIndexOffset = pVectorField.mIndexOffset;
    mInvIndexOffset = pVectorField.mInvIndexOffset;
    mVectors = pVectorField.mVectors;
    mInterpolVectorCount = pVectorField.mInterpolVectorCount;
                        
    return *this;
}
    
template<typename ValueType>
unsigned int
VectorField<ValueType>::fieldDim() const
{
    return mFieldDim;
}
    
template<typename ValueType>
unsigned int
VectorField<ValueType>::vectorDim() const
{
    return mVectorDim;
}
    
template<typename ValueType>
const dab::Array<unsigned int>&
VectorField<ValueType>::size() const
{
    return mSize;
}

template<typename ValueType>
unsigned int
VectorField<ValueType>::calcIndex(const dab::Array<unsigned int>& pIndex) const throw (Exception)
{
    if( pIndex.size() != mFieldDim ) throw Exception("MATH ERROR: index dimension " + std::to_string(pIndex.size()) + " doesn't match field dimension " +  std::to_string(mFieldDim), __FILE__, __FUNCTION__, __LINE__);
    for(int d=0; d<mFieldDim; ++d)
    {
        if(pIndex[d] >= mSize[d]) throw Exception("MATH ERROR: at dimension " + std::to_string(d) + " index " +  std::to_string(pIndex[d]) + " exceeds field size " + std::to_string(mSize[d]), __FILE__, __FUNCTION__, __LINE__);
    }
    
    unsigned int index = 0;
    for(int d = 0; d < mFieldDim; ++d) index += pIndex[d] * mIndexOffset[d];
    return index;
}

template<typename ValueType>
dab::Array<unsigned int>
VectorField<ValueType>::calcIndex(unsigned int pIndex) const throw (Exception)
{
    if(pIndex >= mVectorCount) throw Exception( "MATH ERROR: index " + std::to_string(pIndex) + " exceeds vector count " + std::to_string(mVectorCount), __FILE__, __FUNCTION__, __LINE__ );
    
    dab::Array<unsigned int> index(mFieldDim);
    
    for(int d = mFieldDim - 1; d >= 0; --d)
    {
        index[d] = pIndex * mInvIndexOffset[d];
        pIndex -= index[d] * mIndexOffset[d];
    }
    
    return index;
}
    
template<typename ValueType>
void
VectorField<ValueType>::calcIndex(unsigned int pIndex, dab::Array<unsigned int>& pVecIndex) const throw (Exception)
{
    if(pIndex >= mVectorCount) throw Exception( "MATH ERROR: index " + std::to_string(pIndex) + " exceeds vector count " + std::to_string(mVectorCount), __FILE__, __FUNCTION__, __LINE__ );
    if(pVecIndex.size() != mFieldDim) throw Exception("MATH ERROR: index dimension " + std::to_string(pVecIndex.size()) + " doesn't match field dimension " + std::to_string(mFieldDim), __FILE__, __FUNCTION__, __LINE__);
    
    for(int d = mFieldDim - 1; d >= 0; --d)
    {
        pVecIndex[d] = pIndex * mInvIndexOffset[d];
        pIndex -= pVecIndex[d] * mIndexOffset[d];
    }
}

template<typename ValueType>
unsigned int
VectorField<ValueType>::vectorCount() const
{
    return mVectorCount;
}

template<typename ValueType>
const dab::Array<unsigned int>&
VectorField<ValueType>::indexOffset() const
{
    return mIndexOffset;
}

template<typename ValueType>
std::vector< Eigen::Matrix<ValueType, Eigen::Dynamic, 1> >&
VectorField<ValueType>::vectors()
{
    return mVectors;
}
    
template<typename ValueType>
const std::vector< Eigen::Matrix<ValueType, Eigen::Dynamic, 1> >&
VectorField<ValueType>::vectors() const
{
    return mVectors;
}
    
template<typename ValueType>
Eigen::Matrix<ValueType, Eigen::Dynamic, 1>
VectorField<ValueType>::sum() const
{
        Eigen::Matrix<ValueType, Eigen::Dynamic, 1> tmp(mVectorDim);
        for(int d=0; d<mVectorDim; ++d) tmp[d] = 0.0;
    
        for(unsigned int vI=0; vI<mVectorCount; ++vI)
        {
            tmp += mVectors[vI];
        }
        return tmp;
}
    
template<typename ValueType>
void
VectorField<ValueType>::subField( const dab::Array<unsigned int>& pIndex, VectorField<ValueType>& pSubField ) throw (Exception)
{
    //BUGGY; but in the same way buggy as the original iso version, fix bug
    
    if(pIndex.size() != mFieldDim) throw Exception("MATH ERROR: index dimension " + std::to_string(pIndex.size()) + " doesn't match field dimension " + std::to_string(mFieldDim), __FILE__, __FUNCTION__, __LINE__);
    if(pSubField.mFieldDim != mFieldDim) throw Exception("MATH ERROR: subfield dimension " + std::to_string(pSubField.mFieldDim) + " doesn't match field dimension " + std::to_string(mFieldDim), __FILE__, __FUNCTION__, __LINE__);
    if(pSubField.mVectorDim != mVectorDim) throw Exception("MATH ERROR: subfield vector dimension " + std::to_string(pSubField.mVectorDim) + " doesn't match field vector dimension " + std::to_string(mVectorDim), __FILE__, __FUNCTION__, __LINE__);
    
    try
    {
        unsigned int fDim = mFieldDim;
        unsigned int vDim = mVectorDim;
        const dab::Array<unsigned int>& sfSize = pSubField.size(); // subfield size
        std::vector< Eigen::Matrix<ValueType, Eigen::Dynamic, 1> >& sfVectors = pSubField.vectors(); // subfield vectors
        const dab::Array<unsigned int>& sfIndexOffset = pSubField.mIndexOffset;
        
        for(unsigned int d=0; d<fDim; ++d) if(sfSize[d] > mSize[d]) throw Exception( "MATH ERROR: at dim " + std::to_string(d) + " subfield size: " + std::to_string(sfSize[d]) + " is larger than vector field size " + std::to_string(mSize[d]), __FILE__, __FUNCTION__, __LINE__ );

        dab::Array<unsigned int> fStartIndex(pIndex); // relative field start index
        dab::Array<unsigned int> fEndIndex(fDim); // relative field end index
        for(unsigned int d=0; d<fDim; ++d) fEndIndex[d] = fStartIndex[d] + sfSize[d] - 1;
            
        int fAbsStartIndex = calcIndex( fStartIndex ); // absolute field start index
        int fAbsEndIndex = calcIndex( fEndIndex ); // absolute subfield start index
        dab::Array<unsigned int> fCurRelIndex(fStartIndex); // running relative field index
        dab::Array<unsigned int> sfCurRelIndex(fDim); // running relative subfield index
                        
        int fAbsIndex = fAbsStartIndex; // absolute field index
        int sfAbsIndex = 0; // absolute subfield index
                        
        while( fAbsIndex < fAbsEndIndex )
        {
            // first dimension
            for(fCurRelIndex[0] = fStartIndex[0], sfCurRelIndex[0] = 0; fCurRelIndex[0] <= fEndIndex[0]; fCurRelIndex[0] += 1, sfCurRelIndex[0] += 1, fAbsIndex++, sfAbsIndex++)
            {
                //std::cout << "fCurRelIndex: " << fCurRelIndex << " sfCurRelIndex " << sfCurRelIndex << " fAbsIndex " << fAbsIndex << " sfAbsIndex " << sfAbsIndex << "\n";
                
                sfVectors[sfAbsIndex] = mVectors[fAbsIndex];
            }
            
            if(fAbsIndex >= fAbsEndIndex) break;
            
            // successive dimensions
            for(unsigned int d = 0; d < fDim, fCurRelIndex[d] > fEndIndex[d]; ++d)
            {
                fCurRelIndex[d] = fStartIndex[d];
                fCurRelIndex[d + 1] += 1;
                fAbsIndex -= mIndexOffset[d] * ( fEndIndex[d] - fStartIndex[d] + 1);
                fAbsIndex += mIndexOffset[d + 1];
                
                sfCurRelIndex[d] = 0;
                sfCurRelIndex[d + 1] += 1;
                sfAbsIndex -= sfIndexOffset[d] * ( sfSize[d] );
                sfAbsIndex += sfIndexOffset[d + 1];
            }
        }
    }
    catch(Exception& e)
    {
        e += Exception("MATH ERROR: failed to create vector sub field", __FILE__, __FUNCTION__, __LINE__);
        throw e;
    }
}

template<typename ValueType>
const Eigen::Matrix<ValueType, Eigen::Dynamic, 1>&
VectorField<ValueType>::operator[] (unsigned int pIndex) const throw (Exception)
{
    if(pIndex >= mVectorCount) throw Exception( "MATH ERROR: index out of bounds: max index " + std::to_string(mVectorCount - 1) + " supplied index " + std::to_string(pIndex), __FILE__, __FUNCTION__, __LINE__ );
        
    return mVectors[pIndex];
}
    
template<typename ValueType>
const Eigen::Matrix<ValueType, Eigen::Dynamic, 1>&
VectorField<ValueType>::get(unsigned int pIndex) const throw (Exception)
{
    if(pIndex >= mVectorCount) throw Exception( "MATH ERROR: index out of bounds: max index " + std::to_string(mVectorCount - 1) + " supplied index " + std::to_string(pIndex), __FILE__, __FUNCTION__, __LINE__ );
        
        return mVectors[pIndex];
}

template<typename ValueType>
const Eigen::Matrix<ValueType, Eigen::Dynamic, 1>&
VectorField<ValueType>::get(const dab::Array<unsigned int>& pIndex) const throw (Exception)
{
    try
    {
        return mVectors[calcIndex(pIndex)];
    }
    catch( Exception& e)
    {
        std::string indexString = pIndex;
        e += Exception("MATH ERROR: failed to get vector for index " + indexString, __FILE__, __FUNCTION__, __LINE__);
        throw e;
    }
}

template<typename ValueType>
Eigen::Matrix<ValueType, Eigen::Dynamic, 1>
VectorField<ValueType>::get(const dab::Array<float>& pIndex) const throw (Exception)
{
    if(pIndex.size() != mFieldDim) throw Exception("MATH ERROR: index dimension " + std::to_string(pIndex.size()) + " doesn't match field dimension " + std::to_string(mFieldDim), __FILE__, __FUNCTION__, __LINE__);
    
    Eigen::Matrix<ValueType, Eigen::Dynamic, 1> result(mVectorDim);
    
    try
    {
        get(pIndex,result);
        return result;
    }
    catch(dab::Exception& e)
    {
        throw;
    }
}
    
template<typename ValueType>
void
VectorField<ValueType>::get(const dab::Array<float>& pIndex, Eigen::Matrix<ValueType, Eigen::Dynamic, 1>& pResult) const throw (Exception)
{
    if(pIndex.size() != mFieldDim) throw Exception("MATH ERROR: index dimension " + std::to_string(pIndex.size()) + " doesn't match field dimension " + std::to_string(mFieldDim), __FILE__, __FUNCTION__, __LINE__);
    if(pResult.rows() != mVectorDim) throw Exception("MATH ERROR: result vector dimension " + std::to_string(pResult.rows()) + " doesn't match field vector dimension " + std::to_string(mVectorDim), __FILE__, __FUNCTION__, __LINE__);
    
    unsigned int fDim = mFieldDim;
    unsigned int vDim = mVectorDim;
    
    for(int d=0; d<mFieldDim; ++d)
    {
        if(mSize[d] <= pIndex[d] + 1.0) throw Exception( "MATH ERROR: for dimension " + std::to_string(d) + " index exceeds field size: fieldsize " + std::to_string(mSize[d]) + " supplied " + std::to_string(pIndex[d]), __FILE__, __FUNCTION__, __LINE__ );
        if( pIndex[d] < 0.0 ) throw Exception( "MATH ERROR: for dimension " + std::to_string(d) + " index smaller than zero", __FILE__, __FUNCTION__, __LINE__ );
    }
    
    try
    {
        dab::Array<unsigned int> startIndex(mFieldDim);
		std::vector<unsigned int> localIndex(fDim);
        std::vector<float> ratios(fDim * 2);
        float ratio;
        
        for(unsigned int d=0; d<fDim; ++d)
        {
            startIndex[d] = static_cast<unsigned int>( pIndex[d] );
            ratios[d] = 1.0 - ( pIndex[d] - floor( pIndex[d] ) );
            if(ratios[d] < 0) ratios[d] += 1.0;
            ratios[d + fDim] = 1.0 - ratios[d];
            localIndex[d] = 0;
        }
        
        for(unsigned int d=0; d<vDim; ++d)
        {
            pResult[d] = 0.0;
        }
        
        dab::Array<unsigned int> endIndex(mFieldDim);
        for(int d=0; d<mFieldDim; ++d) endIndex[d] = startIndex[d] + 1;
        
        int absStartIndex = calcIndex( startIndex );
        int absEndIndex = calcIndex( endIndex );
        
        if(absEndIndex >= mVectorCount) return;
        
        int absIndex = absStartIndex;
        
        while( absIndex < absEndIndex )
        {
            // first dimension
            for(localIndex[0] = 0; localIndex[0] < 2; localIndex[0] += 1, absIndex++)
            {
                ratio = 1.0;
                for(unsigned int d=0; d<fDim; ++d) ratio *= ratios[ d + localIndex[d] * fDim ];
                pResult += mVectors[absIndex] * ratio;
            }
            
            if(absIndex >= absEndIndex) break;
            
            // successive dimensions
            for(unsigned int d = 0; d < fDim, localIndex[d] > 1; ++d)
            {
                localIndex[d] = 0;
                localIndex[d + 1] += 1;
                absIndex -= 2 * mIndexOffset[d];
                absIndex += mIndexOffset[d + 1];
            }
        }
    }
    catch(Exception& e)
    {
        std::string indexString = pIndex;
        e += Exception("MATH ERROR: failed to get interpolated vector at index " + indexString, __FILE__, __FUNCTION__, __LINE__);
        throw e;
    }
}
    
template<typename ValueType>
void
VectorField<ValueType>::set(unsigned int pIndex, const Eigen::Matrix<ValueType, Eigen::Dynamic, 1>& pValue) throw (Exception)
{
    if(pValue.rows() != mVectorDim) throw Exception("MATH ERROR: vector dimension " + std::to_string(pValue.rows()) + " doesn't match field vector dimension " + std::to_string(mVectorDim), __FILE__, __FUNCTION__, __LINE__);
    if(pIndex >= mVectorCount) throw Exception( "MATH ERROR: index " + std::to_string(pIndex) + " exceeds vector count " + std::to_string(mVectorCount), __FILE__, __FUNCTION__, __LINE__ );
 
    mVectors[pIndex] = pValue;
}

template<typename ValueType>
void
VectorField<ValueType>::set(const dab::Array<unsigned int>& pIndex, const Eigen::Matrix<ValueType, Eigen::Dynamic, 1>& pValue) throw (Exception)
{
    if(pIndex.size() != mFieldDim) throw Exception("MATH ERROR: index dimension " + std::to_string(pIndex.size()) + " doesn't match field dimension " + std::to_string(mFieldDim), __FILE__, __FUNCTION__, __LINE__);
    if(pValue.rows() != mVectorDim) throw Exception("MATH ERROR: vector dimension " + std::to_string(pValue.rows()) + " doesn't match field vector dimension " + std::to_string(mVectorDim), __FILE__, __FUNCTION__, __LINE__);
    
    try
    {
        mVectors[calcIndex(pIndex)] = pValue;
    }
    catch(Exception& e)
    {
        std::string indexString = pIndex;
        e += Exception("MATH ERROR: failed to set vector at index " + indexString, __FILE__, __FUNCTION__, __LINE__);
        throw e;
    }
}
    
template<typename ValueType>
void
VectorField<ValueType>::set(const dab::Array<float>& pIndex, const Eigen::Matrix<ValueType, Eigen::Dynamic, 1>& pValue) throw (Exception)
{
    if(pIndex.size() != mFieldDim) throw Exception("MATH ERROR: index dimension " + std::to_string(pIndex.size()) + " doesn't match field dimension " + std::to_string(mFieldDim), __FILE__, __FUNCTION__, __LINE__);
    if(pValue.rows() != mVectorDim) throw Exception("MATH ERROR: vector dimension " + std::to_string(pValue.rows()) + " doesn't match field vector dimension " + std::to_string(mVectorDim), __FILE__, __FUNCTION__, __LINE__);
        
    unsigned int fDim = mFieldDim;
    unsigned int vDim = mVectorDim;
    
    for(int d=0; d<fDim; ++d)
    {
        if(mSize[d] <= pIndex[d] ) throw Exception( "MATH ERROR: at dimension " + std::to_string(d) + " index " + std::to_string(mSize[d]) + " exceeds field size " + std::to_string(pIndex[d]), __FILE__, __FUNCTION__, __LINE__ );
    }
                
    try
    {
        dab::Array<unsigned int> startIndex(fDim);
        dab::Array<unsigned int> localIndex(fDim);
        dab::Array<float> ratios(fDim * 2);
        float ratio;
        
        for(unsigned int d=0; d<fDim; ++d)
        {
            localIndex[d] = 0;
            startIndex[d] = static_cast<unsigned int>( pIndex[d] );
            ratios[d] = 1.0 - ( pIndex[d] - floor( pIndex[d] ) );
            if(ratios[d] < 0) ratios[d] += 1.0;
            ratios[d + fDim] = 1.0 - ratios[d];
        }
        
        dab::Array<unsigned int> endIndex(fDim);
        for(unsigned int d=0; d<fDim; ++d) endIndex[d] = startIndex[d] + 1;
        
        int absStartIndex = calcIndex( startIndex );
        int absEndIndex = calcIndex( endIndex );
            
//        std::cout << "start index ";
//        for(int d=0; d<fDim; ++d) std::cout << startIndex[d] << " ";
//        std::cout << "\n";
//        std::cout << "end index ";
//        for(int d=0; d<fDim; ++d) std::cout << endIndex[d] << " ";
//        std::cout << "\n";
//        std::cout << "absStartIndex " << absStartIndex << "\n";
//        std::cout << "absEndIndex " << absEndIndex << "\n";
        
        if(absEndIndex >= mVectorCount) return;
        
        int absIndex = absStartIndex;
        
        float ratioSum = 0;
        
        while( absIndex < absEndIndex )
        {
            // first dimension
            for(localIndex[0] = 0; localIndex[0] < 2; localIndex[0] += 1, absIndex++)
            {
                ratio = 1.0;
                for(unsigned int d=0; d<fDim; ++d) ratio *= ratios[ d + localIndex[d] * fDim ];
                    
//                std::cout << "absIndex " << absIndex;
//                std::cout << " localIndex ";
//                for(int d=0; d<fDim; ++d)
//                {
//                    std::cout << localIndex[d] << " ";
//                }
//                std::cout << "ratio " << ratio << "\n";
                    
                mVectors[absIndex] = pValue * ratio;
                    
                ratioSum += ratio;
            }
                    
            if(absIndex >= absEndIndex) break;
                    
            // successive dimensions
            for(unsigned int d = 0; d < fDim, localIndex[d] > 1; ++d)
            {
                localIndex[d] = 0;
                localIndex[d + 1] += 1;
                absIndex -= 2 * mIndexOffset[d];
                absIndex += mIndexOffset[d + 1];
            }
        }
                    
        //std::cout << ">> ratioSum " << ratioSum << "\n";
    }
    catch(Exception& e)
    {
        std::string indexString = pIndex;
        e += Exception("MATH ERROR: failed to set vector value at index " + indexString, __FILE__, __FUNCTION__, __LINE__);
        throw e;
    }
}

template<typename ValueType>
void
    VectorField<ValueType>::set(const dab::Array<unsigned int>& pStartIndex, const dab::Array<unsigned int>& pEndIndex, const Eigen::Matrix<ValueType, Eigen::Dynamic, 1>& pValue) throw (Exception)
{
    if(pStartIndex.size() != mFieldDim) throw Exception("MATH ERROR: start index dimension " + std::to_string(pStartIndex.size()) + " doesn't match field dimension " + std::to_string(mFieldDim), __FILE__, __FUNCTION__, __LINE__);
    if(pEndIndex.size() != mFieldDim) throw Exception("MATH ERROR: end index dimension " + std::to_string(pEndIndex.size()) + " doesn't match field dimension " + std::to_string(mFieldDim), __FILE__, __FUNCTION__, __LINE__);
    if(pValue.rows() != mVectorDim) throw Exception("MATH ERROR: vector dimension " + std::to_string(pValue.rows()) + " doesn't match field vector dimension " + std::to_string(mVectorDim), __FILE__, __FUNCTION__, __LINE__);
    
    
    unsigned int fDim = mFieldDim;
    unsigned int vDim = mVectorDim;
    
    for(int d=0; d<fDim; ++d)
    {
        if(pStartIndex[d] > pEndIndex[d]) throw Exception( "MATH ERROR: for dimension " + std::to_string(d) + " start index " + std::to_string(pStartIndex[d]) + " is larger than end index " + std::to_string( pEndIndex[d]), __FILE__, __FUNCTION__, __LINE__ );
        if(pEndIndex[d] >= mSize[d]) throw Exception( "MATH ERROR: for dimension " + std::to_string(d) + " end index " + std::to_string(pEndIndex[d]) + " is larger than field size " + std::to_string( mSize[d]), __FILE__, __FUNCTION__, __LINE__ );
    }
    
    try
    {
        int absStartIndex = calcIndex( pStartIndex );
        int absEndIndex = calcIndex( pEndIndex );
        dab::Array<unsigned int> curRelIndex = pStartIndex;
        int absIndex = absStartIndex;
        
        while( absIndex < absEndIndex )
        {
            // first dimension
            for(curRelIndex[0] = pStartIndex[0]; curRelIndex[0] <= pEndIndex[0]; curRelIndex[0] += 1, absIndex++)
            {
//                std::cout << "curRelIndex: [";
//                for(int i=0; i<fDim; ++i) std::cout << curRelIndex[i] << " ";
//                std::cout << " ] ";
//                std::cout << "absIndex: " << absIndex << "\n";
                
                mVectors[absIndex] = pValue;
            }
            
            if(absIndex >= absEndIndex) break;
            
            // successive dimensions
            for(unsigned int d = 0; d < fDim, curRelIndex[d] > pEndIndex[d]; ++d)
            {
                curRelIndex[d] = pStartIndex[d];
                curRelIndex[d + 1] += 1;
                absIndex -= mIndexOffset[d] * ( pEndIndex[d] - pStartIndex[d] + 1);
                absIndex += mIndexOffset[d + 1];
            }
        }
    }
    catch(Exception& e)
    {
        std::string startIndexString = pStartIndex;
        std::string endIndexString = pEndIndex;
        
        std::string valueString;
        for(int d=0; d<vDim; ++d)
        {
            valueString += std::to_string(pValue[d]) + " ";
        }
        
        e += Exception("MATH ERROR: failed to set vectors between start index " + startIndexString + " and end index " + endIndexString + " to value " + valueString, __FILE__, __FUNCTION__, __LINE__);
        throw e;
    }
}

template<typename ValueType>
void
VectorField<ValueType>::set(const Eigen::Matrix<ValueType, Eigen::Dynamic, 1>& pValue) throw (Exception)
{
    if(pValue.rows() != mVectorDim) throw Exception("MATH ERROR: vector dimension " + std::to_string(pValue.rows()) + " doesn't match field vector dimension " + std::to_string(mVectorDim), __FILE__, __FUNCTION__, __LINE__);
    
    for(unsigned int i=0; i<mVectorCount; ++i) mVectors[i] = pValue;
}
                
template<typename ValueType>
void
VectorField<ValueType>::add(unsigned int pIndex, const Eigen::Matrix<ValueType, Eigen::Dynamic, 1>& pValue) throw (Exception)
{

    if(pIndex >= mVectorCount) throw Exception( "MATH ERROR: index " + std::to_string(pIndex) + " out of bounds " + std::to_string(mVectorCount), __FILE__, __FUNCTION__, __LINE__ );
    if(pValue.rows() != mVectorDim) throw Exception("MATH ERROR: vector dimension " + std::to_string(pValue.rows()) + " doesn't match field vector dimension " + std::to_string(mVectorDim), __FILE__, __FUNCTION__, __LINE__);

    for(int d=0; d<mVectorDim; ++d) mVectors[pIndex][d] += pValue[d];
}
    
template<typename ValueType>
void
VectorField<ValueType>::add(const dab::Array<unsigned int>& pIndex, const Eigen::Matrix<ValueType, Eigen::Dynamic, 1>& pValue) throw (Exception)
{
    if(pValue.rows() != mVectorDim) throw Exception("MATH ERROR: vector dimension " + std::to_string(pValue.rows()) + " doesn't match field vector dimension " + std::to_string(mVectorDim), __FILE__, __FUNCTION__, __LINE__);
    
    try
    {
        mVectors[calcIndex(pIndex)] += pValue;
    }
    catch(Exception& e)
    {
        std::string indexString = pIndex;
        std::string valueString;
        for(int d=0; d<mVectorDim; ++d) valueString += std::to_string(pValue[d]) + " ";
        
        e += Exception("MATH ERROR: failed to add value " + valueString + " at index " + indexString, __FILE__, __FUNCTION__, __LINE__);
        throw e;
    }
}
                
template<typename ValueType>
void
VectorField<ValueType>::add(const dab::Array<float>& pIndex, const Eigen::Matrix<ValueType, Eigen::Dynamic, 1>& pValue) throw (Exception)
{
    if(pIndex.size() != mFieldDim) throw Exception("MATH ERROR: index dimension " + std::to_string(pIndex.size()) + " doesn't match field dimension " + std::to_string(mFieldDim), __FILE__, __FUNCTION__, __LINE__);
    if(pValue.rows() != mVectorDim) throw Exception("MATH ERROR: vector dimension " + std::to_string(pValue.rows()) + " doesn't match field vector dimension " + std::to_string(mVectorDim), __FILE__, __FUNCTION__, __LINE__);
        
    unsigned int fDim = mFieldDim;
    unsigned int vDim = mVectorDim;

    try
    {
        dab::Array<unsigned int> startIndex(fDim);
        dab::Array<unsigned int> localIndex(fDim);
        dab::Array<float> ratios(fDim*2);
        float ratio;
        
        for(unsigned int d=0; d<fDim; ++d)
        {
            localIndex[d] = 0;
            startIndex[d] = static_cast<unsigned int>( pIndex[d] );
            ratios[d] = 1.0 - ( pIndex[d] - floor( pIndex[d] ) );
            if(ratios[d] < 0) ratios[d] += 1.0;
            ratios[d + fDim] = 1.0 - ratios[d];
        }
        
        dab::Array<unsigned int> endIndex(fDim);
        for(unsigned int d=0; d<fDim; ++d)
        {
            endIndex[d] = startIndex[d] + 1;
        }
        
        int absStartIndex = calcIndex( startIndex );
        int absEndIndex = calcIndex( endIndex );
        
        if(absEndIndex >= mVectorCount) return;
        
        int absIndex = absStartIndex;
        
        
        while( absIndex < absEndIndex )
        {
            // first dimension
            for(localIndex[0] = 0; localIndex[0] < 2; localIndex[0] += 1, absIndex++)
            {
                ratio = 1.0;
                for(unsigned int d=0; d<fDim; ++d) ratio *= ratios[ d + localIndex[d] * fDim ];
                mVectors[absIndex] += pValue * ratio;
            }
            
            if(absIndex >= absEndIndex) break;
            
            // successive dimensions
            for(unsigned int d = 0; d < fDim, localIndex[d] > 1; ++d)
            {
                localIndex[d] = 0;
                localIndex[d + 1] += 1;
                absIndex -= 2 * mIndexOffset[d];
                absIndex += mIndexOffset[d + 1];
            }
        }
    }
    catch(Exception& e)
    {
        std::string indexString = pIndex;
        std::string valueString;
        for(int d=0; d<mVectorDim; ++d) valueString += std::to_string(pValue[d]) + " ";
        
        e += Exception("MATH ERROR: failed to add value " + valueString + " to index " + indexString, __FILE__, __FUNCTION__, __LINE__);
        throw e;
    }
}
    
template<typename ValueType>
void
VectorField<ValueType>::add(const dab::Array<unsigned int>& pStartIndex, const dab::Array<unsigned int>& pEndIndex, const Eigen::Matrix<ValueType, Eigen::Dynamic, 1>& pValue) throw (Exception)
{
    if(pStartIndex.size() != mFieldDim) throw Exception("MATH ERROR: start index dimension " + std::to_string(pStartIndex.size()) + " doesn't match field dimension " + std::to_string(mFieldDim), __FILE__, __FUNCTION__, __LINE__);
    if(pEndIndex.size() != mFieldDim) throw Exception("MATH ERROR: end index dimension " + std::to_string(pEndIndex.size()) + " doesn't match field dimension " + std::to_string(mFieldDim), __FILE__, __FUNCTION__, __LINE__);
    if(pValue.rows() != mVectorDim) throw Exception("MATH ERROR: vector dimension " + std::to_string(pValue.rows()) + " doesn't match field vector dimension " + std::to_string(mVectorDim), __FILE__, __FUNCTION__, __LINE__);
    
    unsigned int fDim = mFieldDim;
    unsigned int vDim = mVectorDim;
    
    for(int d=0; d<fDim; ++d)
    {
        if(pStartIndex[d] > pEndIndex[d]) throw Exception( "MATH ERROR: for dimension " + std::to_string(d) + " start index " + std::to_string(pStartIndex[d]) + " is larger than end index " + std::to_string( pEndIndex[d]), __FILE__, __FUNCTION__, __LINE__ );
        if(pEndIndex[d] >= mSize[d]) throw Exception( "MATH ERROR: for dimension " + std::to_string(d) + " end index " + std::to_string(pEndIndex[d]) + " is larger than field size " + std::to_string( mSize[d]), __FILE__, __FUNCTION__, __LINE__ );
    }

    try
    {
        int absStartIndex = calcIndex( pStartIndex );
        int absEndIndex = calcIndex( pEndIndex );
        dab::Array<unsigned int> curRelIndex = pStartIndex;
        int absIndex = absStartIndex;
        
        while( absIndex < absEndIndex )
        {
            // first dimension
            for(curRelIndex[0] = pStartIndex[0]; curRelIndex[0] <= pEndIndex[0]; curRelIndex[0] += 1, absIndex++)
            {
                mVectors[absIndex] += pValue;
            }
            
            if(absIndex >= absEndIndex) break;
            
            // successive dimensions
            for(unsigned int d = 0; d < fDim, curRelIndex[d] > pEndIndex[d]; ++d)
            {
                curRelIndex[d] = pStartIndex[d];
                curRelIndex[d + 1] += 1;
                absIndex -= mIndexOffset[d] * ( pEndIndex[d] - pStartIndex[d] + 1);
                absIndex += mIndexOffset[d + 1];
            }
        }
    }
    catch(Exception& e)
    {
        std::string startIndexString = pStartIndex;
        std::string endIndexString = pEndIndex;
        std::string valueString;
        for(int d=0; d<mVectorDim; ++d)
        {
            valueString += std::to_string(pValue[d]) + " ";
        }
        
        e += Exception("MATH ERROR: failed to add value " + valueString + " to vectors between start index " + startIndexString + " and end index " + endIndexString, __FILE__, __FUNCTION__, __LINE__);
        throw e;
    }
}

template<typename ValueType>
void
VectorField<ValueType>::add(const Eigen::Matrix<ValueType, Eigen::Dynamic, 1>& pValue) throw (Exception)
{
    if(pValue.rows() != mVectorDim) throw Exception("MATH ERROR: vector dimension " + std::to_string(pValue.rows()) + " doesn't match field vector dimension " + std::to_string(mVectorDim), __FILE__, __FUNCTION__, __LINE__);
    
    for(unsigned int i=0; i<mVectorCount; ++i) mVectors[i] += pValue;
}

template<typename ValueType>
const VectorField<ValueType>&
    VectorField<ValueType>::operator+=( const Eigen::Matrix<ValueType, Eigen::Dynamic, 1>& pValue) throw (Exception)
{
    if(pValue.rows() != mVectorDim) throw Exception("MATH ERROR: vector dimension " + std::to_string(pValue.rows()) + " doesn't match field vector dimension " + std::to_string(mVectorDim), __FILE__, __FUNCTION__, __LINE__);
    
    for(unsigned int vI=0; vI<mVectorCount; ++vI) mVectors[vI] += pValue;
    return(*this);
}

template<typename ValueType>
const VectorField<ValueType>&
VectorField<ValueType>::operator-=( const Eigen::Matrix<ValueType, Eigen::Dynamic, 1>& pValue) throw (Exception)
{
    if(pValue.rows() != mVectorDim) throw Exception("MATH ERROR: vector dimension " + std::to_string(pValue.rows()) + " doesn't match field vector dimension " + std::to_string(mVectorDim), __FILE__, __FUNCTION__, __LINE__);
    
    for(unsigned int vI=0; vI<mVectorCount; ++vI) mVectors[vI] -= pValue;
    return(*this);
}
    
template<typename ValueType>
const VectorField<ValueType>&
VectorField<ValueType>::operator*=( const Eigen::Matrix<ValueType, Eigen::Dynamic, 1>& pValue ) throw (Exception)
{
    if(pValue.rows() != mVectorDim) throw Exception("MATH ERROR: vector dimension " + std::to_string(pValue.rows()) + " doesn't match field vector dimension " + std::to_string(mVectorDim), __FILE__, __FUNCTION__, __LINE__);
    
    for(unsigned int vI=0; vI<mVectorCount; ++vI)
    {
        for(unsigned int d=0; d<mVectorDim; ++d) mVectors[vI][d] *= pValue[d];
    }
    return(*this);
}
    
template<typename ValueType>
const VectorField<ValueType>&
VectorField<ValueType>::operator/=( const Eigen::Matrix<ValueType, Eigen::Dynamic, 1>& pValue ) throw (Exception)
{
    if(pValue.rows() != mVectorDim) throw Exception("MATH ERROR: vector dimension " + std::to_string(pValue.rows()) + " doesn't match field vector dimension " + std::to_string(mVectorDim), __FILE__, __FUNCTION__, __LINE__);
    
    for(unsigned int vI=0; vI<mVectorCount; ++vI)
    {
        for(unsigned int d=0; d<mVectorDim; ++d) mVectors[vI][d] /= pValue[d];
    }
    return(*this);
}
    
template<typename ValueType>
VectorField<ValueType>
VectorField<ValueType>::operator+( const Eigen::Matrix<ValueType, Eigen::Dynamic, 1>& pValue ) throw (Exception)
{
    if(pValue.rows() != mVectorDim) throw Exception("MATH ERROR: vector dimension " + std::to_string(pValue.rows()) + " doesn't match field vector dimension " + std::to_string(mVectorDim), __FILE__, __FUNCTION__, __LINE__);

    VectorField<ValueType> field(*this);
    field += pValue;
    return field;
}
    
template<typename ValueType>
VectorField<ValueType>
VectorField<ValueType>::operator-( const Eigen::Matrix<ValueType, Eigen::Dynamic, 1>& pValue ) throw (Exception)
{
    if(pValue.rows() != mVectorDim) throw Exception("MATH ERROR: vector dimension " + std::to_string(pValue.rows()) + " doesn't match field vector dimension " + std::to_string(mVectorDim), __FILE__, __FUNCTION__, __LINE__);
    
    VectorField<ValueType> field(*this);
    field -= pValue;
    return field;
}
    
template<typename ValueType>
VectorField<ValueType>
VectorField<ValueType>::operator*( const Eigen::Matrix<ValueType, Eigen::Dynamic, 1>& pValue ) throw (Exception)
{
    if(pValue.rows() != mVectorDim) throw Exception("MATH ERROR: vector dimension " + std::to_string(pValue.rows()) + " doesn't match field vector dimension " + std::to_string(mVectorDim), __FILE__, __FUNCTION__, __LINE__);
    
    VectorField<ValueType> field(*this);
    field *= pValue;
    return field;
}

template<typename ValueType>
VectorField<ValueType>
VectorField<ValueType>::operator/( const Eigen::Matrix<ValueType, Eigen::Dynamic, 1>& pValue ) throw (Exception)
{
    if(pValue.rows() != mVectorDim) throw Exception("MATH ERROR: vector dimension " + std::to_string(pValue.rows()) + " doesn't match field vector dimension " + std::to_string(mVectorDim), __FILE__, __FUNCTION__, __LINE__);
    
    VectorField<ValueType> field(*this);
    field /= pValue;
    return field;
}

template<typename ValueType>
const VectorField<ValueType>&
VectorField<ValueType>::operator+=( const VectorField<ValueType>& pVectorField ) throw (Exception)
{
    if(pVectorField.mFieldDim != mFieldDim) throw Exception("MATH ERROR: supplied field dimension " + std::to_string(pVectorField.mFieldDim) + " doesn't match field dimension " + std::to_string(mFieldDim), __FILE__, __FUNCTION__, __LINE__);
    if(pVectorField.mVectorDim != mVectorDim) throw Exception("MATH ERROR: supplied vector dimension " + std::to_string(pVectorField.mVectorDim) + " doesn't match field vector dimension " + std::to_string(mVectorDim), __FILE__, __FUNCTION__, __LINE__);

    const std::vector< Eigen::Matrix<ValueType, Eigen::Dynamic, 1> >& pVectors = pVectorField.vectors();
    const dab::Array<unsigned int>& pSize = pVectorField.size();
    
    for(int d=0; d<mFieldDim; ++d)
    {
        if( mSize[d] != pSize[d] ) throw Exception("MATH ERROR: field size mismatch for dimension " + std::to_string(d) + " : provided " + std::to_string(pSize[d]) + " expected " + std::to_string(mSize[d]), __FILE__, __FUNCTION__, __LINE__ );
    }
    
    for(unsigned int vI=0; vI<mVectorCount; ++vI) mVectors[vI] += pVectors[vI];
    return(*this);
}

template<typename ValueType>
const VectorField<ValueType>&
VectorField<ValueType>::operator-=( const VectorField<ValueType>& pVectorField ) throw (Exception)
{
    if(pVectorField.mFieldDim != mFieldDim) throw Exception("MATH ERROR: supplied field dimension " + std::to_string(pVectorField.mFieldDim) + " doesn't match field dimension " + std::to_string(mFieldDim), __FILE__, __FUNCTION__, __LINE__);
    if(pVectorField.mVectorDim != mVectorDim) throw Exception("MATH ERROR: supplied vector dimension " + std::to_string(pVectorField.mVectorDim) + " doesn't match field vector dimension " + std::to_string(mVectorDim), __FILE__, __FUNCTION__, __LINE__);
            
    const std::vector< Eigen::Matrix<ValueType, Eigen::Dynamic, 1> >& pVectors = pVectorField.vectors();
    const dab::Array<unsigned int>& pSize = pVectorField.size();
            
    for(int d=0; d<mFieldDim; ++d)
    {
        if( mSize[d] != pSize[d] ) throw Exception("MATH ERROR: field size mismatch for dimension " + std::to_string(d) + " : provided " + std::to_string(pSize[d]) + " expected " + std::to_string(mSize[d]), __FILE__, __FUNCTION__, __LINE__ );
    }
    
    for(unsigned int vI=0; vI<mVectorCount; ++vI) mVectors[vI] -= pVectors[vI];
        return(*this);
}

template<typename ValueType>
const VectorField<ValueType>&
VectorField<ValueType>::operator*=( const VectorField<ValueType>& pVectorField ) throw (Exception)
{
    if(pVectorField.mFieldDim != mFieldDim) throw Exception("MATH ERROR: supplied field dimension " + std::to_string(pVectorField.mFieldDim) + " doesn't match field dimension " + std::to_string(mFieldDim), __FILE__, __FUNCTION__, __LINE__);
    if(pVectorField.mVectorDim != mVectorDim) throw Exception("MATH ERROR: supplied vector dimension " + std::to_string(pVectorField.mVectorDim) + " doesn't match field vector dimension " + std::to_string(mVectorDim), __FILE__, __FUNCTION__, __LINE__);

    
    const std::vector< Eigen::Matrix<ValueType, Eigen::Dynamic, 1> >& pVectors = pVectorField.vectors();
    const dab::Array<unsigned int>& pSize = pVectorField.size();
    
    for(int d=0; d<mFieldDim; ++d)
    {
        if( mSize[d] != pSize[d] ) throw Exception("MATH ERROR: field size mismatch for dimension " + std::to_string(d) + " : provided " + std::to_string(pSize[d]) + " expected " + std::to_string(mSize[d]), __FILE__, __FUNCTION__, __LINE__ );
    }
    
    for(unsigned int vI=0; vI<mVectorCount; ++vI)
    {
        for(int d=0; d<mVectorDim; ++d)
        {
            mVectors[vI][d] *= pVectors[vI][d];
        }
    }
    
    return(*this);
}
    
template<typename ValueType>
const VectorField<ValueType>&
VectorField<ValueType>::operator/=( const VectorField<ValueType>& pVectorField ) throw (Exception)
{
    if(pVectorField.mFieldDim != mFieldDim) throw Exception("MATH ERROR: supplied field dimension " + std::to_string(pVectorField.mFieldDim) + " doesn't match field dimension " + std::to_string(mFieldDim), __FILE__, __FUNCTION__, __LINE__);
    if(pVectorField.mVectorDim != mVectorDim) throw Exception("MATH ERROR: supplied vector dimension " + std::to_string(pVectorField.mVectorDim) + " doesn't match field vector dimension " + std::to_string(mVectorDim), __FILE__, __FUNCTION__, __LINE__);
            
            
    const std::vector< Eigen::Matrix<ValueType, Eigen::Dynamic, 1> >& pVectors = pVectorField.vectors();
    const dab::Array<unsigned int>& pSize = pVectorField.size();
            
    for(int d=0; d<mFieldDim; ++d)
    {
        if( mSize[d] != pSize[d] ) throw Exception("MATH ERROR: field size mismatch for dimension " + std::to_string(d) + " : provided " + std::to_string(pSize[d]) + " expected " + std::to_string(mSize[d]), __FILE__, __FUNCTION__, __LINE__ );
    }
    
    for(unsigned int vI=0; vI<mVectorCount; ++vI)
    {
        for(int d=0; d<mVectorDim; ++d)
        {
            mVectors[vI][d] /= pVectors[vI][d];
        }
    }
    
    return(*this);
}
    
template<typename ValueType>
VectorField<ValueType>
VectorField<ValueType>::operator+( const VectorField<ValueType>& pVectorField ) throw (Exception)
{
    if(pVectorField.mFieldDim != mFieldDim) throw Exception("MATH ERROR: supplied field dimension " + std::to_string(pVectorField.mFieldDim) + " doesn't match field dimension " + std::to_string(mFieldDim), __FILE__, __FUNCTION__, __LINE__);
    if(pVectorField.mVectorDim != mVectorDim) throw Exception("MATH ERROR: supplied vector dimension " + std::to_string(pVectorField.mVectorDim) + " doesn't match field vector dimension " + std::to_string(mVectorDim), __FILE__, __FUNCTION__, __LINE__);
    
    try
    {
        VectorField<ValueType> field(*this);
        field += pVectorField;
        return field;
    }
    catch(Exception& e)
    {
        throw;
    }
}
    
template<typename ValueType>
VectorField<ValueType>
VectorField<ValueType>::operator-( const VectorField<ValueType>& pVectorField ) throw (Exception)
{
    if(pVectorField.mFieldDim != mFieldDim) throw Exception("MATH ERROR: supplied field dimension " + std::to_string(pVectorField.mFieldDim) + " doesn't match field dimension " + std::to_string(mFieldDim), __FILE__, __FUNCTION__, __LINE__);
    if(pVectorField.mVectorDim != mVectorDim) throw Exception("MATH ERROR: supplied vector dimension " + std::to_string(pVectorField.mVectorDim) + " doesn't match field vector dimension " + std::to_string(mVectorDim), __FILE__, __FUNCTION__, __LINE__);
    
    try
    {
        VectorField<ValueType> field(*this);
        field -= pVectorField;
        return field;
    }
    catch(Exception& e)
    {
        throw;
    }
}
    
template<typename ValueType>
VectorField<ValueType>
VectorField<ValueType>::operator*( const VectorField<ValueType>& pVectorField ) throw (Exception)
{
    if(pVectorField.mFieldDim != mFieldDim) throw Exception("MATH ERROR: supplied field dimension " + std::to_string(pVectorField.mFieldDim) + " doesn't match field dimension " + std::to_string(mFieldDim), __FILE__, __FUNCTION__, __LINE__);
    if(pVectorField.mVectorDim != mVectorDim) throw Exception("MATH ERROR: supplied vector dimension " + std::to_string(pVectorField.mVectorDim) + " doesn't match field vector dimension " + std::to_string(mVectorDim), __FILE__, __FUNCTION__, __LINE__);
    
    try
    {
        VectorField<ValueType> field(*this);
        field *= pVectorField;
        return field;
    }
    catch(Exception& e)
    {
        throw;
    }
}

template<typename ValueType>
VectorField<ValueType>
VectorField<ValueType>::operator/( const VectorField<ValueType>& pVectorField ) throw (Exception)
{
    if(pVectorField.mFieldDim != mFieldDim) throw Exception("MATH ERROR: supplied field dimension " + std::to_string(pVectorField.mFieldDim) + " doesn't match field dimension " + std::to_string(mFieldDim), __FILE__, __FUNCTION__, __LINE__);
    if(pVectorField.mVectorDim != mVectorDim) throw Exception("MATH ERROR: supplied vector dimension " + std::to_string(pVectorField.mVectorDim) + " doesn't match field vector dimension " + std::to_string(mVectorDim), __FILE__, __FUNCTION__, __LINE__);
    
    try
    {
        VectorField<ValueType> field(*this);
        field /= pVectorField;
        return field;
    }
    catch(Exception& e)
    {
        throw;
    }
}
    
template<typename ValueType>
void
VectorField<ValueType>::convolve( const VectorField<ValueType>& pKernel ) throw (Exception)
{
    if(pKernel.mFieldDim != mFieldDim) throw Exception("MATH ERROR: kernel field dimension " + std::to_string(pKernel.mFieldDim) + " doesn't match field dimension " + std::to_string(mFieldDim), __FILE__, __FUNCTION__, __LINE__);
    if(pKernel.mVectorDim != mVectorDim) throw Exception("MATH ERROR: kernel vector dimension " + std::to_string(pKernel.mVectorDim) + " doesn't match field vector dimension " + std::to_string(mVectorDim), __FILE__, __FUNCTION__, __LINE__);

    unsigned int fDim = mFieldDim;
    unsigned int vDim = mVectorDim;
    const dab::Array<unsigned int>& kernelSize = pKernel.size();
    

    for(unsigned int d=0; d<fDim; ++d)
    {
        if(kernelSize[d] % 2 != 1) throw Exception( "MATH EXCEPTION: kernel size in dim " + std::to_string(d) + " not odd", __FILE__, __FUNCTION__, __LINE__ );
    }
    try
    {
        dab::Array<unsigned int> kernelCenter(fDim); // center position of kernel
        for(int d=0; d<fDim; ++d) kernelCenter[d] = kernelSize[d] / 2;
        
        VectorField<ValueType> backupField( mSize, mVectors[0] ); // backup of vector field
        VectorField<ValueType> movingField( pKernel.size(), mVectors[0] ); // field section that movies together with kernel
        std::vector< Eigen::Matrix<ValueType, Eigen::Dynamic, 1> >& backupVectors = backupField.vectors();
        dab::Array<unsigned int> startIndex(fDim);
        dab::Array<unsigned int> endIndex(fDim);
        int absIndexKernelCenterOffset = 0;
        for(unsigned int d=0; d<fDim; ++d) 
        {
            startIndex[d] = 0;
            endIndex[d] = static_cast<unsigned int>( static_cast<int>( mSize[d] ) - kernelSize[d] );
            absIndexKernelCenterOffset += mIndexOffset[d] * kernelCenter[d];
        }
        
        int absStartIndex = calcIndex( startIndex );
        int absEndIndex = calcIndex( endIndex );
        dab::Array<unsigned int> curRelIndex(startIndex);
        dab::Array<unsigned int> curRelCenterIndex(startIndex);
        int absIndex = absStartIndex;
        
        while( absIndex < absEndIndex )
        {
            // first dimension
            for(curRelIndex[0] = startIndex[0]; curRelIndex[0] <= endIndex[0]; curRelIndex[0] += 1, absIndex++)
            {
                //std::cout << "curRelIndex: " << curRelIndex << " absIndex: " << absIndex << " absCenterIndex " << absIndex + absIndexKernelCenterOffset << "\n";
                
                subField( curRelIndex, movingField );
                
                //std::cout << "movingField " << movingField << "\n";
                
                movingField *= pKernel;
                
                //std::cout << "movingField * Kernel " << movingField << "\n";
                
                backupVectors[absIndex + absIndexKernelCenterOffset] = movingField.sum();
                
                //std::cout << "sum " << movingField.sum() << "\n";
            }
            
            if(absIndex >= absEndIndex) break;
            
            // successive dimensions
            for(unsigned int d = 0; d < fDim, curRelIndex[d] > endIndex[d]; ++d)
            {
                curRelIndex[d] = startIndex[d];
                curRelIndex[d + 1] += 1;			
                absIndex -= mIndexOffset[d] * ( endIndex[d] - startIndex[d] + 1);
                absIndex += mIndexOffset[d + 1];
            }
        }
        
        for(unsigned int vI=0; vI<mVectorCount; ++vI) mVectors[vI] = backupVectors[vI];
    }
    catch(Exception& e)
    {
        std::string kernelString = pKernel;
        e += Exception("MATH ERROR: failed to comvolve with kernel " + kernelString, __FILE__, __FUNCTION__, __LINE__);
        throw e;
    }
}
    
template<typename ValueType>
VectorField<ValueType>::operator std::string() const
{
    std::stringstream stream;

    stream << "fieldDim: " << mFieldDim << "\n";
    stream << "vectorDim: " << mVectorDim << "\n";
    stream << "size: [ ";
    for(int d=0; d<mFieldDim; ++d) stream << mSize[d] << " ";
    stream << "]\n";
    
    stream << "vectors:\n";
    for(unsigned int i=0; i<mVectorCount; ++i)
    {
        stream << i << " : [ ";
        for(int d=0; d<mVectorDim; ++d) stream << mVectors[i][d] << " ";
        stream << "]\n";
    }

    return stream.str();
}

     

    
};

};