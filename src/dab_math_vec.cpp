/** \file dab_math_vec.cpp
*/

#include "dab_math_vec.h"

using namespace dab;
using namespace dab::math;

glm::vec3
VectorMath::vec3Mat4Mul(const glm::mat4x4& pMatrix, const glm::vec3& pVector)
{
	glm::vec4 tmp = pMatrix * glm::vec4(pVector, 1.0);
	return tmp / tmp[3];
}

glm::quat 
VectorMath::makeQuaternion(const glm::vec3& pFrom, const glm::vec3& pTo)
{
	// Code adapted from ofQuaternion.cpp

	// This routine takes any vector as argument but normalized
	// vectors are necessary, if only for computing the dot product.
	// Too bad the API is that generic, it leads to performance loss.
	// Even in the case the 2 vectors are not normalized but same length,
	// the sqrt could be shared, but we have no way to know beforehand
	// at this point, while the caller may know.
	// So, we have to test... in the hope of saving at least a sqrt

	glm::quat _quat;
	glm::vec3 sourceVector = pFrom;
	glm::vec3 targetVector = pTo;

	float fromLen2 = glm::length2(pFrom);
	float fromLen;
	// normalize only when necessary, epsilon test
	if ((fromLen2 < 1.0 - 1e-7) || (fromLen2 > 1.0 + 1e-7)) {
		fromLen = sqrt(fromLen2);
		sourceVector /= fromLen;
	}
	else fromLen = 1.0;

	float toLen2 = glm::length2(pTo);
	// normalize only when necessary, epsilon test
	if ((toLen2 < 1.0 - 1e-7) || (toLen2 > 1.0 + 1e-7)) {
		float toLen;
		// re-use fromLen for case of mapping 2 vectors of the same length
		if ((toLen2 > fromLen2 - 1e-7) && (toLen2 < fromLen2 + 1e-7)) {
			toLen = fromLen;
		}
		else toLen = sqrt(toLen2);
		targetVector /= toLen;
	}


	// Now let's get into the real stuff
	// Use "dot product plus one" as test as it can be re-used later on
	double dotProdPlus1 = 1.0 + glm::dot(sourceVector, targetVector);

	// Check for degenerate case of full u-turn. Use epsilon for detection
	if (dotProdPlus1 < 1e-7) {

		// Get an orthogonal vector of the given vector
		// in a plane with maximum vector coordinates.
		// Then use it as quaternion axis with pi angle
		// Trick is to realize one value at least is >0.6 for a normalized vector.
		if (fabs(sourceVector.x) < 0.6) {
			const double norm = sqrt(1.0 - sourceVector.x * sourceVector.x);
			_quat.x = 0.0;
			_quat.y = sourceVector.z / norm;
			_quat.z = -sourceVector.y / norm;
			_quat.w = 0.0;
		}
		else if (fabs(sourceVector.y) < 0.6) {
			const double norm = sqrt(1.0 - sourceVector.y * sourceVector.y);
			_quat.x = -sourceVector.z / norm;
			_quat.y = 0.0;
			_quat.z = sourceVector.x / norm;
			_quat.w = 0.0;
		}
		else {
			const double norm = sqrt(1.0 - sourceVector.z * sourceVector.z);
			_quat.x = sourceVector.y / norm;
			_quat.y = -sourceVector.x / norm;
			_quat.z = 0.0;
			_quat.w = 0.0;
		}
	}
	else {
		// Find the shortest angle quaternion that transforms normalized vectors
		// into one other. Formula is still valid when vectors are colinear
		const double s = sqrt(0.5 * dotProdPlus1);
		const glm::vec3 tmp = glm::cross(sourceVector, targetVector) / (2.0 * s);
		_quat.x = tmp.x;
		_quat.y = tmp.y;
		_quat.z = tmp.z;
		_quat.w = s;
	}

	return _quat;
}