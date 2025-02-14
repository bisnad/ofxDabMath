/** \file dab_math_vec.h
*/

#pragma once

#include "dab_singleton.h"
#include "ofVectorMath.h"

namespace dab
{

namespace math
{

	class VectorMath : public Singleton<VectorMath>
	{
	public:
		/** Rotate vec3 by mat4x4 */
		glm::vec3 vec3Mat4Mul(const glm::mat4x4& pMatrix, const glm::vec3& pVector);

		/** Make a rotation Quat which will rotate vec1 to vec2 */
		glm::quat makeQuaternion(const glm::vec3& pFrom, const glm::vec3& pTo);
	};

}

}
