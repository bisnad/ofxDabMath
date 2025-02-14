#include "ofApp.h"
#include "dab_math.h"
#include "dab_math_vector_field.h"
#include "dab_math_roesseler_field_algorithm.h"

//--------------------------------------------------------------
void ofApp::setup()
{
	// create empty vector field
	dab::Array<unsigned int> fieldSize({ 3, 4, 5 });

	std::cout << "fieldSize " << fieldSize << "\n";
	std::cout << "fieldSize size " << fieldSize.size() << "\n";

	Eigen::Vector3f fieldVec;

	// manually set vector values in vector field
	dab::math::VectorField<float> vectorField(fieldSize, fieldVec);

	dab::Array<unsigned int> fieldIndex(3);
	Eigen::Vector3f fieldVector;
	for (int z = 0; z<fieldSize[2]; ++z)
	{
		fieldIndex[2] = z;
		fieldVector[2] = z;

		for (int y = 0; y<fieldSize[1]; ++y)
		{
			fieldIndex[1] = y;
			fieldVector[1] = y;

			for (int x = 0; x<fieldSize[0]; ++x)
			{
				fieldIndex[0] = x;
				fieldVector[0] = x;

				vectorField.set(fieldIndex, fieldVector);
			}
		}
	}
	std::cout << vectorField << "\n";

	// copy vector field
	dab::math::VectorField<float> vectorField2(vectorField);
	std::cout << vectorField2 << "\n";

	// add a vector to all vectors in the vector field
	Eigen::Vector3f mathVec(10.0, 100.0, 1000.0);
	dab::math::VectorField<float> mathField = vectorField + mathVec;

	std::cout << mathField << "\n";

	// get interpolated vector by using a floating point vector index
	std::cout << "intervector " << vectorField.get(dab::Array<float>({ 1.1f, 2.5f, 3.9f })) << "\n";

	// set interpolated vector values by using a floating point vector index
	vectorField.set(dab::Array<float>({ 1.1f, 2.5f, 3.9f }), fieldVector);

	std::cout << vectorField << "\n";

	// sum all vector values in vector field
	Eigen::Vector3f sum = vectorField.sum();
	std::cout << "sum " << sum << "\n";

	// test multidimensional to scalar vector index conversion
	for (int z = 0; z<fieldSize[2]; ++z)
	{
		fieldIndex[2] = z;

		for (int y = 0; y<fieldSize[1]; ++y)
		{
			fieldIndex[1] = y;

			for (int x = 0; x<fieldSize[0]; ++x)
			{
				fieldIndex[0] = x;

				std::cout << "x " << x << " y " << y << " z " << z << " index " << vectorField.calcIndex(fieldIndex) << "\n";
			}
		}
	}

	// test scalar to multidimensional vector index conversion
	for (int i = 0; i<vectorField.vectorCount(); ++i)
	{
		fieldIndex = vectorField.calcIndex(i);

		std::cout << "index " << i;

		for (int d = 0; d<fieldIndex.size(); ++d)
		{
			std::cout << " " << fieldIndex[d];
		}
		std::cout << "\n";
	}

	// create a subfield from an existing vector field
	{
		dab::Array<unsigned int> subFieldStartIndex({ 0, 1, 2 });
		dab::math::VectorField<float> subField(dab::Array<unsigned int>({ 1, 2, 3 }), fieldVec);

		vectorField.subField(subFieldStartIndex, subField);

		std::cout << "subField\n";
		std::cout << subField << "\n";

		dab::math::VectorField<float> anotherField(vectorField);
		std::cout << anotherField << "\n";
	}

	// test application of convolution kernel on vector field
	{
		dab::Array<unsigned int> kernelSize({ 3, 3, 3 });
		Eigen::Matrix<float, 3, 1> kernelVector;
		kernelVector[0] = 1.0;
		kernelVector[1] = 1.0;
		kernelVector[2] = 1.0;
		dab::math::VectorField<float> kernel(kernelSize, kernelVector);

		std::cout << "kernel " << kernel << "\n";

		vectorField.convolve(kernel);

		std::cout << vectorField << "\n";
	}

	// test vector field algorithm
	{
		dab::math::RoesselerFieldAlgorithm* fieldAlgorithm = new dab::math::RoesselerFieldAlgorithm();
		dab::math::VectorField<float>* vectorField = new dab::math::VectorField<float>(dab::Array<unsigned int>({ 10, 10, 10 }), Eigen::Vector3f(0.0, 0.0, 0.0));

		fieldAlgorithm->setVectorField(vectorField);
		fieldAlgorithm->update();

		std::cout << *vectorField << "\n";
	}
}


//--------------------------------------------------------------
void ofApp::update() {

}

//--------------------------------------------------------------
void ofApp::draw() {

}

//--------------------------------------------------------------
void ofApp::keyPressed(int key) {

}

//--------------------------------------------------------------
void ofApp::keyReleased(int key) {

}

//--------------------------------------------------------------
void ofApp::mouseMoved(int x, int y) {

}

//--------------------------------------------------------------
void ofApp::mouseDragged(int x, int y, int button) {

}

//--------------------------------------------------------------
void ofApp::mousePressed(int x, int y, int button) {

}

//--------------------------------------------------------------
void ofApp::mouseReleased(int x, int y, int button) {

}

//--------------------------------------------------------------
void ofApp::mouseEntered(int x, int y) {

}

//--------------------------------------------------------------
void ofApp::mouseExited(int x, int y) {

}

//--------------------------------------------------------------
void ofApp::windowResized(int w, int h) {

}

//--------------------------------------------------------------
void ofApp::gotMessage(ofMessage msg) {

}

//--------------------------------------------------------------
void ofApp::dragEvent(ofDragInfo dragInfo) {

}
