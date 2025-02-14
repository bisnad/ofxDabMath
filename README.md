## ofxDabMath

**Author**: Daniel Bisig - Coventry University, UK - [ad5041@coventry.ac.uk](ad5041@coventry.ac.uk) - Zurich University of the Arts, CH - [daniel.bisig@zhdk.ch](daniel.bisig@zhdk.ch)

**Dependencies**: [ofxDabBase](https://bitbucket.org/dbisig/ofxdabbase_011/src/master/), [Eigen](https://eigen.tuxfamily.org/) (included)

---

## Summary

ofxDabMath provides a variety of functions for conducting mathematical operations. In addition, it also a vector field class and methods for calculating vector fields.  Most of the vector calculations are done using the Eigen linear algebra library. The code is compatible with OpenFrameworks 0.11 and has been tested on Windows and MacOS. The following classes are available.

### Mathematical Operations

**Math**: provides functions for calculating random values, factorials, window envelopes, multidimensional array indices, and regular or normalised sigmoid values.

**VectorMath**: provides functions for rotation a three dimensional vector by a 4x4 matrix and to calculate a quaternion that rotates one three dimensional vector into another. 

### Vector Fields

**VectorField**: a vector field of arbitrary dimensions.

**FieldAlgorithm**: base class for calculating vector values for a vector field.

**RoesselerFieldAlgorithm**: example field algorithm based on the Roesseler Attractor.