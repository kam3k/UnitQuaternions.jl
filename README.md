# UnitQuaternions.jl
A Julia implementation of unit quaternions for representing rotations. Created by [Marc Gallant](http://kam3k.github.io), originally for use in the [Mining Systems Laboratory](https://msl.engineering.queensu.ca).

## Installation
From the Julia REPL, simply clone this package via

    julia> Pkg.clone("https://github.com/kam3k/UnitQuaternions.jl.git")
    
and then add it to your current session with

    julia> using UnitQuaternions

## Description
Unit quaternions in this package are represented by the 4 &times; 1 column 

**q** = [**&epsilon;** &eta;]<sup>T</sup>

where **&epsilon;** is a 3-vector (the *vector* part of the unit quaternion) and &eta; is a scalar (the *scalar* part of the unit quaternion). Because **q**&equiv;-**q**, &eta; &ge; 0 is always maintained to make every **q** unique.

Given a rotation parameterized by an axis of rotation **a** and an angle &theta;, the same rotation parameterized by a unit quaternion is

**q** = [**a**sin(&theta;/2) cos(&theta;/2)]<sup>T</sup>.

This package supports several constructors, operators, and methods on unit quaternions.

## API

The API is given in the form of examples. When necessary, mathematical background or references are included to supplement the examples.

### Constructors

    julia> UnitQuaternion()
	ϵ = [0.000, 0.000, 0.000], η = 1.000
An empty argument list constructs the identity unit quaternion.

    julia> q = UnitQuaternion(1,2,3,4)
    ϵ = [0.183, 0.365, 0.548], η = 0.730
An argument list of four real numbers assigns the first three to the vector part and the fourth to the scalar part. The entries are normalized.

    julia> UnitQuaternion([1,2,3,4])
    ϵ = [0.183, 0.365, 0.548], η = 0.730
A 4-vector of real numbers given as the only argument is equivalent to `UnitQuaternion(1,2,3,4)`.

    julia> UnitQuaternion([1,2,3])
    ϵ = [-0.255, -0.511, -0.766], η = 0.296
A 3-vector of real numbers as the only argument is assumed to be a rotation vector. A rotation vector **r** = &theta;**a** is the product of an angle of rotation &theta; and an axis of rotation **a**, where **a** is a 3-vector. Put differently, the magnitude of **r** is the angle of rotation, and the direction of **r** is the axis of rotation.

### Operators

    julia> q = UnitQuaternion(1,2,3,4)
    ϵ = [0.183, 0.365, 0.548], η = 0.730

    julia> +(q)
	4x4 Array{Float64,2}:
	  0.730297   0.547723  -0.365148  0.182574
	 -0.547723   0.730297   0.182574  0.365148
	  0.365148  -0.182574   0.730297  0.547723
	 -0.182574  -0.365148  -0.547723  0.730297
The `+` operator with a single unit quaternion as its argument applies the left-hand compound operator to **q** (see [Pose estimation using linearized rotations and quaternion algebra](http://www.sciencedirect.com/science/article/pii/S0094576510002407) by Barfoot *et al*.), returning a 4 &times; 4 matrix.

	julia> q = UnitQuaternion(1,2,3,4)
	ϵ = [0.183, 0.365, 0.548], η = 0.730
	
	julia> ⊕(q)
	4x4 Array{Float64,2}:
	  0.730297  -0.547723   0.365148  0.182574
	  0.547723   0.730297  -0.182574  0.365148
	 -0.365148   0.182574   0.730297  0.547723
	 -0.182574  -0.365148  -0.547723  0.730297
The `⊕` operator with a single unit quaternion as its argument applies the right-hand compound operator to **q** (see [Pose estimation using linearized rotations and quaternion algebra](http://www.sciencedirect.com/science/article/pii/S0094576510002407) by Barfoot *et al*.), returning a 4 &times; 4 matrix.

## Help
Use the Julia help system (i.e., enter a question mark in the REPL) to query information about the methods. For example,

    help?> rotatevector
    search: rotatevector

      rotatevector(q, v)

      Given a unit quaternion q and a 3-vector v, performs a vector rotation and returns the rotated vector v.