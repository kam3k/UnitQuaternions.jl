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

The API is given in the form of examples. When necessary, mathematical background or references are included to supplement the examples. The examples described in this section use the following quaternions:

	julia> q = UnitQuaternion(1,2,3,4)
	ϵ = [0.183, 0.365, 0.548], η = 0.730
	
	julia> p = UnitQuaternion(-5,3,-2,4)
	ϵ = [-0.680, 0.408, -0.272], η = 0.544

There is nothing special about these unit quaternions, they are simply used to demonstrate various operations.

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

#### Quaternion product
	julia> q + p
	ϵ = [-0.075, 0.820, -0.224], η = 0.522
Returns the quaternion product **q** &otimes; **p**. Equivalent to `p ⊕ q`.

	julia> q ⊕ p
	ϵ = [-0.721, 0.174, 0.422], η = 0.522
Returns the quaternion product **p** &otimes; **q**. Equivalent to `p + q`.

#### Compound operators
For background information on compound operators, see [Pose estimation using linearized rotations and quaternion algebra](http://www.sciencedirect.com/science/article/pii/S0094576510002407) by Barfoot *et al*.
 
    julia> +(q)
	4x4 Array{Float64,2}:
	  0.730297   0.547723  -0.365148  0.182574
	 -0.547723   0.730297   0.182574  0.365148
	  0.365148  -0.182574   0.730297  0.547723
	 -0.182574  -0.365148  -0.547723  0.730297
The `+` operator with a single unit quaternion as its argument applies the left-hand compound operator to **q**, returning a 4 &times; 4 matrix.

	julia> q = UnitQuaternion(1,2,3,4)
	ϵ = [0.183, 0.365, 0.548], η = 0.730
	
	julia> ⊕(q)
	4x4 Array{Float64,2}:
	  0.730297  -0.547723   0.365148  0.182574
	  0.547723   0.730297  -0.182574  0.365148
	 -0.365148   0.182574   0.730297  0.547723
	 -0.182574  -0.365148  -0.547723  0.730297
The `⊕` operator with a single unit quaternion as its argument applies the right-hand compound operator to **q**, returning a 4 &times; 4 matrix.

#### Manifold operations
For background information on manifold operations, see [Integrating generic sensor fusion algorithms with sound state representations through encapsulation of manifolds](http://www.sciencedirect.com/science/article/pii/S1566253511000571) by Hertzberg *et al*.

	julia> q ⊞ [0.1, 0.2, -0.1]
	ϵ = [0.290, 0.399, 0.507], η = 0.707
Perturbs **q** by the rotation vector `[0.1, 0.2, -0.1]`.

	julia> q ⊟ p
	3-element Array{Float64,1}:
	 2.47319
	 0.601587
	 0.467901
Calculates the rotational difference between **q** and **p** and returns the result as a rotation vector.

#### Miscellaneous

	julia> q^1.273
	ϵ = [0.219, 0.437, 0.656], η = 0.576
The `^n` operator raises **q** to the power of the real number `n`.

	julia> q == p
	false
	
	julia> q == UnitQuaternion(1,2,3,4)
	true
The `==` operator checks for equivalency between two unit quaternions.

	julia> q[3]
	0.5477225575051661
`q[n]` returns the `n`-th index of **q**, where `n` = 1, 2, 3, or 4. `q[4]` is the scalar part of **q**.

### Methods

The methods are listed in alphabetical order.

	julia> angle(q)
	1.5040801783846713
Returns the angle of rotation (in radians) of the rotation parameterized by **q**.

	julia> axis(q)
	3-element Array{Float64,1}:
	 0.267261
	 0.534522
	 0.801784
Returns the axis of rotation (unit vector) of the rotation parameterized by **q**.

	julia> covariance([p,q])
	3x3 Array{Float64,2}:
	 1.52917   0.37196    0.289302
	 0.37196   0.0904767  0.0703708
	 0.289302  0.0703708  0.0547328
	 
	julia> covariance([p,q], [1,2])
	3x3 Array{Float64,2}:
	 1.53098   0.3724     0.289645
	 0.3724    0.0905839  0.0704541
	 0.289645  0.0704541  0.0547977
Returns the 3 &times; 3 covariance matrix of a vector of unit quaternions (`[p,q]` in this case). A second argument of weights (`[1,2]` in this case) can be provided.

	julia> inv(q)
	ϵ = [-0.183, -0.365, -0.548], η = 0.730
Returns the inverse of **q**. Recall that for unit quaternions, the inverse is simply **q**<sup>-1</sup> = [-**&epsilon;** &eta;]<sup>T</sup>.

	julia> mean([q,p])
	ϵ = [-0.312, 0.485, 0.173], η = 0.799
	
	julia> mean([q,p],[1,2])
	ϵ = [-0.583, 0.455, -0.128], η = 0.661
Returns the mean unit quaternion of a vector of unit quaternions (`[p,q]` in this case). A second argument of weights (`[1,2]` in this case) can be provided.

	julia> rotateframe(q,[1,0,0])
	3-element Array{Float64,1}:
	  0.133333
	 -0.666667
	  0.733333
Given a unit quaternion that represents the rotation of coordinate frame B with respect to coordinate frame A, and a vector in coordinate frame A (`[1,0,0]` in this case), returns the vector in coordinate frame B. Equivalent to `rotatevector(inv(q),[1,0,0])`.

	julia> rotatevector(q,[1,0,0])
	3-element Array{Float64,1}:
	  0.133333
	  0.933333
	 -0.333333
Rotates a vector (`[1,0,0]` in this case) by the rotation parameterized by **q**. Equivalent to `rotateframe(inv(q),[1,0,0])`.

	julia> rotationmatrix(q)
	3x3 Array{Float64,2}:
	  0.133333  0.933333  -0.333333
	 -0.666667  0.333333   0.666667
	  0.733333  0.133333   0.666667
Returns the rotation matrix (sometimes called a direction cosine matrix) parameterization of the rotation parameterized by **q**. Uses the natural order of rotations. In other words, `rotationmatrix(q) * [1,0,0]` is equivalent to `rotateframe(q,[1,0,0])`.

	julia> log(q)
	3-element Array{Float64,1}:
	 0.200991
	 0.401982
	 0.602974
Returns the log of unit quaternion `q`.

	julia> slerp(q,p,0.2)
	ϵ = [-0.018, 0.435, 0.417], η = 0.798
Performs spherical linear interpolation from **q** to **p**. The interpolation parameter `t` (`0.2` in this case) specifies the fraction of the arc to traverse. The returned unit quaternion parameterizes a rotation from **q** to the end point specified by `t`.

	julia> vector(q)
	4-element Array{Float64,1}:
	 0.182574
	 0.365148
	 0.547723
	 0.730297
Returns the elements of **q** as a 4-vector.

## Help
Use the Julia help system (i.e., enter a question mark in the REPL) to query information about the methods. For example,

	help?> slerp
	search: slerp islower selectperm selectperm! sylvester Serializer
	
	  slerp(p, q, t)
	
	  Performs spherical linear interpolation from unit quaternion p to unit quaternion q. The interpolation parameter t specifies the fraction of arc to traverse. The returned unit quaternion parameterizes a rotation from p to the end point specified by t.
