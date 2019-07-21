#include <math.h>
#include "glm/vec4.hpp"

//***********
//QUATERNIONS
//***********
struct quat 
{
    float x = 0; 
    float y = 0; 
    float z = 0;
    float w = 1;
};

//Vector Adding 
glm::vec4 quatnerion_add(glm::vec4 vector1, glm::vec4 vector2)
{
    return {vector1.x + vector2.x, vector1.y + vector2.y, vector1.z + vector2.z, 1};
};

//Vector Subtracting 
glm::vec4 vector_subtract(glm::vec4 vector1, glm::vec4 vector2)
{
    return {vector1.x - vector2.x, vector1.y - vector2.y, vector1.z - vector2.z, 1};
};

//Vector Multiplication with Constant
glm::vec4 vector_multiply(glm::vec4 vector, float constant)
{
    return {vector.x * constant, vector.y * constant, vector.z * constant, 1};
};

//Multiplying two quaternions using Hamilton Product  
quat quaternion_divide(quat quaternion, float constant)
{   
    return {quaternion.x / constant, quaternion.y / constant, quaternion.z / constant, 1};
};

//Multiplying two quaternions using Hamilton Product  
quat quaternion_multiply(quat quaternion1, quat quaternion2)
{   
    quat result;

    result.x = (quaternion1.w*quaternion2.x + quaternion1.x*quaternion2.w + quaternion1.y*quaternion2.z - quaternion1.z*quaternion2.y);
    result.y = (quaternion1.w*quaternion2.y - quaternion1.x*quaternion2.z + quaternion1.y*quaternion2.w + quaternion1.z*quaternion2.x);
    result.z = (quaternion1.w*quaternion2.z + quaternion1.x*quaternion2.y - quaternion1.y*quaternion2.x + quaternion1.z*quaternion2.w);
    result.w = (quaternion1.w*quaternion2.w - quaternion1.x*quaternion2.x - quaternion1.y*quaternion2.y - quaternion1.z*quaternion2.z);

    return result;
};

//Normalization
void quaternion_normalize(quat vector)
{   
    float length = pow(pow(vector.x, 2) + pow(vector.y, 2) + pow(vector.z, 2) + pow(vector.w, 2), 0.5);

    vector = quaternion_divide(vector, length);
};

//The basic quaternion structure 
quat quaternion_structure(glm::vec3 axis, float angle)
{   
    quat quaternion;

    quaternion.x = axis.x * sinf( angle/2 );
    quaternion.y = axis.y * sinf( angle/2 );
    quaternion.z = axis.z * sinf( angle/2 );
    quaternion.w = cosf( angle/2 );

    quaternion_normalize(quaternion);

    return quaternion;
};

//Seting up the new global rotation based on input axis, change in angles and total quaternion
quat quaternion_setup(quat total_quaternion, glm::vec3 angle, glm::vec3 x_axis, glm::vec3 y_axis, glm::vec3 z_axis)
{   
    quat quaternion_x = quaternion_structure(x_axis, angle.x);
    quat quaternion_y = quaternion_structure(y_axis, angle.y);
    quat quaternion_z = quaternion_structure(z_axis, angle.z);

    //Multiplying change in quaternion by universal quaternion then rotating point
    quat quaternion = quaternion_multiply(quaternion_multiply(quaternion_multiply(quaternion_z, quaternion_y), quaternion_x), total_quaternion);

    return quaternion;

    /* This Way Does it all at once and more efficiently and should be reimplented (tried and didnt work)
    //precompute to save on processing time
    float cosX = cos( angle[1] / 2 );
    float cosY = cos( angle[2] / 2 );
    float cosZ = cos( angle[0] / 2 );
    float sinX = sin( angle[1] / 2 );
    float sinY = sin( angle[2] / 2 );
    float sinZ = sin( angle[0] / 2 );

    //this is apparaently equavalent to rotation around all 3 axis
    result[0] = cosX * cosY * cosZ - sinX * sinY * sinZ;
    result[1] = sinX * cosY * cosZ - cosX * sinY * sinZ;
    result[2] = cosX * sinY * cosZ + sinX * cosY * sinZ;
    result[3] = cosX * cosY * sinZ + sinX * sinY * cosZ;
    */

};

//Provides the conjugate of a quaternion
quat quaternion_conjugatation(quat quaternion)
{
    return {-quaternion.x, -quaternion.y, -quaternion.z, quaternion.w};
};

//Converting a Vector into a quaternion to match Data Types 
quat quaternion_from_vector(glm::vec4 vector)
{
    return{vector.x, vector.y, vector.z, vector.w};
};

//Converting a Vector into a quaternion to match Data Types 
glm::vec4 vector_from_quaternion(quat quaternion)
{
    return{quaternion.x, quaternion.y, quaternion.z, quaternion.w};
};

//Rotates a point or vector based on: R = P*Q*P^-1
glm::vec4 quaternion_rotation(quat quaternion, glm::vec4 position)
{   
    quat conjugate = quaternion_conjugatation(quaternion);

    return(vector_from_quaternion(quaternion_multiply(quaternion_multiply(quaternion, quaternion_from_vector(position)), conjugate)));
}; 

// Angles between Vectors
glm::vec4 quaternion_to_euler(quat quaternion)
{		
	glm::vec4 euler;

	// if the input quaternion is normalized, this is exactly one. Otherwise, this acts as a correction factor for the quaternion's not-normalizedness
	float unit = (quaternion.x * quaternion.x) + (quaternion.y * quaternion.y) + (quaternion.z * quaternion.z) + (quaternion.w * quaternion.w);

	// this will have a magnitude of 0.5 or greater if and only if this is a singularity case
	float test = quaternion.x * quaternion.w - quaternion.y * quaternion.z;

	if (test > 0.4995f * unit) // singularity at north pole
	{
		euler.x = 3.141593 / 2;
		euler.y = 2.0f * atan2(quaternion.y, quaternion.x);
		euler.z = 0;
	}
	else if (test < -0.4995f * unit) // singularity at south pole
	{
		euler.x = -3.141593 / 2;
		euler.y = -2.0f * atan2(quaternion.y, quaternion.x);
		euler.z = 0;
	}
	else // no singularity - this is the majority of cases
	{
		euler.x = asin(2.0f * (quaternion.w * quaternion.x - quaternion.y * quaternion.z));
		euler.y = atan2(2.0f * quaternion.w * quaternion.y + 2.0f * quaternion.z * quaternion.x, 1 - 2.0f * (quaternion.x * quaternion.x + quaternion.y * quaternion.y));
		euler.z = atan2(2.0f * quaternion.w * quaternion.z + 2.0f * quaternion.x * quaternion.y, 1 - 2.0f * (quaternion.z * quaternion.z + quaternion.x * quaternion.x));
	}

	return euler;
};

//*********************
//Collision Detetection AABB (voxel) and Polygon (Theory by Tomas Akenine-Moller, code from someelse)
//*********************

//This checks for intersection between polygon and voxel
#pragma once
#include <cmath>

inline void findMinMax(float x0, float x1, float x2, float &min, float &max) 
{
	min = max = x0;
	if (x1 < min)
		min = x1;
	if (x1 > max)
		max = x1;
	if (x2 < min)
		min = x2;
	if (x2 > max)
		max = x2;
}

inline bool planeBoxOverlap(glm::vec3 normal, glm::vec3 vert, glm::vec3 maxbox) 
{
	glm::vec3 vmin, vmax;
	float v;

    v = vert.x;
    if (normal.x > 0.0f) {
        vmin.x = -maxbox.x - v;
        vmax.x = maxbox.x - v;
    } else {
        vmin.x = maxbox.x - v;
        vmax.x = -maxbox.x - v;
    }

    v = vert.y;
    if (normal.y > 0.0f) {
        vmin.y = -maxbox.y - v;
        vmax.y = maxbox.y - v;
    } else {
        vmin.y = maxbox.y - v;
        vmax.y = -maxbox.y - v;
    }

    v = vert.z;
    if (normal.z > 0.0f) {
        vmin.z = -maxbox.z - v;
        vmax.z = maxbox.z - v;
    } else {
        vmin.z = maxbox.z - v;
        vmax.z = -maxbox.z - v;
    }
	
	if (glm::dot(normal, vmin) > 0.0f)
		return false;
	if (glm::dot(normal, vmax) >= 0.0f)
		return true;

	return false;
}

/*======================== X-tests ========================*/

inline bool axisTestX01(float a, float b, float fa, float fb, const glm::vec3 &v0, const glm::vec3 &v2, const glm::vec3 &boxhalfsize, float &rad, float &min, float &max, float &p0, float &p2) 
{
	p0 = a * v0.y - b * v0.z;
	p2 = a * v2.y - b * v2.z;
	if (p0 < p2) {
		min = p0;
		max = p2;
	} else {
		min = p2;
		max = p0;
	}
	rad = fa * boxhalfsize.y + fb * boxhalfsize.z;
	if (min > rad || max < -rad)
		return false;
	return true;
}
inline bool axisTestX2(float a, float b, float fa, float fb, const glm::vec3 &v0, const glm::vec3 &v1, const glm::vec3 &boxhalfsize, float &rad, float &min, float &max, float &p0, float &p1) 
{
	p0 = a * v0.y - b * v0.z;
	p1 = a * v1.y - b * v1.z;
	if (p0 < p1) {
		min = p0;
		max = p1;
	} else {
		min = p1;
		max = p0;
	}
	rad = fa * boxhalfsize.y + fb * boxhalfsize.z;
	if (min > rad || max < -rad)
		return false;
	return true;
}

/*======================== Y-tests ========================*/

inline bool axisTestY02(float a, float b, float fa, float fb, const glm::vec3 &v0, const glm::vec3 &v2, const glm::vec3 &boxhalfsize, float &rad, float &min, float &max, float &p0, float &p2) 
{
	p0 = -a * v0.x + b * v0.z;
	p2 = -a * v2.x + b * v2.z;
	if (p0 < p2) {
		min = p0;
		max = p2;
	} else {
		min = p2;
		max = p0;
	}
	rad = fa * boxhalfsize.x + fb * boxhalfsize.z;
	if (min > rad || max < -rad)
		return false;
	return true;
}

inline bool axisTestY1(float a, float b, float fa, float fb, const glm::vec3 &v0, const glm::vec3 &v1, const glm::vec3 &boxhalfsize, float &rad, float &min, float &max, float &p0, float &p1) 
{
	p0 = -a * v0.x + b * v0.z;
	p1 = -a * v1.x + b * v1.z;
	if (p0 < p1) {
		min = p0;
		max = p1;
	} else {
		min = p1;
		max = p0;
	}
	rad = fa * boxhalfsize.x + fb * boxhalfsize.z;
	if (min > rad || max < -rad)
		return false;
	return true;
}

/*======================== Z-tests ========================*/
inline bool axisTestZ12(float a, float b, float fa, float fb, const glm::vec3 &v1, const glm::vec3 &v2, const glm::vec3 &boxhalfsize, float &rad, float &min, float &max, float &p1, float &p2) 
{
	p1 = a * v1.x - b * v1.y;
	p2 = a * v2.x - b * v2.y;
	if (p1 < p2) {
		min = p1;
		max = p2;
	} else {
		min = p2;
		max = p1;
	}
	rad = fa * boxhalfsize.x + fb * boxhalfsize.y;
	if (min > rad || max < -rad)
		return false;
	return true;
}

inline bool axisTestZ0(float a, float b, float fa, float fb, const glm::vec3 &v0, const glm::vec3 &v1, const glm::vec3 &boxhalfsize, float &rad, float &min, float &max, float &p0, float &p1) 
{
	p0 = a * v0.x - b * v0.y;
	p1 = a * v1.x - b * v1.y;
	if (p0 < p1) {
		min = p0;
		max = p1;
	} else {
		min = p1;
		max = p0;
	}
	rad = fa * boxhalfsize.x + fb * boxhalfsize.y;
	if (min > rad || max < -rad)
		return false;
	return true;
}

bool voxel_mesh_intersection(glm::vec3 boxcenter, glm::vec3 boxhalfsize, glm::vec3 tv0, glm::vec3 tv1, glm::vec3 tv2) {
	glm::vec3 v0, v1, v2;
	float min, max, p0, p1, p2, rad, fex, fey, fez;
	glm::vec3 normal, e0, e1, e2;

	v0 = tv0 - boxcenter;
	v1 = tv1 - boxcenter;
	v2 = tv2 - boxcenter;

	e0 = v1 - v0;
	e1 = v2 - v1;
	e2 = v2 - v1;

	fex = fabsf(e0.x);
	fey = fabsf(e0.y);
	fez = fabsf(e0.z);

	if (!axisTestX01(e0.z, e0.y, fez, fey, v0, v2, boxhalfsize, rad, min, max, p0, p2))
		return 0;
	if (!axisTestY02(e0.z, e0.x, fez, fex, v0, v2, boxhalfsize, rad, min, max, p0, p2))
		return 0;
	if (!axisTestZ12(e0.y, e0.x, fey, fex, v1, v2, boxhalfsize, rad, min, max, p1, p2))
		return 0;

	fex = fabsf(e1.x);
	fey = fabsf(e1.y);
	fez = fabsf(e1.z);

	if (!axisTestX01(e1.z, e1.y, fez, fey, v0, v2, boxhalfsize, rad, min, max, p0, p2))
		return 0;
	if (!axisTestY02(e1.z, e1.x, fez, fex, v0, v2, boxhalfsize, rad, min, max, p0, p2))
		return 0;
	if (!axisTestZ0(e1.y, e1.x, fey, fex, v0, v1, boxhalfsize, rad, min, max, p0, p1))
		return 0;

	fex = fabsf(e2.x);
	fey = fabsf(e2.y);
	fez = fabsf(e2.z);
	if (!axisTestX2(e2.z, e2.y, fez, fey, v0, v1, boxhalfsize, rad, min, max, p0, p1))
		return 0;
	if (!axisTestY1(e2.z, e2.x, fez, fex, v0, v1, boxhalfsize, rad, min, max, p0, p1))
		return 0;
	if (!axisTestZ12(e2.y, e2.x, fey, fex, v1, v2, boxhalfsize, rad, min, max, p1, p2))
		return 0;

	findMinMax(v0.x, v1.x, v2.x, min, max);
	if (min > boxhalfsize.x || max < -boxhalfsize.x)
		return 0;

	findMinMax(v0.y, v1.y, v2.y, min, max);
	if (min > boxhalfsize.y || max < -boxhalfsize.y)
		return 0;

	findMinMax(v0.z, v1.z, v2.z, min, max);
	if (min > boxhalfsize.z || max < -boxhalfsize.z)
		return 0;

	normal = glm::cross(e0, e1);
	if (!planeBoxOverlap(normal, v0, boxhalfsize))
		return 0;

	return 1;
}