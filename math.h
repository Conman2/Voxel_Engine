#include <math.h>

//*******
//VECTORS
//*******

struct vect 
{
    float i = 0; 
    float j = 0; 
    float k = 0;
    float w = 1;
};

//Vector Adding 
vect vector_add(vect vector1, vect vector2)
{
    return {vector1.i + vector2.i, vector1.j + vector2.j, vector1.k + vector2.k, 1};
}

//Vector Subtracting 
vect vector_subtract(vect vector1, vect vector2)
{
    return {vector1.i - vector2.i, vector1.j - vector2.j, vector1.k - vector2.k, 1};
}

//Vector Multiplication with Constant
vect vector_multiply(vect vector, float constant)
{
    return {vector.i * constant, vector.j * constant, vector.k * constant, 1};
}

vect vector_divide(vect vector, float constant)
{
    return {vector.i / constant, vector.j / constant, vector.k / constant, 1};
}

//Dot Product 
float vector_dot_product(vect vector1, vect vector2)
{
    return {vector1.i * vector2.i + vector1.j * vector2.j + vector1.k * vector2.k}; 
}

//Cross Product
vect vector_cross_product(vect vector1, vect vector2)
{   
    return {vector1.j * vector2.k - vector1.k * vector2.j, vector1.k * vector2.i - vector1.i * vector2.k, vector1.i * vector2.j - vector1.j * vector2.i, 1};
}

//Magnitude of a vector
float vector_magnitude(vect vector)
{
	return pow(pow(vector.i, 2) + pow(vector.j, 2) + pow(vector.k, 2), 0.5);
};

//Normalization
vect vector_normalize(vect vector)
{  
    return vector_divide(vector, vector_magnitude(vector));
}

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
}

//The basic quaternion structure 
quat quaternion_structure(vect axis, float angle)
{   
    quat quaternion;

    quaternion.x = axis.i * sinf( angle/2 );
    quaternion.y = axis.j * sinf( angle/2 );
    quaternion.z = axis.k * sinf( angle/2 );
    quaternion.w = cosf( angle/2 );

    quaternion_normalize(quaternion);

    return quaternion;
};

//Seting up the new global rotation based on input axis, change in angles and total quaternion
quat quaternion_setup(quat total_quaternion, vect angle, vect x_axis, vect y_axis, vect z_axis)
{   
    quat quaternion_x = quaternion_structure(x_axis, angle.i);
    quat quaternion_y = quaternion_structure(y_axis, angle.j);
    quat quaternion_z = quaternion_structure(z_axis, angle.k);

    //Multiplying change in quaternion by universal quaternion then rotating point
    quat quaternion = quaternion_multiply(quaternion_multiply(quaternion_multiply(quaternion_z, quaternion_y), quaternion_x), total_quaternion);

    return quaternion;

    /* This Way Does it all at once and more efficiently and should be reimplented 
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
quat quaternion_from_vector(vect vector)
{
    return{vector.i, vector.j, vector.k, vector.w};
};

//Converting a Vector into a quaternion to match Data Types 
vect vector_from_quaternion(quat quaternion)
{
    return{quaternion.x, quaternion.y, quaternion.z, quaternion.w};
};

//Rotates a point or vector based on: R = P*Q*P^-1
vect quaternion_rotation(quat quaternion, vect position)
{   
    quat conjugate = quaternion_conjugatation(quaternion);

    return(vector_from_quaternion(quaternion_multiply(quaternion_multiply(quaternion, quaternion_from_vector(position)), conjugate)));
}; 

//******** 
//MATRICES
//********

//Matrix Vector Multiplication (note that with this function input and output must be different variables)
vect matrix_vector_multiplication(vect vector, float matrix[4][4])
{   
    vect result;

    result.i = vector.i* matrix[0][0] + vector.j * matrix[0][1] + vector.k * matrix[0][2] + vector.w * matrix[0][3];
    result.j = vector.i* matrix[1][0] + vector.j * matrix[1][1] + vector.k * matrix[1][2] + vector.w * matrix[1][3];
    result.k = vector.i* matrix[2][0] + vector.j * matrix[2][1] + vector.k * matrix[2][2] + vector.w * matrix[2][3];
    result.w = vector.i* matrix[3][0] + vector.j * matrix[3][1] + vector.k * matrix[3][2] + vector.w * matrix[3][3];

    return result;
}

//Look at Matrix
void matrix_lookat(float matrix[4][4], vect camera_position, vect camera_direction, vect world_up, vect camera_right)
{
    // //Calculating the Camera Looking Direction
    //vect camera_direction = vector_subtract(camera_position, camera_target); 
    //vector_normalize(camera_direction);

    //Calculating new up
    vect camera_up = vector_subtract(vector_multiply(camera_direction, vector_dot_product(world_up, camera_direction)), world_up);//Might need t0 be reversed
    vector_normalize(camera_up);

    //Calculating Camera Right Direction
    //vect camera_right = cross_product(camera_up, camera_direction);

    //Look At Matrix 
    matrix[0][0] = camera_right.i;     matrix[0][1] = camera_right.j;     matrix[0][2] = camera_right.k;     matrix[0][3] = -vector_dot_product(camera_right, camera_position);
    matrix[1][0] = camera_up.i;        matrix[1][1] = camera_up.j;        matrix[1][2] = camera_up.k;        matrix[1][3] = -vector_dot_product(camera_up, camera_position);
    matrix[2][0] = camera_direction.i; matrix[2][1] = camera_direction.j; matrix[2][2] = camera_direction.k; matrix[2][3] = -vector_dot_product(camera_direction, camera_position);
    matrix[3][0] = 0;                  matrix[3][1] = 0;                  matrix[3][2] = 0;                  matrix[3][3] = 1; 
}

void matrix_projection(float matrix[4][4], float camera_view_angle, float screen_height, float screen_width, float z_max_distance, float z_min_distance)
{
    //Projection Matrix Set-up
    float view_angle =    (float) camera_view_angle*(3.14159 / 180);
    float aspect_ratio =  (float) screen_height/screen_width;
    float z_normalize =   (float) z_max_distance/(z_max_distance - z_min_distance);
    float z_other =       (float) (-z_max_distance * z_min_distance)/(z_max_distance - z_min_distance);
    float feild_of_view = (float) 1/tan(view_angle / 2);

    //Assigning values to the correct position in the matrix
    matrix[0][0] = aspect_ratio*feild_of_view; 
    matrix[1][1] = feild_of_view; 
    matrix[2][2] = z_normalize; 
    matrix[2][3] = z_other; 
    matrix[3][2] = 1; 
};

//Clears a 4x4 matrix and sets the content to 0
void matrix_clear(float matrix[4][4])
{
	for (int i = 0; i < 5; i ++)
	{     
		for (int j = 0; j < 5; j ++)
		{
            matrix[i][j] = 0.0;
        }
    }
}  

//Clears a 4x4 matrix and sets the content to 0
void matrix_equalize_bool4x6(bool matrix1[6][4], const bool matrix2[6][4])
{
	for (int i = 0; i < 6; i ++)
	{     
		for (int j = 0; j < 4; j ++)
		{
            matrix1[i][j] = matrix2[i][j];
        }
    }
}  

//*********************
//Collision Detetection 
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

inline bool planeBoxOverlap(vect normal, vect vert, vect maxbox) 
{
	vect vmin, vmax;
	float v;

    v = vert.i;
    if (normal.i > 0.0f) {
        vmin.i = -maxbox.i - v;
        vmax.i = maxbox.i - v;
    } else {
        vmin.i = maxbox.i - v;
        vmax.i = -maxbox.i - v;
    }

    v = vert.j;
    if (normal.j > 0.0f) {
        vmin.j = -maxbox.j - v;
        vmax.j = maxbox.j - v;
    } else {
        vmin.j = maxbox.j - v;
        vmax.j = -maxbox.j - v;
    }

    v = vert.k;
    if (normal.k > 0.0f) {
        vmin.k = -maxbox.k - v;
        vmax.k = maxbox.k - v;
    } else {
        vmin.k = maxbox.k - v;
        vmax.k = -maxbox.k - v;
    }
	
	if (vector_dot_product(normal, vmin) > 0.0f)
		return false;
	if (vector_dot_product(normal, vmax) >= 0.0f)
		return true;

	return false;
}

/*======================== X-tests ========================*/

inline bool axisTestX01(float a, float b, float fa, float fb, const vect &v0, const vect &v2, const vect &boxhalfsize, float &rad, float &min, float &max, float &p0, float &p2) 
{
	p0 = a * v0.j - b * v0.k;
	p2 = a * v2.j - b * v2.k;
	if (p0 < p2) {
		min = p0;
		max = p2;
	} else {
		min = p2;
		max = p0;
	}
	rad = fa * boxhalfsize.j + fb * boxhalfsize.k;
	if (min > rad || max < -rad)
		return false;
	return true;
}
inline bool axisTestX2(float a, float b, float fa, float fb, const vect &v0, const vect &v1, const vect &boxhalfsize, float &rad, float &min, float &max, float &p0, float &p1) 
{
	p0 = a * v0.j - b * v0.k;
	p1 = a * v1.j - b * v1.k;
	if (p0 < p1) {
		min = p0;
		max = p1;
	} else {
		min = p1;
		max = p0;
	}
	rad = fa * boxhalfsize.j + fb * boxhalfsize.k;
	if (min > rad || max < -rad)
		return false;
	return true;
}

/*======================== Y-tests ========================*/

inline bool axisTestY02(float a, float b, float fa, float fb, const vect &v0, const vect &v2, const vect &boxhalfsize, float &rad, float &min, float &max, float &p0, float &p2) 
{
	p0 = -a * v0.i + b * v0.k;
	p2 = -a * v2.i + b * v2.k;
	if (p0 < p2) {
		min = p0;
		max = p2;
	} else {
		min = p2;
		max = p0;
	}
	rad = fa * boxhalfsize.i + fb * boxhalfsize.k;
	if (min > rad || max < -rad)
		return false;
	return true;
}

inline bool axisTestY1(float a, float b, float fa, float fb, const vect &v0, const vect &v1, const vect &boxhalfsize, float &rad, float &min, float &max, float &p0, float &p1) 
{
	p0 = -a * v0.i + b * v0.k;
	p1 = -a * v1.i + b * v1.k;
	if (p0 < p1) {
		min = p0;
		max = p1;
	} else {
		min = p1;
		max = p0;
	}
	rad = fa * boxhalfsize.i + fb * boxhalfsize.k;
	if (min > rad || max < -rad)
		return false;
	return true;
}

/*======================== Z-tests ========================*/
inline bool axisTestZ12(float a, float b, float fa, float fb, const vect &v1, const vect &v2, const vect &boxhalfsize, float &rad, float &min, float &max, float &p1, float &p2) 
{
	p1 = a * v1.i - b * v1.j;
	p2 = a * v2.i - b * v2.j;
	if (p1 < p2) {
		min = p1;
		max = p2;
	} else {
		min = p2;
		max = p1;
	}
	rad = fa * boxhalfsize.i + fb * boxhalfsize.j;
	if (min > rad || max < -rad)
		return false;
	return true;
}

inline bool axisTestZ0(float a, float b, float fa, float fb, const vect &v0, const vect &v1, const vect &boxhalfsize, float &rad, float &min, float &max, float &p0, float &p1) 
{
	p0 = a * v0.i - b * v0.j;
	p1 = a * v1.i - b * v1.j;
	if (p0 < p1) {
		min = p0;
		max = p1;
	} else {
		min = p1;
		max = p0;
	}
	rad = fa * boxhalfsize.i + fb * boxhalfsize.j;
	if (min > rad || max < -rad)
		return false;
	return true;
}

bool voxel_mesh_intersection(vect boxcenter, vect boxhalfsize, vect tv0, vect tv1, vect tv2) {
	vect v0, v1, v2;
	float min, max, p0, p1, p2, rad, fex, fey, fez;
	vect normal, e0, e1, e2;

	v0 = vector_subtract(tv0, boxcenter);
	v1 = vector_subtract(tv1, boxcenter);
	v2 = vector_subtract(tv2, boxcenter);

	e0 = vector_subtract(v1, v0);
	e1 = vector_subtract(v2, v1);
	e2 = vector_subtract(v2, v1);

	fex = fabsf(e0.i);
	fey = fabsf(e0.j);
	fez = fabsf(e0.k);

	if (!axisTestX01(e0.k, e0.j, fez, fey, v0, v2, boxhalfsize, rad, min, max, p0, p2))
		return 0;
	if (!axisTestY02(e0.k, e0.i, fez, fex, v0, v2, boxhalfsize, rad, min, max, p0, p2))
		return 0;
	if (!axisTestZ12(e0.j, e0.i, fey, fex, v1, v2, boxhalfsize, rad, min, max, p1, p2))
		return 0;

	fex = fabsf(e1.i);
	fey = fabsf(e1.j);
	fez = fabsf(e1.k);

	if (!axisTestX01(e1.k, e1.j, fez, fey, v0, v2, boxhalfsize, rad, min, max, p0, p2))
		return 0;
	if (!axisTestY02(e1.k, e1.i, fez, fex, v0, v2, boxhalfsize, rad, min, max, p0, p2))
		return 0;
	if (!axisTestZ0(e1.j, e1.i, fey, fex, v0, v1, boxhalfsize, rad, min, max, p0, p1))
		return 0;

	fex = fabsf(e2.i);
	fey = fabsf(e2.j);
	fez = fabsf(e2.k);
	if (!axisTestX2(e2.k, e2.j, fez, fey, v0, v1, boxhalfsize, rad, min, max, p0, p1))
		return 0;
	if (!axisTestY1(e2.k, e2.i, fez, fex, v0, v1, boxhalfsize, rad, min, max, p0, p1))
		return 0;
	if (!axisTestZ12(e2.j, e2.i, fey, fex, v1, v2, boxhalfsize, rad, min, max, p1, p2))
		return 0;

	findMinMax(v0.i, v1.i, v2.i, min, max);
	if (min > boxhalfsize.i || max < -boxhalfsize.i)
		return 0;

	findMinMax(v0.j, v1.j, v2.j, min, max);
	if (min > boxhalfsize.j || max < -boxhalfsize.j)
		return 0;

	findMinMax(v0.k, v1.k, v2.k, min, max);
	if (min > boxhalfsize.k || max < -boxhalfsize.k)
		return 0;

	normal = vector_cross_product(e0, e1);
	if (!planeBoxOverlap(normal, v0, boxhalfsize))
		return 0;

	return 1;
}