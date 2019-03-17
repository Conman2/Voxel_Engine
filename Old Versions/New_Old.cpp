#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <SDL2/SDL.h>

//Vector Data Structure
struct vec 
{
    float x = 0; 
    float y = 0; 
    float z = 0;
    float w = 1;
};

//The Vector Class
class vector
{   
    public:
    
        //Vector Adding 
        vec add(vec vector1, vec vector2)
        {
            return {vector1.x + vector2.x, vector1.y + vector2.y, vector1.z + vector2.z, 1};
        }

        //Vector Subtracting 
        vec subtract(vec vector1, vec vector2)
        {
            return {vector1.x - vector2.x, vector1.y - vector2.y, vector1.z - vector2.z, 1};
        }

        //Vector Multiplication with Constant
        vec multiply(vec vector, float constant)
        {
            return {vector.x * constant, vector.y * constant, vector.z * constant, 1};
        }

        vec divide(vec vector, float constant)
        {
            return {vector.x / constant, vector.y / constant, vector.z / constant, 1};
        }

        //Dot Product 
        float dot_product(vec vector1, vec vector2)
        {
            return {vector1.x * vector2.x + vector1.y * vector2.y + vector1.z * vector2.z}; 
        }

        //Cross Product
        vec cross_product(vec vector1, vec vector2)
        {   
            return {vector1.y * vector2.z - vector1.z * vector2.y, vector1.z * vector2.x - vector1.x * vector2.z, vector1.x * vector2.y - vector1.y * vector2.x, 1};
        }

        //Normalization
        void normalize(vec vector)
        {   
            float length = pow(pow(vector.x, 2) + pow(vector.y, 2) + pow(vector.z, 2), 0.5);

            vector = divide(vector, length);
        }
};

//The Vector Class
class quaternion
{
    public:
        //Multiplying two quaternions using Hamilton Product  
        vec divide(vec quaternion, float constant)
        {   
            return {quaternion.x / constant, quaternion.y / constant, quaternion.z / constant, 1};
        };

        //Multiplying two quaternions using Hamilton Product  
        vec multiply(vec quaternion1, vec quaternion2)
        {   
            vec result;

            result.x = (quaternion1.w*quaternion2.x + quaternion1.x*quaternion2.w + quaternion1.y*quaternion2.z - quaternion1.z*quaternion2.y);
            result.y = (quaternion1.w*quaternion2.y - quaternion1.x*quaternion2.z + quaternion1.y*quaternion2.w + quaternion1.z*quaternion2.x);
            result.z = (quaternion1.w*quaternion2.z + quaternion1.x*quaternion2.y - quaternion1.y*quaternion2.x + quaternion1.z*quaternion2.w);
            result.w = (quaternion1.w*quaternion2.w - quaternion1.x*quaternion2.x - quaternion1.y*quaternion2.y - quaternion1.z*quaternion2.z);

            return result;
        };

        //Normalization
        void normalize(vec vector)
        {   
            float length = pow(pow(vector.x, 2) + pow(vector.y, 2) + pow(vector.z, 2) + pow(vector.w, 2), 0.5);

            vector = divide(vector, length);
        }

        //The basic quaternion structure 
        vec structure(vec axis, float angle)
        {   
            vec quaternion;

            quaternion.x = axis.x * sinf( angle/2 );
            quaternion.y = axis.y * sinf( angle/2 );
            quaternion.z = axis.z * sinf( angle/2 );
            quaternion.w = cosf( angle/2);

            normalize(quaternion);

            return quaternion;
        };

        //Seting up the new global rotation based on input axis, change in angles and total quaternion
        vec setup(vec total_quaternion, vec angle, vec x_axis, vec y_axis, vec z_axis)
        {   
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

            vec quaternion_x = structure(x_axis, angle.x);
            vec quaternion_y = structure(y_axis, angle.y);
            vec quaternion_z = structure(z_axis, angle.z);

            //Multiplying change in quaternion by universal quaternion then rotating point
            vec quaternion = multiply(multiply(multiply(quaternion_z, quaternion_y), quaternion_x), total_quaternion);

            return quaternion;
        };

        //Provides the conjugate of a quaternion
        vec conjugatation(vec quaternion)
        {
            return {-quaternion.x, -quaternion.y, -quaternion.z, quaternion.w};
        };

        //Rotates a point or vector based on: R = P*Q*P^-1
        vec rotation(vec quaternion, vec position)
        {
            vec conjugate = conjugatation(quaternion);
            vec rotated = multiply(multiply(quaternion, position), conjugate);

            return(rotated);
        }; 
};

//The Matrix Class
class matrix 
{
    public:
        //Declaring the vector class
        vector vect;

        //Matrix Vector Multiplication (note that with this function input and output must be different variables)
        vec vector_multiplication(vec vector, float matrix[4][4])
        {   
            vec result;

            result.x = vector.x* matrix[0][0] + vector.y * matrix[0][1] + vector.z * matrix[0][2] + vector.w * matrix[0][3];
            result.y = vector.x* matrix[1][0] + vector.y * matrix[1][1] + vector.z * matrix[1][2] + vector.w * matrix[1][3];
            result.z = vector.x* matrix[2][0] + vector.y * matrix[2][1] + vector.z * matrix[2][2] + vector.w * matrix[2][3];
            result.w = vector.x* matrix[3][0] + vector.y * matrix[3][1] + vector.z * matrix[3][2] + vector.w * matrix[3][3];

            return result;
        }

        //Look at Matrix
        void lookat(float matrix[4][4], vec camera_position, vec camera_direction, vec world_up, vec camera_right)
        {
            // //Calculating the Camera Looking Direction
            //vec camera_direction = vector_subtract(camera_position, camera_target); 
            //vector_normalize(camera_direction);

            //Calculating new up
            vec camera_up = vect.subtract(vect.multiply(camera_direction, vect.dot_product(world_up, camera_direction)), world_up);//Might need t0 be reversed
            vect.normalize(camera_up);

            //Calculating Camera Right Direction
            //vec camera_right = cross_product(camera_up, camera_direction);

            //Look At Matrix 
            matrix[0][0] = camera_right.x;     matrix[0][1] = camera_right.y;     matrix[0][2] = camera_right.z;     matrix[0][3] = -vect.dot_product(camera_right, camera_position);
            matrix[1][0] = camera_up.x;        matrix[1][1] = camera_up.y;        matrix[1][2] = camera_up.z;        matrix[1][3] = -vect.dot_product(camera_up, camera_position);
            matrix[2][0] = camera_direction.x; matrix[2][1] = camera_direction.y; matrix[2][2] = camera_direction.z; matrix[2][3] = -vect.dot_product(camera_direction, camera_position);
            matrix[3][0] = 0;                  matrix[3][1] = 0;                  matrix[3][2] = 0;                  matrix[3][3] = 1; 
        }

        void projection(float matrix[4][4], float camera_view_angle, float screen_height, float screen_width, float z_max_distance, float z_min_distance)
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
        void clear(float matrix[4][4])
        {
            for(int i = 0; i++; i < 5)
            {
                for(int j = 0; j++; j < 5)
                {
                    matrix[i][j] = 0.0;
                }
            }
        }  
};

//Propulsion Class
class propulsion
{
    public: 
        //
        //Piston Engine 
        //

        //Engine Levers (range 0 - 1)
        float lever_throttle = 0.0f;
        float lever_mixture = 0.0f;
        float lever_propellor = 0.0f;

        float max_rpm = 2700; 
        float max_hoursepower = 160;

        float rpm; 

        float fuel_weight_total = 500; 
        float fuel_mixture; 
        float fuel_air_ratio; 
        float fuel_flow_rate;

        float power_factor; //To account for mixture
        float power_static; 
        float power_loss; 
        float power; 

        float manifold_pressure; 
        float manifold_temperature; //Assumed to be ambient 

        //Need to be brought over from the Aircraft
        float pressure;
        float density; 

        //
        //Propellor 
        //
        float propellor_diameter;
        float propellor_thickness; 
        float propellor_cord_root;
        float propellor_cord_tip;
        float propellor_pitch_root;
        float propellor_pitch_tip; 
        float section_amount; 

        float coefficient_thrust; 
        float coefficient_power; 
        float advance_ratio; 

        //Predicting engine behvior 
        //TODO convert to SI if possible
        void piston_engine()
        {   
            rpm = max_rpm*pow(lever_throttle, 2);
            manifold_pressure = pressure + (0.04635f*lever_throttle - 0.0469f)*rpm; 
            power_static = manifold_pressure*(0.0039f*rpm - 1);
            power_loss = 0.0413f*pow(rpm, 2)/max_rpm; 

            fuel_mixture = 2.0f - lever_mixture; 
            fuel_air_ratio = fuel_mixture*density*0.1f; 
            fuel_flow_rate = (0.235f*power_static + 0.0125f*rpm - 9.69f)*lever_mixture*0.000225f;
            power = power_static*power_factor - power_loss; 
        };

        //Predicting Propellor Performance
        void propellor()
        {
            for (int i = 0; i < section_amount; i++)
            {
                float pitch = propellor_pitch_root - propellor_pitch_tip*(i/section_amount);
                float cord  = propellor_cord_root - propellor_cord_tip*(i/section_amount);
            };
        };


};

//The Object Class Used to track a object in the world space 
class object
{
    private: 
        //The Universal Up, Foward and Right
        vec world_right  = {1, 0, 0, 0};
        vec world_up     = {0, 1, 0, 0};
        vec world_foward = {0, 0, 1, 0};

        vec world_origin = {0, 0, 0, 0};

    public:
        //To acess vector and Quaternion Functions 
        vector vect; 
        quaternion quat;

        //Object Behavior 
        vec position;
        vec quaternion;
        vec euler; 

        vec velocity; 
        vec angular_velocity; 

        vec acceleration;
        vec angular_acceleration;

        //Objects local Axis 
        vec up;
        vec right;
        vec foward;

        //Updates the objects position and angle
        void update(vec delta_angle, vec delta_position)
        {      
            //Updating Position 
            position.x += delta_position.x;
            position.y += delta_position.y;
            position.z += delta_position.z;

            //Updating Euler Angle 
            euler.x += delta_angle.x;
            euler.y += delta_angle.y;
            euler.z += delta_angle.z;
            
            //Updating Quaternion Angle
            quaternion = quat.setup(quaternion, delta_angle, right, up, foward);

            //Updating Objects Local Axis from rotation
            foward = quat.rotation(quaternion, world_foward);
            vect.normalize(foward);

            up = quat.rotation(quaternion, world_up); 
            vect.normalize(up);

            right = vect.cross_product(foward, up);
            vect.normalize(right);
        };

        //Converting quaternion to euler angles (To make sure euler and quaternion representations are equal)
        void quaternion_to_euler()
        {
            double test = quaternion.x*quaternion.y + quaternion.z*quaternion.w;
            
            if (test > 0.499) { // singularity at north pole
                euler.x = 2 * atan2(quaternion.x,quaternion.w);
                euler.y = 3.1412/2;
                euler.z = 0;
                return;
            }
            if (test < -0.499) { // singularity at south pole
                euler.x = -2 * atan2(quaternion.x,quaternion.w);
                euler.y = - 3.1412/2;
                euler.z = 0;
                return;
            }

            double sqx = quaternion.x*quaternion.x;
            double sqy = quaternion.y*quaternion.y;
            double sqz = quaternion.z*quaternion.z;

            euler.x = atan2(2*quaternion.y*quaternion.w-2*quaternion.x*quaternion.z , 1 - 2*sqy - 2*sqz);
            euler.y = asin(2*test);
            euler.z = atan2(2*quaternion.x*quaternion.w-2*quaternion.y*quaternion.z , 1 - 2*sqx - 2*sqz);
        };
};

/*
class aircraft: public object 
{      
    //Note that world space is z into screen, x horizontal and y vertical 
    //For the aircraft body space, x is foward, y points down wing and z upwards 
    private:
        //Atmosheric Constants (SI Units)
        float universal_gas_constant = 286.97f;
        float troposhere_temperature_constant = 0.00651f; 
        float temperature_kelvin = 273.15f; 
        float ratio_specfic_heat = 1.4f; //Assumed constant 
        float troposhere_altitude = 11000.0f; 
        float gravity = 9.81f; 

        //Sea Level Conditions 
        float sea_level_temperature = 288.16f;
        float sea_level_pressure = 101325.0f;
        float sea_level_density = 1.225f; 

        //
        //Mass and Interial Properties 
        //
        float mass = 5000;
        float Ixx = 21000;
        float Iyy = 50000;
        float Izz = 60000;
        float Ixz = -1700;

        //Inertia Constants for Ron
        float c0 = Ixx*Izz - pow(Ixz, 2);
        float c1 = Izz/c0; 
        float c2 = Ixz/c0; 
        float c3 = c2*(Ixx - Iyy + Izz);
        float c4 = c1*(Iyy - Izz) - c2*Ixz;
        float c5 = 1/Iyy; 
        float c6 = c5*Ixz;
        float c7 = c5*(Izz - Ixx);
        float c8 = Ixx/c0; 
        float c9 = c8*(Ixx - Iyy) + c2*Ixz; 


    public: 

        // 
        //Axes 
        //
        vec wind; //?

        vec body_axis; 
        vec wind_axis;
        vec earth_axis;  

        vec aircraft_velocity;
        vec aircraft_forces; 
        vec aircraft_moments; 

        //Angle between body and wind
        vec stability_angles; 

        //
        //Geometery
        //
        float wing_area = 20.0f;
        float wing_span = 10.0f; 
        float wing_sweep = 0.0f; 
        float taper_ratio = 0.25f; 

        float lift_curve_infinite = 2.0f*3.1416f;
        float zero_lift_drag = 0.02f; 
        float zero_lift_alpha = 0.0f; //Needs to be calculated based on stuff
        float oswald_factor = 0.8f;
        float induce_drag; 

        float lift_coefficient;
        float drag_coefficient; 
        float pitching_moment_coefficient = -0.1; //Assumed constant across Alpha 

        //Calculated values 
        float lift_curve = (3.1416f*aspect_ratio)/(1.0f + pow(1 + (3.1416f*aspect_ratio)/(lift_curve_infinite*cos(wing_sweep)), 0.5)); 
        float aspect_ratio = pow(wing_span, 2)/wing_area; 
 
        //
        //Control Surfaces 
        //
        vec control_surface_angle;
        vec trim_tap_angle; 
        
        //
        //Atmoshere (can add extension into troposhere)
        //
        float density;
        float pressure;
        float dynamic_pressure; 
        float temperature;
        float speed_of_sound;

        ///
        ///Forces and Moments 
        ///

        vec propulsive_forces_body;
        vec forces_stability;
        vec aerodynamic_forces; 
        vec aerodynamic_moments;

        vec angular_acceleration_body;

        //Atmoshperic changes with altitude 
        void atmosphere(float altitude)
        {       
            if(altitude < troposhere_altitude)
            {
                temperature = sea_level_temperature - troposhere_temperature_constant*altitude;
                pressure = sea_level_pressure*pow((temperature/sea_level_temperature), 5.256); 
                density = sea_level_density*pow((temperature/sea_level_temperature), 4.256);
                speed_of_sound = pow(ratio_specfic_heat*universal_gas_constant*temperature, 0.5); 
                dynamic_pressure = 0.5f*density*pow(aircraft_velocity.z, 2); //Assuming z velocity is majority of velocity (should pythag maybe?)
            }; 

        };
        
        //Wind behavior (should return a vector)
        void wind_currents()
        {      
        }; 

        //Predicting plane behavior
        void equations_of_motion()
        {
            //Calculating Angular Acceleration Body Axis


        };

        //Currently Assuming that twist = 0 in CL calc;
        void aerodynamics()
        {   
            

            lift_coefficient = euler.y*lift_curve;
            induce_drag = pow(lift_coefficient, 2)/(3.1416f*oswald_factor*aspect_ratio);
            drag_coefficient = zero_lift_drag + induce_drag; 

        }; 

        void moments_and_forces()
        {   
            //Aerodynamic Forces
            aerodynamic_force.x = ;
            aerodynamic_force.y = ;
            aer0dynamic_force.z = ;

            //Calculating the total forces in the stability axis 
            forces_stability.y = propulsive_forces_body.y + aerodynamic_forces.y;
            forces_stability.x = propulsive_forces_body.x*cos(euler.y) - propulsive_forces_body.y*sin(euler.x) + aerodynamic_forces.z;
            forces_stability.z = propulsive_forces_body.y*cos(euler.x) + propulsive_forces_body.z*sin(euler.x) + aerodynamic_forces.y;

            //Calculating the total moments in the stability axis 

        };

        //Transforming the body to earth axes 
        void body_to_earth_axes()
        {

        };
};
*/

class aircraft: public object 
{      
    private:
        float mass = 2000; 

        float Ixx = 21000;
        float Iyy = 50000;
        float Izz = 60000;
        float Ixz = -1700;

        //Inertia Constants for Ron
        float c0 = Ixx*Izz - pow(Ixz, 2);
        float c1 = Izz/c0; 
        float c2 = Ixz/c0; 
        float c3 = c2*(Ixx - Iyy + Izz);
        float c4 = c1*(Iyy - Izz) - c2*Ixz;
        float c5 = 1/Iyy; 
        float c6 = c5*Ixz;
        float c7 = c5*(Izz - Ixx);
        float c8 = Ixx/c0; 
        float c9 = c8*(Ixx - Iyy) + c2*Ixz; 

    public: 
        //Axes and Stuff F = Force, M = Moment, V = Velocity, A = Angle
        //Assuming X is foward, Y is along right wing, Z is up though aircraft
        //Earth position and angles are stored in object class
        vec F_stab;
        vec F_erth; 
        
        vec F_aero; 
        vec F_prop;

        vec M_body;
        vec M_aero;
        vec M_prop; 

        float true_airspeed; 
        vec body_speed; 
        vec airspeed;
        vec wind;

        vec A_wind; //Stores Alpha and Beta

        //Used to convert to Earth Axis
        vec L;
        vec M; 
        vec N; 

        //Storing Time
        float old_time;
        float delta_time;
        
        void forces_and_moments(float current_time)
        {   
            //Change in time since last loop 
            delta_time = old_time - current_time; 
            old_time = current_time;

            //Quaternion Components used for Transformations
            L.x = 2*(pow(quaternion.x, 2) + pow(quaternion.y, 2)) - 1;
            L.y = 2*(quaternion.y*quaternion.z + quaternion.x*quaternion.w);
            L.z = 2*(quaternion.y*quaternion.w - quaternion.x*quaternion.z);

            M.x = 2*(quaternion.y*quaternion.z - quaternion.x*quaternion.w);
            M.y = 2*(pow(quaternion.x, 2) + pow(quaternion.z, 2)) - 1;
            M.z = 2*(quaternion.z*quaternion.w + quaternion.x*quaternion.y);

            N.x = 2*(quaternion.y*quaternion.w + quaternion.x*quaternion.z);
            N.y = 2*(quaternion.z*quaternion.w - quaternion.x*quaternion.y);
            N.z = 2*(pow(quaternion.x, 2) + pow(quaternion.w, 2)) - 1;

            //Forces from Propulsion and Aerodynamic in Stability axis
            F_stab.x = F_prop.x*cos(A_wind.y) + F_prop.z*sin(A_wind.y) + F_aero.x;  
            F_stab.y = F_prop.y + F_aero.y;
            F_stab.z = F_prop.z*cos(A_wind.y) - F_prop.x*sin(A_wind.y) + F_aero.z;  

            //Forces Transformed into Earth Axis
            F_erth.x = (L.x*cos(A_wind.y) + N.x*sin(A_wind.x))*F_stab.x + M.x*F_stab.y + (N.x*cos(A_wind.y) - L.x*sin(A_wind.y))*F_stab.z;
            F_erth.y = (L.y*cos(A_wind.y) + N.y*sin(A_wind.x))*F_stab.x + M.y*F_stab.y + (N.y*cos(A_wind.y) - L.y*sin(A_wind.y))*F_stab.z;
            F_erth.z = (L.z*cos(A_wind.y) + N.z*sin(A_wind.x))*F_stab.x + M.z*F_stab.y + (N.z*cos(A_wind.y) - L.z*sin(A_wind.y))*F_stab.z;

            //Moments around Body Axis
            M_body.x = M_aero.x*cos(A_wind.y) - M_aero.y*sin(A_wind.y) + M_prop.x; 
            M_body.y = M_aero.y + M_prop.y; 
            M_body.z = M_aero.x*cos(A_wind.y) + M_aero.y*sin(A_wind.y) + M_prop.y; 

            //Angular Accelerations
            angular_acceleration.x = M_body.x*c1 + M_body.z*c2 + (angular_velocity.x*c3 + angular_velocity.z*c4)*angular_velocity.y;
            angular_acceleration.y = M_body.y*c5 + (pow(angular_velocity.z, 2) - pow(angular_velocity.x, 2))*c6 + angular_velocity.x*angular_velocity.z*c7;
            angular_acceleration.z = M_body.z*c8 + M_body.x*c2 + (angular_velocity.x*c9 - angular_velocity.z*c3)*angular_velocity.y;

            //Earth Accelerations 
            acceleration.x = F_erth.x/mass;
            acceleration.y = F_erth.y/mass;
            acceleration.z = F_erth.z/mass;

            //Earth Velocities 
            velocity.x += acceleration.x*delta_time;
            velocity.y += acceleration.y*delta_time;
            velocity.z += acceleration.z*delta_time;

            //Angular Velocities 
            angular_velocity.x = angular_acceleration.x*delta_time;
            angular_velocity.y = angular_acceleration.y*delta_time;
            angular_velocity.z = angular_acceleration.z*delta_time;

            //Earth Positions
            position.x += velocity.x*delta_time;
            position.y += velocity.y*delta_time;
            position.z += velocity.z*delta_time;            

            //Airspeed
            airspeed.x = velocity.x - wind.x; 
            airspeed.y = velocity.y - wind.y;
            airspeed.z = velocity.z - wind.z;

            //Aircraft speed
            body_speed.x = L.x*airspeed.x + L.y*airspeed.y + L.z*airspeed.z;
            body_speed.y = M.x*airspeed.x + M.y*airspeed.y + M.z*airspeed.z;
            body_speed.z = N.x*airspeed.x + N.y*airspeed.y + N.z*airspeed.z;

            //Angles 
            true_airspeed = pow(pow(body_speed.x, 2) + pow(body_speed.y, 2) + pow(body_speed.z, 2), 0.5);
            A_wind.y = atan(body_speed.x / body_speed.z); //Alpha
            A_wind.z = atan(body_speed.y / true_airspeed); //Beta
        };
};

//Global Variables 
const int screen_width = 1000;
const int screen_height = 1000;
const float camera_view_angle = 90; 
const float z_max_distance = 1000;
const float z_min_distance = 0.1;

//Game Run Condition
int run = 1; 

//Point Definition
float points[8][4] = 
{
    {0, 0 ,0, 1},
    {1, 0 ,0, 1},
    {1, 1, 0, 1},
    {0, 1, 0, 1},
    {0, 0, 1, 1},
    {1, 0, 1, 1},
    {1, 1, 1, 1},
    {0, 1, 1, 1}
};

float updated[8][2];
float look_at[4][4];
float projection[4][4];

//
//Main Function
//

int main(int argc, char *argv[]) 
{  
    //Starting Setup
    printf("Starting Setup\n");

    //Setting up SDL
    SDL_Init(SDL_INIT_VIDEO);
    SDL_Window*window  = SDL_CreateWindow("Testing", SDL_WINDOWPOS_UNDEFINED, SDL_WINDOWPOS_UNDEFINED, screen_width, screen_height, SDL_WINDOW_RESIZABLE | SDL_WINDOW_OPENGL);
    SDL_Renderer*renderer = SDL_CreateRenderer(window, -1, SDL_RENDERER_ACCELERATED | SDL_RENDERER_PRESENTVSYNC);
    SDL_SetHint(SDL_HINT_RENDER_SCALE_QUALITY, "nearest");
    SDL_SetHint(SDL_HINT_RENDER_VSYNC, "1");
    SDL_Event window_event;

    //Declaring Objects 
    object camera; 
    vector vect;
    matrix matx; 
    quaternion quat;
    
    //Defining the contenet of the matrix as 0
    matx.clear(look_at);
    matx.clear(projection);

    //Creating Projection Matrix 
    matx.projection(projection, camera_view_angle,  screen_height,  screen_width,  z_max_distance,  z_min_distance);

    //Velocities for the Camera
    vec point;
    vec result;
    vec camera_view;
    vec delta_position;
    vec delta_angle;

    float left_speed = 0;     float yaw_left_speed = 0;
    float right_speed = 0;    float yaw_right_speed = 0;
    float down_speed = 0;     float pitch_down_speed = 0;
    float up_speed = 0;       float pitch_up_speed = 0;
    float foward_speed = 0;   float roll_right_speed = 0;
    float backward_speed = 0; float roll_left_speed = 0;

    float velocity = 8;       float angular_velocity = 2;

    float game_time = SDL_GetTicks(); 
    float time_elapsed = 0; 

    //Main Loop 
    while (run == 1)
    {
        while (SDL_PollEvent(&window_event)){
            switch( window_event.type ){
                case SDL_QUIT:
                    run = 0;
                    break;
                    // Look for a keypress
                case SDL_KEYDOWN:
                    // Check the SDLKey values and move change the coords
                    switch( window_event.key.keysym.sym ){
                        //Translation Keys 
                        case SDLK_LEFT:
                            left_speed = velocity;
                            break; 
                        case SDLK_RIGHT:
                            right_speed = velocity;                      
                            break; 
                        case SDLK_UP:
                            up_speed = velocity;
                            break; 
                        case SDLK_DOWN:
                            down_speed = velocity;
                            break;
                        case SDLK_r:
                            foward_speed = velocity;
                            break; 
                        case SDLK_f:
                            backward_speed = velocity;
                            break;                            break; 

                        //Rotation Keys  
                        case SDLK_w:
                            yaw_left_speed = angular_velocity;
                            break; 
                        case SDLK_s:
                            yaw_right_speed = angular_velocity;                      
                            break; 
                        case SDLK_a:
                            pitch_up_speed = angular_velocity;
                            break; 
                        case SDLK_d:
                            pitch_down_speed = angular_velocity;
                            break;  
                        case SDLK_e:
                            roll_right_speed = angular_velocity;
                            break; 
                        case SDLK_q:
                            roll_left_speed = angular_velocity;
                            break;                              

                        //Defualt      
                        default:    
                            break;
                    }
                    break;

                case SDL_KEYUP:
                    // Check the SDLKey values and move change the coords
                    switch( window_event.key.keysym.sym ){
                        //Translation Keys 
                        case SDLK_LEFT:
                            left_speed = 0;
                            break; 
                        case SDLK_RIGHT:
                            right_speed = 0;                      
                            break; 
                        case SDLK_UP:
                            up_speed = 0;
                            break; 
                        case SDLK_DOWN:
                            down_speed = 0;
                            break; 
                        case SDLK_r:
                            foward_speed = 0;
                            break; 
                        case SDLK_f:
                            backward_speed = 0;
                            break;   

                        //Rotation Keys  
                        case SDLK_w:
                            yaw_left_speed = 0;
                            break; 
                        case SDLK_s:
                            yaw_right_speed = 0;                      
                            break; 
                        case SDLK_a:
                            pitch_up_speed = 0;
                            break; 
                        case SDLK_d:
                            pitch_down_speed = 0;
                            break;  
                        case SDLK_e:
                            roll_right_speed = 0;
                            break; 
                        case SDLK_q:
                            roll_left_speed = 0;
                            break;                                 

                        //Defualt      
                        default:    
                            break;
                    }
                    break;
            }
        }

        //Updating Time (in seconds)
        time_elapsed = (game_time - SDL_GetTicks())/1000;
        game_time = SDL_GetTicks();

        //Updating Positions and Angles
        delta_angle.x = (yaw_right_speed - yaw_left_speed)*time_elapsed;
        delta_angle.y = (roll_right_speed - roll_left_speed)*time_elapsed;
        delta_angle.z = (pitch_up_speed - pitch_down_speed)*time_elapsed; 

        delta_position.x = (left_speed - right_speed)*time_elapsed;
        delta_position.y = (up_speed - down_speed)*time_elapsed; 
        delta_position.z = (foward_speed - backward_speed)*time_elapsed;

        //Updating the Camera Oobject 
        camera.update(delta_angle, delta_position);

        //Generating the LookAt matrix
        matx.lookat(look_at, camera.position, camera.foward, camera.up, camera.right);

        //Loop for each point
        for(int l = 0; l < 50; l++)
        {
            for(int k = 0; k < 50; k++)
            {
                for (int i = 0; i < 8; i++)
                {   
                    //Moving data
                    point.x = points[i][0] + l;
                    point.y = points[i][1];
                    point.z = points[i][2] + k; 

                    //Camera Manipulation
                    camera_view = matx.vector_multiplication(point, look_at);

                    //Projectring from 3d into 2d then normalizing  
                    result = matx.vector_multiplication(camera_view, projection);
                    result = vect.divide(result, result.w);

                    //Scaling into screen (because Projection space is from -1 to 1)
                    updated[i][0] = (result.x + 1.0)*screen_width/2;
                    updated[i][1] = (result.y + 1.0)*screen_width/2; 
                };

                //Drawing the Sides of the Cubes 
                SDL_SetRenderDrawColor(renderer, 255, 255, 255, 255);

                //Front 
                SDL_RenderDrawLine(renderer, updated[0][0], updated[0][1], updated[1][0], updated[1][1]);
                SDL_RenderDrawLine(renderer, updated[1][0], updated[1][1], updated[2][0], updated[2][1]);
                SDL_RenderDrawLine(renderer, updated[2][0], updated[2][1], updated[3][0], updated[3][1]);
                SDL_RenderDrawLine(renderer, updated[3][0], updated[3][1], updated[0][0], updated[0][1]);

                //Back
                SDL_RenderDrawLine(renderer, updated[4][0], updated[4][1], updated[5][0], updated[5][1]);
                SDL_RenderDrawLine(renderer, updated[5][0], updated[5][1], updated[6][0], updated[6][1]);
                SDL_RenderDrawLine(renderer, updated[6][0], updated[6][1], updated[7][0], updated[7][1]);
                SDL_RenderDrawLine(renderer, updated[7][0], updated[7][1], updated[4][0], updated[4][1]);

                //Left
                SDL_RenderDrawLine(renderer, updated[1][0], updated[1][1], updated[2][0], updated[2][1]);
                SDL_RenderDrawLine(renderer, updated[2][0], updated[2][1], updated[6][0], updated[6][1]);
                SDL_RenderDrawLine(renderer, updated[6][0], updated[6][1], updated[5][0], updated[5][1]);
                SDL_RenderDrawLine(renderer, updated[5][0], updated[5][1], updated[1][0], updated[1][1]);
                
                ////Right
                SDL_RenderDrawLine(renderer, updated[0][0], updated[0][1], updated[3][0], updated[3][1]);
                SDL_RenderDrawLine(renderer, updated[3][0], updated[3][1], updated[7][0], updated[7][1]);
                SDL_RenderDrawLine(renderer, updated[7][0], updated[7][1], updated[4][0], updated[4][1]);
                SDL_RenderDrawLine(renderer, updated[4][0], updated[4][1], updated[0][0], updated[0][1]);

            };
        };

        //Render the screen
        SDL_RenderPresent(renderer); 

        // Clear the current renderer
        SDL_SetRenderDrawColor(renderer, 0, 0, 0, 255);
        SDL_RenderClear(renderer);
    }

    //Destroy and Clode SDL after program finishes
    SDL_DestroyRenderer(renderer);
    SDL_DestroyWindow(window);

    SDL_Quit();

    return 0;
}