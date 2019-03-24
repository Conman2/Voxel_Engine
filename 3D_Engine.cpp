//TODO 
//Remove Data outside Camera Volume 

#include <stdlib.h>
#include <stdio.h>
#include <math.h> 
#include <algorithm>
#include <vector>
#include <SDL2/SDL.h>

//Global Constants (SI Units)
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

//World State  Variables 
const int screen_width = 1000;
const int screen_height = 1000;
const float camera_view_angle = 90; 
const float z_max_distance = 1000;
const float z_min_distance = 0.1;

//Game Run Condition
int run = 1; 

//Acessing Map Generation
int generate_terrain(int , double , double , double , double, float **);

//These Classes are used to store Data Structures and Functions 
//The Vector Class
class vect
{   
    public:
        float i = 0; 
        float j = 0; 
        float k = 0;
        float w = 1;
        
        //Vector Adding 
        vect add(vect vector1, vect vector2)
        {
            return {vector1.i + vector2.i, vector1.j + vector2.j, vector1.k + vector2.k, 1};
        }

        //Vector Subtracting 
        vect subtract(vect vector1, vect vector2)
        {
            return {vector1.i - vector2.i, vector1.j - vector2.j, vector1.k - vector2.k, 1};
        }

        //Vector Multiplication with Constant
        vect multiply(vect vector, float constant)
        {
            return {vector.i * constant, vector.j * constant, vector.k * constant, 1};
        }

        vect divide(vect vector, float constant)
        {
            return {vector.i / constant, vector.j / constant, vector.k / constant, 1};
        }

        //Dot Product 
        float dot_product(vect vector1, vect vector2)
        {
            return {vector1.i * vector2.i + vector1.j * vector2.j + vector1.k * vector2.k}; 
        }

        //Cross Product
        vect cross_product(vect vector1, vect vector2)
        {   
            return {vector1.j * vector2.k - vector1.k * vector2.j, vector1.k * vector2.i - vector1.i * vector2.k, vector1.i * vector2.j - vector1.j * vector2.i, 1};
        }

        //Normalization
        void normalize(vect vector)
        {   
            float length = pow(pow(vector.i, 2) + pow(vector.j, 2) + pow(vector.k, 2), 0.5);

            vector = divide(vector, length);
        }
};

//The Vector Class
class quat
{
    public:
        float x = 0; 
        float y = 0; 
        float z = 0;
        float w = 1;

        //Multiplying two quaternions using Hamilton Product  
        quat divide(quat quaternion, float constant)
        {   
            return {quaternion.x / constant, quaternion.y / constant, quaternion.z / constant, 1};
        };

        //Multiplying two quaternions using Hamilton Product  
        quat multiply(quat quaternion1, quat quaternion2)
        {   
            quat result;

            result.x = (quaternion1.w*quaternion2.x + quaternion1.x*quaternion2.w + quaternion1.y*quaternion2.z - quaternion1.z*quaternion2.y);
            result.y = (quaternion1.w*quaternion2.y - quaternion1.x*quaternion2.z + quaternion1.y*quaternion2.w + quaternion1.z*quaternion2.x);
            result.z = (quaternion1.w*quaternion2.z + quaternion1.x*quaternion2.y - quaternion1.y*quaternion2.x + quaternion1.z*quaternion2.w);
            result.w = (quaternion1.w*quaternion2.w - quaternion1.x*quaternion2.x - quaternion1.y*quaternion2.y - quaternion1.z*quaternion2.z);

            return result;
        };

        //Normalization
        void normalize(quat vector)
        {   
            float length = pow(pow(vector.x, 2) + pow(vector.y, 2) + pow(vector.z, 2) + pow(vector.w, 2), 0.5);

            vector = divide(vector, length);
        }

        //The basic quaternion structure 
        quat structure(vect axis, float angle)
        {   
            quat quaternion;

            quaternion.x = axis.i * sinf( angle/2 );
            quaternion.y = axis.j * sinf( angle/2 );
            quaternion.z = axis.k * sinf( angle/2 );
            quaternion.w = cosf( angle/2 );

            normalize(quaternion);

            return quaternion;
        };

        //Seting up the new global rotation based on input axis, change in angles and total quaternion
        quat setup(quat total_quaternion, vect angle, vect x_axis, vect y_axis, vect z_axis)
        {   
            quat quaternion_x = structure(x_axis, angle.i);
            quat quaternion_y = structure(y_axis, angle.j);
            quat quaternion_z = structure(z_axis, angle.k);

            //Multiplying change in quaternion by universal quaternion then rotating point
            quat quaternion = multiply(multiply(multiply(quaternion_z, quaternion_y), quaternion_x), total_quaternion);

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
        quat conjugatation(quat quaternion)
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
        vect rotation(quat quaternion, vect position)
        {   
            quat conjugate = conjugatation(quaternion);

            return(vector_from_quaternion(multiply(multiply(quaternion, quaternion_from_vector(position)), conjugate)));
        }; 
};

//The Matrix Class
class mat 
{
    public:
        //Declaring the vector class
        vect vector;

        //Matrix Vector Multiplication (note that with this function input and output must be different variables)
        vect vector_multiplication(vect vector, float matrix[4][4])
        {   
            vect result;

            result.i = vector.i* matrix[0][0] + vector.j * matrix[0][1] + vector.k * matrix[0][2] + vector.w * matrix[0][3];
            result.j = vector.i* matrix[1][0] + vector.j * matrix[1][1] + vector.k * matrix[1][2] + vector.w * matrix[1][3];
            result.k = vector.i* matrix[2][0] + vector.j * matrix[2][1] + vector.k * matrix[2][2] + vector.w * matrix[2][3];
            result.w = vector.i* matrix[3][0] + vector.j * matrix[3][1] + vector.k * matrix[3][2] + vector.w * matrix[3][3];

            return result;
        }

        //Look at Matrix
        void lookat(float matrix[4][4], vect camera_position, vect camera_direction, vect world_up, vect camera_right)
        {
            // //Calculating the Camera Looking Direction
            //vect camera_direction = vector_subtract(camera_position, camera_target); 
            //vector_normalize(camera_direction);

            //Calculating new up
            vect camera_up = vector.subtract(vector.multiply(camera_direction, vector.dot_product(world_up, camera_direction)), world_up);//Might need t0 be reversed
            vector.normalize(camera_up);

            //Calculating Camera Right Direction
            //vect camera_right = cross_product(camera_up, camera_direction);

            //Look At Matrix 
            matrix[0][0] = camera_right.i;     matrix[0][1] = camera_right.j;     matrix[0][2] = camera_right.k;     matrix[0][3] = -vector.dot_product(camera_right, camera_position);
            matrix[1][0] = camera_up.i;        matrix[1][1] = camera_up.j;        matrix[1][2] = camera_up.k;        matrix[1][3] = -vector.dot_product(camera_up, camera_position);
            matrix[2][0] = camera_direction.i; matrix[2][1] = camera_direction.j; matrix[2][2] = camera_direction.k; matrix[2][3] = -vector.dot_product(camera_direction, camera_position);
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
        vect world_right  = {1, 0, 0, 0};
        vect world_up     = {0, 1, 0, 0};
        vect world_foward = {0, 0, 1, 0};

        vect world_origin = {0, 0, 0, 0};

    public:
        //To acess vector and Quaternion Functions 
        vect vector; 
        quat quater;

        //Object Behavior 
        vect position;
        quat quaternion;
        vect euler; 

        vect velocity; 
        vect angular_velocity; 

        vect acceleration;
        vect angular_acceleration;

        //Objects local Axis 
        vect up;
        vect right;
        vect foward;

        //Updates the objects position and angle
        void update(vect delta_angle, vect delta_position)
        {      
            //Updating Position 
            position.i += delta_position.i*right.i + delta_position.j*up.i + delta_position.k*foward.i;
            position.j += delta_position.i*right.j + delta_position.j*up.j + delta_position.k*foward.j;
            position.k += delta_position.i*right.k + delta_position.j*up.k + delta_position.k*foward.k;

            //Updating Euler Angle 
            euler.i += delta_angle.i;
            euler.j += delta_angle.j;
            euler.k += delta_angle.k;
            
            //Updating Quaternion Angle
            quaternion = quater.setup(quaternion, delta_angle, right, up, foward);

            //Updating Objects Local Axis from rotation
            foward = quater.rotation(quaternion, world_foward);
            vector.normalize(foward);

            up = quater.rotation(quaternion, world_up); 
            vector.normalize(up);

            right = vector.cross_product(foward, up);
            vector.normalize(right);
        };
};

class aircraft: public object 
{      
    private:
        //
        //Mass Properties (SI Units)
        //
        float weight; 
        float mass = 2000; 
        float Ixx = 21000;
        float Iyy = 50000;
        float Izz = 60000;
        float Ixz = -1700;

        //Inertia Constants
        float c0 = Ixx*Izz - pow(Ixz, 2);
        float c1 = Izz/c0; 
        float c2 = Ixz/c0; 
        float c3 = c2*(Ixx - Iyy + Izz);
        float c4 = c1*(Iyy - Izz) - c2*Ixz;
        float c5 = 1.0f/Iyy; 
        float c6 = c5*Ixz;
        float c7 = c5*(Izz - Ixx);
        float c8 = Ixx/c0; 
        float c9 = c8*(Ixx - Iyy) + c2*Ixz; 

        //
        //Geometery 
        //
        float wing_area = 20.0f;
        float wing_span = 10.0f; 
        float wing_mean_cord = 5.0f; 
        float wing_sweep = 0.0f; 
        float aspect_ratio; 
        float taper_ratio = 0.25f; 

        float dist_cg_ac = 1.0f;

        float tail_area = 5.0f; 
        float tail_dist_ac = 15.0f;  
        float tail_volume; 
        float tail_efficieny; 

        float fuselage_volume = 50.0f;
        float fuselage_moment_curve; 
        float fuselage_moment; 

        //User Inputs
        vect control_surface_angle;
        vect trim_tab_angle; 

        //Aerodynamics Properties
        float lift_curve_infinite = 2.0f*3.1416f;
        float zero_lift_drag = 0.02f; 
        float zero_lift_alpha = 0.0f; //Needs to be calculated based on stuff
        float zero_alpha_lift = 0.1f; 
        float zero_alpha_moment = 0.1f; 
        float oswald_factor = 0.8f;
        float induce_drag; 

        float lift_coefficient;
        float drag_coefficient; 
        float lift_coefficient_wing;
        float drag_coefficient_wing; 
        float lift_coefficient_tail;
        float drag_coefficient_tail; 
        float pitching_moment_coefficient = -0.1; //Assumed constant across Alpha        

        float lift_curve; 
        float pitching_moment_cofficient; 

        //
        //Atmoshperic Variables 
        //
        float density;
        float pressure;
        float dynamic_pressure; 
        float temperature;
        float speed_of_sound;

        //
        //Axes and Stuff F = Force, M = Moment, V = Velocity, A = Angle
        //
        //Assuming X is foward, Y is along right wing, Z is up though aircraft
        //Not i = x axis, j = y_axis, k = z axis
        //Earth position and angles are stored in object class
        vect F_stab;
        vect F_erth; 
        
        vect F_aero; 
        vect F_prop;

        vect M_body;
        vect M_aero;
        vect M_prop; 

        //Used to convert to Earth Axis
        vect L;
        vect M; 
        vect N; 

        //Rotation Crap
        quat quater;
        quat temp_quaternion; 

        float true_airspeed; 
        vect body_speed; 
        vect airspeed;
        vect wind;

        vect A_wind; //Stores Alpha and Beta

        //Storing Time
        float old_time;
        float delta_time;

    public: 
        //Intial calculations 
        void setup ()
        {   
            weight = mass*gravity; 
            lift_curve = (3.1416f*aspect_ratio)/(1.0f + pow(1 + (3.1416f*aspect_ratio)/(lift_curve_infinite*cos(wing_sweep)), 0.5)); 
            aspect_ratio = pow(wing_span, 2)/wing_area; 
            tail_volume = (tail_area*tail_dist_ac)/(wing_area*wing_mean_cord);
            fuselage_moment_curve = (2*fuselage_volume)/(wing_area*wing_mean_cord);
        };

        //Atmoshperic changes with altitude 
        //TODO Add Stratosphere calculations 
        void atmosphere()
        {       
            if(position.k < troposhere_altitude)
            {
                temperature = sea_level_temperature - troposhere_temperature_constant*position.k;
                pressure = sea_level_pressure*pow((temperature/sea_level_temperature), 5.256); 
                density = sea_level_density*pow((temperature/sea_level_temperature), 4.256);
                speed_of_sound = pow(ratio_specfic_heat*universal_gas_constant*temperature, 0.5); 
                dynamic_pressure = 0.5f*density*pow(true_airspeed, 2);
            }; 
        };

        //Calculating Aerodynamic Forces and Moments
        void aerodynamics(vect delta_control_surfaces, vect delta_trim_tabs)
        {   
            //Updating Control Surface Positions
            control_surface_angle.i += delta_control_surfaces.i;
            control_surface_angle.j += delta_control_surfaces.j;
            control_surface_angle.k += delta_control_surfaces.k;

            //Updating Trim Tab Positions
            trim_tab_angle.i += delta_trim_tabs.i;
            trim_tab_angle.j += delta_trim_tabs.j;
            trim_tab_angle.k += delta_trim_tabs.k;

            //Coefficients
            lift_coefficient_wing = A_wind.j*lift_curve + zero_alpha_lift;
            induce_drag = pow(lift_coefficient, 2)/(3.1416f*oswald_factor*aspect_ratio);
            drag_coefficient_wing = zero_lift_drag + induce_drag; 

            //Forces (Might need to add negitives???)
            F_aero.k = 0.5f * lift_coefficient * wing_area * dynamic_pressure; 
            F_aero.i = 0.5f * drag_coefficient * wing_area * dynamic_pressure;

            //Moments 
            fuselage_moment = fuselage_moment_curve*A_wind.j; 

            //Total Moments
            pitching_moment_cofficient = zero_alpha_moment + (dist_cg_ac/wing_mean_cord)*(lift_coefficient_wing + tail_efficieny*(tail_area/wing_area)*lift_coefficient_tail) - tail_efficieny*tail_volume*lift_coefficient_tail + fuselage_moment_curve;
            M_aero.j = pitching_moment_coefficient * wing_area * wing_mean_cord * dynamic_pressure;
        };

        //Calculating Propulsive Forces and Moments
        void propulsion()
        {

        };
        
        //Resolving the Forces and Moments into Earth Axis 
        void resolving_forces_and_moments(float current_time)
        {   
            //Change in time since last loop 
            //This might be in ticks so need to convert to Seconds
            delta_time = old_time - current_time; 
            old_time = current_time;

            //Quaternion Components used for Transformations
            L.i = 2*(pow(quaternion.x, 2) + pow(quaternion.y, 2)) - 1;
            L.j = 2*(quaternion.y*quaternion.z + quaternion.x*quaternion.w);
            L.k = 2*(quaternion.y*quaternion.w - quaternion.x*quaternion.z);

            M.i = 2*(quaternion.y*quaternion.z - quaternion.x*quaternion.w);
            M.j = 2*(pow(quaternion.x, 2) + pow(quaternion.z, 2)) - 1;
            M.k = 2*(quaternion.z*quaternion.w + quaternion.x*quaternion.y);

            N.i = 2*(quaternion.y*quaternion.w + quaternion.x*quaternion.z);
            N.j = 2*(quaternion.z*quaternion.w - quaternion.x*quaternion.y);
            N.k = 2*(pow(quaternion.x, 2) + pow(quaternion.w, 2)) - 1;

            //Forces from Propulsion and Aerodynamic in Stability axis
            F_stab.i = F_prop.i*cos(A_wind.j) + F_prop.k*sin(A_wind.j) + F_aero.i;  
            F_stab.j = F_prop.j + F_aero.j;
            F_stab.k = F_prop.k*cos(A_wind.j) - F_prop.i*sin(A_wind.j) + F_aero.k;  

            //Forces Transformed into Earth Axis
            F_erth.i = (L.i*cos(A_wind.j) + N.i*sin(A_wind.j))*F_stab.i + M.i*F_stab.j + (N.i*cos(A_wind.j) - L.i*sin(A_wind.j))*F_stab.k;
            F_erth.j = (L.j*cos(A_wind.j) + N.j*sin(A_wind.j))*F_stab.i + M.j*F_stab.j + (N.j*cos(A_wind.j) - L.j*sin(A_wind.j))*F_stab.k;
            F_erth.k = (L.k*cos(A_wind.j) + N.k*sin(A_wind.j))*F_stab.i + M.k*F_stab.j + (N.k*cos(A_wind.j) - L.k*sin(A_wind.j))*F_stab.k + weight;

            //Moments around Body Axis
            M_body.i = M_aero.i*cos(A_wind.j) - M_aero.j*sin(A_wind.j) + M_prop.i; 
            M_body.j = M_aero.j + M_prop.j; 
            M_body.k = M_aero.i*cos(A_wind.j) + M_aero.j*sin(A_wind.j) + M_prop.j; 

            //Angular Accelerations
            angular_acceleration.i = M_body.i*c1 + M_body.k*c2 + (angular_velocity.i*c3 + angular_velocity.k*c4)*angular_velocity.j;
            angular_acceleration.j = M_body.j*c5 + (pow(angular_velocity.k, 2) - pow(angular_velocity.i, 2))*c6 + angular_velocity.i*angular_velocity.k*c7;
            angular_acceleration.k = M_body.k*c8 + M_body.i*c2 + (angular_velocity.i*c9 - angular_velocity.k*c3)*angular_velocity.j;

            //Earth Accelerations 
            acceleration.i = F_erth.i/mass;
            acceleration.j = F_erth.j/mass;
            acceleration.k = F_erth.k/mass;

            //Earth Velocities 
            velocity.i += acceleration.i*delta_time;
            velocity.j += acceleration.j*delta_time;
            velocity.k += acceleration.k*delta_time;

            //Angular Velocities 
            angular_velocity.i += angular_acceleration.i*delta_time;
            angular_velocity.j += angular_acceleration.j*delta_time;
            angular_velocity.k += angular_acceleration.k*delta_time;

            //Angular Velocity back into Quaternion and Normalize
            temp_quaternion.x = -0.5f * (angular_velocity.i*quaternion.y + angular_velocity.j*quaternion.z + angular_velocity.k*quaternion.w); 
            temp_quaternion.y =  0.5f * (angular_velocity.i*quaternion.x - angular_velocity.j*quaternion.w + angular_velocity.k*quaternion.z);
            temp_quaternion.z =  0.5f * (angular_velocity.i*quaternion.w + angular_velocity.j*quaternion.x - angular_velocity.k*quaternion.y);
            temp_quaternion.w = -0.5f * (angular_velocity.i*quaternion.z - angular_velocity.j*quaternion.y - angular_velocity.k*quaternion.x);

            quater.normalize(temp_quaternion);

            //Update Global Quaternion
            quaternion = quater.multiply(temp_quaternion, quaternion);

            //Earth Positions
            position.i += velocity.i*delta_time;
            position.j += velocity.j*delta_time;
            position.k += velocity.k*delta_time;            

            //Airspeed
            airspeed.i = velocity.i - wind.i; 
            airspeed.j = velocity.j - wind.j;
            airspeed.k = velocity.k - wind.k;

            //Aircraft speed
            body_speed.i = L.i*airspeed.i + L.j*airspeed.j + L.k*airspeed.k;
            body_speed.j = M.i*airspeed.i + M.j*airspeed.j + M.k*airspeed.k;
            body_speed.k = N.i*airspeed.i + N.j*airspeed.j + N.k*airspeed.k;

            //Angles 
            true_airspeed = pow(pow(body_speed.i, 2) + pow(body_speed.j, 2) + pow(body_speed.k, 2), 0.5);
            A_wind.j = atan(body_speed.i / body_speed.k); //Alpha
            A_wind.k = atan(body_speed.j / true_airspeed); //Beta
        };
};

//Colour Struct
struct colour
{   
    //Colour Struct (Default is Solid White)
    int r = 0; 
    int g = 0;
    int b = 0; 
    int a = 255; 
};

//Voxel Struct 
struct voxel
{
    vect position; 
    float size; 
    colour colour;
};

//Making some colours for ron
colour none = {0, 0, 0, 0};
colour red = {255, 0, 0};
colour blue = {0, 0, 255};
colour green = {0, 255, 0};

colour vox[3][5][5] = 
{
    {
        {red, red, red, red, red},
        {red, red, red, red, red},
        {red, red, red, red, red},
        {red, red, red, red, red},  
        {red, red, red, red, red},     
    },

    {
        {none, none, none, none, none},
        {none, green, green, green, none},
        {none, green, green, green, none},
        {none, green, green, green, none}, 
        {none, none, none, none, none},
    },

    {
        {none, none, none, none, none},
        {none, none, none, none, none},
        {none, none, blue, none, none},
        {none, none, none, none, none}, 
        {none, none, none, none, none},       
    }
};

float look_at[4][4];
float projection[4][4];

//
//Main Function
//
int main(int argc, char *argv[]) 
{  
    //Starting Setup
    printf("Starting Setup\n");

    //Setting up SDL (Credit to Jack)
    SDL_Init(SDL_INIT_VIDEO);
    SDL_Window*window  = SDL_CreateWindow("Testing", SDL_WINDOWPOS_UNDEFINED, SDL_WINDOWPOS_UNDEFINED, screen_width, screen_height, SDL_WINDOW_RESIZABLE | SDL_WINDOW_OPENGL);
    SDL_Renderer*renderer = SDL_CreateRenderer(window, -1, SDL_RENDERER_ACCELERATED | SDL_RENDERER_PRESENTVSYNC);
    SDL_SetHint(SDL_HINT_RENDER_SCALE_QUALITY, "nearest");
    SDL_SetHint(SDL_HINT_RENDER_VSYNC, "1");
    SDL_Event window_event;

    //Generating Terrain 
    //generate_terrain(terrain_size, terrain_size, terrain_size, 3.0, 1, height);

    //Declaring Objects (Used for there functions)
    object camera; 
    vect vectr;
    mat  matrix; 
    quat quater;
    
    //Defining the contenet of the matrix as 0
    matrix.clear(look_at);
    matrix.clear(projection);

    //Creating Projection Matrix 
    matrix.projection(projection, camera_view_angle,  screen_height,  screen_width,  z_max_distance,  z_min_distance);

    //Velocities for the Camera
    vect point2; 
    vect point;
    vect test2;
    vect result;
    vect camera_view;
    vect delta_position;
    vect delta_angle;
    voxel temp; 

    //Initializing Speeds 
    float left_speed = 0;     float yaw_left_speed = 0;
    float right_speed = 0;    float yaw_right_speed = 0;
    float down_speed = 0;     float pitch_down_speed = 0;
    float up_speed = 0;       float pitch_up_speed = 0;
    float foward_speed = 0;   float roll_right_speed = 0;
    float backward_speed = 0; float roll_left_speed = 0;

    //Setting Angular and Translational Velocities
    float velocity = 8;       float angular_velocity = 2;

    //Initializing Game Time 
    float game_time = SDL_GetTicks(); 
    float time_elapsed = 0; 

    //Defining an SDL Rectangle 
    SDL_Rect rect; 

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
        delta_angle.i = (yaw_right_speed - yaw_left_speed)*time_elapsed;
        delta_angle.j = (roll_right_speed - roll_left_speed)*time_elapsed;
        delta_angle.k = (pitch_up_speed - pitch_down_speed)*time_elapsed; 

        delta_position.i = (left_speed - right_speed)*time_elapsed;
        delta_position.j = (down_speed - up_speed)*time_elapsed; 
        delta_position.k = (backward_speed - foward_speed)*time_elapsed;

        //Updating the Camera Oobject 
        camera.update(delta_angle, delta_position);

        //Generating the Look At matrix
        matrix.lookat(look_at, camera.position, camera.foward, camera.up, camera.right);

        //Stored projected x, y, z, size and colour (Reset each loop)
        std::vector<voxel> vox_proj;

        //Loop for each point in Sprite
        for (int k = 0; k < 3; k++)
        {   
            for (int j = 0; j < 5; j++)
            {
                for (int i = 0; i < 5; i++)
                {   
                    if (vox[k][j][i].a != 0)
                    {
                        //Moving data
                        point.i = i;
                        point.j = j;    
                        point.k = k; 

                        //This is used to move the point 1 unit perpindicular to the camera so we can find size
                        point2.i = point.i + camera.right.i;
                        point2.j = point.j + camera.right.j;
                        point2.k = point.k + camera.right.k;

                        //Camera Manipulation
                        camera_view = matrix.vector_multiplication(point, look_at);
                        vect test   = matrix.vector_multiplication(point2, look_at);

                        //Projectring from 3d into 2d then normalizing  
                        result = matrix.vector_multiplication(camera_view, projection);
                        result = vectr.divide(result, result.w);

                        test2 = matrix.vector_multiplication(test, projection);
                        test2 = vectr.divide(test2, test2.w);

                        //Storing the Projected Positions and other Voxel Characteristics 
                        temp.position.i = (result.i + 1.0)*screen_width/2; 
                        temp.position.j = (result.j + 1.0)*screen_width/2; 
                        temp.position.k = (result.k + 1.0)*screen_width/2; 

                        //Calculating square Size
                        temp.size = (pow(pow(result.i - test2.i ,2) + pow(result.j - test2.j ,2), 0.5))*screen_width/2.0f;

                        //Looking Up the Colour from the Structure
                        temp.colour.r = vox[k][j][i].r; 
                        temp.colour.g = vox[k][j][i].g;
                        temp.colour.b = vox[k][j][i].b;

                        //Stores Each Voxel in Projection space
                        vox_proj.push_back(temp);
                    };
                }; 
            };
        };
        
        //Sorting Voxels from back to front
        std::sort(vox_proj.begin(), vox_proj.end(), [](voxel vox1, voxel vox2)
        {
            return vox1.position.k > vox2.position.k;
        });

        //Render Using Painter Algorthim (Some fancy vector loop thingo thats uses Range C++11)
        for(auto temp : vox_proj)
        {   
            //Defining the Rectangle 
            rect.x = temp.position.i; 
            rect.y = temp.position.j;
            rect.w = temp.size;
            rect.h = temp.size;

            //Drawing the Rectangle 
            SDL_SetRenderDrawColor(renderer, temp.colour.r, temp.colour.g, temp.colour.b, 255);
            SDL_RenderFillRect(renderer, &rect);
        }

        //Render the screen
        SDL_RenderPresent(renderer); 

        // Clear the current renderer
        SDL_SetRenderDrawColor(renderer, 0, 0, 0, 255);
        SDL_RenderClear(renderer);
    };

    //Destroy and Clode SDL after program finishes
    SDL_DestroyRenderer(renderer);
    SDL_DestroyWindow(window);

    SDL_Quit();

    return 0;
}

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
        vect wind; //?

        vect body_axis; 
        vect wind_axis;
        vect earth_axis;  

        vect aircraft_velocity;
        vect aircraft_forces; 
        vect aircraft_moments; 

        //Angle between body and wind
        vect stability_angles; 

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
        vect control_surface_angle;
        vect trim_tap_angle; 
        
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

        vect propulsive_forces_body;
        vect forces_stability;
        vect aerodynamic_forces; 
        vect aerodynamic_moments;

        vect angular_acceleration_body;

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