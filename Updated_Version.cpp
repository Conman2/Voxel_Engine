#include "vector.h"

#include <vector>
#include <fstream>
#include <strstream>
#include <algorithm>
#include <SDL2/SDL.h>

//TODO
//Delete the mesh after voxelization 
//Move the final rendering point out of the class 

//World State Variables 
const float voxel_size = 0.002f; 
const int screen_width = 1000;
const int screen_height = 1000;
const float camera_view_angle = 90.0f; 
const float z_max_distance = 1000.0f;
const float z_min_distance = 0.1f;

//The Universal Vectors
const vect world_origin = {0, 0, 0, 0};
const vect world_right  = {1, 0, 0, 0};
const vect world_up     = {0, 1, 0, 0};
const vect world_foward = {0, 0, 1, 0};

//Colour Struct (Default is Solid White)
struct colr
{   
    int r = 0; 
    int g = 0;
    int b = 0; 
    int a = 255; 
};

struct triangle
{
	vect p[3];
};

struct mesh
{
	std::vector<triangle> tris;

    //This stores the AABB for voxelization
    float min_x = 0;
    float max_x = 0; 
    float min_y = 0; 
    float max_y = 0; 
    float min_z = 0;
    float max_z = 0; 

	bool load_mesh(std::string sFilename)
	{
		std::ifstream f(sFilename);
		if (!f.is_open())
			return false;

		// Local cache of verts
		std::vector<vect> verts;

		while (!f.eof())
		{
			char line[128];
			f.getline(line, 128);

			std::strstream s;
			s << line;

			char junk;

			if (line[0] == 'v')
			{
				vect v;
				s >> junk >> v.i >> v.j >> v.k;
				verts.push_back(v);

                //Record the largest and smallest vertix for the Axis Aligned Bounding Box
                if(v.i > max_x)
                {
                    max_x = v.i;
                } 
                else if (v.i < min_x)
                {
                    min_x = v.i;
                }

                if(v.j > max_y)
                {
                    max_y = v.j;
                } 
                else if (v.j < min_y)
                {
                    min_y = v.j;
                }

                if(v.k > max_z)
                {
                    max_z = v.k; 
                } 
                else if (v.k < min_z)
                {
                    min_z = v.k; 
                };
			}

			if (line[0] == 'f')
			{
				int f[3];
				s >> junk >> f[0] >> f[1] >> f[2];
				tris.push_back({ verts[f[0] - 1], verts[f[1] - 1], verts[f[2] - 1] });
			}
		}

		return true;
	}

};

//Voxel Data Structure 
struct voxel
{
    vect position; 
    vect normal; 
    float size; 
    colr colour;
};

//This dynamic array stores all of the worlds voxels 
std::vector<voxel> voxel_projected;

//Stores Position and Attitude and object file (Can be used for cameras as well)
class object
{
    private: 
        //This objects Look At matrix (if it is a camera)
        float look_at[4][4];

        //Defining an SDL Rectangle 
        SDL_Rect rect; 

        //This stores the projection voxels
        vect voxel_volume;  

    public:
        //Mesh Object 
        mesh polygon_mesh;  

        //This stores this objects voxels 
        std::vector<voxel> voxels; 

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
        void update_position_attitude(vect delta_angle, vect delta_position)
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
            quaternion = quaternion_setup(quaternion, delta_angle, right, up, foward);

            //Updating Objects Local Axis from rotation
            foward = quaternion_rotation(quaternion, world_foward);
            vector_normalize(foward);

            up = quaternion_rotation(quaternion, world_up); 
            vector_normalize(up);

            right = vector_cross_product(foward, up);
            vector_normalize(right);
        };

        //This updates the lookat Matrix (if this object is a camera)
        void update_look_at()
        {
            matrix_lookat(look_at, position, foward, up, right);
        };

        //This initializes a Voxel Object 
        void voxelization()
        {    
            //Find the bouding Box 
            int counter = 0; 

            //Precalculating the bounding box in terms of 0 - size instead of Min - Max
            float boxi = polygon_mesh.max_x - polygon_mesh.min_x; 
            float boxj = polygon_mesh.max_y - polygon_mesh.min_y; 
            float boxk = polygon_mesh.max_z - polygon_mesh.min_z; 

            float largets = std::max(ceil(boxi/voxel_size), std::max(ceil(boxj/voxel_size), ceil(boxk/voxel_size))); 
            voxel_volume = {largets, largets, largets};
            vect voxel_half_size = {voxel_size/2, voxel_size/2, voxel_size/2};

            //This loop goes through each voxel in the grid 
            for (int i = 0; i < voxel_volume.i; i++)
            {
                for (int j = 0; j < voxel_volume.j; j++)
                {
                    for (int k = 0; k < voxel_volume.k; k++)
                    {   
                        //Finding voxel half size and voxel middle position
                        vect voxel_position = vector_add({(float) i*voxel_size, (float) j*voxel_size, (float) k*voxel_size}, voxel_half_size);

                        //For every Polygon
                        voxel new_voxel;
                        new_voxel.colour.a = 0; 
                        for(auto poly: polygon_mesh.tris)
                        {   
                            //calculate the current polygon
                            vect tri1 = {poly.p[0].i - polygon_mesh.min_x, poly.p[0].j - polygon_mesh.min_y, poly.p[0].k - polygon_mesh.min_z};
                            vect tri2 = {poly.p[1].i - polygon_mesh.min_x, poly.p[1].j - polygon_mesh.min_y, poly.p[1].k - polygon_mesh.min_z}; 
                            vect tri3 = {poly.p[2].i - polygon_mesh.min_x, poly.p[2].j - polygon_mesh.min_y, poly.p[2].k - polygon_mesh.min_z};

                            //Check for intersection
                            if (voxel_mesh_intersection(voxel_position, voxel_half_size, tri1, tri2, tri3) == 1){
                                //tell it that its solid 
                                new_voxel.colour.a = 1;

                                //Calculate the normal
                                new_voxel.normal = vector_normalize(vector_cross_product(vector_subtract(tri2, tri1), vector_subtract(tri3, tri1))); 

                                break; 
                            };
                        }; 

                        //Add result to the voxel array
                        voxels.push_back(new_voxel); 
                    }; 
                }; 
            };
        };

        //Right now this renders this object (But this wont work when rendering multiple objects so will need to update )
        void render_voxels(SDL_Renderer *renderer, object camera, float projection[4][4], std::vector<voxel> voxel_projected)
        {
            //Light Direction 
            vect light_direction = { 0.0f, 1.0f, 1.0f };
            vector_normalize(light_direction);

            //Loop through the voxels and render the ones you want
            int i = 0; int j = 0; int k = 0; 
            for (auto voxs: voxels)
            {   
                //Should remove voxels with normals facing away from camera
                if (voxs.colour.a != 0)
                {  
                    float dp = std::min(1.0f, std::max(0.2f, vector_dot_product(light_direction, voxs.normal)));

                    //Voxel Position in space 
                    vect voxel_position = vector_add(position, {(float) i*voxel_size, (float) j*voxel_size, (float) k*voxel_size});

                    //Predeclaring Variables 
                    vect point2; 
                    vect camera_view;
                    vect result; 
                    vect sizer;
                    voxel temp; 

                    //This is used to move the point 1 unit perpindicular to the camera so we can find size at its distance from the screen
                    point2.i = voxel_position.i + camera.up.i*voxel_size;
                    point2.j = voxel_position.j + camera.up.j*voxel_size;
                    point2.k = voxel_position.k + camera.up.k*voxel_size;

                    //Camera Manipulation
                    camera_view = matrix_vector_multiplication(voxel_position, camera.look_at);
                    vect thing  = matrix_vector_multiplication(point2, camera.look_at);

                    //Projectring from 3d into 2d then normalizing  
                    result = matrix_vector_multiplication(camera_view, projection);
                    result = vector_divide(result, result.w);

                    sizer = matrix_vector_multiplication(thing, projection);
                    sizer = vector_divide(sizer, sizer.w);

                    //Storing the Projected Positions and other Voxel Characteristics 
                    temp.position.i = (result.i + 1.0)*screen_width/2; 
                    temp.position.j = (result.j + 1.0)*screen_height/2; 
                    temp.position.k = result.k; 

                    //Calculating square Size
                    temp.size = (pow(pow(result.i - sizer.i, 2) + pow(result.j - sizer.j, 2), 0.5))*screen_width*1.05;

                    //only create a object to calculate if in Camera View Space
                    if (temp.position.k > 0 && temp.position.i + temp.size > 0 && temp.position.i < screen_width && temp.position.j + temp.size > 0 && temp.position.j < screen_height)
                    {
                        //Looking Up the Colour from the Structure
                        temp.colour.r = 255*dp; 
                        temp.colour.g = 255*dp;
                        temp.colour.b = 255*dp;   

                        //Stores Each Voxel in Projection space
                        voxel_projected.push_back(temp);
                    };
                }; 

                //Deriving the new position 
                i ++; 
                if (i == voxel_volume.i)
                {
                    i = 0; 
                    j ++; 

                    if (j == voxel_volume.j)
                    {
                        j = 0; 
                        k ++; 
                    };
                }; 
            }; 

            //***********************************
            //This Should be moved somewhere else
            //***********************************
            
            //Sort using painter algortithim from close to far
            std::sort(voxel_projected.begin(), voxel_projected.end(), [](voxel vox1, voxel vox2)
            {
                return vox1.position.k < vox2.position.k;
            });

            //Render Using Painter Algorthim (Some fancy vector loop thingo thats uses Range C++11)
            for(auto temp : voxel_projected)
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

            //Clear the stored voxels (otherwise would get infinitely big)
            voxel_projected.clear(); 
        };                                                                                    
}; 

//Handels the behavior of an aircraft 
class aircraft: public object 
{
    public: 

    private:
}; 

//Handels the behavior of an jet engine
class jet_engine {
    // Effectively a sideways helicopter
    public:
        // The ratio of exhuast speed/rotational speed for the "propeller"
        float prop_pitch = 0.06;
        // Resistance from bearings, etc. in torque
        float friction_torque_max = 19000; 
        // The engines mass or whatever, resistance to acceleration
        float rotational_inertia = 10;
        
        // Core speed, changes constantly
        float rotational_speed = 0;
        // Run once per frame
        vect update(float throttle_percentage, float inlet_speed, float timestep) {
            // Calculates generated power for this update, assume linear increase from zero
            float generated_torque = 7*this->rotational_speed*throttle_percentage;
            // Force on the engine, from the air
            float air_force = pow(inlet_speed - this->rotational_speed * this->prop_pitch, 2);
            // The torque on the core caused by the acceleration of air from inlet to outlet
            float air_torque = air_force * this->prop_pitch;
            // Mechanical friction acts to bring the rpm to zero, with a maximum effect
            float mechanical_friction = std::min(friction_torque_max, rotational_speed * rotational_inertia / timestep);
            
            float rotational_change = ((generated_torque - air_torque - mechanical_friction) * timestep) / rotational_inertia;
            this->rotational_speed = this->rotational_speed + rotational_change;

            return(vect {air_force,0,0});
        }
};

//The main function 
int main(int argc, char *argv[]) 
{  
    //Setting up SDL (Credit to Jack)
    SDL_Init(SDL_INIT_VIDEO);
    SDL_Window*window  = SDL_CreateWindow("Testing", SDL_WINDOWPOS_UNDEFINED, SDL_WINDOWPOS_UNDEFINED, screen_width, screen_height, SDL_WINDOW_RESIZABLE | SDL_WINDOW_OPENGL);
    SDL_Renderer*renderer = SDL_CreateRenderer(window, -1, SDL_RENDERER_ACCELERATED | SDL_RENDERER_PRESENTVSYNC);
    SDL_SetHint(SDL_HINT_RENDER_SCALE_QUALITY, "nearest");
    SDL_SetHint(SDL_HINT_RENDER_VSYNC, "1");
    SDL_Event window_event;

    //Creating Universal rojection Matrix 
    float projection[4][4];
    matrix_clear(projection);
    matrix_projection(projection, camera_view_angle,  screen_height,  screen_width,  z_max_distance,  z_min_distance);
    
    //Creating the Camera object
    object camera;
    object test; 
    test.polygon_mesh.load_mesh("bunny.obj");
    test.voxelization(); 
    jet_engine big_engine; 

    //Initializing Speeds 
    float left_speed = 0;     float yaw_left_speed = 0;
    float right_speed = 0;    float yaw_right_speed = 0;
    float down_speed = 0;     float pitch_down_speed = 0;
    float up_speed = 0;       float pitch_up_speed = 0;
    float foward_speed = 0;   float roll_right_speed = 0;
    float backward_speed = 0; float roll_left_speed = 0;

    //Setting Angular and Translational Velocities
    float velocity = 1;       float angular_velocity = 2;

    //Stores Change in position/angle
    vect delta_position;       vect delta_angle;

    //Initializing Game Time 
    float game_time = SDL_GetTicks(); float time_elapsed = 0; 

    //Main Loop
    bool run = 1;  
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
                            break;                          

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
        };

        //Updating Time (in seconds)
        time_elapsed = (SDL_GetTicks() - game_time)/1000;
        game_time = SDL_GetTicks();

        //Updating Positions and Angles
        delta_angle.i = (yaw_right_speed - yaw_left_speed)*time_elapsed;
        delta_angle.j = (roll_right_speed - roll_left_speed)*time_elapsed;
        delta_angle.k = (pitch_up_speed - pitch_down_speed)*time_elapsed; 

        delta_position.i = (left_speed - right_speed)*time_elapsed;
        delta_position.j = (down_speed - up_speed)*time_elapsed; 
        delta_position.k = (backward_speed - foward_speed)*time_elapsed;

        //Updating the Camera Object 
        camera.update_position_attitude(delta_angle, delta_position);
        camera.update_look_at(); 

        test.render_voxels(renderer, camera, projection, voxel_projected);

        //Jacks Broken Engine
        vect engine_force = big_engine.update(1.0,0.0, time_elapsed);
    
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