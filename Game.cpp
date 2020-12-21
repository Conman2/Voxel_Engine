
//Standard Libaries 
#include <ctime>
#include <math.h>
#include <vector>
#include <fstream>
#include <sstream>
#include <glm/glm.hpp>

//Graphics Linaries
#include <GL/glew.h>
#include <SDL2/SDL.h>

//Home-Made Libaries
#include "headers/math.hpp"
#include "headers/font.hpp"
#include "headers/map_generation.hpp"

//
//TODO
//

//Draw less voxels the further away you are (LOD) (Voxel Octree)
//Normalize mesh size (rabbit reall small) also add scaling input to the mesh. Done?
//Make voxelization more efficent (Also clean-up AABB-Tri intersection function)
//Mess with aspect ratio stuff (voxel size varying with aspect ratio)
//Holes open up when the object in corner of screen because you are essentially looking at the non-existant face of the cube (the one parrel with the view vector) To fixe this you must take into account the 3d shape of a voxel
//Rotate around an objets center (currently around the edge of the object)
//Check if anypart of the Object AABB is in camera before cehcking each voxel (More efficient for lots of objects?) 
//Move all font data into a single matrix (making initializing it easier)
//Create a Initialize funtion and a cleanup function which handles OpenGl and SDL2 attributes ect.
//Move the shader programs into a class
//Data streaming to increase amount of voxels which can be rendered (http://voidptr.io/blog/2016/04/28/ldEngine-Part-1.html)
//Remove the stupid screen scale crap from the font scale thing 
//Make Character/HUD arrays textures
//Colour to models (Maybe write function to allow colouring)
//Aspect ratio being streamed to GPU in translations. should be made a uniform variable 
//Shadows (may need to store voxels in bool array)

//Global Variables 
const int terrain_size = 600; //Currently need to change this value in header file as well cus i am a retard shut up. 
const float voxel_size = 0.01f; //This is world voxel size
const int world_voxel_limit = 500000; //This the limit to voxels rendered at one time
int temp_voxel_limit; //This is used to store the amount of voxels if it is less than the world limit
const int screen_width = 1500;
const int screen_height = 1000;
float aspect_ratio = (float) screen_height/screen_width; 

//World Space Defintions 
const vect world_origin = {0, 0, 0, 1};
const vect world_right  = {1, 0, 0, 0};
const vect world_up     = {0, 1, 0, 0};
const vect world_foward = {0, 0, 1, 0};

//Font Storage
int font_size = 6; //It seems important that this number isn't odd
//bool character_array[7][5][95];
bool small_character_array[5][4][95];

//Physics Variables 
const float velocity = 1;
const float angular_velocity = 0.6;

//Default Variables 
float left_speed = 0,     yaw_left_speed = 0;
float right_speed = 0,    yaw_right_speed = 0;
float down_speed = 0,     pitch_down_speed = 0;
float up_speed = 0,       pitch_up_speed = 0;
float foward_speed = 0,   roll_right_speed = 0;
float backward_speed = 0, roll_left_speed = 0;

//Stores Change in position/angle
vect delta_position;      vect delta_angle;

//Fps Counter Variable
float fps[5]; int fps_counter = 0;  int avg_fps; float fps_timer; int avg_vox; 

//Used to store a mesh
struct mesh
{
    struct triangle
    {
        vect p[3];
    };

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
			char line[136];
			f.getline(line, 136);

			std::stringstream s;
			s << line;

			char junk;

			if (line[0] == 'v')
			{
				vect v;
				s >> junk >> v.x >> v.y >> v.z;
				verts.push_back(v);

                //Record the largest and smallest vertix for the Axis Aligned Bounding Box
                if      (v.x > max_x){ max_x = v.x;} 
                else if (v.x < min_x){ min_x = v.x;}
                if      (v.y > max_y){ max_y = v.y;} 
                else if (v.y < min_y){ min_y = v.y;}
                if      (v.z > max_z){ max_z = v.z;} 
                else if (v.z < min_z){ min_z = v.z;};
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

//This stores all seen voxels 
std::vector<voxel> voxel_projected;

////////////////
//Camera Class//
////////////////
class camera
{
    private: 
        //Camera data 
        float view_angle = 90.0f; 
        float z_max_distance = 10.0f;
        float z_min_distance = 0.1f;

    public:
        //Object Behavior 
        vect position;
        quat quaternion;
        vect euler; 

        //Objects Local Axes 
        vect up;
        vect right;
        vect foward;

        //This stores this objects voxels 
        std::vector<voxel> voxels; 

        //This objects Look At matrix (if it is a camera)
        float projection[4][4];
        float look_at[4][4]; 



        //Sets up projection matrix
        void intialize_projection_matrix()
        {
            matrix_projection(projection, view_angle, screen_height, screen_width, z_max_distance, z_min_distance);
        };
        
        //Updates the objects position and angle
        void update(vect delta_angle, vect delta_position)
        {     
            //Attitude Update
            quaternion = quaternion_setup(quaternion, delta_angle, right, up, foward) ;

            //Updating Objects Local Axis from rotation
            up     = vector_normalize(quaternion_rotation(quaternion, world_up));
            right  = vector_normalize(quaternion_rotation(quaternion, world_right)); 
            foward = vector_normalize(quaternion_rotation(quaternion, world_foward));

            //Position Update
            position.x += vector_dot_product(delta_position, right); 
            position.y += vector_dot_product(delta_position, up); 
            position.z += vector_dot_product(delta_position, foward); 

            //Converting Quaternion Angle into Euler (so our puny minds can comprehend whats happening)
            euler = vector_add(vector_multiply(quaternion_to_euler(quaternion), 180/3.141593), {180, 180, 180}); 

            //Update the LookAt Matrix 
            matrix_lookat(look_at, position, foward, up, right);
        };

        //HUD System 
        void intialize_HUD()
        {   
            //Adding Full Size HUD Texture
            for(int i = 0; i < 4; i++)
            {   
                for (int j = 0; j < 36*42; j++)
                {
                   //Assign colour
                   movingHUDfull[i][j] = movingHUDsection[i][j%36];
                };
            };               
        }; 
        
        //Draw the HUD 
        std::vector<voxel> render_HUD(std::vector<voxel> voxel_projected, int position_x, int position_y, float font_size, colr colour)
        {
            //Cropping location of scrolling texture
            int index = (euler.y/10)*36; 

            for(int i = 0; i < 8; i++)
            {   
                for (int j = 0; j < 216; j++)
                {         
                    if ( i < 4 )
                    {
                        if (movingHUDfull[i][j + index] == 1)
                        {
                            //Creating the voxel 
                            voxel temp;
                            
                            //Character Positions
                            temp.position.x = (position_x + j*font_size)/screen_width; 
                            temp.position.y = (position_y + i*font_size)/screen_height; 
                            temp.position.z = 10000; 
                            temp.size = font_size; 
                            temp.colour = colour; 

                            //Add the voxel to the list 
                            voxel_projected.push_back(temp); 
                        }
                    }
                    else 
                    {
                        if (stationaryHUD[i - 4][j] == 1)
                        {
                            //Creating the voxel 
                            voxel temp;
                            
                            //Character Positions
                            temp.position.x = (position_x + j*font_size)/screen_width; 
                            temp.position.y = (position_y + i*font_size)/screen_height; 
                            temp.position.z = 10000; 
                            temp.size = font_size; 
                            temp.colour = colour; 

                            //Add the voxel to the list 
                            voxel_projected.push_back(temp); 
                        };
                    };
                };
            };

            return voxel_projected; 
        }; 
}; 

////////////////
//Object Class//
////////////////
class object
{
    private: 
        //This stores the voxel infromation
        vect voxel_volume;  
        vect AABB;  

    public:
        //Mesh Object 
        mesh polygon_mesh;  
        vect half_size; 

        //This stores this objects voxels 
        std::vector<voxel> voxels; 

        //Object Behavior 
        vect position;
        quat quaternion;

        //Objects local Axis 
        vect up;
        vect right;
        vect foward;

        //Updates the objects position and angle
        void update_position_attitude(vect delta_angle, vect delta_position)
        {      
            //Updating Position 
            position = vector_add(position, delta_position);
            
            //Updating Quaternion Angle
            quaternion = quaternion_setup(quaternion, delta_angle, right, up, foward);

            //Updating Objects Local Axis from rotation
            up     = vector_normalize(quaternion_rotation(quaternion, world_up)); 
            right  = vector_normalize(quaternion_rotation(quaternion, world_right));
            foward = vector_normalize(quaternion_rotation(quaternion, world_foward));
        };

        //This initializes a Voxel Object 
        void voxelization(std::string sprite_name, float scale)
        {    
            //Make the files adress
            std::string filepath = "voxelized/" + sprite_name + "_scale=" + std::to_string(scale) + "_voxelsize=" + std::to_string(voxel_size) + ".txt"; 

            //Check if file exists
            if ((bool) std::ifstream(filepath) == 1)
            {   
                //Open the text file in read mode and transfer the data 
                std::ifstream file(filepath); 

                //This temporaly stores the extracted voxel
                voxel temporary; 
                
                int counter = 0; 
                while (!file.eof())
                {
                    //Not sure what this does
                    char line[136];
                    file.getline(line, 136);

                    //Temporaly storing the line data
                    std::stringstream thing;
                    thing << line;

                    if (counter == 5)
                    {
                        thing >> half_size.x >> half_size.y >> half_size.z;
                        half_size = vector_multiply(half_size, voxel_size);
                    }
                    else if (counter > 7)
                    {
                        thing >> temporary.colour.r >> temporary.colour.g >> temporary.colour.b >> temporary.colour.a >> temporary.normal.x >> temporary.normal.y >> temporary.normal.z >> temporary.position.x >> temporary.position.y >> temporary.position.z; 
                        
                        voxels.push_back(temporary); 
                    }

                    counter ++; 
                }; 
            }

            //If it doesn't exist create it and fill it with voxel information
            else 
            {
                //Create the text file in edit mode
                std::ofstream file(filepath); 

                //Find the bouding Box 
                int counter = 0; 

                //Precalculating the bounding box in terms of 0 - size instead of Min - Max
                AABB.x = polygon_mesh.max_x - polygon_mesh.min_x; 
                AABB.y = polygon_mesh.max_y - polygon_mesh.min_y; 
                AABB.z = polygon_mesh.max_z - polygon_mesh.min_z; 

                float largest_axis = std::max(std::max(AABB.x, AABB.y), AABB.z);
                vect normalized_AABB = vector_divide(AABB, largest_axis);

                //Object Half Size 
                half_size = vector_divide(normalized_AABB, 2); 

                //Creates the box shape for storing the voxels 
                voxel_volume = {ceil(normalized_AABB.x/voxel_size), ceil(normalized_AABB.y/voxel_size), ceil(normalized_AABB.z/voxel_size)};
                vect voxel_half_size = {voxel_size/2, voxel_size/2, voxel_size/2};

                //Fill in the voxel_volume data
                time_t current_time = time(0); 
                char* date_time = ctime(&current_time); 

                //This is data for the user 
                file << "//Date and Time Voxelized: \n";
                file << date_time; 
                file << "//World Voxel Size: \n";
                file << voxel_size << "\n"; 
                file << "//Voxel Amount in each direction (x, y, z): \n";
                file << voxel_volume.x << "\t" << voxel_volume.y << "\t" << voxel_volume.z << "\n"; 
                file << "//The Voxel Information (First 4 numbers rgba, next 3 numbers voxel position (x, y, z), remaining normal vector (x, y, z)): \n";

                //This loop goes through each voxel in the grid 
                for (int i = 0; i < voxel_volume.x; i++)
                {   
                    //Calculating i position
                    float position_i = i*voxel_size; 
                    for (int j = 0; j < voxel_volume.y; j++)
                    {
                        //Calculating j position
                        float position_j = j*voxel_size; 
                        for (int k = 0; k < voxel_volume.z; k++)
                        {   
                            //Finding voxel half size and voxel middle position
                            vect left_corner_position = {position_i, position_j, k*voxel_size};
                            vect voxel_position = vector_add(left_corner_position, voxel_half_size);

                            //For every Polygon
                            for(auto poly: polygon_mesh.tris)
                            {   
                                //calculate the current polygon
                                vect tri1 = {(poly.p[0].x - polygon_mesh.min_x)/AABB.x*normalized_AABB.x, (poly.p[0].y - polygon_mesh.min_y)/AABB.y*normalized_AABB.y, (poly.p[0].z - polygon_mesh.min_z)/AABB.z*normalized_AABB.z};
                                vect tri2 = {(poly.p[1].x - polygon_mesh.min_x)/AABB.x*normalized_AABB.x, (poly.p[1].y - polygon_mesh.min_y)/AABB.y*normalized_AABB.y, (poly.p[1].z - polygon_mesh.min_z)/AABB.z*normalized_AABB.z}; 
                                vect tri3 = {(poly.p[2].x - polygon_mesh.min_x)/AABB.x*normalized_AABB.x, (poly.p[2].y - polygon_mesh.min_y)/AABB.y*normalized_AABB.y, (poly.p[2].z - polygon_mesh.min_z)/AABB.z*normalized_AABB.z};
    
                                //Check for intersection
                                if (triBoxOverlap({voxel_position.x, voxel_position.y, voxel_position.z}, {voxel_half_size.x, voxel_half_size.y, voxel_half_size.z}, {tri1.x, tri1.y, tri1.z}, {tri2.x, tri2.y, tri2.z}, {tri3.x, tri3.y, tri3.z}) == 1){
                                //if (voxel_mesh_intersection(voxel_position, voxel_half_size, tri1, tri2, tri3) == 1){
                                    
                                    //Create the voxel 
                                    voxel new_voxel;

                                    new_voxel.position = left_corner_position; 

                                    //Calculate the normal
                                    new_voxel.normal = vector_normalize(vector_cross_product(vector_subtract(tri2, tri1), vector_subtract(tri3, tri1))); 

                                    //Adding the voxel information to the file 
                                    file << new_voxel.colour.r << "\t" << new_voxel.colour.g << "\t" << new_voxel.colour.b << "\t" << new_voxel.colour.a << "\t" << new_voxel.normal.x << "\t" << new_voxel.normal.y << "\t" << new_voxel.normal.z << "\t" << left_corner_position.x << "\t" << left_corner_position.y << "\t" << left_corner_position.z << "\n"; 
                                    
                                    //Add result to the voxel array
                                    voxels.push_back(new_voxel); 

                                    break; 
                                };
                            }; 
                        }; 
                    }; 
                };
            }; 
        };

        //Add terrain voxels
        void terrain_voxelization(float height_map[terrain_size][terrain_size])
        {
            for(int i = 0; i < terrain_size; i++)
            {
                for(int j = 0; j < terrain_size; j++)
                { 
                    //Create a voxel 
                    voxel new_voxel; 

                    //Define its Attrabutes
                    new_voxel.position.x = i*voxel_size; 
                    new_voxel.position.y = j*voxel_size; 
                    new_voxel.size = voxel_size; 

                    //Sea level 
                    if (height_map[i][j] < 0.45)
                    {
                        new_voxel.position.z = 0.45;
                        new_voxel.colour = {0.0, 0.41176, 0.58039, 1.0};

                        new_voxel.normal = {0.0, 0.0, 1.0, 1.0};
                    }

                    //Sand Dunes 
                    else if (height_map[i][j] < 0.47)  
                    {
                        new_voxel.colour = {1.0, 0.662745, 0.372549, 1.0};
                        new_voxel.position.z = height_map[i][j];

                        //Calculating Normal
                        vect point1 = new_voxel.position;
                        vect point2 = {(i + 1)*voxel_size, j*voxel_size, height_map[i + 1][j], 1};
                        vect point3 = {i*voxel_size, (j + 1)*voxel_size, height_map[i][j + 1], 1};
                        new_voxel.normal = vector_normalize(vector_cross_product(vector_subtract(point1, point2), vector_subtract(point1, point3))); 
                    }

                    //Grass
                    else 
                    {
                        new_voxel.colour = {0.008, 0.392, 0.251, 1.0};
                        new_voxel.position.z = height_map[i][j];

                        //Calculating Normal
                        vect point1 = new_voxel.position;
                        vect point2 = {(i + 1)*voxel_size, j*voxel_size, height_map[i + 1][j], 1};
                        vect point3 = {i*voxel_size, (j + 1)*voxel_size, height_map[i][j + 1], 1};
                        new_voxel.normal = vector_normalize(vector_cross_product(vector_subtract(point1, point2), vector_subtract(point1, point3))); 
                    };

                    //Add voxel to object 
                    voxels.push_back(new_voxel); 

                };
            };
        }; 

        //Right now this renders this object (But this wont work when rendering multiple objects so will need to update )
        std::vector<voxel> project_voxels(camera camera, float projection[4][4], std::vector<voxel> voxel_projected)
        {
            //Loop through the voxels and render the ones you want
            for (auto voxs: voxels)
            {   
                //
                //Moving into World Space 
                //

                vect normal_direction = quaternion_rotation(quaternion, voxs.normal); //Rotating Normal
                vect voxel_position = vector_add(quaternion_rotation(quaternion, voxs.position), position); //Rotating/Translation the Voxel Position

                //Removing unessecery voxels 
                if (vector_dot_product(vector_normalize(vector_subtract(camera.position, voxel_position)), normal_direction) > -1)  
                {
                    //
                    //Moving Into View Space 
                    //

                    //vect camera_view = quaternion_rotation(camera.quaternion, vector_subtract(voxel_position, camera.position));
                    vect camera_view = matrix_vector_multiplication(vector_subtract(voxel_position, camera.position), camera.look_at);

                    //This removes Voxels behind the camera (Mirror Realm Rabbit)
                    if (camera_view.z < 0.0f)
                    {
                        //
                        //Moving Into Projection Space 
                        // 

                        vect result = matrix_vector_multiplication(camera_view, projection);
                        result = vector_divide(result, result.w);

                        //Storing the Projected Positions and other Voxel Characteristics 
                        voxel temp; 

                        //This is just to test
                        temp.position.x = result.x; 
                        temp.position.y = result.y; 
                        temp.position.z = result.z; 

                        //Calculating Voxel Size (I think I made this up: 1/(distance from camera)*screen_width*voxelsize) (Needs improvement)
                        temp.size = 1/vector_magnitude(vector_subtract(camera.position, voxel_position))*screen_width*voxel_size*2; 

                        //only create a object to calculate if in Camera View Space
                        if (temp.position.x > -1 && temp.position.x < 1 && temp.position.y > -1 && temp.position.y < 1)
                        {
                            //Normalized Light Direction 
                            vect light_direction =  vector_normalize({0.4f, 0.5f, 0.0f});
                            vect light_colour = {1.0f, 1.0f, 1.0f};

                            ///Ambient Light 
                            float ambient_light_strength = 0.2; 
                            vect ambient_light = vector_multiply(light_colour, ambient_light_strength);

                            //Difused Light 
                            vect difused_light = vector_multiply(light_colour, std::max(0.0f, vector_dot_product(light_direction, normal_direction)));

                            //Specular Light 
                            float specular_strength = 0.0; 
                            vect specular_light = vector_multiply(light_colour, pow(std::max(0.0f, vector_dot_product(vector_normalize(vector_subtract(camera.position, temp.position)), vector_reflect(light_direction, voxs.normal))), 32) * specular_strength); 

                            //Total Light
                            vect total_light = vector_add(vector_add(ambient_light, difused_light), specular_light); 

                            //Looking Up the Colour from the Structure
                            temp.colour = vector_colour_multiply(total_light, voxs.colour);

                            //Stores Each Voxel in Projection space
                            voxel_projected.push_back(temp);
                        };
                    };
                }; 
            }; 
            
            // //Sort using painter algortithim from close to far
            // std::sort(voxel_projected.begin(), voxel_projected.end(), [](voxel vox1, voxel vox2)
            // {
            //     return vox1.position.z < vox2.position.z;
            // });

            return voxel_projected;
        };                                                                            
}; 

////////////////////
//Shader Functions//
////////////////////
std::string read_file(std::string file_name)
{
    //Open the text file in read mode and transfer the data 
    std::ifstream file(file_name);   

    //Return data  (and convert std::ifstream to std::string)
    return {std::istreambuf_iterator<char>(file), std::istreambuf_iterator<char>()};
};

static unsigned int compile_shader(std::string file_name, unsigned int shader_type) 
{
    //Create the OpenGL Shader
    unsigned int shader_id = glCreateShader(shader_type); 

    //Read The Souce File 
    std::string shader_data = read_file(file_name); 

    //Transfer it into the datatype OpenGL wants 
    const char* source = shader_data.c_str(); 

    //Link Source code to Shader Id
    glShaderSource(shader_id, 1, &source, nullptr); 

    //Compile the Shader
    glCompileShader(shader_id);  

    //Error Checking and Outputting
    int compile_status; 
    glGetShaderiv(shader_id, GL_COMPILE_STATUS, &compile_status);
    if (compile_status == GL_FALSE)
    {   
        //If there is an error print error message
        int length; 
        glGetShaderiv(shader_id, GL_INFO_LOG_LENGTH, &length);

        //Does something wack
        char* message = (char*)alloca(length* sizeof(char)); 
        glGetShaderInfoLog(shader_id, length, &length, message); 
        
        //Write error message to file
        std::string error_log = "error_logs/"+file_name+".txt"; 

        //Get time and date of the error message 
        time_t current_time = time(0); 
        char* date_time = ctime(&current_time); 

        //Write the error log
        std::ofstream file;
        file.open(error_log, std::ofstream::out | std::ofstream::trunc);
        file << date_time << "/n" << message; 

        //Clean Up 
        glDeleteShader(shader_id);

        return 0; 
    };

    //Return the compiled shader
    return shader_id; 
};

static unsigned int create_shaders(std::string vertex_shader, std::string fragment_shader)
{   
    //Create the Overall Program
    unsigned int shader_program = glCreateProgram(); 

    //Create and Compile the Shaders 
    static unsigned int vertex_shader_id = compile_shader(vertex_shader, GL_VERTEX_SHADER); 
    static unsigned int fragment_shader_id = compile_shader(fragment_shader, GL_FRAGMENT_SHADER); 

    //Link the shader programs 
    glAttachShader(shader_program, vertex_shader_id);
    glAttachShader(shader_program, fragment_shader_id);
    glLinkProgram(shader_program); 
    glValidateProgram(shader_program);

    //Delete the Unneeded shader data
    glDeleteShader(vertex_shader_id);
    glDeleteShader(fragment_shader_id);

    //Return the final program
    return shader_program; 
}; 

/////////////////////
//Writing to Screen//
/////////////////////
// void create_font_array(bool character_array[7][5][95], int index, const bool font_layout[7][5])
// {
//     for (int i = 0; i < 7; i++)
//     {
//         for (int j = 0; j < 5; j++)
//         {
//             character_array[i][j][index] = font_layout[i][j]; 
//         };
//     };
// }; 

void create_font_array_small(bool character_array[5][4][95], int index, const bool font_layout[5][4])
{
    for (int i = 0; i < 5; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            character_array[i][j][index] = font_layout[i][j]; 
        };
    };
}; 

void initialize_font()
{
    //Move the characters into a big array 
    // create_font_array(character_array, ' ' - 32,  blk); 
    // create_font_array(character_array, '0' - 32,  zer); 
    // create_font_array(character_array, '1' - 32,  one); 
    // create_font_array(character_array, '2' - 32,  two); 
    // create_font_array(character_array, '3' - 32,  thr); 
    // create_font_array(character_array, '4' - 32,  fur); 
    // create_font_array(character_array, '5' - 32,  fiv); 
    // create_font_array(character_array, '6' - 32,  six); 
    // create_font_array(character_array, '7' - 32,  sev); 
    // create_font_array(character_array, '8' - 32,  eig);
    // create_font_array(character_array, '9' - 32,  nin);
    // create_font_array(character_array, 'a' - 32,  lca);
    // create_font_array(character_array, 'Q' - 32,  cpq);
    // create_font_array(character_array, 'W' - 32,  cpw);
    // create_font_array(character_array, 'X' - 32,  cpx);
    // create_font_array(character_array, 'Y' - 32,  cpy);
    // create_font_array(character_array, 'Z' - 32,  cpz);
    // create_font_array(character_array, ':' - 32,  col);
    // create_font_array(character_array, '.' - 32,  stp);
    // create_font_array(character_array, '-' - 32,  neg);    

    //Move the characters into a big array  
    create_font_array_small(small_character_array, ' ' - 32,  sml_blk); 
    create_font_array_small(small_character_array, '0' - 32,  sml_zer); 
    create_font_array_small(small_character_array, '1' - 32,  sml_one); 
    create_font_array_small(small_character_array, '2' - 32,  sml_two); 
    create_font_array_small(small_character_array, '3' - 32,  sml_thr); 
    create_font_array_small(small_character_array, '4' - 32,  sml_fur); 
    create_font_array_small(small_character_array, '5' - 32,  sml_fiv); 
    create_font_array_small(small_character_array, '6' - 32,  sml_six); 
    create_font_array_small(small_character_array, '7' - 32,  sml_sev); 
    create_font_array_small(small_character_array, '8' - 32,  sml_eig);
    create_font_array_small(small_character_array, '9' - 32,  sml_nin);
};

std::vector<voxel> write_to_screen(std::vector<voxel> voxel_projected, bool character_array[5][4][95], std::string message, float position_x, float position_y, float font_size, colr colour)
{
    //Moving from screen co-oridinates to OpenGl co-oridinates
    position_x = position_x - screen_width; 
    position_y = position_y - screen_height; 

    //For each Character in the String
    for(char& character: message)
    {   
        //Selecting Character
        int index = character - 32; 

        //Add the voxel_projected 
        for (int i = 0; i < 5; i++)
        {
            for (int j = 0; j < 4; j++)
            {   
                if (character_array[i][j][index] == 1)
                {
                    //Creating the voxel 
                    voxel temp;
                    
                    //Character Positions
                    temp.position.x = (position_x + j*font_size)/screen_width; 
                    temp.position.y = (position_y + i*font_size)/screen_height; 
                    temp.position.z = 10000; 
                    temp.size = font_size; 
                    temp.colour = colour; 

                    //Add the voxel to the list 
                    voxel_projected.push_back(temp); 
                }; 
            };
        };

        //Update x position for each letter
        position_x += 6*font_size; 
    };  

    return voxel_projected; 
}

//The main function 
int main(int argc, char *argv[]) 
{  
    printf("\nStarting Game... \n"); 

    ///////////////
    //Asset Setup//
    ///////////////
    
    //Create Terrrain
    printf("Creating Terrain...  \n"); 

    // float height_map[terrain_size][terrain_size];
    // generate_terrain(terrain_size, terrain_size, terrain_size, 3.0, 1, height_map, 1);
    // object terrain; 
    // terrain.terrain_voxelization(height_map);

    //Initialize font 
    initialize_font(); 

    //This creates a camera for the world 
    camera camera;
    camera.intialize_projection_matrix(); 
    camera.intialize_HUD(); 

    //This creates and Object asigns a mesh to it and then voxelizes that mesh
    printf("Creating Game Assets...  \n"); 

    object test1; 
    test1.polygon_mesh.load_mesh("meshes/bunny.obj");
    test1.voxelization("bunny_", 1);  
    test1.quaternion = {0, 0.7068252, 0, 0.7073883}; 
    test1.position =  {-test1.half_size.y, 0, 0};

    /////////////////////////
    //Setting up SDL/OpenGL//
    /////////////////////////
    printf("Initializing SDL2 and OpenGL...  \n"); 

    SDL_Init(SDL_INIT_VIDEO);
    SDL_Window* window = SDL_CreateWindow("Game Window", SDL_WINDOWPOS_UNDEFINED, SDL_WINDOWPOS_UNDEFINED, screen_width, screen_height, SDL_WINDOW_OPENGL);
    SDL_GLContext context = SDL_GL_CreateContext(window);  
    SDL_Event window_event; 
	SDL_GL_SetAttribute(SDL_GL_CONTEXT_PROFILE_MASK, SDL_GL_CONTEXT_PROFILE_CORE); 
	SDL_GL_SetAttribute(SDL_GL_CONTEXT_MAJOR_VERSION, 4); 
	SDL_GL_SetAttribute(SDL_GL_CONTEXT_MINOR_VERSION, 4); 
	SDL_GL_SetAttribute(SDL_GL_DOUBLEBUFFER, 1); 
    SDL_GL_SetSwapInterval(1);  

    ///////////////////
    //Setting up GLEW//
    ///////////////////
    printf("Intializing GLEW...  \n"); 

    glewInit(); 
    glEnable(GL_DEPTH_TEST); 
    glClearColor(0.0, 0.0, 0.0, 1.0); 
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT); 
    SDL_GL_SwapWindow(window); 

    //////////////////////
    //Setting Up Shaders//
    //////////////////////
    printf("Compiling Shader Programs...  \n"); 

    static unsigned int shader_program = create_shaders("shaders/vertex.glsl", "shaders/fragment.glsl");
    glUseProgram(shader_program);

    ////////////////////////
    //Intializing Textures//
    ////////////////////////
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_NEAREST);
    
    //////////////////////////////
    //Intializing Voxel VBOs VAO//
    //////////////////////////////

    printf("Creating VAO and VBOs...  \n"); 

    //This stores the vertex data for the mesh 
    const GLfloat voxel[4][3] =
    {
        { -0.5, 0.5,  0.5 },
        { -0.5, -0.5, 0.5 }, 
        {  0.5, 0.5,  0.5 }, 
        { 0.5,  -0.5, 0.5 },
    };

    //The adress for the buffer objects 
    unsigned int VBO[4], VAO[1]; 
    
    //Create 4 VBOs and specify the array which stores there adresses
    glGenBuffers(4, VBO);  

    //Create the VAO
    glGenVertexArrays(1, VAO); 

    //Bind the VAO as our currently used object 
    glBindVertexArray(VAO[0]);

    /////////////////////////////////////////
    //Generate the VBO to store vertex data//
    /////////////////////////////////////////

    // Bind our first VBO as being the active buffer and storing vertex attributes (coordinates)
    glBindBuffer(GL_ARRAY_BUFFER, VBO[0]);
    
    // Copy the vertex data from dice to our buffer (4 and 3 are size of mesh array)
    uint32_t buffer_size =  (4 * 3) * sizeof(GLfloat);
    glBufferData(GL_ARRAY_BUFFER, buffer_size, voxel, GL_STATIC_DRAW);

    //This tells OpenGL the data layout
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(float) * 3, 0);

    //Enable the VBO within the VAO 
    glEnableVertexAttribArray(0);
    
    //////////////////////////////////////////////
    //Generate the VBO to store translation data//
    //////////////////////////////////////////////

    //Generate the buffer which will be used to store the translations 
    typedef std::vector<vect> trans;
    trans translations(world_voxel_limit);

    // Bind our first VBO as being the active buffer and storing vertex translations
    glBindBuffer(GL_ARRAY_BUFFER, VBO[1]);
    
    // Copy the vertex data from dice to our buffer 
    buffer_size =  world_voxel_limit * sizeof(vect);
    glBufferData(GL_ARRAY_BUFFER, buffer_size, &translations[0], GL_DYNAMIC_DRAW);

    //This tells OpenGL the data layout
    glVertexAttribPointer(1, 4, GL_FLOAT, GL_FALSE, sizeof(vect), 0); 
    
    //Enable the VBO withing the VAO 
    glEnableVertexAttribArray(1);

    //Important for instanced rendering 
    glVertexAttribDivisor(1, 1);

    //////////////////////////////////////////
    //Generate the VBO to store scaling data//
    //////////////////////////////////////////

    //Generate the buffer which will be used to store the translations 
    typedef std::vector<float> scale;
    scale scaling(world_voxel_limit); 

    // Bind our first VBO as being the active buffer and storing vertex translations
    glBindBuffer(GL_ARRAY_BUFFER, VBO[2]);
    
    // Copy the vertex data from dice to our buffer 
    buffer_size =  world_voxel_limit * sizeof(float);
    glBufferData(GL_ARRAY_BUFFER, buffer_size, &scaling[0], GL_DYNAMIC_DRAW);

    //This tells OpenGL the data layout
    glVertexAttribPointer(2, 1, GL_FLOAT, GL_FALSE, sizeof(float), 0); 
    
    //Enable the VBO withing the VAO 
    glEnableVertexAttribArray(2);

    //Important for instanced rendering 
    glVertexAttribDivisor(2, 1);

    /////////////////////////////////////////
    //Generate the VBO to store colour data//
    /////////////////////////////////////////

    //Generate the buffer which will be used to store the translations 
    typedef std::vector<colr> colour;
    colour colours(world_voxel_limit);

    // Bind our first VBO as being the active buffer and storing vertex translations
    glBindBuffer(GL_ARRAY_BUFFER, VBO[3]);
    
    // Copy the vertex data from dice to our buffer 
    buffer_size =  world_voxel_limit * sizeof(colr);
    glBufferData(GL_ARRAY_BUFFER, buffer_size, &colours[0], GL_DYNAMIC_DRAW);

    //This tells OpenGL the data layout
    glVertexAttribPointer(3, 4, GL_FLOAT, GL_FALSE, sizeof(colr), 0); 
    
    //Enable the VBO withing the VAO 
    glEnableVertexAttribArray(3);

    //Important for instanced rendering 
    glVertexAttribDivisor(3, 1);

    //Unbind the VAO
    glBindVertexArray(0);

    /////////////
    //Main Loop//
    /////////////

    //Initializing Game Time 
    float game_time = SDL_GetTicks(); float time_elapsed = 0; 

    //Defining Exit condition
    bool run = 1;  
    int counter = 0;

    printf("Entering Game Loop...  \n"); 

    while (run == 1)
    {
        //SDL Loop Conditions (Key Presses)
        while (SDL_PollEvent(&window_event)){
            switch( window_event.type ){
                case SDL_QUIT:
                    printf("\nQuiting Game...  \n"); 
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
                        case SDLK_q:
                            yaw_left_speed = angular_velocity;
                            break; 
                        case SDLK_e:
                            yaw_right_speed = angular_velocity;                      
                            break; 
                        case SDLK_w:
                            pitch_up_speed = angular_velocity;
                            break; 
                        case SDLK_s:
                            pitch_down_speed = angular_velocity;
                            break;  
                        case SDLK_a:
                            roll_right_speed = angular_velocity;
                            break; 
                        case SDLK_d:
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
                        case SDLK_q:
                            yaw_left_speed = 0;
                            break; 
                        case SDLK_e:
                            yaw_right_speed = 0;                      
                            break; 
                        case SDLK_w:
                            pitch_up_speed = 0;
                            break; 
                        case SDLK_s:
                            pitch_down_speed = 0;
                            break;  
                        case SDLK_a:
                            roll_right_speed = 0;
                            break; 
                        case SDLK_d:
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
        delta_angle.x = (pitch_up_speed - pitch_down_speed)*time_elapsed; 
        delta_angle.y = (yaw_right_speed - yaw_left_speed)*time_elapsed;
        delta_angle.z = -(roll_right_speed - roll_left_speed)*time_elapsed;

        delta_position.x = (left_speed - right_speed)*time_elapsed;
        delta_position.y = -(down_speed - up_speed)*time_elapsed; 
        delta_position.z = (backward_speed - foward_speed)*time_elapsed;
    
        //Updating the Camera Object 
        camera.update(delta_angle, delta_position);

        //Write to screen average FPS
        fps_counter ++; fps_timer += time_elapsed; 
        if(fps_timer > 0.3) 
        {
            avg_fps = (int) fps_counter/fps_timer; 
            fps_timer = 0; fps_counter = 0; 
            avg_vox = temp_voxel_limit; 
        }
        
        //Render Writing to Screen 
        voxel_projected = write_to_screen(voxel_projected, small_character_array, std::to_string((int) avg_fps), 10.0f,                10.0f, font_size, {1, 1, 1, 1});
        voxel_projected = write_to_screen(voxel_projected, small_character_array, std::to_string((int) avg_vox), 10.0f, screen_height*2 - 50, font_size, {1, 1, 1, 1});

        //Render HUD to Screen
        voxel_projected = camera.render_HUD(voxel_projected, -screen_width/2 - 150, -900, font_size, {1, 1, 1, 1}); 

        //Project the terrain voxels 
        //voxel_projected = terrain.project_voxels(camera, camera.projection, voxel_projected); 

        //Render and update attitude
        voxel_projected = test1.project_voxels(camera, camera.projection, voxel_projected);

        //Render Setup
        counter = 0;
        for(auto voxels : voxel_projected)
        {   
            if (counter < world_voxel_limit)
            {
                //Not sure why but that "1/" is needed. Probable becuase it cant be a value larger than 1
                translations[counter] = (vect) {voxels.position.x, voxels.position.y, 1/voxels.position.z, aspect_ratio};
                scaling[counter] = voxels.size/((screen_width + screen_height)/2);
                colours[counter++] = voxels.colour; 
            }
            else 
            {
                break; 
            }
        }

        //Record how many voxels were used
        temp_voxel_limit = counter; 

        //Clear the stored voxels 
        voxel_projected.clear(); 

        ///////////////
        //Update Data//
        ///////////////

        //Bind the VAO
        glBindVertexArray(VAO[0]);

        //Update the translation buffer storage
        glBindBuffer(GL_ARRAY_BUFFER, VBO[1]);
        glBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(vect)*temp_voxel_limit, &translations[0]); 

        //Update the scaling buffer storage
        glBindBuffer(GL_ARRAY_BUFFER, VBO[2]);
        glBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(float)*temp_voxel_limit, &scaling[0]); 

        //Update the colour buffer storage
        glBindBuffer(GL_ARRAY_BUFFER, VBO[3]);
        glBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(colr)*temp_voxel_limit, &colours[0]); 

        ///////////////
        //Render Call//
        ///////////////

        //OpenGL Clear Render
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        //Instanced Array
        glDrawArraysInstanced(GL_TRIANGLE_STRIP, 0, 4, temp_voxel_limit);

        //Unbind the VAO
        glBindVertexArray(0);

        //Swap Rendering Buffers
        SDL_GL_SwapWindow(window);
    };
    
    ///////////////////
    //Program Cleanup//
    ///////////////////

    printf("Deleting Graphics Assets... \n"); 

    SDL_GL_DeleteContext(context); 
    SDL_DestroyWindow(window);
    SDL_Quit();

    printf("Closing... \n"); 

    return 0;
};