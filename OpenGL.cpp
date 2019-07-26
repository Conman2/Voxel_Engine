
//Standard Libaries 
#include <ctime>
#include <vector>
#include <fstream>
#include <strstream>
#include <algorithm>
#include <GL/glew.h>
#include <SDL2/SDL.h>

//Home Made Libaries
#include "headers/math.hpp"
#include "headers/font.hpp"

//Global Variables 
const float voxel_size = 0.003f; //This is world voxel size
const int world_voxel_limit = 150000; //This the limit to voxels rendered at one time
int temp_voxel_limit; //This is used to store the amount of voxels if it is less than the world limit
const int screen_width = 1000;
const int screen_height = 1000;

//World Space Defintions 
const vect world_origin = {0, 0, 0, 1};
const vect world_right  = {1, 0, 0, 0};
const vect world_up     = {0, 1, 0, 0};
const vect world_foward = {0, 0, 1, 0};

//Font Storage
int font_size = 3; 
bool character_array[7][5][95];

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
float fps[5]; int fps_counter = 0;  int avg_fps; float fps_timer;

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

			std::strstream s;
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

//
//Camera Class 
//
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

        //This objects Look At matrix (if it is a camera)
        float projection[4][4];

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
        };

        //Sets up projection matrix
        void intialize_projection_matrix()
        {
            matrix_projection(projection, view_angle, screen_height, screen_width, z_max_distance, z_min_distance);
        };
}; 

//
//Object Class
//
class object
{
    private: 
        //This stores the voxel infromation
        vect voxel_volume;  
        voxel AABB_corners; 

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
        void voxelization(std::string sprite_name)
        {    
            //Make the files adress
            std::string filepath = "voxelized/voxelized_"+sprite_name+std::to_string(voxel_size)+".txt"; 

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
                    //not sure what this does
                    char line[136];
                    file.getline(line, 136);

                    //Temporaly storing the line data
                    std::strstream thing;
                    thing << line;

                    if (counter > 5)
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
                float boxi = polygon_mesh.max_x - polygon_mesh.min_x; 
                float boxj = polygon_mesh.max_y - polygon_mesh.min_y; 
                float boxk = polygon_mesh.max_z - polygon_mesh.min_z; 

                //Object Half Size 
                half_size = {boxi/2, boxj/2, boxk/2}; 

                //Creates the box shape for storing the voxels 
                voxel_volume = {ceil(boxi/voxel_size), ceil(boxj/voxel_size), ceil(boxk/voxel_size)};
                vect voxel_half_size = {voxel_size/2, voxel_size/2, voxel_size/2};

                //Fill in the voxel_volume data
                file << "//World Voxel Size: \n";
                file << voxel_size << "\n"; 
                file << "//Voxel Amount in each direction (x, y, z): \n";
                file << voxel_volume.x << "\t" << voxel_volume.y << "\t" << voxel_volume.z << "\n"; 
                file << "//The Voxel Information (First 4 numbers rgba, remaining 3 numbers normal vector): \n";

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
                            voxel new_voxel; new_voxel.colour.a = 0; 

                            for(auto poly: polygon_mesh.tris)
                            {   
                                //calculate the current polygon
                                vect tri1 = {poly.p[0].x - polygon_mesh.min_x, poly.p[0].y - polygon_mesh.min_y, poly.p[0].z - polygon_mesh.min_z};
                                vect tri2 = {poly.p[1].x - polygon_mesh.min_x, poly.p[1].y - polygon_mesh.min_y, poly.p[1].z - polygon_mesh.min_z}; 
                                vect tri3 = {poly.p[2].x - polygon_mesh.min_x, poly.p[2].y - polygon_mesh.min_y, poly.p[2].z - polygon_mesh.min_z};
    
                                //Check for intersection
                                if (voxel_mesh_intersection(voxel_position, voxel_half_size, tri1, tri2, tri3) == 1){
                                    //tell it that its solid 
                                    new_voxel.colour.a = 1;
                                    new_voxel.position = left_corner_position; 

                                    //Calculate the normal
                                    new_voxel.normal = vector_normalize(vector_cross_product(vector_subtract(tri2, tri1), vector_subtract(tri3, tri1))); 

                                    //Adding the voxel information to the file 
                                    file << new_voxel.colour.r << "\t" << new_voxel.colour.g << "\t" << new_voxel.colour.b << "\t" << new_voxel.colour.a << "\t" << new_voxel.normal.x << "\t" << new_voxel.normal.y << "\t" << new_voxel.normal.z << "\t" << left_corner_position.x << "\t" << left_corner_position.y << "\t" << left_corner_position.z << "\n"; 
                                    
                                    break; 
                                };
                            }; 

                            //Add result to the voxel array
                            voxels.push_back(new_voxel); 
                        }; 
                    }; 
                };
            }; 
        };

        //Right now this renders this object (But this wont work when rendering multiple objects so will need to update )
        std::vector<voxel> project_voxels(camera camera, float projection[4][4], std::vector<voxel> voxel_projected)
        {
            //Normalized Light Direction 
            vect light_direction =  vector_normalize({0.0f, 1.0f, 1.0f});

            //Loop through the voxels and render the ones you want
            for (auto voxs: voxels)
            {   
                //
                //Moving into World Space 
                //
                vect normal_direction = quaternion_rotation(quaternion, voxs.normal); //Rotating Normal
                vect voxel_position = vector_add(quaternion_rotation(quaternion, voxs.position), position); //Rotating/Translation the Voxel Position

                //Removing unessecery voxels 
                if (vector_dot_product(vector_normalize(vector_subtract(camera.position, voxel_position)), normal_direction) > -0.35f)  
                {
                    //
                    //Moving Into View Space 
                    //
                    vect camera_view = quaternion_rotation(camera.quaternion, vector_subtract(voxel_position, camera.position));

                    //This removes Voxels behind the camera (Mirror Realm Rabbit)
                    if (camera_view.z < 0.0f)
                    {
                        //Basic Lighting 
                        float dot_product = std::min(1.0f, std::max(0.2f, vector_dot_product(light_direction, normal_direction)));
                        
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
                        temp.size = 1/vector_magnitude(vector_subtract(camera.position, voxel_position))*screen_width*voxel_size; 

                        //only create a object to calculate if in Camera View Space
                        if (temp.position.x > -1 && temp.position.x < 1 && temp.position.y > -1 && temp.position.y < 1)
                        {
                            //Looking Up the Colour from the Structure
                            temp.colour.r = dot_product; 
                            temp.colour.g = dot_product;
                            temp.colour.b = dot_product;   
                            temp.colour.a = 1.0; 

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

//
//Shader Functions 
//
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

//
//Writing to Screen
//
void create_font_array(bool character_array[7][5][95], int index, const bool font_layout[7][5])
{
    for (int i = 0; i < 7; i++)
    {
        for (int j = 0; j < 5; j++)
        {
            character_array[i][j][index] = font_layout[i][j]; 
        };
    };
}; 

void initialize_font(bool character_array[7][5][95])
{
    //Create all the textures 
    create_font_array(character_array, ' ' - 32,  blk); 
    create_font_array(character_array, '0' - 32,  zer); 
    create_font_array(character_array, '1' - 32,  one); 
    create_font_array(character_array, '2' - 32,  two); 
    create_font_array(character_array, '3' - 32,  thr); 
    create_font_array(character_array, '4' - 32,  fur); 
    create_font_array(character_array, '5' - 32,  fiv); 
    create_font_array(character_array, '6' - 32,  six); 
    create_font_array(character_array, '7' - 32,  sev); 
    create_font_array(character_array, '8' - 32,  eig);
    create_font_array(character_array, '9' - 32,  nin);
    create_font_array(character_array, 'a' - 32,  lca);
    create_font_array(character_array, 'Q' - 32,  cpq);
    create_font_array(character_array, 'W' - 32,  cpw);
    create_font_array(character_array, 'X' - 32,  cpx);
    create_font_array(character_array, 'Y' - 32,  cpy);
    create_font_array(character_array, 'Z' - 32,  cpz);
    create_font_array(character_array, ':' - 32,  col);
    create_font_array(character_array, '.' - 32,  stp);
    create_font_array(character_array, '-' - 32,  neg);    
};

std::vector<voxel> write_to_screen(std::vector<voxel> voxel_projected, bool character_array[7][5][95], std::string message, int position_x, int position_y, float font_size, colr colour)
{
    //For each Character in the String
    for(char& character: message)
    {   
        //Selecting Character
        int index = character - 32; 

        //Add the voxel_projected 
        for (int i = 0; i < 7; i++)
        {
            for (int j = 0; j < 5; j++)
            {   
                if (character_array[i][j][index] == 1)
                {
                    //Creating the voxel 
                    voxel temp;
                    
                    temp.position.x = position_x/screen_width - 1;
                    temp.position.y = position_y/screen_height - 1;
                    temp.position.z = 1;
                    temp.size = font_size; 
                    temp.colour = colour; 

                    //Add the voxel to the list 
                    voxel_projected.push_back(temp); 
                }; 
                
                //Update x position for each letter
                position_x += (5+1)*font_size; 
            };
        };
    };  

    return voxel_projected; 
}

//The main function 
int main(int argc, char *argv[]) 
{  
    printf("\nStarting Up Game... \n"); 

    ///////////////
    //Asset Setup//
    ///////////////

    //Initialize font 
    initialize_font(character_array); 

    //This creates a camera for the world 
    camera camera;
    camera.intialize_projection_matrix(); 

    //This creates and Object asigns a mesh to it and then voxelizes that mesh
    printf("Creating Game Assets...  \n"); 

    object test1; 
    test1.polygon_mesh.load_mesh("meshes/bunny.obj");
    test1.voxelization("rabbit_"); 
    test1.position = {0.0, 0, 0}; 

    object test2;
    test2.polygon_mesh.load_mesh("meshes/bunny.obj");
    test2.voxelization("rabbit_");
    test2.position = {0.5, 0, 0}; 

    object test3;
    test3.polygon_mesh.load_mesh("meshes/bunny.obj");
    test3.voxelization("rabbit_");
    test3.position = {1.0, 0, 0}; 

    object test4;
    test4.polygon_mesh.load_mesh("meshes/bunny.obj");
    test4.voxelization("rabbit_");
    test4.position = {1.5, 0, 0.0}; 

     object test5; 
    test5.polygon_mesh.load_mesh("meshes/bunny.obj");
    test5.voxelization("rabbit_"); 
    test5.position = {0.0, 0.5, 0}; 

    object test6;
    test6.polygon_mesh.load_mesh("meshes/bunny.obj");
    test6.voxelization("rabbit_");
    test6.position = {0.5, 0.5, 0}; 

    object test7;
    test7.polygon_mesh.load_mesh("meshes/bunny.obj");
    test7.voxelization("rabbit_");
    test7.position = {1.0, 0.5, 0}; 

    object test8;
    test8.polygon_mesh.load_mesh("meshes/bunny.obj");
    test8.voxelization("rabbit_");
    test8.position = {1.5, 0.5, 0.0}; 

    object test9; 
    test9.polygon_mesh.load_mesh("meshes/bunny.obj");
    test9.voxelization("rabbit_"); 
    test9.position = {0.0, 1, 0}; 

    object test10;
    test10.polygon_mesh.load_mesh("meshes/bunny.obj");
    test10.voxelization("rabbit_");
    test10.position = {0.5, 1, 0}; 

    object test11;
    test11.polygon_mesh.load_mesh("meshes/bunny.obj");
    test11.voxelization("rabbit_");
    test11.position = {1.0, 1, 0}; 

    object test12;
    test12.polygon_mesh.load_mesh("meshes/bunny.obj");
    test12.voxelization("rabbit_");
    test12.position = {1.5, 1, 0.0}; 

    object test13; 
    test13.polygon_mesh.load_mesh("meshes/bunny.obj");
    test13.voxelization("rabbit_"); 
    test13.position = {0.0, 1.5, 0}; 

    object test14;
    test14.polygon_mesh.load_mesh("meshes/bunny.obj");
    test14.voxelization("rabbit_");
    test14.position = {0.5, 1.5, 0}; 

    object test15;
    test15.polygon_mesh.load_mesh("meshes/bunny.obj");
    test15.voxelization("rabbit_");
    test15.position = {1.0, 1.5, 0}; 

    object test16;
    test16.polygon_mesh.load_mesh("meshes/bunny.obj");
    test16.voxelization("rabbit_");
    test16.position = {1.5, 1.5, 0.0}; 

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
    
    //////////////////////////////
    //Intializing Voxel VBOs VAO//
    //////////////////////////////

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
    printf("Entering Game Loop...  \n"); 

    bool run = 1;  
    int counter = 0;
    while (run == 1)
    {
        //SDL Loop Conditions (Key Presses)
        while (SDL_PollEvent(&window_event)){
            switch( window_event.type ){
                case SDL_QUIT:
                    printf("Exiting Game...  \n"); 

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

        //Dont run the loop if game is closed
        if (run != 1){break;}; 

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
            printf("FPS: %i \n", avg_fps);
            printf("Voxel Number: %i \n", temp_voxel_limit);
        }
        
       voxel_projected = write_to_screen(voxel_projected, character_array, std::to_string((int) avg_fps), 0.0f, 0.0f, font_size, {1, 1, 1, 1});

        //Render and update attitude
        voxel_projected = test1.project_voxels(camera, camera.projection, voxel_projected);
        voxel_projected = test2.project_voxels(camera, camera.projection, voxel_projected);
        voxel_projected = test3.project_voxels(camera, camera.projection, voxel_projected);
        voxel_projected = test4.project_voxels(camera, camera.projection, voxel_projected);
        voxel_projected = test5.project_voxels(camera, camera.projection, voxel_projected);
        voxel_projected = test6.project_voxels(camera, camera.projection, voxel_projected);
        voxel_projected = test7.project_voxels(camera, camera.projection, voxel_projected);
        voxel_projected = test8.project_voxels(camera, camera.projection, voxel_projected);
        voxel_projected = test9.project_voxels(camera, camera.projection, voxel_projected);
        voxel_projected = test10.project_voxels(camera, camera.projection, voxel_projected);
        voxel_projected = test11.project_voxels(camera, camera.projection, voxel_projected);
        voxel_projected = test12.project_voxels(camera, camera.projection, voxel_projected);
        voxel_projected = test13.project_voxels(camera, camera.projection, voxel_projected);
        voxel_projected = test14.project_voxels(camera, camera.projection, voxel_projected);
        voxel_projected = test15.project_voxels(camera, camera.projection, voxel_projected);
        voxel_projected = test16.project_voxels(camera, camera.projection, voxel_projected);

        //Render All of the voxels 
        int counter = 0;
        for(auto voxels : voxel_projected)
        {   
            if (counter < world_voxel_limit)
            {
                //Not sure why but that "1/" is needed. Probable becuase it cant be a value larger than 1
                translations[counter] = (vect) {voxels.position.x, voxels.position.y, 1/voxels.position.z, 1.0};
                scaling[counter] = voxels.size/screen_width*1.5; 
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

        ///////
        //FPS//
        ///////
        // //Averaging FPS
        // fps_counter ++; fps_timer += time_elapsed; 
        // if(fps_timer > 0.2) 
        // {
        //     avg_fps = fps_counter/fps_timer; 
        //     fps_timer = 0; fps_counter = 0; 

        //     printf("%f \n", avg_fps);
        // };
    };
        
    //Program Cleanup 
    printf("Deleting Graphics Assets... \n"); 

    SDL_GL_DeleteContext(context); 
    SDL_DestroyWindow(window);
    SDL_Quit();

    return 0;
};