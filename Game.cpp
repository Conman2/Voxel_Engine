//Standard Libaries 
#include <ctime>
#include <math.h>
#include <list>
#include <vector>
#include <fstream>
#include <sstream>
#include <glm/glm.hpp>

//Graphics Libaries
#include <GL/glew.h>
#include <SDL2/SDL.h>

//Home-Made Libary
#include "headers/math.hpp"

////////
//TODO//
////////

//Proper Terrain
//Draw less voxels the further away you are (LOD) (Voxel Octree)
//Rotate around an objets center (currently around the edge of the object)
//Check if anypart of the Object AABB is in camera before cehcking each voxel (More efficient for lots of objects???) 
//Create a Initialize funtion and a cleanup function which handles OpenGl and SDL2 attributes ect.
//Make Character/HUD arrays textures
//Colour to models (Maybe write function to allow colouring)
//Shadows? (may need to store voxels in Octree)
//Improve Lighting 
//Reduce the colour resolution so that it feels more like pixel art maybe.
//Fix the Error logs for the shader progermas
//Add a render check on CPU to see if obeject is in camera view or not?
//Move keypress checks into function, maybe

//Global Variables
const float voxel_size = 0.015f; //This is world voxel size
const int screen_width = 1500;
const int screen_height = 1000;
float aspect_ratio = (float) screen_height/screen_width; 
float scaling_factor = 0.02; //screen_width * voxel_size * 2;

//World Space Defintions 
const vect world_origin = {0, 0, 0, 1};
const vect world_right  = {1, 0, 0, 0};
const vect world_up     = {0, 1, 0, 0};
const vect world_foward = {0, 0, 1, 0};

//Camera Variables 
const float velocity = 1.5;
const float angular_velocity = 1.25;

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
        float z_max_distance = 100.0f;
        float z_min_distance = 0.1f;

    public:
        //Object Behavior 
        vect position;
        quat quaternion;
        quat rotation_quaternion;
        vect euler; 

        //Objects Local Axes 
        vect up = world_up;
        vect right = world_right;
        vect foward = world_foward;

        //This stores this objects voxels 
        std::vector<voxel> voxels; 

        //This objects Look At matrix (if it is a camera)
        float projection[4][4];
        float look_at[4][4]; 

        //Sets up projection matrix
        void intialize_projection_matrix()
        {   
            //Initialize projection matrix
            matrix_projection(projection, view_angle, screen_height, screen_width, z_max_distance, z_min_distance);
        };
        
        //Updates the camera position and angle
        void update(vect delta_angle, vect delta_position)
        {     
            //Attitude Update
            quat quaternion_x = quaternion_structure(right, delta_angle.x);
            quat quaternion_y = quaternion_structure(up, delta_angle.y);
            quat quaternion_z = quaternion_structure(foward, delta_angle.z);

            //Multiplying change in quaternion by universal quaternion then rotating point
            rotation_quaternion = quaternion_multiply(quaternion_multiply(quaternion_x, quaternion_y), quaternion_z);
            quaternion = quaternion_multiply(rotation_quaternion, quaternion);

            //Updating Objects Local Axis from rotation
            up     = vector_normalize(quaternion_rotation(quaternion, world_up));
            right  = vector_normalize(quaternion_rotation(quaternion, world_right)); 
            foward = vector_normalize(quaternion_rotation(quaternion, world_foward));

            position = vector_add(vector_add(vector_add(vector_multiply(right, delta_position.x), vector_multiply(up, delta_position.y)), vector_multiply(foward, delta_position.z)), position);

            //Converting Quaternion Angle into Euler
            euler = vector_add(vector_multiply(quaternion_to_euler(quaternion), 180/3.141593), {180, 180, 180}); 

            //Update the LookAt Matrix 
            matrix_lookat(look_at, position, foward, up, right);
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
        std::vector<voxel> voxel_object; 
        int voxel_amount;

        //Object Behavior 
        vect position;
        quat rotation_quaternion;
        quat quaternion;
        float model[4][4];
        float rotation[4][4];

        //Objects local Axis 
        vect up = world_up;
        vect right = world_right;
        vect foward = world_foward;

        //Render Buffers
        unsigned int VAO[1]; 
        unsigned int VBO[3];

        //Voxel Data
        const GLfloat mesh[4][3] =
        {
            {(GLfloat) -0.5 * scaling_factor, (GLfloat) 0.5 / aspect_ratio  * scaling_factor, 0.5 },
            {(GLfloat) 0.5  * scaling_factor, (GLfloat) 0.5 / aspect_ratio  * scaling_factor, 0.5 }, 
            {(GLfloat) 0.5  * scaling_factor, (GLfloat) -0.5 / aspect_ratio * scaling_factor, 0.5 },
            {(GLfloat) -0.5 * scaling_factor, (GLfloat) -0.5 / aspect_ratio * scaling_factor, 0.5 }, 
        };

        //Initialise object 
        void initialise()
        {
            //Prepare Data
            std::vector<vect> translations(voxel_amount);
            std::vector<vect> normals(voxel_amount);

            //Populate Data
            int counter = 0;
            for(auto voxels : voxel_object)
            {
                translations[counter] = voxels.position;
                normals[counter] = {voxels.normal.x, voxels.normal.y, voxels.normal.z, 0.0}; //Zero cancels out translation in model matrix 
                counter++;
            };

            //Create the VAO and VBO
            glGenBuffers(3, VBO);
            glGenVertexArrays(1, VAO); 

            //Set the VAO to the current use item
            glBindVertexArray(VAO[0]); 

            /////////////////////////////////////////
            //Generate the VBO to store vertex data//
            /////////////////////////////////////////

            // Bind our first VBO as being the active buffer and storing vertex attributes (coordinates)
            glBindBuffer(GL_ARRAY_BUFFER, VBO[0]);
            
            // Copy the vertex data from dice to our buffer (4 and 3 are size of mesh array)
            uint32_t buffer_size =  (4 * 3) * sizeof(GLfloat);
            glBufferData(GL_ARRAY_BUFFER, buffer_size, mesh, GL_STATIC_DRAW);

            //This tells OpenGL the data layout
            glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(float) * 3, 0);

            //Enable the VBO within the VAO 
            glEnableVertexAttribArray(0);
            
            //////////////////////////////////////////////
            //Generate the VBO to store translation data//
            //////////////////////////////////////////////

            // Bind our first VBO as being the active buffer and storing vertex translations
            glBindBuffer(GL_ARRAY_BUFFER, VBO[1]);
            
            // Copy the vertex data from dice to our buffer 
            buffer_size =  voxel_amount * sizeof(vect);
            glBufferData(GL_ARRAY_BUFFER, buffer_size, &translations[0], GL_STATIC_DRAW);

            //This tells OpenGL the data layout
            glVertexAttribPointer(1, 4, GL_FLOAT, GL_FALSE, sizeof(vect), 0); 
            
            //Enable the VBO within the VAO 
            glEnableVertexAttribArray(1);

            //Important for instanced rendering 
            glVertexAttribDivisor(1, 1);

            /////////////////////////////////////////
            //Generate the VBO to store noraml data//
            /////////////////////////////////////////

            // Bind our first VBO as being the active buffer and storing vertex translations
            glBindBuffer(GL_ARRAY_BUFFER, VBO[2]);
            
            // Copy the vertex data from dice to our buffer 
            buffer_size =  voxel_amount * sizeof(vect);
            glBufferData(GL_ARRAY_BUFFER, buffer_size, &normals[0], GL_STATIC_DRAW);

            //This tells OpenGL the data layout
            glVertexAttribPointer(2, 4, GL_FLOAT, GL_FALSE, sizeof(vect), 0); 
            
            //Enable the VBO within the VAO 
            glEnableVertexAttribArray(2);

            //Important for instanced rendering 
            glVertexAttribDivisor(2, 1);

            //Unbind the VAO
            glBindVertexArray(0);
        };
    
        //Updates the objects position and angle
        void update_position_attitude(vect delta_angle, vect delta_position)
        {     
            //Attitude Update
            quat quaternion_x = quaternion_structure(right, delta_angle.x);
            quat quaternion_y = quaternion_structure(up, delta_angle.y);
            quat quaternion_z = quaternion_structure(foward, delta_angle.z);

            //Multiplying change in quaternion by universal quaternion then rotating point
            rotation_quaternion = quaternion_multiply(quaternion_multiply(quaternion_x, quaternion_y), quaternion_z);
            quaternion = quaternion_multiply(rotation_quaternion, quaternion);

            //Updating Objects Local Axis from rotation
            up     = vector_normalize(quaternion_rotation(quaternion, world_up));
            right  = vector_normalize(quaternion_rotation(quaternion, world_right)); 
            foward = vector_normalize(quaternion_rotation(quaternion, world_foward));

            //Updating Position 
            position = vector_add(position, delta_position);

            //Creating the model matrix
            matrix_model(model, position, foward, up, right); 
            //matrix_rotation(rotation, foward, up, right);
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
                        
                        voxel_object.push_back(temporary); 
                    }

                    counter ++; 
                }; 

                voxel_amount = counter;
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

                //Creates the box shape for storing the voxel_object 
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
                                    voxel_object.push_back(new_voxel); 

                                    break; 
                                };
                            }; 
                        }; 
                    }; 
                };

                //Record total amount of voxels
                voxel_amount = voxel_object.size(); 
            }; 
        };

        //Render function
        void render(GLuint ModelMatrixID, GLuint ViewMatrixID, float look_at[4][4])
        {   
            //Bind the Array Object
            glBindVertexArray(VAO[0]);

            //Update the Model and View Matrices
            glUniformMatrix4fv(ModelMatrixID, 1, NULL, &model[0][0]);
            glUniformMatrix4fv(ViewMatrixID, 1, NULL, &look_at[0][0]);

            //Draw The Voxels
            glDrawArraysInstanced(GL_QUADS, 0, 4, voxel_amount); 

            //Unbind the Array Object
            glBindVertexArray(0);
        }; 

        //Clean up function
        void clean_up()
        {
            glDeleteBuffers(3, VBO);
            glDeleteVertexArrays(1, VAO);

            voxel_object.clear();
        }
}; 

////////////////
//Shader Class//
////////////////
class shader
{  
    public:
        unsigned int shader_program;
        GLuint ModelMatrixID;
        GLuint ViewMatrixID;
        GLuint ProjectionMatrixID;
        GLuint ScalingFactorID;

        std::string read_file(std::string file_name)
        {
            //Open the text file in read mode and transfer the data 
            std::ifstream file(file_name);   

            //Return data  (and convert std::ifstream to std::string)
            return {std::istreambuf_iterator<char>(file), std::istreambuf_iterator<char>()};
        };

        unsigned int compile_shader(std::string file_name, unsigned int shader_type) 
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

                printf(message);

                return 0; 
            };

            //Return the compiled shader
            return shader_id; 
        };

        unsigned int create_shaders(std::string vertex_shader, std::string fragment_shader)
        {   
            //Create the Overall Program
            shader_program = glCreateProgram(); 

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

        void intialize_shader(float projection[4][4])
        {
            //Setup the ID for the transformation matrices in the shader
            ModelMatrixID = glGetUniformLocation(shader_program, "ModelMatrix");
            ViewMatrixID = glGetUniformLocation(shader_program, "ViewMatrix");
            ProjectionMatrixID = glGetUniformLocation(shader_program, "ProjectionMatrix"); //Not really necessery unless camaera parameters are to be changed mid-game

            //Assign the projection matrix as it is a constant at runtime
            glUniformMatrix4fv(ProjectionMatrixID, 1, GL_FALSE, &projection[0][0]);

        };
};

//The main function 
int main(int argc, char *argv[]) 
{  
    printf("\nStarting Game... \n"); 

    ///////////////
    //Asset Setup//
    ///////////////

    //This creates a camera for the world 
    camera camera;
    camera.position = {0.0, 0.0, 0}; 
    camera.intialize_projection_matrix(); 

    //This creates and Object asigns a mesh to it and then voxelizes that mesh
    printf("Creating Game Assets...  \n"); 

    int object_number = 4;
    object game_objects[(int) pow(object_number, 3)];

    //Create Objects
    int iter = 0;
    for(int i = 0; i < object_number; i++)
    {
        for(int j = 0; j < object_number; j++)
        {
            for(int k = 0; k < object_number; k++)
            {
                game_objects[iter].polygon_mesh.load_mesh("meshes/bunny.obj");
                game_objects[iter].voxelization("bunny_", 1); 
                game_objects[iter].quaternion = {0, 0.7068252, 0, 0.7073883}; 
                game_objects[iter].position =  {(float) i * 2, (float) j * 2, (float) k * 2};

                iter++;
            };
        };
    };

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
    glEnable(GL_DEPTH_TEST);;
     
    glClearColor(0.0, 0.0, 0.0, 1.0); 
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT); 
    SDL_GL_SwapWindow(window); 

    //////////////////////
    //Setting Up Shaders//
    //////////////////////
    printf("Compiling Shader Programs...  \n"); 
    
    shader shader;
    shader.create_shaders("shaders/vertex.glsl", "shaders/fragment.glsl");
    glUseProgram(shader.shader_program);
    shader.intialize_shader(camera.projection);

    //////////////////////
    //Intializing Assets//
    //////////////////////
    printf("Initialising Assets...  \n"); 

    for(int i = 0; i < pow(object_number, 3); i++)
    {
        game_objects[i].initialise();
    };

    /////////////
    //Main Loop//
    /////////////

    //Initializing Game Time 
    float game_time = SDL_GetTicks(); 
    float time_elapsed = 0; 

    //Defining Exit condition
    bool run = true;  
    int counter = 0;

    printf("Entering Game Loop...  \n"); 

    while (run == true)
    {
        //SDL Loop Conditions (Key Presses)
        while (SDL_PollEvent(&window_event)){
            switch( window_event.type ){
                case SDL_QUIT:
                    printf("\nQuiting Game...  \n"); 
                    run = false;
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

        ///////////////
        //Update Call//
        ///////////////

        //Updating Time (in seconds)
        time_elapsed = (SDL_GetTicks() - game_time)/1000;
        game_time = SDL_GetTicks();

        //Updating Positions and Angles
        delta_angle.x = -(pitch_up_speed - pitch_down_speed)*time_elapsed; 
        delta_angle.y = (yaw_right_speed - yaw_left_speed)*time_elapsed;
        delta_angle.z = (roll_right_speed - roll_left_speed)*time_elapsed;

        delta_position.x = -(left_speed - right_speed)*time_elapsed;
        delta_position.y = -(down_speed - up_speed)*time_elapsed; 
        delta_position.z = -(backward_speed - foward_speed)*time_elapsed;
    
        //Updating the Camera Object 
        camera.update(delta_angle, delta_position);

        //Update attitude
        for(int i = 0; i <  pow(object_number, 3); i++)
        {
            game_objects[i].update_position_attitude({0.00, 0.00, 0.00}, {0.0, 0.0, 0});
        };

        ///////////////
        //Render Call//
        ///////////////

        //OpenGL Clear Render
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        //Render Object
        for(int i = 0; i <  pow(object_number, 3); i++)
        {
            game_objects[i].render(shader.ModelMatrixID, shader.ViewMatrixID, camera.look_at);
        }

        //Swap Rendering Buffers
        SDL_GL_SwapWindow(window);
    };
    
    ///////////////////
    //Program Cleanup//
    ///////////////////

    printf("Deleting Assets... \n"); 

    for(int i = 0; i < pow(object_number, 3); i++)
    {
        game_objects[i].clean_up();
    }

    //Cleanup OpenGL Assets
    glDeleteProgram(shader.shader_program); 

    //Cleanup SDL Assets
    SDL_GL_DeleteContext(context); 
    SDL_DestroyWindow(window);
    SDL_Quit();

    printf("Closing... \n"); 

    return 0;
};