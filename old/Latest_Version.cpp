//Standard Libaries 
#include <ctime>
#include <vector>
#include <fstream>
#include <strstream>
#include <algorithm>
#include <GL/glew.h>
#include <SDL2/SDL.h>

//Quaternion Libary
#include "quaternion.hpp"

//GLM Libaries
#include "glm/geometric.hpp"
#include "glm/vec3.hpp"
#include "glm/vec4.hpp"
#include "glm/mat4x4.hpp"
#include "glm/gtc/matrix_transform.hpp"

//
//TODO
//
//Add a Z-buffer??? (To replace the painter algorithim) (Through OpenGl) (Maybe painter through opengl))
//Look Into Vertix Buffing (Instancing actually)
//Add a refernce plane (so its not a easy to get lost)
//Look into terrain (what method of rendering? Using voxels same as sprites or using perlin-noise and line voxel thing)
//Draw less voxels the further away you are (LOD)
//Make a HUD compass
//Move font functions into header file 
//Normalize mesh size (rabbit reall small) also add scaling input to the mesh
//Fix the voxelization algorithim (make more efficent)
//Swithch to OpenGL projection and stuff 
//Make change in angle/position irrelevant to time
//Flickering voxels when pitch is 0. Maybe add if pitch==0 pitch = 0.01 or something
//Mess with aspect ratio stuff (voxel size varying with aspect ratio) (Currently limited to square to work propely) (maybe diagonal screen distance)
//Behind Font Texture is black instead of alpha 
//Holes open up when the object in corner of screen because you are essentially looking at the non-existant face of the cube (the one parrel with the view vector) To fixe this you must take into account the 3d shape of a voxel
//Setup shaders to do stuff (projection/lighting)
//Simplty camera class functions to an intialize, update, close
//Rotate around an objets center (currently around the edge of the object)
//Remove LookAt Matrix infastructure (Not being used anymore)
//Check if anypart of the Object AABB is in camera before cehcking each voxel
//Scale Function (Zoom) by multiplying the diagnonal values of the view matrix by a scale factor (values that would normally be one)
//Move all font data into a single matrix (making initializing it easier)
//Fix FPS counter
//Create a Initialize funtion and a cleanup function which handles OpenGl and SDL2 attributes ect.
//Move the shader programs into a class

//World State Variables 
typedef Uint32 colour_format; 
const Uint32 pixel_format_id = SDL_PIXELFORMAT_RGBA32;

//This stores all font textures 
SDL_Texture* character_textures[95];

//Random Constants 
const int fontsize = 3;
const float voxel_size = 0.003f; //This is world voxel size 
const float screen_width = 1000;
const float screen_height = 1000;
const int world_voxel_limit = 1000; 

//The Universal Origin and Axes 
const glm::vec4 world_origin = {0, 0, 0, 1};
const glm::vec4 world_right  = {1, 0, 0, 0};
const glm::vec4 world_up     = {0, 1, 0, 0};
const glm::vec4 world_foward = {0, 0, 1, 0};

//Colour Struct (Default is Solid White)
struct colr
{   
    unsigned char r = 0; 
    unsigned char g = 0;
    unsigned char b = 0; 
    unsigned char a = 255; 
};

//Used to store a mesh
struct mesh
{
    struct triangle
    {
        glm::vec4 p[3];
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
		std::vector<glm::vec4> verts;

		while (!f.eof())
		{
			char line[136];
			f.getline(line, 136);

			std::strstream s;
			s << line;

			char junk;

			if (line[0] == 'v')
			{
				glm::vec4 v;
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
    glm::vec4 position; 
    glm::vec4 normal; 
    float size; 
    colr colour;
};

//This stores all seen voxels 
std::vector<voxel> voxel_projected;

//Camera Class 
class camera
{
    private: 
        //Camera data 
        float view_angle = 100.0f; 
        float z_max_distance = 10.0f;
        float z_min_distance = 0.1f;

    public:
        //Object Behavior 
        glm::vec4 position = world_origin;
        quat quaternion;
        glm::vec4 euler; 

        //Objects Local Axes 
        glm::vec4 up;
        glm::vec4 right;
        glm::vec4 foward;

        //This objects Look At matrix (if it is a camera)
        glm::mat4 look_at;
        glm::mat4 projection;

        void initialize()
        {
            projection = glm::perspective(glm::radians(view_angle), screen_width/screen_height, z_min_distance, z_max_distance);  
        }

        //Updates the objects position and angle
        void update(glm::vec4 delta_angle, glm::vec4 delta_position)
        {     
            //Attitude Update
            //Updating Quaternion Angles (Axis Angle to Quaternion)
            quaternion = quaternion_setup(quaternion, delta_angle, right, up, foward) ;

            //Updating Objects Local Axis from rotation
            up     = glm::normalize(quaternion_rotation(quaternion, world_up)); 
            right  = glm::normalize(quaternion_rotation(quaternion, world_right));
            foward = glm::normalize(quaternion_rotation(quaternion, world_foward));

            //Converting Quaternion Angle into Euler (so our puny minds can comprehend whats happening)
            //euler = quaternion_add(quaternion_multiply(quaternion_to_euler(quaternion), 180/3.141593), {180, 180, 180}); 

            //Update Look At matrix 
            look_at = glm::lookAt((glm::vec3) position, (glm::vec3) (position + foward), (glm::vec3) up);

            //Position Update
            position.x += glm::dot(delta_position, right); 
            position.y += glm::dot(delta_position, up); 
            position.z += glm::dot(delta_position, foward); 
        };
}; 

//Stores Position and Attitude and object file
class object
{
    private: 
        //This stores the voxel infromation
        glm::vec4 voxel_volume;  
        voxel AABB_corners; 

    public:
        //Mesh Object 
        mesh polygon_mesh;  
        glm::vec4 half_size; 

        //This stores this objects voxels 
        std::vector<voxel> voxels; 

        //Object Behavior 
        glm::vec4 position;
        quat quaternion;
        //glm::vec4 euler; 

        glm::vec4 velocity; 
        glm::vec4 angular_velocity; 

        glm::vec4 acceleration;
        glm::vec4 angular_acceleration;

        //Objects local Axis 
        glm::vec4 up;
        glm::vec4 right;
        glm::vec4 foward;

        //Updates the objects position and angle
        void update_position_attitude(glm::vec4 delta_angle, glm::vec4 delta_position)
        {      
            //Updating Position 
            position = position + delta_position;
            
            //Updating Quaternion Angle
            quaternion = quaternion_setup(quaternion, delta_angle, right, up, foward);

            //Updating Objects Local Axis from rotation
            up     = glm::normalize(quaternion_rotation(quaternion, world_up)); 
            right  = glm::normalize(quaternion_rotation(quaternion, world_right));
            foward = glm::normalize(quaternion_rotation(quaternion, world_foward));
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
                half_size = {boxi/2, boxj/2, boxk/2, 1}; 

                //Creates the box shape for storing the voxels 
                voxel_volume = {ceil(boxi/voxel_size), ceil(boxj/voxel_size), ceil(boxk/voxel_size), 1};
                glm::vec4 voxel_half_size = {voxel_size/2, voxel_size/2, voxel_size/2, 1};

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
                            glm::vec4 left_corner_position = {position_i, position_j, k*voxel_size, 1};
                            glm::vec4 voxel_position = left_corner_position + voxel_half_size;

                            //For every Polygon
                            voxel new_voxel; new_voxel.colour.a = 0; 

                            for(auto poly: polygon_mesh.tris)
                            {   
                                //calculate the current polygon
                                glm::vec3 tri1 = {poly.p[0].x - polygon_mesh.min_x, poly.p[0].y - polygon_mesh.min_y, poly.p[0].z - polygon_mesh.min_z};
                                glm::vec3 tri2 = {poly.p[1].x - polygon_mesh.min_x, poly.p[1].y - polygon_mesh.min_y, poly.p[1].z - polygon_mesh.min_z}; 
                                glm::vec3 tri3 = {poly.p[2].x - polygon_mesh.min_x, poly.p[2].y - polygon_mesh.min_y, poly.p[2].z - polygon_mesh.min_z};
    
                                //Check for intersection
                                if (voxel_mesh_intersection(voxel_position, voxel_half_size, tri1, tri2, tri3) == 1){
                                    //tell it that its solid 
                                    new_voxel.colour.a = 1;
                                    new_voxel.position = left_corner_position; 

                                    //Calculate the normal
                                    glm::vec3 vector1 = tri2 - tri1; 
                                    glm::vec3 vector2 = tri2 - tri1; 
                                    new_voxel.normal = {glm::cross(vector1, vector2), 1};
                                    glm::normalize(new_voxel.normal); 

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
        std::vector<voxel> project_voxels(camera camera, glm::mat4 projection, std::vector<voxel> voxel_projected)
        {
            //Normalized Light Direction 

            glm::vec4 light_direction = {0.0f, 1.0f, 1.0f, 1};
            glm::normalize(light_direction);

            //Loop through the voxels and render the ones you want
            for (auto voxs: voxels)
            {   
                //
                //Moving into World Space 
                //
                glm::vec4 normal_direction = quaternion_rotation(quaternion, voxs.normal); //Rotating Normal
                glm::vec4 voxel_position = quaternion_rotation(quaternion, voxs.position) + position; //Rotating/Translation the Voxel Position

                //Removing unessecery voxels 
                if (glm::dot(glm::normalize(camera.position - voxel_position), normal_direction) > -0.35f)  
                {
                    //
                    //Moving Into View Space 
                    //
                    //glm::vec4 camera_view = quaternion_rotation(camera.quaternion, voxel_position - camera.position);
                    glm::vec4 camera_view = camera.look_at*voxel_position;

                    //This removes Voxels behind the camera (Mirror Realm Rabbit)
                    if (camera_view.z < 0.0f)
                    {
                        //Basic Lighting 
                        float dp = std::min(1.0f, std::max(0.2f, glm::dot(light_direction, normal_direction)));
                        
                        //
                        //Moving Into Projection Space 
                        // 
                        glm::vec4 result = camera_view * projection;
                        result = result/result.w;

                        //Storing the Projected Positions and other Voxel Characteristics 
                        voxel temp; 
                        // temp.position.x = (result.x + 1.0)*screen_width/2; 
                        // temp.position.y = (result.y + 1.0)*screen_height/2; 

                        //This is just to test
                        temp.position.x = result.x; 
                        temp.position.y = result.y; 
                        temp.position.z = result.z; 

                        //Calculating Voxel Size (I think I made this up: 1/(distance from camera)*screen_width*voxelsize) (Needs improvement)
                        temp.size = 1/glm::length(camera.position - voxel_position)*screen_width*voxel_size; 

                        //only create a object to calculate if in Camera View Space
                        if (temp.position.x + temp.size > 0 && temp.position.x < screen_width && temp.position.y + temp.size > 0 && temp.position.y < screen_height)
                        {
                            //Looking Up the Colour from the Structure
                            temp.colour.r = 255*dp; 
                            temp.colour.g = 255*dp;
                            temp.colour.b = 255*dp;   

                            //Stores Each Voxel in Projection space
                            voxel_projected.push_back(temp);
                        };
                    };
                }; 
            }; 
            
            //Sort using painter algortithim from close to far
            std::sort(voxel_projected.begin(), voxel_projected.end(), [](voxel vox1, voxel vox2)
            {
                return vox1.position.z < vox2.position.z;
            });

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

//The main function 
int main(int argc, char *argv[]) 
{  
    //
    //Setting up SDL
    //
    SDL_Init(SDL_INIT_VIDEO);
    SDL_Window* window  = SDL_CreateWindow("Testing", SDL_WINDOWPOS_UNDEFINED, SDL_WINDOWPOS_UNDEFINED, screen_width, screen_height, SDL_WINDOW_OPENGL);
    SDL_Renderer* renderer = SDL_CreateRenderer(window, -1, SDL_RENDERER_ACCELERATED | SDL_RENDERER_PRESENTVSYNC);
    SDL_Event window_event;
    SDL_SetHint(SDL_HINT_RENDER_SCALE_QUALITY, "nearest");

    //Defining Pixel Formating 
    SDL_PixelFormat *pixel_format = SDL_AllocFormat(pixel_format_id);

    //Defining an SDL Rectangle 
    SDL_Rect rectangle;

    //
    //Setting up OpenGL
    //
	SDL_GLContext render_context = SDL_GL_CreateContext(window); // Links SDL window to render context
    SDL_GL_SetAttribute(SDL_GL_CONTEXT_PROFILE_MASK, SDL_GL_CONTEXT_PROFILE_CORE);
    SDL_GL_SetAttribute(SDL_GL_CONTEXT_MAJOR_VERSION, 3);
    SDL_GL_SetAttribute(SDL_GL_CONTEXT_MINOR_VERSION, 2);	
    SDL_GL_SetAttribute(SDL_GL_DOUBLEBUFFER, 1);
    //glEnable(GL_DEPTH_TEST);

    //
    //Initializing GLEW
    //
	glewExperimental = GL_TRUE;
	glewInit();

    //Clear back buffer
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    //Set Clear Colour
    glClearColor(0.0, 0.0, 0.0, 1.0);

    //
    //Setting Up VBO/VBA
    //

    // Our object has 4 points, each point has 3 values
    const uint32_t points = 4;
    const uint32_t floatsPerPoint = 3;
    
    // This is the object we'll draw ( a simple square
    const GLfloat square[points][floatsPerPoint] = 
    {
        { -0.5,  0.5,  0.5 }, // Top left
        {  0.5,  0.5,  0.5 }, // Top right
        {  0.5, -0.5,  0.5 }, // Bottom right 
        { -0.5, -0.5,  0.5 }, // Bottom left
    };
    
    // Create variables for storing the ID of our VAO and VBO
    unsigned int VBO[4], VAO[1]; 
    
    // Allocate and assign two Vertex Buffer Objects to our handle
    glGenBuffers(4, VBO);
    
    // Allocate and assign a Vertex Array Object to our handle
    glGenVertexArrays(1, VAO);
    
    // Bind our first VBO as being the active buffer and storing vertex attributes (coordinates)
    glBindBuffer(GL_ARRAY_BUFFER, VBO[0]);
    
    // Copy the vertex data from dice to our buffer
    uint32_t buffer_size =  ( points * floatsPerPoint) * sizeof(GLfloat);
    glBufferData(GL_ARRAY_BUFFER, buffer_size, square, GL_STATIC_DRAW);
    
    // Bind our Vertex Array Object as the current used object
    glBindVertexArray(VAO[0]);
    
    // Specify that our coordinate data is going into attribute index 0, and contains three floats per vertex
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(float) * 3, 0);
    
    // Enable our attribute within the current VAO
    glEnableVertexAttribArray(0);  

    //
    // VBO for translation
    //
    //This is to store the buffer data
    glm::vec2 translations[world_voxel_limit]; 

    //Bind the instance VAO 
    glBindVertexArray(VAO[0]); 

    //Selecting the buffer we want to work on
    glBindBuffer(GL_ARRAY_BUFFER, VBO[1]);

    //Transfer data into the VBO
    //Can try out different data types (Gl_STREAM_DRAW..) to see if it changes performance
    glBufferData(GL_ARRAY_BUFFER, sizeof(glm::vec2)*world_voxel_limit, &translations[0], GL_DYNAMIC_DRAW);


    //This specifies the layout of the data
    glVertexAttribPointer(2, 2, GL_FLOAT, GL_FALSE, sizeof(glm::vec2), 0); 

    //This enables the attribute (links it with shader)
    glEnableVertexAttribArray(2);

    //glVertexAttribDivisor(2, 1); 
  
    //
    // VOB for scale 
    //
    //This is to store the buffer data
    float scaling[world_voxel_limit]; 

    //Bind the instance VAO 
    glBindVertexArray(VAO[0]); 

    //Selecting the buffer we want to work on
    glBindBuffer(GL_ARRAY_BUFFER, VBO[2]);

    //Transfer data into the VBO
    //Can try out different data types (Gl_STREAM_DRAW..) to see if it changes performance
    glBufferData(GL_ARRAY_BUFFER, sizeof(float)*world_voxel_limit, translations, GL_DYNAMIC_DRAW);

    //This specifies the layout of the data
    glVertexAttribPointer(3, 1, GL_FLOAT, GL_FALSE, sizeof(float), 0); 

    //This enables the attribute (links it with shader)
    glEnableVertexAttribArray(3);

    //
    // VBO for colour 
    //


    //
    //Shader Setup 
    //
    std::string vertex_shader_location = "vertex.glsl";
    std::string fragment_shader_location = "fragment.glsl";
    static unsigned int shader_program = create_shaders(vertex_shader_location, fragment_shader_location);
    glUseProgram(shader_program);
    
    //Creating the Camera object
    camera camera;
    camera.initialize(); 

    //
    //Creating Game Assets 
    //

    //This creates and Object asigns a mesh to it and then voxelizes that mesh
    object test1; 
    test1.polygon_mesh.load_mesh("meshes/bunny.obj");
    test1.voxelization("rabbit_"); 

    //Initializing Speeds 
    float fps[5]; int fps_counter = 0;  float avg_fps; float fps_timer;
    float left_speed = 0;     float yaw_left_speed = 0;
    float right_speed = 0;    float yaw_right_speed = 0;
    float down_speed = 0;     float pitch_down_speed = 0;
    float up_speed = 0;       float pitch_up_speed = 0;
    float foward_speed = 0;   float roll_right_speed = 0;
    float backward_speed = 0; float roll_left_speed = 0;

    //Setting Angular and Translational Velocities
    float velocity = 1;       float angular_velocity =  0.6;

    //Stores Change in position/angle
    glm::vec4 delta_position;       glm::vec4 delta_angle;

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

        //Render and update attitude
        // test1.update_position_attitude({0.005,0.005,0.005},{0.0,0,0});
        voxel_projected = test1.project_voxels(camera, camera.projection, voxel_projected);

        //Render All of the voxels 
        int counter = 0;
        for(auto voxels : voxel_projected)
        {   
            if (counter < world_voxel_limit)
            {
                translations[counter] = {voxels.position.x, voxels.position.y};
                scaling[counter++] = voxels.size; 
            }
            else 
            {
                break; 
            }
        }

        //Clear the stored voxels (otherwise would get infinitely big)
        voxel_projected.clear(); 

        //Update the buffer storage (Maybe try glMapBuffer instead)
        glBindBuffer(GL_ARRAY_BUFFER, VBO[1]);
        glBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(glm::vec2)*world_voxel_limit, &translations); 

        //Update the buffer storage (Maybe try glMapBuffer instead)
        glBindBuffer(GL_ARRAY_BUFFER, VBO[2]);
        glBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(float)*world_voxel_limit, &scaling); 


        //OpenGL Clear Render
        glClear(GL_COLOR_BUFFER_BIT);


        glBindVertexArray(VAO[0]);

        //OpenGL Rendering 
        //This should become glDrawInstanceArrays
        glDrawArraysInstanced(GL_TRIANGLE_FAN, 0, 4, world_voxel_limit); 

        glBindVertexArray(0);

        //Swap Rendering Buffers
        SDL_GL_SwapWindow(window);
    };

    //Destroy and Close SDL after program finishes
    SDL_DestroyRenderer(renderer);
    SDL_DestroyWindow(window);

    SDL_Quit();

    return 0;
};

// #version 330 core
// layout (location = 0) in vec3 aPos;
// layout (location = 2) in vec2 aTexCoords;
// layout (location = 3) in mat4 instanceMatrix;

// out vec2 TexCoords;

// uniform mat4 projection;
// uniform mat4 view;

// void main()
// {
//     gl_Position = projection * view * instanceMatrix * vec4(aPos, 1.0); 
//     TexCoords = aTexCoords;
// }