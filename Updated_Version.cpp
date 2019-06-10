#include "math.h"
#include "font.h"
#include <vector>
#include <fstream>
#include <strstream>
#include <algorithm>
#include <SDL2/SDL.h>
// #include <GL/gl.h>
// #include <GL/glu.h>

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
//Fix the voxelization algorithim 
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

//World State Variables 
typedef Uint32 colour_format; 
const Uint32 pixel_format_id = SDL_PIXELFORMAT_RGBA32;

//This stores all font textures 
SDL_Texture* character_textures[95];

//Random Constants 
const int fontsize = 3;
const float voxel_size = 0.003f; //This is world voxel size 
const int screen_width = 1000;
const int screen_height = 1000;

//The Universal Origin and Axes 
const vect world_origin = {0, 0, 0, 0};
const vect world_right  = {1, 0, 0, 1};
const vect world_up     = {0, 1, 0, 1};
const vect world_foward = {0, 0, 1, 1};

//Defining the mesh of a square
// static const GLfloat vertex_data[] =
// {
//     -0.5f, -0.5f, 0.0f,
//     0.5f,  -0.5f, 0.0f,
//     -0.5f, 0.5f,  0.0f,
//     0.5f,  0.5f,  0.0f,    
// };

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
				s >> junk >> v.i >> v.j >> v.k;
				verts.push_back(v);

                //Record the largest and smallest vertix for the Axis Aligned Bounding Box
                if      (v.i > max_x){ max_x = v.i;} 
                else if (v.i < min_x){ min_x = v.i;}
                if      (v.j > max_y){ max_y = v.j;} 
                else if (v.j < min_y){ min_y = v.j;}
                if      (v.k > max_z){ max_z = v.k;} 
                else if (v.k < min_z){ min_z = v.k;};
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

//Camera Class 
class camera
{
    private: 
        //HUD Information
        int HUDsize = fontsize; //same as fontsize looks best (consistant pixel size)
        SDL_Texture* stationaryHUD_texture; 
        SDL_Texture* movingHUD_texture; 

        //Camera data 
        float view_angle = 270.0f; 
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
        float look_at[4][4];
        float projection[4][4];

        //Updates the objects position and angle
        void update(vect delta_angle, vect delta_position)
        {     
            //Updating Quaternion Angles (Axis Angle to Quaternion)
            quaternion = quaternion_setup(quaternion, delta_angle, right, up, foward) ;

            //Updating Objects Local Axis from rotation
            up     = vector_normalize(quaternion_rotation(quaternion, world_up)); 
            right  = vector_normalize(quaternion_rotation(quaternion, world_right));
            foward = vector_normalize(quaternion_rotation(quaternion, world_foward));

            //Updating Position 
            position.i += delta_position.i*right.i + delta_position.j*up.i + delta_position.k*foward.i;
            position.j += delta_position.i*right.j + delta_position.j*up.j + delta_position.k*foward.j;
            position.k += delta_position.i*right.k + delta_position.j*up.k + delta_position.k*foward.k;

            //Converting Quaternion Angle into Euler (so our puny minds can comprehend whats happening)
            euler = vector_add(vector_multiply(quaternion_to_euler(quaternion), 180/3.142), {180, 180, 180}); 
        
            //Update Look At Matrix
            matrix_lookat(look_at, position, foward, up, right);
        };

        //Sets up projection matrix
        void intialize_projection_matrix()
        {
            matrix_projection(projection, view_angle, screen_height, screen_width, z_max_distance, z_min_distance);
        };

        //HUD System 
        void intialize_HUD(SDL_Renderer *renderer, colr colour, SDL_PixelFormat* pixel_format)
        {   
            //
            //Stationary Section 
            //

            //Defining Texture 
            bool stationaryHUD[4][216] = 
            {
                1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,
                1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,
                1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,
            };     

            //Creating texture
            stationaryHUD_texture = SDL_CreateTexture(renderer, pixel_format_id, SDL_TEXTUREACCESS_STREAMING, 216, 4); 

            // Create the texture
            colour_format texture_data[4*216];

            int counter = 0; 
            for(int i = 0; i < 4; i++)
            {   
                for (int j = 0; j < 216; j++)
                {
                    //Assign texture colour
                    if (stationaryHUD[i][j] == 1) 
                    {
                        texture_data[counter] = SDL_MapRGBA(pixel_format, colour.r,  colour.g, colour.b, colour.a); 
                    }
                    else
                    {
                        texture_data[counter] = SDL_MapRGBA(pixel_format, 0, 0, 0, 0);
                    }; 

                    //Update Total Positions 
                    counter++; 
                };
            }; 

            //Assign the texture
            SDL_UpdateTexture(stationaryHUD_texture, NULL, texture_data, sizeof(colour_format)*216);
            
            //
            //Moving Texture Section  
            //
            bool movingHUD[4][36] = 
            {
                1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
                
            }; 

            // Create the texture
            colour_format moving_texture_data[4*36*42];

            movingHUD_texture = SDL_CreateTexture(renderer, pixel_format_id, SDL_TEXTUREACCESS_STREAMING, 36*42, 4);

            //For each degree in a compass
            counter = 0;
            for(int i = 0; i < 4; i++)
            {   
                for (int j = 0; j < 36*42; j++)
                {
                    //Assign texture colour
                    if (movingHUD[i][j%36] == 1) 
                    {
                        moving_texture_data[counter] = SDL_MapRGBA(pixel_format, colour.r,  colour.g, colour.b, colour.a); 
                    }
                    else
                    {
                        moving_texture_data[counter] = SDL_MapRGBA(pixel_format, 0, 0, 0, 0);
                    };

                    //Update Total Positions 
                    counter++; 
                };
            };                
            
            SDL_UpdateTexture(movingHUD_texture, NULL, moving_texture_data, sizeof(colour_format)*36*42);
        }; 
        
        //Draw the HUD 
        void render_HUD(SDL_Renderer *renderer, int position_x, int position_y )
        {
            //Rendering the stationary texture
            SDL_Rect HUD_position = {position_x - 216*HUDsize/2, position_y, 216*HUDsize, 4*HUDsize}; 
            SDL_RenderCopy(renderer, stationaryHUD_texture, NULL, &HUD_position);

            //Preperation for the scrolling texture 
            SDL_Rect HUD_position2 = {position_x - 216*HUDsize/2, position_y - 4*HUDsize, 216*HUDsize, 4*HUDsize}; 

            //Cropping location of scrolling texture
            int index =  (euler.j/10)*36; 
            SDL_Rect cropping = {index, 0, 216, 4}; 
            
            //Cutting out the desired section of the scrolling texture
            SDL_Texture* cropped = SDL_CreateTexture(renderer, pixel_format_id, SDL_TEXTUREACCESS_TARGET, cropping.w, cropping.h);          
            SDL_SetRenderTarget(renderer, cropped);
            SDL_RenderCopy(renderer, movingHUD_texture, &cropping, NULL);
            SDL_SetRenderTarget(renderer, NULL);  

            //Rendering the moving texture
            SDL_RenderCopy(renderer, cropped, NULL, &HUD_position2);

            //Destroy the texture now it has been used
            SDL_DestroyTexture(cropped); 
        }; 
}; 

//Stores Position and Attitude and object file
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
        //vect euler; 

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
                        thing >> temporary.colour.r >> temporary.colour.g >> temporary.colour.b >> temporary.colour.a >> temporary.normal.i >> temporary.normal.j >> temporary.normal.k >> temporary.position.i >> temporary.position.j >> temporary.position.k; 
                        
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
                file << voxel_volume.i << "\t" << voxel_volume.j << "\t" << voxel_volume.k << "\n"; 
                file << "//The Voxel Information (First 4 numbers rgba, remaining 3 numbers normal vector): \n";

                //This loop goes through each voxel in the grid 
                for (int i = 0; i < voxel_volume.i; i++)
                {   
                    //Calculating i position
                    float position_i = i*voxel_size; 
                    for (int j = 0; j < voxel_volume.j; j++)
                    {
                        //Calculating j position
                        float position_j = j*voxel_size; 
                        for (int k = 0; k < voxel_volume.k; k++)
                        {   
                            //Finding voxel half size and voxel middle position
                            vect left_corner_position = {position_i, position_j, k*voxel_size};
                            vect voxel_position = vector_add(left_corner_position, voxel_half_size);

                            //For every Polygon
                            voxel new_voxel; new_voxel.colour.a = 0; 

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
                                    new_voxel.position = left_corner_position; 

                                    //Calculate the normal
                                    new_voxel.normal = vector_normalize(vector_cross_product(vector_subtract(tri2, tri1), vector_subtract(tri3, tri1))); 

                                    //Adding the voxel information to the file 
                                    file << new_voxel.colour.r << "\t" << new_voxel.colour.g << "\t" << new_voxel.colour.b << "\t" << new_voxel.colour.a << "\t" << new_voxel.normal.i << "\t" << new_voxel.normal.j << "\t" << new_voxel.normal.k << "\t" << left_corner_position.i << "\t" << left_corner_position.j << "\t" << left_corner_position.k << "\n"; 
                                    
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
            //Light Direction 
            vect light_direction = {0.0f, 1.0f, 1.0f};
            vector_normalize(light_direction);

            //Loop through the voxels and render the ones you want
            int i = 0; int j = 0; int k = 0; 
            for (auto voxs: voxels)
            {   
                //Rotating Normal
                vect normal_direction = quaternion_rotation(quaternion, voxs.normal);

                //Rotating the voxel
                vect voxel_position = vector_add(quaternion_rotation(quaternion, voxs.position), position);

                //Should remove voxels with normals facing away from camera
                //vector_dot_product(camera.foward, vector_subtract(camera.position, voxs.position)) > 0.0f <- This kills the dreaded mirror realm rabbit
                if (vector_dot_product(camera.foward, normal_direction) > -0.35f && vector_dot_product(camera.foward, vector_subtract(camera.position, voxel_position)) > 0.0f)
                {  
                    float dp = std::min(1.0f, std::max(0.2f, vector_dot_product(light_direction, normal_direction)));

                    //Camera Manipulation (LookAt Matrix Was Not Working!)
                    vect camera_view = quaternion_rotation(camera.quaternion, vector_subtract(voxel_position, camera.position));

                    //Projectring from 3d into 2d then normalizing  
                    vect result = matrix_vector_multiplication(camera_view, projection);
                    result = vector_divide(result, result.w);

                    //Storing the Projected Positions and other Voxel Characteristics 
                    voxel temp; 
                    temp.position.i = (result.i + 1.0)*screen_width/2; 
                    temp.position.j = (result.j + 1.0)*screen_height/2; 
                    temp.position.k = result.k; 

                    //Calculating Voxel Size (I think I made this up: 1/(distance from camera)*screen_width*voxelsize) (Needs improvement)
                    temp.size = 1/vector_magnitude(vector_subtract(camera.position, voxel_position))*screen_width*voxel_size; 

                    //only create a object to calculate if in Camera View Space
                    if (temp.position.i + temp.size > 0 && temp.position.i < screen_width && temp.position.j + temp.size > 0 && temp.position.j < screen_height)
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
            
            //Sort using painter algortithim from close to far
            std::sort(voxel_projected.begin(), voxel_projected.end(), [](voxel vox1, voxel vox2)
            {
                return vox1.position.k < vox2.position.k;
            });

            return voxel_projected;
        };                                                                                    
}; 

//Writes to screen 
void write_to_screen(SDL_Renderer* renderer, SDL_Texture** character_textures, std::string string, int position_x, int position_y)
{  
    //For each Character in the String
    for(char& character: string)
    {   
        //Calculating layer
        int index = character - 32; 

        SDL_Rect font_position = {position_x, position_y, 5*fontsize, 7*fontsize}; 
        SDL_RenderCopy(renderer, character_textures[index], NULL, &font_position);

        position_x += (5+1)*fontsize; 
    };  
};

//Writes the Charactes to a Surface or Texture 
void create_character_texture(SDL_Renderer *renderer, int index, const bool character[7][5], colr colour, SDL_PixelFormat* pixel_format)
{
    //Creating the texture
    character_textures[index] = SDL_CreateTexture(renderer, pixel_format_id, SDL_TEXTUREACCESS_STREAMING, 5, 7); 

    // Create the texture
    colour_format texture_data[7*5];

    int counter = 0; 
    for(int i = 0; i < 7; i++)
    {   
        for (int j = 0; j < 5; j++)
        {
            //Assign texture colour
            if (character[i][j] == 1) 
            {
                texture_data[counter] = SDL_MapRGBA(pixel_format, colour.r,  colour.g, colour.b, colour.a); 
            }
            else
            {
                texture_data[counter] = SDL_MapRGBA(pixel_format, 0, 0, 0, 0);
            }; 

            //Update Total Positions 
            counter++; 
        };
    }; 

    //Assign the texture
    SDL_UpdateTexture(character_textures[index], NULL, texture_data, sizeof(colour_format)*5);
}; 

void initialize_characters(SDL_Renderer *renderer, SDL_PixelFormat* pixel_format, colr colour)
{
    //Create all the textures 
    create_character_texture(renderer, ' ' - 32,  blk, colour, pixel_format); 
    create_character_texture(renderer, '0' - 32,  zer, colour, pixel_format); 
    create_character_texture(renderer, '1' - 32,  one, colour, pixel_format); 
    create_character_texture(renderer, '2' - 32,  two, colour, pixel_format); 
    create_character_texture(renderer, '3' - 32,  thr, colour, pixel_format); 
    create_character_texture(renderer, '4' - 32,  fur, colour, pixel_format); 
    create_character_texture(renderer, '5' - 32,  fiv, colour, pixel_format); 
    create_character_texture(renderer, '6' - 32,  six, colour, pixel_format); 
    create_character_texture(renderer, '7' - 32,  sev, colour, pixel_format); 
    create_character_texture(renderer, '8' - 32,  eig, colour, pixel_format);
    create_character_texture(renderer, '9' - 32,  nin, colour, pixel_format);
    create_character_texture(renderer, 'a' - 32,  lca, colour, pixel_format);
    create_character_texture(renderer, 'Q' - 32,  cpq, colour, pixel_format);
    create_character_texture(renderer, 'W' - 32,  cpw, colour, pixel_format);
    create_character_texture(renderer, 'X' - 32,  cpx, colour, pixel_format);
    create_character_texture(renderer, 'Y' - 32,  cpy, colour, pixel_format);
    create_character_texture(renderer, 'Z' - 32,  cpz, colour, pixel_format);
    create_character_texture(renderer, ':' - 32,  col, colour, pixel_format);
    create_character_texture(renderer, '.' - 32,  stp, colour, pixel_format);
    create_character_texture(renderer, '-' - 32,  neg, colour, pixel_format);    
};

//The main function 
int main(int argc, char *argv[]) 
{  
    //
    //Setting up SDL
    //
    SDL_Init(SDL_INIT_VIDEO);
    SDL_Window*window  = SDL_CreateWindow("Testing", SDL_WINDOWPOS_UNDEFINED, SDL_WINDOWPOS_UNDEFINED, screen_width, screen_height, SDL_WINDOW_OPENGL);
    SDL_Renderer*renderer = SDL_CreateRenderer(window, -1, SDL_RENDERER_ACCELERATED | SDL_RENDERER_PRESENTVSYNC);
    SDL_Event window_event;
    SDL_SetHint(SDL_HINT_RENDER_SCALE_QUALITY, "nearest");

    //Defining Pixel Formating 
    SDL_PixelFormat *pixel_format = SDL_AllocFormat(pixel_format_id);

    //Defining an SDL Rectangle 
    SDL_Rect rectangle;

    //
    //Setting up OpenGL
    //
    // SDL_GL_SetAttribute(SDL_GL_CONTEXT_PROFILE_MASK, SDL_GL_CONTEXT_PROFILE_CORE);
    // SDL_GL_SetAttribute(SDL_GL_CONTEXT_MINOR_VERSION, 2);	
    // SDL_GL_SetAttribute(SDL_GL_DOUBLEBUFFER, 1);
    // SDL_GLContext SDL_GL_CreateContext(window); 
    //glEnable(GL_DEPTH_TEST);

    // Clear back buffer
    //glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    
    //Creating the Camera object
    camera camera;
    camera.intialize_projection_matrix(); 

    //This creates and Object asigns a mesh to it and then voxelizes that mesh
    object test1; 
    test1.polygon_mesh.load_mesh("meshes/bunny.obj");
    test1.voxelization("rabbit_"); 

    object test2;
    test2.polygon_mesh.load_mesh("meshes/bunny.obj");
    test2.voxelization("rabbit_");
    test2.position = {0.5, 0, 0}; 

    object test3;
    test3.polygon_mesh.load_mesh("meshes/bunny.obj");
    test3.voxelization("rabbit_");
    test3.position = {0, 0, 0.5}; 

    object test4;
    test4.polygon_mesh.load_mesh("meshes/bunny.obj");
    test4.voxelization("rabbit_");
    test4.position = {0.5, 0, 0.5}; 

    //Setup font 
    initialize_characters(renderer, pixel_format, {255, 255, 255, 255}); 
    
    //Setup HUD 
    camera.intialize_HUD(renderer, {255, 255, 255, 255}, pixel_format); 

    //Initializing Speeds 
    float fps[5];  int fps_counter = 0;  float avg_fps; float fps_timer;
    float left_speed = 0;     float yaw_left_speed = 0;
    float right_speed = 0;    float yaw_right_speed = 0;
    float down_speed = 0;     float pitch_down_speed = 0;
    float up_speed = 0;       float pitch_up_speed = 0;
    float foward_speed = 0;   float roll_right_speed = 0;
    float backward_speed = 0; float roll_left_speed = 0;

    //Setting Angular and Translational Velocities
    float velocity = 1;       float angular_velocity =  0.6;

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
        delta_angle.i = (pitch_up_speed - pitch_down_speed)*time_elapsed; 
        delta_angle.j = (yaw_right_speed - yaw_left_speed)*time_elapsed;
        delta_angle.k = -(roll_right_speed - roll_left_speed)*time_elapsed;

        delta_position.i = (left_speed - right_speed)*time_elapsed;
        delta_position.j = -(down_speed - up_speed)*time_elapsed; 
        delta_position.k = (backward_speed - foward_speed)*time_elapsed;
    
        //Updating the Camera Object 
        camera.update(delta_angle, delta_position);

        //Render and update attitude
        // test1.update_position_attitude({0.005,0.005,0.005},{0.0,0,0});
        // test2.update_position_attitude({0.005,0.005,0.005},{0.0,0,0});
        // test3.update_position_attitude({0.005,0.005,0.005},{0.0,0,0});
        // test4.update_position_attitude({0.005,0.005,0.005},{0.0,0,0});
        voxel_projected = test1.project_voxels(camera, camera.projection, voxel_projected);
        voxel_projected = test2.project_voxels(camera, camera.projection, voxel_projected);
        voxel_projected = test3.project_voxels(camera, camera.projection, voxel_projected);
        voxel_projected = test4.project_voxels(camera, camera.projection, voxel_projected);

        //Render All of the voxels 
        for(auto voxels : voxel_projected)
        {   
            //Set Colour  
            SDL_SetRenderDrawColor(renderer, voxels.colour.r, voxels.colour.g, voxels.colour.b, 255);

            //If the rectangle is larger than a pixel
            if(voxels.size > 1)
            {
                //Defining the Rectangle 
                rectangle.x = voxels.position.i; 
                rectangle.y = voxels.position.j;
                rectangle.w = voxels.size;
                rectangle.h = voxels.size;

                //Drawing the Rectangle 
                SDL_RenderFillRect(renderer, &rectangle);
            } 

            //Else if the rectangle is smaller than a pixel 
            else 
            {
                SDL_RenderDrawPoint(renderer, voxels.position.i, voxels.position.j);
            }; 
        }

        //Clear the stored voxels (otherwise would get infinitely big)
        voxel_projected.clear(); 

        //Draw HUD 
        camera.render_HUD(renderer, screen_width/2, 60.0f);

        //Averaging FPS
        fps_counter ++; fps_timer += time_elapsed; 
        if(fps_timer > 0.2) 
        {
            avg_fps = fps_counter/fps_timer; 
            fps_timer = 0; fps_counter = 0; 
        }
        
        //Writing FPS to screen
        write_to_screen(renderer, character_textures, std::to_string((int) avg_fps), 5.0f, 5.0f);

        //Write Position
        std::string positions = "X:" + std::to_string(camera.position.i) + " Y:" + std::to_string(camera.position.j) + " Z:" + std::to_string(camera.position.k);
        write_to_screen(renderer, character_textures, positions, 5.0f, screen_height - 100);

        //Write Angles 
        std::string angles = "aX:" + std::to_string(camera.euler.i) + " aY:" + std::to_string(camera.euler.j) + " aZ:" + std::to_string(camera.euler.k);
        write_to_screen(renderer, character_textures, angles, 5.0f, screen_height - 50);        

        //
        //OpenGL Rendering
        //
        //glClearColor(1.0, 1.0, 1.0, 1.0);
        //glClear(GL_COLOR_BUFFER_BIT);
        //SDL_GL_SwapWindow(window);

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
};