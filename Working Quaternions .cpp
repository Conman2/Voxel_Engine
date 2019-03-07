#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <SDL2/SDL.h>

//
// Initial Set Up
//

//Data Structures
struct vect 
{
    float x = 0; 
    float y = 0; 
    float z = 0;
    float w = 1;
};

// Screen Properties
int screen_width = 1000;
int screen_height = 1000;

//Variable Declerations
vect camera_position  = {0, 0, 0, 1};
vect camera_direction = {0, 0, 1, 0};
vect camera_angle     = {0, 0, 0, 0};
vect world_up         = {0, 1, 0, 0};
vect quaternion       = {1, 0, 0, 0};

float camera_view_angle = 90; 
float z_max_distance = 1000;
float z_min_distance = 0.1;

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

vect point;
vect result;
vect camera;
vect camera_up;
vect camera_target; 
vect camera_right;
vect look_direction;


//
//Vector Functions
//

//Vector Adding 
vect vector_add(vect vector1, vect vector2)
{
    return {vector1.x + vector2.x, vector1.y + vector2.y, vector1.z + vector2.z, 1};
}

//Vector Subtracting 
vect vector_subtract(vect vector1, vect vector2)
{
    return {vector1.x - vector2.x, vector1.y - vector2.y, vector1.z - vector2.z, 1};
}

//Vector Multiplication with Constant
vect vector_multiplication(vect vector, float constant)
{
    return {vector.x * constant, vector.y * constant, vector.z * constant, 1};
}

vect vector_division(vect vector, float constant)
{
    return {vector.x / constant, vector.y / constant, vector.z / constant, 1};
}

//Dot Product 
float dot_product(vect vector1, vect vector2)
{
    return {vector1.x * vector2.x + vector1.y * vector2.y + vector1.z * vector2.z}; 
}

//Cross Product
vect cross_product(vect vector1, vect vector2)
{   
    return {vector1.y * vector2.z - vector1.z * vector2.y, vector1.z * vector2.x - vector1.x * vector2.z, vector1.x * vector2.y - vector1.y * vector2.x, 1};
}

//Normalization
void vector_normalize(vect vector)
{   
    float length = pow(pow(vector.x, 2) + pow(vector.y, 2) + pow(vector.z, 2), 0.5);

    vector = vector_division(vector, length);
}

//
// Quaternion Functions
//

//Normalization
void quaternion_normalize(vect vector)
{   
    float length = pow(pow(vector.x, 2) + pow(vector.y, 2) + pow(vector.z, 2) + pow(vector.w, 2), 0.5);

    vector = vector_division(vector, length);
}

//Multiplying two quaternions using Hamilton Product  
vect quaternion_multiplication(vect quaternion1, vect quaternion2)
{   
    vect result;

    result.x = (quaternion1.w*quaternion2.x + quaternion1.x*quaternion2.w + quaternion1.y*quaternion2.z - quaternion1.z*quaternion2.y);
    result.y = (quaternion1.w*quaternion2.y - quaternion1.x*quaternion2.z + quaternion1.y*quaternion2.w + quaternion1.z*quaternion2.x);
    result.z = (quaternion1.w*quaternion2.z + quaternion1.x*quaternion2.y - quaternion1.y*quaternion2.x + quaternion1.z*quaternion2.w);
    result.w = (quaternion1.w*quaternion2.w - quaternion1.x*quaternion2.x - quaternion1.y*quaternion2.y - quaternion1.z*quaternion2.z);

    return result;
};

//The basic quaternion structure 
vect quaternion_structure(vect axis, float angle)
{   
    vect quaternion;

    quaternion.x = axis.x * sinf( angle/2 );
    quaternion.y = axis.y * sinf( angle/2 );
    quaternion.z = axis.z * sinf( angle/2 );
    quaternion.w = cosf( angle/2);

    quaternion_normalize(quaternion);

    return quaternion;
};

//Seting up the new global rotation based on input axis, change in angles and total quaternion
vect quaternion_setup(vect total_quaternion, vect angle, vect x_axis, vect y_axis, vect z_axis)
{   
    vect quaternion_x = quaternion_structure(x_axis, angle.x);
    vect quaternion_y = quaternion_structure(y_axis, angle.y);
    vect quaternion_z = quaternion_structure(z_axis, angle.z);

    //Multiplying change in quaternion by universal quaternion then rotating point
    vect quaternion = quaternion_multiplication(quaternion_multiplication(quaternion_multiplication(quaternion_z, quaternion_y), quaternion_x), total_quaternion);

    return quaternion;
};

//Provides the conjugate of a quaternion
vect quaternion_conjugate(vect quaternion)
{
    return {-quaternion.x, -quaternion.y, -quaternion.z, quaternion.w};
};

//Rotates a point or vector based on: R = P*Q*P^-1
vect quaternion_rotation(vect quaternion, vect position)
{
    vect conjugate = quaternion_conjugate(quaternion);
    vect rotated = quaternion_multiplication(quaternion_multiplication(quaternion, position), conjugate);

    return(rotated);
};

//
//Matrix Functions 
//

//Matrix Vector Multiplication (note that with this function input and output must be different variables)
vect matrix_vector_multiplication(vect vector, float matrix[4][4])
{   
    vect result;

    result.x = vector.x* matrix[0][0] + vector.y * matrix[0][1] + vector.z * matrix[0][2] + vector.w * matrix[0][3];
    result.y = vector.x* matrix[1][0] + vector.y * matrix[1][1] + vector.z * matrix[1][2] + vector.w * matrix[1][3];
    result.z = vector.x* matrix[2][0] + vector.y * matrix[2][1] + vector.z * matrix[2][2] + vector.w * matrix[2][3];
    result.w = vector.x* matrix[3][0] + vector.y * matrix[3][1] + vector.z * matrix[3][2] + vector.w * matrix[3][3];

    return result;
}

//Look at Matrix
vect matrix_lookat(float matrix[4][4], vect camera_position, vect camera_target, vect world_up)
{
    //Calculating the Camera Looking Direction
    vect camera_direction = vector_subtract(camera_position, camera_target); 
    vector_normalize(camera_direction);

    //Calculating new up
    vect camera_up = vector_subtract(vector_multiplication(camera_direction, dot_product(world_up, camera_direction)), world_up);//Might need t0 be reversed
    vector_normalize(camera_up);

    //Calculating Camera Right Direction
    vect camera_right = cross_product(camera_up, camera_direction);

    //Look At Matrix 
    matrix[0][0] = camera_right.x;     matrix[0][1] = camera_right.y;     matrix[0][2] = camera_right.z;     matrix[0][3] = -dot_product(camera_right, camera_position);
    matrix[1][0] = camera_up.x;        matrix[1][1] = camera_up.y;        matrix[1][2] = camera_up.z;        matrix[1][3] = -dot_product(camera_up, camera_position);
    matrix[2][0] = camera_direction.x; matrix[2][1] = camera_direction.y; matrix[2][2] = camera_direction.z; matrix[2][3] = -dot_product(camera_direction, camera_position);
    matrix[3][0] = 0;                  matrix[3][1] = 0;                  matrix[3][2] = 0;                  matrix[3][3] = 1; 

    return camera_right;
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
void clear_matrix(float matrix[4][4])
{
    for(int i = 0; i++; i < 5)
    {
        for(int j = 0; j++; j < 5)
        {
            matrix[i][j] = 0.0;
        }
    }
}

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
    
    //Defining the contenet of the matrix as 0
    clear_matrix(look_at);
    clear_matrix(projection);

    matrix_projection(projection, camera_view_angle,  screen_height,  screen_width,  z_max_distance,  z_min_distance);

    float left_speed = 0;
    float right_speed = 0;
    float down_speed = 0;
    float up_speed = 0;
    float yaw_left_speed = 0;
    float yaw_right_speed = 0;
    float pitch_down_speed = 0;
    float pitch_up_speed = 0;
    float roll_right_speed = 0;
    float roll_left_speed = 0;

    float angular_velocity = 2;
    float velocity = 8; 

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
        camera_angle.x = (yaw_right_speed - yaw_left_speed)*time_elapsed;
        camera_angle.y = (roll_right_speed - roll_left_speed)*time_elapsed;
        camera_angle.z = (pitch_up_speed - pitch_down_speed)*time_elapsed; 

        camera_position.x += (left_speed - right_speed)*time_elapsed;
        camera_position.y += (up_speed - down_speed)*time_elapsed; 
        camera_position.z = 5;

        //Rotation Matrix. note that it is rotating around the cameras local axis from last frame, which is already a unit direction
        quaternion = quaternion_setup(quaternion, camera_angle, camera_right, camera_up, look_direction);

        look_direction = quaternion_rotation(quaternion, camera_direction);
        vector_normalize(look_direction);

        //Changing Up Vector Based on camera angle
        camera_up = quaternion_rotation(quaternion, world_up); 
        vector_normalize(camera_up);

        //Camera Control Setup
        camera_target = vector_add(camera_position, look_direction);
        camera_right = matrix_lookat(look_at, camera_position, camera_target, camera_up);

        //Loop for each point
        for(int l = 0; l < 50; l++)
        {
            for(int k = 0; k < 50; k++)
            {
                for (int i = 0; i < 8; i++)
                {   
                    //Moving data
                    point.x = points[i][0] + l*2;
                    point.y = points[i][1];
                    point.z = points[i][2] + k*2; 

                    //Camera Manipulation
                    camera = matrix_vector_multiplication(point, look_at);

                    //Projectring from 3d into 2d then normalizing
                    result = matrix_vector_multiplication(camera, projection);
                    result = vector_division(result, result.w);

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
