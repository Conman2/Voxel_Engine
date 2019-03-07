#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <SDL2/SDL.h>

//
// Initial Set Up
//

// Screen Properties
int screen_width = 1000;
int screen_height = 1000;

//result[0]iable Declerations
float camera_position[4]  = {0, 0, 0, 1};
float camera_direction[4] = {0, 0, 1, 0};
float camera_angle[4]     = {0, 0, 0, 0};
float world_up[4]         = {0, 1, 0, 0};

float camera_view_angle = 90; 
float z_max_distance = 1000;
float z_min_distance = 0.1;

//Run Condition
int run = 1; 

//Point Definition (resultory)
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
float point[4];
float result[4];

float transformed[4];
float translated[4];
float camera[4];
float camera_up[4];

float look_at[4][4];
float rotation_x[4][4];
float rotation_y[4][4];
float rotation_z[4][4];
float translation[4][4];
float rotation_matrix[4][4];

float camera_target[4]; 
float rotating_x[4];
float rotating_y[4];
float rotating_z[4];
float translating[4];

//
//Function Definition 
//

//Matrix Vector Multiplication (note that with this function input and output must be different variables)
int matrix_vector_multiplication(float result[4], float vector[4], float matrix[4][4])
{   
    for(int i = 0; i < 5; i++)
    {
        result[i] = vector[0]* matrix[i][0] + vector[1] * matrix[i][1] + vector[2] * matrix[i][2] + vector[3]*matrix[i][3];
    }

    return 0;
}

//Clears a 4x4 matrix and sets the content to 0
int clear_matrix(float matrix[4][4])
{
    for(int i = 0; i++; i < 5)
    {
        for(int j = 0; j++; j < 5)
        {
            matrix[i][j] = 0.0;
        }
    }
}

//Rotation Matrix X
int matrix_rotation_x(float matrix[4][4], float angle)
{   
    matrix[0][0] = 1.0;
    matrix[1][1] = cos(angle);
    matrix[1][2] = -sin(angle);
    matrix[2][1] = sin(angle);
    matrix[2][2] = cos(angle);
    matrix[3][3] = 1.0;

    return 0;
}

//Rotation Matrix Y
int matrix_rotation_y(float matrix[4][4], float angle)
{   
    matrix[0][0] = cos(angle);
    matrix[0][2] = -sin(angle);
    matrix[2][0] = sin(angle);
    matrix[1][1] = 1.0;
    matrix[2][2] = cos(angle);
    matrix[3][3] = 1.0;

    return 0;
}

//Rotation Matrix Z
int matrix_rotation_z(float matrix[4][4], float angle)
{   
    matrix[0][0] = cos(angle);
    matrix[0][1] = -sin(angle);
    matrix[1][0] = sin(angle);
    matrix[1][1] = cos(angle);
    matrix[2][2] = 1.0;
    matrix[3][3] = 1.0;

    return 0;
}

// Translating
int matrix_translation(float matrix[4][4], float position[4])
{  
    matrix[0][0] = 1.0;
    matrix[1][1] = 1.0;
    matrix[2][2] = 1.0;
    matrix[3][3] = 1.0;
    matrix[0][3] = position[0];
    matrix[1][3] = position[1];
    matrix[2][3] = position[2];

    return 0;
}

//Vector Adding 
int vector_add(float result[4], float vector1[4], float vector2[4])
{
    for(int i = 0; i < 4; i++)
    {
        result[i] = vector1[i] + vector2[i];
    }

    return 0;
}

//Vector Subtracting 
int vector_subtract(float result[4], float vector1[4], float vector2[4])
{
    for(int i = 0; i < 4; i++)
    {
        result[i] = vector1[i] - vector2[i];
    }

    return 0;
}

//Vector Multiplication with Constant
int vector_multiplication(float result[4], float vector[4], float constant)
{
    for(int i = 0; i < 4; i++)
    {
        result[i] = vector[i]*constant;
    }

    return 0;
}

int vector_division(float result[4], float vector[4], float constant)
{
    for(int i = 0; i < 4; i++)
    {
        result[i] = vector[i]/constant;
    }

    return 0;
}

//Dot Product 
float dot_product(float vector1[4], float vector2[4])
{
    float result = 0.0;

    for(int i = 0; i < 4; i++)
    {
        result += vector1[i]*vector2[i]; 
    }

    return result;
}

//Cross Product
int cross_product(float result[4], float vector1[4], float vector2[4])
{   

    result[0] = vector1[1]*vector2[2] - vector1[2]*vector2[1];
    result[1] = vector1[2]*vector2[0] - vector1[0]*vector2[2];
    result[2] = vector1[0]*vector2[1] - vector1[1]*vector2[0];

    return 0;
}

//Normalization
int normalize(float result[4], float vector[4])
{   
    float length = pow(pow(result[0],2) + pow(result[1],2) + pow(result[2],2) , 0.5);

    vector_division(result, vector, length);

    return 0;
}

//Point at Matrix
int matrix_lookat(float matrix[4][4], float camera_position[4], float camera_target[4], float world_up[4])
{
    //Calculating the Camera Looking Direction
    float camera_direction[4];
    vector_subtract(camera_direction, camera_position, camera_target); 
    normalize(camera_direction, camera_direction);

    //Calculating new up
    float camera_up[4];
    vector_multiplication(camera_up, camera_direction, dot_product(world_up, camera_direction));
    vector_subtract(camera_up, world_up, camera_up);
    normalize(camera_up, camera_up);

    //Calculating Camera Right Direction
    float camera_right[4];
    cross_product(camera_right, camera_up, camera_direction);

    //Look At Matrix 
    matrix[0][0] = camera_right[0];     matrix[0][1] = camera_right[1];     matrix[0][2] = camera_right[2];     matrix[0][3] = 0;
    matrix[1][0] = camera_up[0];        matrix[1][1] = camera_up[1];        matrix[1][2] = camera_up[2];        matrix[1][3] = 0;
    matrix[2][0] = camera_direction[0]; matrix[2][1] = camera_direction[1]; matrix[2][2] = camera_direction[2]; matrix[2][3] = 0;
    matrix[3][0] = 0;                   matrix[3][1] = 0;                   matrix[3][2] = 0;                   matrix[3][3] = 1; 

    return 0;
}

//Rotating 
float quaternion(float result[4], float axis[4], float angle)
{
    result[0] = axis[0] * sin(angle / 2);
    result[1] = axis[1] * sin(angle / 2);
    result[2] = axis[2] * sin(angle / 2);
    result[3] = cos(angle / 2);

    normalize(result, result);
};

void quaternion_multiplication(float result[4], float quaternion1[4], float quaternion2[4])
{   
    result[0] = (quaternion1[3]*quaternion2[0] + quaternion1[0]*quaternion2[3] + quaternion1[1]*quaternion2[3] - quaternion1[2]*quaternion2[1]);
    result[1] = (quaternion1[3]*quaternion2[1] - quaternion1[0]*quaternion2[3] + quaternion1[1]*quaternion2[3] + quaternion1[2]*quaternion2[0]);
    result[2] = (quaternion1[3]*quaternion2[3] + quaternion1[0]*quaternion2[1] - quaternion1[1]*quaternion2[0] + quaternion1[2]*quaternion2[3]);
    result[3] = (quaternion1[3]*quaternion2[3] - quaternion1[0]*quaternion2[0] - quaternion1[1]*quaternion2[1] - quaternion1[2]*quaternion2[3]);
};

void rotation(float rotation_matrix[4][4], float angle[4], float position[4])
{
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

    normalize(result, result);

    rotation_matrix[0][0] = 1 - 2*pow(result[1], 2)- 2*pow(result[2], 2);
    rotation_matrix[0][1] = 2*result[0]*result[1] - 2*result[3]*result[2];
    rotation_matrix[0][2] = 2*result[0]*result[2] + 2*result[3]*result[1];
    rotation_matrix[1][0] = 2*result[0]*result[1] + 2*result[3]*result[2];
    rotation_matrix[1][1] = 1 - 2*pow(result[0], 2)- 2*pow(result[2], 2);
    rotation_matrix[1][2] = 2*result[1]*result[2] + 2*result[3]*result[0];
    rotation_matrix[2][0] = 2*result[0]*result[2] - 2*result[3]*result[1];
    rotation_matrix[2][1] = 2*result[1]*result[2] - 2*result[3]*result[0];
    rotation_matrix[2][2] = 1 - 2*pow(result[0], 2)- 2*pow(result[1], 2);
    rotation_matrix[3][3] = 1;

    /*
    float x_axis[4] = {1, 0, 0, 1};
    float y_axis[4] = {0, 1, 0, 1};
    float z_axis[4] = {0, 0, 1, 1};

    float quaternion_x_axis[4];
    float quaternion_y_axis[4];
    float quaternion_z_axis[4];  

    quaternion(quaternion_x_axis, position, angle[0]); 
    quaternion(quaternion_y_axis, position, angle[1]); 
    quaternion(quaternion_z_axis, position, angle[2]); 

    float quaternion_xy[4];
    float quat[4];

    quaternion_multiplication(quaternion_xy, quaternion_z_axis, quaternion_y_axis);
    quaternion_multiplication(quat, quaternion_xy,     quaternion_x_axis);
    
    rotation_matrix[0][0] = 1 - 2*pow(quat[1], 2)- 2*pow(quat[2], 2);
    rotation_matrix[0][1] = 2*quat[0]*quat[1] - 2*quat[3]*quat[2];
    rotation_matrix[0][2] = 2*quat[0]*quat[2] + 2*quat[3]*quat[1];
    rotation_matrix[1][0] = 2*quat[0]*quat[1] + 2*quat[3]*quat[2];
    rotation_matrix[1][1] = 1 - 2*pow(quat[0], 2)- 2*pow(quat[2], 2);
    rotation_matrix[1][2] = 2*quat[1]*quat[2] + 2*quat[3]*quat[0];
    rotation_matrix[2][0] = 2*quat[0]*quat[2] - 2*quat[3]*quat[1];
    rotation_matrix[2][1] = 2*quat[1]*quat[2] - 2*quat[3]*quat[0];
    rotation_matrix[2][2] = 1 - 2*pow(quat[0], 2)- 2*pow(quat[1], 2);
    rotation_matrix[3][3] = 1;
    */
};

//
//Main Function
//

int main(int argc, char *argv[]) 
{  
    //Starting Setup
    printf("Starting Setup\n");
    
    //Setting up SDL
    SDL_Init(SDL_INIT_VIDEO);
    SDL_Window*window  = SDL_CreateWindow("Hopefully this is 3d", SDL_WINDOWPOS_UNDEFINED, SDL_WINDOWPOS_UNDEFINED, screen_width, screen_height, SDL_WINDOW_RESIZABLE | SDL_WINDOW_OPENGL);
    SDL_Renderer*renderer = SDL_CreateRenderer(window, -1, SDL_RENDERER_ACCELERATED | SDL_RENDERER_PRESENTVSYNC);
    SDL_SetHint(SDL_HINT_RENDER_SCALE_QUALITY, "nearest");
    SDL_SetHint(SDL_HINT_RENDER_VSYNC, "1");
    SDL_Event window_event;

    //Projection Matrix Set-up
    float view_angle =    (float) camera_view_angle*(3.14159/180); //Turning Degrees to Radians 
    float aspect_ratio =  (float) screen_height/screen_width;
    float z_normalize =   (float) z_max_distance/(z_max_distance - z_min_distance);
    float z_other =       (float) (-z_max_distance*z_min_distance)/(z_max_distance - z_min_distance);
    float feild_of_view = (float) 1/tan(view_angle/2);

    float projection_matrix[4][4] = 
    {
        {aspect_ratio*feild_of_view, 0,             0,           0      },
        {0,                          feild_of_view, 0,           0      },
        {0,                          0,             z_normalize, z_other},
        {0,                          0,             1,           0      }
    };
    
    //Defining the contenet of the matrix as 0
    clear_matrix(rotation_matrix);
    clear_matrix(translation);
    clear_matrix(rotation_x);
    clear_matrix(rotation_y);
    clear_matrix(rotation_z);
    clear_matrix(look_at);

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
        camera_angle[0] -= (yaw_right_speed - yaw_left_speed)*time_elapsed;
        camera_angle[1] -= (pitch_up_speed - pitch_down_speed)*time_elapsed; 
        camera_angle[2] -= (roll_right_speed - roll_left_speed)*time_elapsed;

        camera_position[0] += (right_speed - left_speed)*time_elapsed;
        camera_position[1] += (up_speed - down_speed)*time_elapsed; 
        camera_position[2] = 5;

        //Rotation Matrices
        matrix_rotation_x(rotation_x, camera_angle[0]);
        matrix_rotation_y(rotation_y, camera_angle[1]);
        matrix_rotation_z(rotation_z, camera_angle[2]); 

        rotation(rotation_matrix, camera_angle, camera_position);

        //Changing Target Vector accoarding to camera angles
        matrix_vector_multiplication(rotating_x, camera_direction, rotation_matrix);
        //matrix_vector_multiplication(rotating_y, rotating_x,       rotation_y);
        normalize(rotating_x, rotating_x);

        //Changing Up Vector Based on camera angle
        matrix_vector_multiplication(camera_up, world_up, rotation_matrix); 
        normalize(camera_up, camera_up);

        //Translation Matrix  
        matrix_translation(translation, camera_position);

        //Camera Control Setup
        vector_add(camera_target, camera_position, rotating_x);
        matrix_lookat(look_at, camera_position, camera_target, camera_up);
        
        //Draw More Cubes 

                //Loop for each point
                for (int i = 0; i < 8; i++)
                {   
                    //Moving data
                    for(int j = 0; j < 4; j++)
                    {
                        point[j] = points[i][j];
                    }; 

                    //Translated
                    matrix_vector_multiplication(translating, point, translation);

                    //Camera Manipulation
                    matrix_vector_multiplication(camera, translating, look_at);

                    //Projectring from 3d into 2d then normalizing
                    matrix_vector_multiplication(result, camera, projection_matrix);
                    vector_division(result, result, result[3]);

                    //Scaling into screen 
                    result[0] = (result[0] + 1.0)*screen_width/2;
                    result[1] = (result[1] + 1.0)*screen_width/2; 

                    updated[i][0] = result[0];
                    updated[i][1] = result[1];
                };

                //Drawing the Sides of the Cubes 
                SDL_SetRenderDrawColor(renderer, 255, 255, 255, 255);

                //Front 
                SDL_RenderDrawLine(renderer,updated[0][0],updated[0][1], updated[1][0],updated[1][1]);
                SDL_RenderDrawLine(renderer,updated[1][0],updated[1][1], updated[2][0],updated[2][1]);
                SDL_RenderDrawLine(renderer,updated[2][0],updated[2][1], updated[3][0],updated[3][1]);
                SDL_RenderDrawLine(renderer,updated[3][0],updated[3][1], updated[0][0],updated[0][1]);

                //Back
                SDL_RenderDrawLine(renderer,updated[4][0],updated[4][1], updated[5][0],updated[5][1]);
                SDL_RenderDrawLine(renderer,updated[5][0],updated[5][1], updated[6][0],updated[6][1]);
                SDL_RenderDrawLine(renderer,updated[6][0],updated[6][1], updated[7][0],updated[7][1]);
                SDL_RenderDrawLine(renderer,updated[7][0],updated[7][1], updated[4][0],updated[4][1]);

                //Left
                SDL_RenderDrawLine(renderer,updated[1][0],updated[1][1], updated[2][0],updated[2][1]);
                SDL_RenderDrawLine(renderer,updated[2][0],updated[2][1], updated[6][0],updated[6][1]);
                SDL_RenderDrawLine(renderer,updated[6][0],updated[6][1], updated[5][0],updated[5][1]);
                SDL_RenderDrawLine(renderer,updated[5][0],updated[5][1], updated[1][0],updated[1][1]);
                
                ////Right
                SDL_RenderDrawLine(renderer,updated[0][0],updated[0][1], updated[3][0],updated[3][1]);
                SDL_RenderDrawLine(renderer,updated[3][0],updated[3][1], updated[7][0],updated[7][1]);
                SDL_RenderDrawLine(renderer,updated[7][0],updated[7][1], updated[4][0],updated[4][1]);
                SDL_RenderDrawLine(renderer,updated[4][0],updated[4][1], updated[0][0],updated[0][1]);

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
