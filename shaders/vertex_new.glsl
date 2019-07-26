#version 440 core

//Inputs
layout (location = 0) in vec3 aPosition;
layout (location = 0) in vec4 aNormal;
layout (location = 1) in vec4 aObjectPosition; 
layout (location = 2) in vec4 aLightDirection;

uniform mat4 aProjection; 

//Outputs
out vec4 aColor;
out vec4 bPosition
out float aSize

//Functions 
vec4 quaternion_rotation( vec4 q, vec3 v )
{ 
	return (v + 2.0*cross(cross(v, q.xyz ) + q.w*v, q.xyz), 1.0);
};

void main()
{
    //Moving into World Space 
    vect4 normal_direction = quaternion_rotation(quaternion, aNormal);
    vect4 voxel_position = quaternion_rotation(quaternion, aPosition) + aObjectPosition;

    //Removing unessecery voxels 
    if (dot(normalize(camera.position - voxel_position), normal_direction) > -0.35f)  
    {
        //Moving Into View Space 
        vect4 camera_view = quaternion_rotation(camera.quaternion, voxel_position - camera.position);

        //This removes Voxels behind the camera (Mirror Realm Rabbit)
        if (camera_view.k < 0.0f)
        {
            //Lighting
            float dp = min(1.0f, max(0.2f, dot(light_direction, normal_direction)));
            
            //Projecting
            vect4 result = camera_view*aProjection;
            result = result / result.w; 

            //Recording Position
            gl_Position = result; 

            //Recording Size
            aSize = 1/length(camera.position - voxel_position) * screen_width * voxel_size; 

            //only create a object to calculate if in Camera View Space
            if (result.x + aSize > 0 && result.x < screen_width && result.y + aSize > 0 && result.y < screen_height)
            {
                //Looking Up the Colour from the Structure
                aColour = vec4(dp, dp, dp, 1); 
            };
        };
    }; 
}
