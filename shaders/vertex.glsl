#version 440 core

//Matrices
uniform mat4 ModelMatrix;
uniform mat4 ViewMatrix;
uniform mat4 ProjectionMatrix;

//Inputs (Note these are in static buffers)
layout (location = 0) in vec3 aMesh;
layout (location = 1) in vec4 aPosition; //Can be a vec3?
layout (location = 2) in vec4 aNormal;

//Outputs
out vec4 aColour;

void main()
{    
    //Rotate the Normal to match object rotation
    vec4 bNormal = normalize(ModelMatrix * aNormal);

    //Apply the model transformations 
    vec4 aTranslation = ModelMatrix * aPosition;

    //Apply the camera transformations
    aTranslation = ViewMatrix * aTranslation;

    //Calculate the voxel size (Important this goes between View and Projection)
    float aScale = 1 / aTranslation.z;

    //Apply the screen projection 
    aTranslation = ProjectionMatrix * aTranslation;
    aTranslation = aTranslation / aTranslation.w;

    //Submit the final position
    gl_Position = vec4((aMesh.x * aScale + aTranslation.x), (aMesh.y * aScale + aTranslation.y), aTranslation.z, 1.0);

    //Basic Lighting
    aColour = min(1.0, max(0.1, dot(vec4(0.707, 0.707, 0.0, 1.0), bNormal))) * vec4(1.0, 1.0, 1.0, 1.0);
};

