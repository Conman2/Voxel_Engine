#version 440 core

//Matrices
uniform mat4 ModelMatrix;
uniform mat4 ViewMatrix;
uniform mat4 ProjectionMatrix;

//Inputs
layout (location = 0) in vec4 aPosition;
layout (location = 1) in vec4 aNormal;
layout (location = 1) in vec4 aColour; 

//Outputs
out vec4 fColor;
out vec4 fNormal;

void main()
{   
    fNormal = ModelMatrix * aNormal; 

    gl_Position = ModelMatrix * ViewMatrix * ProjectionMatrix * aPosition; 

    fColor = aColour;
}

