#version 440 core

//Inputs
layout (location = 0) in vec3 aPosition;
layout (location = 1) in vec4 aTranslation; 
layout (location = 2) in float aScale;
layout (location = 3) in vec4 aColour; 

//Outputs
out vec4 fColor;

void main()
{
    gl_Position = vec4((aPosition.x*aScale + aTranslation.x), -(aPosition.y*aScale/aTranslation.w + aTranslation.y), aTranslation.z, 1.0);
    fColor = aColour;
}

