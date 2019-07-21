#version 440 core

layout (location = 0) in vec3 aPosition;
layout (location = 1) in vec4 aTranslation; 
layout (location = 2) in float aScale;
layout (location = 3) in vec4 aColour; 

out vec4 fColor;

void main()
{
    gl_Position = vec4((aPosition.x*aScale + aTranslation.x), (aPosition.y*aScale + aTranslation.y), aPosition.z, 1.0);
    fColor = aColour;
}
