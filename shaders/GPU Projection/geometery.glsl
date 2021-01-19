//This geometery shader inputs a point and size then outputs a quad to be drawn via a triangle_strip
#version 440 core

//Vertex Definitions
layout(points) in; 
layout(triangle_strip, max_vertices = 4) out; 

//Inputs 
in vec4 aColour; 
in float size; 

//Outputs
out vec4 bColour; 

void main()
{
    //Vertex 1
    gl_Position = gl_in[0]; 
    EmitVertex(); 

    //Vertex 2
    gl_Position = gl_in[0] + vec4(-size, 0.0, 0.0, 0.0); 
    EmitVertex(); 

    //Vertex 3
    gl_Position = gl_in[0] + vec4(0.0, size, 0.0, 0.0); 
    EmitVertex(); 

    //Vertex 4
    gl_Position = gl_in[0] + vec4(size, -size, 0.0, 0.0); 
    EmitVertex(); 

    //Outputs 
    bColour = aColour; 
}