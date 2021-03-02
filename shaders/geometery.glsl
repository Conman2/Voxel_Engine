#version 440 core
 
layout(points) in;
layout(points, max_vertices = 1) out;
 
in vec4 aPosition;
in int aRender;
 
out vec4 CulledPosition;
 
void main() 
{
    if (aRender == 1)
    {
        CulledPosition = aPositition;
        EmitVertex();
        EndPrimitive();
    };
};