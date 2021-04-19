#version 330 core

layout(location = 0) in vec3 position; // Position of the vertex
layout(location = 1) in vec3 normal;   // Normal of the vertex

out VS_OUT
{
    vec2 positionSS;        // screen space location
} vs_out;

void main(void)
{
    vs_out.positionSS = clamp((position.xy + 1.0) * 0.5, 0.0, 1.0);
    gl_Position = vec4(position, 1.0);
}
