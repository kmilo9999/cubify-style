#version 330 core

layout(location = 0) in vec3 position; // Position of the vertex
layout(location = 1) in vec3 normal;   // Normal of the vertex

uniform mat4 m;
uniform mat4 vp;

out VS_OUT
{
    vec3 positionWS;
    vec3 normalWS;
} vs_out;


void main(void)
{
    vec4 positionWS = m * vec4(position, 1.0);
    vs_out.positionWS = positionWS.xyz / positionWS.w;
    vec4 normalWS = vec4(normalize(mat3(transpose(inverse(m))) * normal), 0);
    vs_out.normalWS = normalWS.xyz;
    gl_Position = vp * (positionWS + normalWS * 0.03);
}
