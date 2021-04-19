#version 420 core

layout(location = 0) in vec3 position; // Position of the vertex
layout(location = 1) in vec3 normal;   // Normal of the vertex
layout(location = 2) in vec2 uv;        // uv coordinate

uniform mat4 m;
uniform mat4 vp;

out VS_OUT
{
    vec3 positionWS;
    vec4 positionCS;        // clip space position
    vec3 normalWS;
    vec2 uv;
} vs_out;

void main(void)
{
    vec4 positionWS = m * vec4(position, 1.0);
    vs_out.positionWS = positionWS.xyz / positionWS.w;
    vec4 normalWS = vec4(normalize(mat3(transpose(inverse(m))) * normal), 0);
    vs_out.normalWS = normalWS.xyz;
    vs_out.uv = uv;
    vs_out.positionCS = vp * positionWS;
    vs_out.positionCS = vs_out.positionCS / vs_out.positionCS.w;
    gl_Position = vs_out.positionCS;
}
