#version 330 core
out vec4 FragColor;

in VS_OUT
{
    vec2 positionSS;
} fs_in;

uniform sampler2D GBuffer;

vec3 ACESFilm(vec3 x)
{
    float a = 2.51f;
    float b = 0.03f;
    float c = 2.43f;
    float d = 0.59f;
    float e = 0.14f;
    return clamp((x*(a*x+b))/(x*(c*x+d)+e), 0.0, 1.0);
}

void main(void)
{
    vec4 color = texture(GBuffer, fs_in.positionSS);
    // aces
    color.xyz = ACESFilm(color.xyz);
    FragColor = color;
}
