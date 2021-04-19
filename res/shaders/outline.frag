#version 330 core
out vec4 FragColor;

in VS_OUT
{
    vec3 positionWS;
    vec3 normalWS;
} fs_in;

void main(void)
{
    FragColor = vec4(0.0, 0.0, 0.0, 1.0);
}
