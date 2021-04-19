#version 420 core
out vec4 FragColor;

in VS_OUT
{
    vec3 positionWS;
    vec4 positionCS;
    vec3 normalWS;
    vec2 uv;
} fs_in;

uniform mat4 vp;

// Don't change the layout and order.
layout(std140, binding=1) uniform Lighting
{
    vec4 lightColor[4];
    vec4 lightDirection[4];
    vec4 lightPos[4];
    int lightNumber;
};

// material uniforms
uniform vec3 diffuseColor;
uniform vec3 specularColor;
uniform vec3 ambientColor;
uniform vec3 sssColor;          // subsurface scattering, hacked via warped lighting

uniform float specularPower;

uniform float sssPower = 0;
uniform float sssRadius = 0;

// actually, color ramp is better using a texture, but here I code it.
uniform vec3 colorLevels;
uniform vec4 colorPowers;
uniform vec3 colorRadius;

// textures
uniform sampler2D diffuseTex;

// Don't worry to clamp color or gamma correction, it's HDR support and gamma correction in post processing

vec3 computeDiffuse(vec3 L, vec3 N)
{
    float LdN = dot(L, N);
    float retPower = 0;
    for (int i = 0; i < 3; ++i) {
        float a = smoothstep(colorLevels[0] - colorRadius[0], colorLevels[0] + colorRadius[0], LdN);
        retPower = mix(colorPowers[i+1], colorPowers[i], a);
    }
    return vec3(retPower, retPower, retPower);
}

vec3 computeSpecular(vec3 L, vec3 N)
{
    // use clip space position to compute direction
    vec4 dirCS = fs_in.positionCS - vec4(0,0,-1,0);
    dirCS.w = 0;            // because it's direction
    vec4 viewDir = inverse(vp) * dirCS;
    vec3 V = normalize(viewDir.xyz);

    // half vector
    vec3 H = normalize(V + L);
    float NdH = dot(N, H);
    float intensity = pow(clamp(NdH, 0, 1), specularPower);
    return vec3(intensity, intensity, intensity);
}

void main(void)
{
    // diffuse color
    vec3 texColor = texture(diffuseTex, fs_in.uv).xyz;
    vec3 ambient = ambientColor;
    vec3 diffuse = vec3(0);
    vec3 specular = vec3(0);

    // lighting
    vec3 N = normalize(fs_in.normalWS);
    for (int i = 0; i < lightNumber; ++i) {
        vec3 L = -normalize(lightDirection[i].xyz);
        vec3 LDiff = computeDiffuse(L, N) * lightColor[i].xyz;
        diffuse += LDiff;

        vec3 LSpec = computeSpecular(L, N) * lightColor[i].xyz;
        specular += LSpec;
    }
    diffuse = diffuse * diffuseColor * texColor;
    specular = specular * specularColor * texColor;
    vec3 finalColor = ambient + diffuse + specular;
    FragColor = vec4(finalColor, 1.0);
}
