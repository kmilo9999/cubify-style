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

vec3 computeDiffuse(float LdN)
{
    float retPower = 0;
    float a1 = smoothstep(colorLevels[0] - colorRadius[0], colorLevels[0] + colorRadius[0], LdN);
    retPower = mix(colorPowers[0], colorPowers[1], a1);
    float a2 = smoothstep(colorLevels[1] - colorRadius[1], colorLevels[1] + colorRadius[1], LdN);
    retPower = mix(retPower, colorPowers[2], a2);
    float a3 = smoothstep(colorLevels[2] - colorRadius[2], colorLevels[2] + colorRadius[2], LdN);
    retPower = mix(retPower, colorPowers[3], a3);
    // four color spaces:
    // a1 = 0 | a1=1 a2=0 | a2=1 a3=0 | a3=1
    return vec3(retPower, retPower, retPower);
}

vec3 computeSpecular(vec3 L, vec3 N)
{
    // use clip space position to compute direction
    vec4 p1 = fs_in.positionCS;
    p1.w = 1;
    vec4 p2 = p1;
    p2.z = -1;
    mat4 inv_vp = inverse(vp);
    vec4 world_p1 = inv_vp * p1;
    world_p1 /= world_p1.w;
    vec4 world_p2 = inv_vp * p2;
    world_p2 /= world_p2.w;

    vec3 V = normalize(world_p1.xyz - world_p2.xyz);

    // half vector
    vec3 H = normalize(V + L);
    float NdH = dot(N, H);
    float intensity = pow(clamp(NdH, 0, 1), specularPower);
    return vec3(intensity, intensity, intensity);
}

vec3 computeSSS(float LdN)
{
    float warp = max(0, (LdN + sssRadius) / (1 + sssRadius)) * sssPower;
    return warp * sssColor;
}

void main(void)
{
    // diffuse color
    vec3 texColor = texture(diffuseTex, fs_in.uv).xyz;
    vec3 ambient = ambientColor;
    vec3 diffuse = vec3(0);
    vec3 specular = vec3(0);
    vec3 sss = vec3(0);

    // lighting
    vec3 N = normalize(fs_in.normalWS);
    for (int i = 0; i < lightNumber; ++i) {
        vec3 L = -normalize(lightDirection[i].xyz);
        float LdN = dot(L, N);
        vec3 LDiff = computeDiffuse(LdN);
        diffuse += LDiff * lightColor[i].xyz;

        vec3 LSpec = computeSpecular(L, N) * lightColor[i].xyz;
        specular += LSpec;

        sss += computeSSS(LdN);
    }
    diffuse = diffuse * diffuseColor * texColor;
    specular = specular * specularColor * texColor;
    // vec3 finalColor = ambient * 0 + diffuse + specular + sss;
    vec3 finalColor = ambient + diffuse + specular + sss;
    FragColor = vec4(finalColor, 1.0);
}
