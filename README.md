# cubify-style

Download [libigl](https://github.com/libigl/libigl) and place it in the libs folder. you should be able to build and compile

Rendering Instruction:
The scene is defined in DemoScene.cpp under graphics/scenes/DemoScene.cpp. I suggest you only modify this file.

**The cubify process will automatically start and you can press P to pause or resume it.**

You can modify the constructor of DemoScene to add geometries.

The general process of adding geometries is that:

```
std::shared_ptr<ToonMaterial> toon = std::make_shared<ToonMaterial>();  // create a toon material
// Modify the material property like diffuse color.
// Or you can load a diffuse texture.
std::shared_ptr<Material> mat = std::dynamic_pointer_cast<Material>(toon); // convert to a Material Pointer
addPrimitive(mesh_file_path, mat)        // add primitive.
```

`By default, all added geometries will be send to the cubifyprocessor later.`

`If you want to add some geometries that isn't cubified. You can add it to the m_noCubifySet.`

Here is the property of ToonMaterial:
You can also only modify the colors, these attributes all have default values.
|Property|Data Type|Meaning|
|---|---|---|
|diffuseColor|Vector3f|The diffuse color|
|specularColor|Vector3f|The specular color|
|ambiendColor|Vector3f|The ambient color|
|sssColor|Vector3f|The subsurface scattering (SSS) color|
|specularPower|float|The power of blinn-phong specular term|
|sssPower|float|The power of SSS, it's the weight used to multiply the SSS shading|
|sssRadius|float [0,1]|The radius of SSS|
|colorLevels|Vector3f|Explained below|
|colorPowers|Vector4f|Explained below|
|colorRadius|Vector3f|Explained below|

In toon shading, I used a color ramp. A more conventional way is to use a 2D texture as lookup table, but I just compute it by code.

The dot(L, N) is mapped to several levels in toon shading. **colorLevels** defines the 3 boundaries of the levels which split the [-1,1] into 4 subspaces, so please make sure for **colorLevels**, its 3 elements are in increasing order.

The dot(L, N) will be mapped to one of **colorPowers** or the mix of two of them.

|colorPowers.x|colorPowers.y|colorPowers.z|colorPowers.w|

These five boundaries are **-1, colorLevels.x, colorLevelx.y, colorLevel.z, 1** from left to right.

**colorRadius** is used to soften the change of colorPowers.

You can refer to res/shaders/toon.frag, func: computeDiffuse for the implementation.

But I think you can just leave these 3 attributes untouched.

## Textures

If you want to use textures, first you need to make sure the UV coordinates are bijective, which means each vertex only has one UV coordinates, because when read in the obj file, the UV coordinates will be truncated to fit the number of vertices.

To use the diffuse texture, you can simply call:

```
ToonMaterial::loadDiffuseTexture(const std::string& texture_path)
```

Then the shader will use this texture.

