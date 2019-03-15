# INF585-Computer_Animation
Physical-based fluid simulation using SPH.

The base code is the vcl/ folder from the TD sessions. All the project work is in the directory "src/exercises/04_physics_deformable_models/".

To change the parameters for fluid simulation, please change the variables in the method "void scene_exercise::initialize_sph()" in "fluid_sph.cpp".

At the beginning of the file "fluid_sph.hpp", one can change cube_size and display_method. Choose 1 to visualize voxels, and 2 for marching cubes.

![Fluid rendering with voxel](https://github.com/leonzheng2/INF585-Computer_Animation/blob/master/fluid_voxel2.png "Fluid rendering with voxel")

![Fluid rendering with marching cubes](https://github.com/leonzheng2/INF585-Computer_Animation/blob/master/fluid_mc.png "Fluid rendering with voxel")
