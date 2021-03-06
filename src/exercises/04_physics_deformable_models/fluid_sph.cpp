
#include "fluid_sph.hpp"
// #include "transvoxel.hpp"
#include "mctable.h"

#include <random>
#include <math.h>
#include <string.h>

#ifdef EXERCISE_FLUID_SPH

using namespace vcl;
using index3 = std::array<unsigned int, 3>;

// Random value generator
std::default_random_engine generator;
std::uniform_real_distribution<float> distrib(0.0,1.0);

// Counter used to save image on hard drive
int counter_image = 0;

/*
 *  Kernel methods
 */
double smoothKernel(vcl::vec3 p, float h) {
    double dist = sqrt(double(dot(p,p)));
    if(dist > h)
        return 0;
    return 315/(64*M_PI*pow(double(h),3))*pow(1-pow(dist/double(h), 2), 3);
}

vcl::vec3 gradKernel(vcl::vec3 p, float h){
    double dist = sqrt(double(dot(p,p)));
    if(dist > h)
        return vcl::vec3(0, 0, 0);
    return -6*315/(64*M_PI*pow(double(h),5))*pow(1-pow(dist/double(h), 2), 2)*p;
}

/*
 *  Initialization parameters: choose parameters here
 */

void scene_exercise::initialize_sph()
{
    // Influence distance of a particle (size of the kernel)
    const float h = 0.05f;

    // Rest density (consider 1000 Kg/m^3)
    const float rho0 = 1000.0f;

    // Stiffness (consider ~2000 - used in tait equation)
    const float stiffness = 1750.0f;

    // Viscosity parameter
    const float nu = 0.05f;

    // Total mass of a particle (consider rho0 h^2)
    const float m = rho0*h*h;

    // Rotation speed for the centrifuge force
    const float omega = 0.0f;

    // Artificial rotation acceleration so that the particles can rotate
    const float rotation_acceleration = 4.f;
    // const float rotation_acceleration = 0.f;

    // Initial particle spacing (relative to h)
    const float c = 0.8;
//     const float c = 0.85f;

    // Scale the size of the particle cube
    const float scale_factor = 0.7;

    // Fill a square with particles

    const int start_from_corners = 0;
    const float epsilon = 1e-3f;

    sph_param.h = h;
    sph_param.rho0 = rho0;
    sph_param.nu = nu;
    sph_param.stiffness = stiffness;
    sph_param.m = m;
    sph_param.scale_factor = scale_factor;
    sph_param.omega = omega;
    sph_param.c = c;
    sph_param.epsilon = epsilon;
    sph_param.rotation_acceleration = rotation_acceleration;

    if (start_from_corners == 1) {

        define_form(0, (vec3){-cube_size+2*h, cube_size-scale_factor*cube_size-h, -cube_size+2*h});
        define_form(0, (vec3){cube_size-scale_factor*cube_size-h, cube_size-scale_factor*cube_size-h, -cube_size+2*h});
        define_form(0, (vec3){cube_size-scale_factor*cube_size-h, cube_size-scale_factor*cube_size-h, cube_size-scale_factor*cube_size-h});
        define_form(0, (vec3){-cube_size+2*h, cube_size-scale_factor*cube_size-h, cube_size-scale_factor*cube_size-h});

    } else {

        define_form(0, (vec3){-scale_factor*cube_size/2+h, -h, -scale_factor*cube_size/2+h});

    }
}

void scene_exercise::define_form(int form, vec3 base_position) {
    // succession of particles.push_back instructions
    float bpx = base_position.x;
    float bpy = base_position.y;
    float bpz = base_position.z;
    if (form == 0) {
        // square
        for(float x=bpx; x<bpx+sph_param.scale_factor*cube_size-sph_param.h; x=x+sph_param.c*sph_param.h)
        {
            for(float y=bpy; y<bpy+sph_param.scale_factor*cube_size-sph_param.h; y=y+sph_param.c*sph_param.h)
            {
                for (float z=bpz; z < bpz+sph_param.scale_factor*cube_size-sph_param.h; z=z+sph_param.c*sph_param.h) {
                    particle_element particle;
                    particle.p = {x+sph_param.epsilon*distrib(generator),y,z+sph_param.epsilon * distrib(generator)}; // a zero value in z position will lead to a 2D simulation
                    particles.push_back(particle);
                }
            }
        }
    }
}

void scene_exercise::setup_data(std::map<std::string,GLuint>& shaders, scene_structure& , gui_structure& gui)
{
    gui.show_frame_camera = false;
    shaders["segment_immediate_mode"] = create_shader_program("shaders/segment_immediate_mode/segment_immediate_mode.vert.glsl","shaders/segment_immediate_mode/segment_immediate_mode.frag.glsl");

    sphere = mesh_drawable( mesh_primitive_sphere(1.0f));
    voxel = mesh_drawable( mesh_primitive_parallelepiped({0., 0., 0.}));
    std::vector<vec3> borders_segments = {{-cube_size,-cube_size,-cube_size},{cube_size,-cube_size,-cube_size}, {cube_size,-cube_size,-cube_size},{cube_size,cube_size,-cube_size}, {cube_size,cube_size,-cube_size},{-cube_size,cube_size,-cube_size}, {-cube_size,cube_size,-cube_size},{-cube_size,-cube_size,-cube_size},
                                          {-cube_size,-cube_size,cube_size} ,{cube_size,-cube_size,cube_size},  {cube_size,-cube_size,cube_size}, {cube_size,cube_size,cube_size},  {cube_size,cube_size,cube_size}, {-cube_size,cube_size,cube_size},  {-cube_size,cube_size,cube_size}, {-cube_size,-cube_size,cube_size},
                                          {-cube_size,-cube_size,-cube_size},{-cube_size,-cube_size,cube_size}, {cube_size,-cube_size,-cube_size},{cube_size,-cube_size,cube_size}, {cube_size,cube_size,-cube_size},{cube_size,cube_size,cube_size},   {-cube_size,cube_size,-cube_size},{-cube_size,cube_size,cube_size}};
    borders = segments_gpu(borders_segments);
    borders.uniform_parameter.color = {0,0,0};


    timer.stop();
    initialize_sph();

    field_image.voxel_influence = 0.2f;

    gui_param.display_field = true;
    gui_param.display_particles = true;
    gui_param.save_field = false;


}

void scene_exercise::update_acceleration()
{
    // Weight
    const size_t N = particles.size();
    for(size_t i=0; i<N; ++i)
    {
        particles[i].a = vec3{0,-9.81f,0};
    }

    // Add contribution of SPH forces
    // ...
    update_density();
    update_pression();

    // Add pression force
    for(particle_element& part_i: particles){
        vec3 f_pression(0., 0., 0.);
        const std::vector<particle_element*> neighbors = sph_grid.findPotentialNeighbors(part_i.p);
        for(size_t j=0; j<neighbors.size(); j++){
            const particle_element part_j = *neighbors[j];
            const vec3 v = - sph_param.m * (part_i.pression/pow(part_i.rho, 2) + part_j.pression/pow(part_j.rho, 2)) * gradKernel(part_i.p - part_j.p, sph_param.h);
            f_pression += v;
        }
        part_i.a += f_pression;
    }

    // Add viscosity force
    for(particle_element& part_i: particles){
        vec3 f_viscosity(0., 0., 0.);
        const std::vector<particle_element*> neighbors = sph_grid.findPotentialNeighbors(part_i.p);
        for(size_t j=0; j<neighbors.size(); j++){
            const particle_element part_j = *neighbors[j];
            const vec3 v = sph_param.m/part_j.rho * dot(part_i.p-part_j.p, part_i.v-part_j.v)/(dot(part_i.p-part_j.p, part_i.p-part_j.p)+0.01*pow(sph_param.h, 2)) * gradKernel(part_i.p-part_j.p, sph_param.h);
            f_viscosity += v;
        }
        part_i.a += 2*sph_param.nu*f_viscosity;
    }

    // Add centrifuge force
    for(particle_element& part_i: particles)
    {
        const vec3 r = vec3(part_i.p.x, 0, part_i.p.z);
        const vec3 f_centrifuge = pow(sph_param.omega, 2) * r;
        part_i.a += f_centrifuge;
    }

    // Add a rotation force
    for(particle_element& part_i: particles)
    {
        if(part_i.p.y < -cube_size + 2*sph_param.h){
            const vec3 r = vec3(part_i.p.x, 0, part_i.p.z) / sqrt(pow(part_i.p.x, 2) + pow(part_i.p.z, 2));
            part_i.a += sph_param.rotation_acceleration * vec3(-r.z, 0, r.x);
        }
    }


}


void scene_exercise::frame_draw(std::map<std::string,GLuint>& shaders, scene_structure& scene, gui_structure& gui)
{
    const float dt = timer.update();
    set_gui();

    // Force constant time step
    float h = dt<=1e-6f? 0.0f : timer.scale*0.003f;

    // Create uniform grid
    int n = floor(cube_size * 2.2/sph_param.h);
    std::vector<float> boundingBox;
    for(int i=0; i<3; i++){
        boundingBox.push_back(-1.1*cube_size);
        boundingBox.push_back(1.1*cube_size);
    }
    sph_grid = uniform_grid(n, boundingBox, particles);

    // Update acceleration
    update_acceleration();

    // Numerical integration
    const float damping = 0.5f;
    const size_t N = particles.size();
    for(size_t k=0; k<N; ++k)
    {
        vec3& p = particles[k].p;
        vec3& v = particles[k].v;
        vec3& a = particles[k].a;

        v = (1-h*damping)*v + h*a;
        p = p + h*v;
    }

    // Collision with border
    for(size_t k=0; k<N; ++k)
    {
        vec3& p = particles[k].p;
        vec3& v = particles[k].v;

        if( p.y<-cube_size ) {p.y = -cube_size; v.y *= -0.5f;}
        if( p.x<-cube_size ) {p.x = -cube_size; v.x *= -0.5f;}
        if( p.x>cube_size )  {p.x = cube_size;  v.x *= -0.5f;}
        if( p.z<-cube_size ) {p.z = -cube_size; v.z *= -0.5f;}
        if( p.z>cube_size )  {p.z = cube_size;  v.z *= -0.5f;}
    }

    display(shaders, scene, gui);

}


void scene_exercise::display(std::map<std::string,GLuint>& shaders, scene_structure& scene, gui_structure& )
{
    borders.draw(shaders["curve"], scene.camera);

    // Display particles
    if(gui_param.display_particles)
    {
        const size_t N = particles.size();
        sphere.uniform_parameter.scaling = sph_param.h/5.0f;
        sphere.uniform_parameter.color = {0,0.5,1};
        for(size_t k=0; k<N; ++k)
        {
            sphere.uniform_parameter.translation = particles[k].p;
            sphere.draw(shaders["mesh"], scene.camera);
        }
    }

    // Grid structure for accelerate 3D vizualization
    int n = floor(cube_size * 2.2/(field_image.voxel_influence*1.1));
    std::vector<float> boundingBox;
    for(int i=0; i<3; i++){
        boundingBox.push_back(-1.1*cube_size);
        boundingBox.push_back(1.1*cube_size);
    }
    voxel_grid = uniform_grid(n, boundingBox, particles);

    // Update field image
    if(gui_param.display_field)
    {
        const float resolution = 30;
        float voxel_size = 2 * cube_size / resolution;
        voxel.uniform_parameter.scaling = voxel_size;
        float threshold = 3;

        if (display_method == 1) {
            // voxel
            threshold = 0.75;
            float x = -cube_size;
            float y = -cube_size;
            float z = -cube_size;
            for(int i = 0; i < resolution; i++) {
                for(int j = 0; j < resolution; j++) {
                    for(int k = 0; k < resolution; k++) {
                        const float f = evaluate_display_field({x -voxel_size / 2,y -voxel_size / 2,z -voxel_size / 2});
                        const float value = 0.5f * f;
                        if (value > threshold) {
                            float r = 1-value;
                            float g = 1-value;
                            float b = 1;
                            voxel.uniform_parameter.color = {r, g, b};
                            voxel.uniform_parameter.translation = {x, y, z};
                            voxel.draw(shaders["mesh"], scene.camera);
                        }
                        z += voxel_size;
                    }
                    z = -cube_size;
                    y += voxel_size;
                }
                z = -cube_size;
                y = -cube_size;
                x += voxel_size;
            }

        } else if (display_method == 2) {  // Marching cubes

            vec3 corners[8];
            vec3 vertList[12];
            float densities[8];
            std::vector<vec3> new_position;
            std::vector<index3> new_connectivity;
            
            float x = -cube_size;
            float y = -cube_size;
            float z = -cube_size;
            for (int i = 0; i < resolution; ++i) {
                for (int j = 0; j < resolution; ++j) {
                    for (int k = 0; k < resolution; ++k) {
                        int cubeindex = 0;
                        corners[0] = (vec3){x, y, z};
                        corners[1] = (vec3){x + voxel_size, y, z};
                        corners[2] = (vec3){x + voxel_size, y, z + voxel_size};
                        corners[3] = (vec3){x, y, z + voxel_size};
                        corners[4] = (vec3){x, y + voxel_size, z};
                        corners[5] = (vec3){x + voxel_size, y + voxel_size, z};
                        corners[6] = (vec3){x + voxel_size, y + voxel_size, z + voxel_size};
                        corners[7] = (vec3){x, y + voxel_size, z + voxel_size};

                        for (int l = 0; l < 8 ; ++l) {
                            densities[l] = evaluate_display_field(corners[l]);
                        }

                        if (densities[0] < threshold) {cubeindex |= 1;};
                        if (densities[1] < threshold) {cubeindex |= 2;};
                        if (densities[2] < threshold) {cubeindex |= 4;};
                        if (densities[3] < threshold) {cubeindex |= 8;};
                        if (densities[4] < threshold) {cubeindex |= 16;};
                        if (densities[5] < threshold) {cubeindex |= 32;};
                        if (densities[6] < threshold) {cubeindex |= 64;};
                        if (densities[7] < threshold) {cubeindex |= 128;};

                        if (edgeTable[cubeindex] != 0 && edgeTable[cubeindex] != 255) {

                            if ((edgeTable[cubeindex] & 1) == 1) {
                                vertList[0] = vertexInterpolation(threshold, corners[0], corners[1], densities[0], densities[1]);
                            }
                            if ((edgeTable[cubeindex] & 2) == 2) {
                                vertList[1] = vertexInterpolation(threshold, corners[1], corners[2], densities[1], densities[2]);
                            }
                            if ((edgeTable[cubeindex] & 4) == 4) {
                                vertList[2] = vertexInterpolation(threshold, corners[2], corners[3], densities[2], densities[3]);
                            }
                            if ((edgeTable[cubeindex] & 8) == 8) {
                                vertList[3] = vertexInterpolation(threshold, corners[3], corners[0], densities[3], densities[0]);
                            }
                            if ((edgeTable[cubeindex] & 16) == 16) {
                                vertList[4] = vertexInterpolation(threshold, corners[4], corners[5], densities[4], densities[5]);
                            }
                            if ((edgeTable[cubeindex] & 32) == 32) {
                                vertList[5] = vertexInterpolation(threshold, corners[5], corners[6], densities[5], densities[6]);
                            }
                            if ((edgeTable[cubeindex] & 64) == 64) {
                                vertList[6] = vertexInterpolation(threshold, corners[6], corners[7], densities[6], densities[7]);
                            }
                            if ((edgeTable[cubeindex] & 128) == 128) {
                                vertList[7] = vertexInterpolation(threshold, corners[7], corners[4], densities[7], densities[4]);
                            }
                            if ((edgeTable[cubeindex] & 256) == 256) {
                                vertList[8] = vertexInterpolation(threshold, corners[0], corners[4], densities[0], densities[4]);
                            }
                            if ((edgeTable[cubeindex] & 512) == 512) {
                                vertList[9] = vertexInterpolation(threshold, corners[1], corners[5], densities[1], densities[5]);
                            }
                            if ((edgeTable[cubeindex] & 1024) == 1024) {
                                vertList[10] = vertexInterpolation(threshold, corners[2], corners[6], densities[2], densities[6]);
                            }
                            if ((edgeTable[cubeindex] & 2048) == 2048) {
                                vertList[11] = vertexInterpolation(threshold, corners[3], corners[7], densities[3], densities[7]);
                            }

                            for (int l = 0; triTable[cubeindex][l] != -1; l += 3) {
                                int index = (int)(new_position.size());
                                new_position.push_back(vertList[triTable[cubeindex][l + 2]]);
                                new_position.push_back(vertList[triTable[cubeindex][l + 1]]);
                                new_position.push_back(vertList[triTable[cubeindex][l]]);
                                index3 new_triangle{ {index, index + 1, index + 2} };
                                new_connectivity.push_back(new_triangle);
                            }
                        }
                        z += voxel_size;
                    }

                    z = -cube_size;
                    y += voxel_size;

                }

                z = -cube_size;
                y = -cube_size;
                x += voxel_size;


            liquid.position = new_position;
            liquid.normal.clear();
            liquid.color.clear();
            liquid.texture_uv.clear();
            liquid.connectivity = new_connectivity;
            liquid.fill_color_uniform({0, 0, 1});
            liquid.fill_empty_fields();
            normal(liquid.position, liquid.connectivity, liquid.normal);
            mesh_drawable liquid_drawable = mesh_drawable(liquid);
            liquid_drawable.draw(shaders["mesh"], scene.camera);
            }

        }

        // Save texture on hard drive
        if( gui_param.save_field )
        {
            const std::string filename = vcl::zero_fill(std::to_string(counter_image),3);
            image_save_png("output/fluid/file_"+filename+".png",field_image.im);
            ++counter_image;
        }
    }

}

vec3 vertexInterpolation(float threshold, vec3 p1, vec3 p2, float val1, float val2) {
    float mu;
    vec3 p = {0, 0, 0};
    mu = (threshold - val1) / (val2 - val1);
    p.x = p1.x + mu * (p2.x - p1.x);
    p.y = p1.y + mu * (p2.y - p1.y);
    p.z = p1.z + mu * (p2.z - p1.z);
    return p;
}

void scene_exercise::set_gui()
{
    // Can set the speed of the animation
    float scale_min = 0.05f;
    float scale_max = 2.0f;
    ImGui::SliderScalar("Time scale", ImGuiDataType_Float, &timer.scale, &scale_min, &scale_max, "%.2f s");

    ImGui::Checkbox("Display field", &gui_param.display_field);
    ImGui::Checkbox("Display particles", &gui_param.display_particles);
    ImGui::Checkbox("Save field on disk", &gui_param.save_field);

    // Start and stop animation
    if (ImGui::Button("Stop"))
        timer.stop();
    if (ImGui::Button("Start"))
        timer.start();

}

// Fill an image with field computed as a distance function to the particles
float scene_exercise::evaluate_display_field(const vcl::vec3& p)
{
    // change this value to get wider/narrower visible particle influence on the texture
    const float d = 0.1f;

    float field = 0.0f;

    const std::vector<particle_element*> neighbors = voxel_grid.findPotentialNeighbors(p);

    size_t count = 0;
    for(size_t j=0; j<neighbors.size(); j++){
        const particle_element part_j = *neighbors[j];
        const vec3& pi = part_j.p;
        const float r  = norm(p-pi);
        const float u = r/d;
        if(u < field_image.voxel_influence/d){
            count += 1;
            field += std::exp(-u*u);
        }
    }
    return field;
}

// Initialize an image where local density is displayed
void scene_exercise::initialize_field_image()
{
    size_t N = 50; // Image dimension (adapt this value to set the texture precision)

    field_image.quad = mesh_primitive_quad({-cube_size,-cube_size,0},{cube_size,-cube_size,0},{-cube_size,cube_size,0});
    field_image.im.width = N;
    field_image.im.height = N;
    field_image.im.data.resize(4*field_image.im.width*field_image.im.height);
    field_image.texture_id = texture_gpu(field_image.im);

    field_image.quad.uniform_parameter.shading.ambiant = 1.0f;
    field_image.quad.uniform_parameter.shading.diffuse = 0.0f;
    field_image.quad.uniform_parameter.shading.specular = 0.0f;
}

/*
 * SPH methods
 */

void scene_exercise::update_density(){
    for(particle_element& part_i: particles){
        // Update density of particle i
        float density = 0;
        const std::vector<particle_element*> neighbors = sph_grid.findPotentialNeighbors(part_i.p);
        for(size_t j=0; j<neighbors.size(); j++){
            const particle_element part_j = *neighbors[j];
            density += sph_param.m * smoothKernel(part_i.p - part_j.p, sph_param.h);
        }
        part_i.rho = density;
    }
}

void scene_exercise::update_pression(){
    for(particle_element& part_i: particles){
        // Update pression of particle i
        part_i.pression = sph_param.stiffness * pow(part_i.rho/sph_param.rho0 - 1, 1.7);
    }
}

/*
 * Uniform grid class
 */

std::vector<size_t> uniform_grid::findCellIndices(const vec3& p) const{
    std::vector<size_t> indices;
    indices.push_back(floor((p.x-xMin)/(xMax-xMin) * n));
    indices.push_back(floor((p.y-yMin)/(yMax-yMin) * n));
    indices.push_back(floor((p.z-zMin)/(zMax-zMin) * n));
    return indices;
}

uniform_grid::uniform_grid(size_t n, const std::vector<float>& boundingBox, std::vector<particle_element>& particles){
    this->n = n;
    // Creating the bounding box
    xMin = boundingBox[0];
    xMax = boundingBox[1];
    yMin = boundingBox[2];
    yMax = boundingBox[3];
    zMin = boundingBox[4];
    zMax = boundingBox[5];

    // Filling the cells with particles
    cells.resize(n*n*n); // cell (i, j, k) is accessed with the index k+N*j+N^2*i
    for(particle_element& part: particles){
        // Add the pointer of the particle in the cell
        int index = findCellIndex(findCellIndices(part.p), n);
        cells[index].push_back(&part);
    }
}

std::vector<particle_element*> uniform_grid::findPotentialNeighbors(const vcl::vec3& p){
    std::vector<size_t> indices = findCellIndices(p);

    const int N = n;
    const int i = indices[0];
    const int j = indices[1];
    const int k = indices[2];

    std::vector<particle_element*> neighborsPointers = cells[findCellIndex(indices, n)];
    // Corners
    if(i-1>=0 && j-1>=0 && k-1>=0){
        for(particle_element* neighbor: cells[(k-1)+n*(j-1)+n*n*(i-1)]){
            neighborsPointers.push_back(neighbor);
        }
    }
    if(i-1>=0 && j-1>=0 && k+1<N){
        for(particle_element* neighbor: cells[(k+1)+n*(j-1)+n*n*(i-1)]){
            neighborsPointers.push_back(neighbor);
        }
    }
    if(i-1>=0 && j+1<N && k-1>=0){
        for(particle_element* neighbor: cells[(k-1)+n*(j+1)+n*n*(i-1)]){
            neighborsPointers.push_back(neighbor);
        }
    }
    if(i-1>=0 && j+1<N && k+1<N){
        for(particle_element* neighbor: cells[(k+1)+n*(j+1)+n*n*(i-1)]){
            neighborsPointers.push_back(neighbor);
        }
    }
    if(i+1<N && j-1>=0 && k-1>=0){
        for(particle_element* neighbor: cells[(k-1)+n*(j-1)+n*n*(i+1)]){
            neighborsPointers.push_back(neighbor);
        }
    }
    if(i+1<N && j-1>=0 && k+1<N){
        for(particle_element* neighbor: cells[(k+1)+n*(j-1)+n*n*(i+1)]){
            neighborsPointers.push_back(neighbor);
        }
    }
    if(i+1<N && j+1<N && k-1>=0){
        for(particle_element* neighbor: cells[(k-1)+n*(j+1)+n*n*(i+1)]){
            neighborsPointers.push_back(neighbor);
        }
    }
    if(i+1<N && j+1<N && k+1<N){
        for(particle_element* neighbor: cells[(k+1)+n*(j+1)+n*n*(i+1)]){
            neighborsPointers.push_back(neighbor);
        }
    }

    // Edges, k fixed
    if(i-1>=0 && j-1>=0){
        for(particle_element* neighbor: cells[k+n*(j-1)+n*n*(i-1)]){
            neighborsPointers.push_back(neighbor);
        }
    }
    if(i-1>=0 && j+1<N){
        for(particle_element* neighbor: cells[k+n*(j+1)+n*n*(i-1)]){
            neighborsPointers.push_back(neighbor);
        }
    }
    if(i+1<N && j-1>=0){
        for(particle_element* neighbor: cells[k+n*(j-1)+n*n*(i+1)]){
            neighborsPointers.push_back(neighbor);
        }
    }
    if(i+1<N && j+1<N){
        for(particle_element* neighbor: cells[k+n*(j+1)+n*n*(i+1)]){
            neighborsPointers.push_back(neighbor);
        }
    }

    // Edges, i fixed
    if(k-1>=0 && j-1>=0){
        for(particle_element* neighbor: cells[(k-1)+n*(j-1)+n*n*i]){
            neighborsPointers.push_back(neighbor);
        }
    }
    if(k-1>=0 && j+1<N){
        for(particle_element* neighbor: cells[(k-1)+n*(j+1)+n*n*i]){
            neighborsPointers.push_back(neighbor);
        }
    }
    if(k+1<N && j-1>=0){
        for(particle_element* neighbor: cells[(k+1)+n*(j-1)+n*n*i]){
            neighborsPointers.push_back(neighbor);
        }
    }
    if(k+1<N && j+1<N){
        for(particle_element* neighbor: cells[(k+1)+n*(j+1)+n*n*i]){
            neighborsPointers.push_back(neighbor);
        }
    }

    // Edges, j fixed
    if(k-1>=0 && i-1>=0){
        for(particle_element* neighbor: cells[(k-1)+n*j+n*n*(i-1)]){
            neighborsPointers.push_back(neighbor);
        }
    }
    if(k-1>=0 && i+1<N){
        for(particle_element* neighbor: cells[(k-1)+n*j+n*n*(i+1)]){
            neighborsPointers.push_back(neighbor);
        }
    }
    if(k+1<N && i-1>=0){
        for(particle_element* neighbor: cells[(k+1)+n*j+n*n*(i-1)]){
            neighborsPointers.push_back(neighbor);
        }
    }
    if(k+1<N && i+1<N){
        for(particle_element* neighbor: cells[(k+1)+n*j+n*n*(i+1)]){
            neighborsPointers.push_back(neighbor);
        }
    }

    // Centers
    if(i-1>=0){
        for(particle_element* neighbor: cells[k+n*j+n*n*(i-1)]){
            neighborsPointers.push_back(neighbor);
        }
    }
    if(i+1<N){
        for(particle_element* neighbor: cells[k+n*j+n*n*(i+1)]){
            neighborsPointers.push_back(neighbor);
        }
    }
    if(j-1>=0){
        for(particle_element* neighbor: cells[k+n*(j-1)+n*n*i]){
            neighborsPointers.push_back(neighbor);
        }
    }
    if(j+1<N){
        for(particle_element* neighbor: cells[k+n*(j+1)+n*n*i]){
            neighborsPointers.push_back(neighbor);
        }
    }
    if(k-1>=0){
        for(particle_element* neighbor: cells[(k-1)+n*j+n*n*i]){
            neighborsPointers.push_back(neighbor);
        }
    }
    if(k+1<N){
        for(particle_element* neighbor: cells[(k+1)+n*j+n*n*i]){
            neighborsPointers.push_back(neighbor);
        }
    }
    return neighborsPointers;
}

#endif
