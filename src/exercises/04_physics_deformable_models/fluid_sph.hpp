#pragma once

#include "../../exercises/base_exercise/base_exercise.hpp"
#include<cstring>

#ifdef EXERCISE_FLUID_SPH

const float cube_size = 0.5;
const int display_method = 2;

// Kernel

double smoothKernel(vcl::vec3 p, float h);
vcl::vec3 gradKernel(vcl::vec3 p, float h);
vcl::vec3 vertexInterpolation(float threshold, vcl::vec3 p1, vcl::vec3 p2, float val1, float val2);

// SPH Particle
struct particle_element
{
    vcl::vec3 p; // Position
    vcl::vec3 v; // Speed
    vcl::vec3 a; // Acceleration

    // local density and pression
    float rho;
    float pression;

    particle_element() : p{0,0,0},v{0,0,0},a{0,0,0},rho(0),pression(0) {}

};

// SPH simulation parameters
struct sph_parameters
{
    float h;     // influence distance of a particle
    float rho0;  // rest density
    float m;     // total mass of a particle
    float stiffness; // constant of tait equation (relation density / pression)
    float nu;    // viscosity parameter
    float omega; // Rotation speed
    float scale_factor;
    float c;
    float epsilon;
    float rotation_acceleration;
};

// Image used to display the water appearance
struct field_display
{
    vcl::image im;           // Image storage on CPU
    GLuint texture_id;       // Texture stored on GPU
    vcl::mesh_drawable quad; // Mesh used to display the texture
    float voxel_influence;
};

// User parameters available in the GUI
struct gui_parameters
{
    bool display_field;
    bool display_particles;
    bool save_field;
};

class uniform_grid
{
public:
    uniform_grid(){}
    uniform_grid(size_t n, const std::vector<float>& boundingBox, std::vector<particle_element>& particles); // Construct a grid in which each cell contains the pointers of the particles contained in the cell
    std::vector<particle_element*> findPotentialNeighbors(const vcl::vec3& p); // Using the acceleration structure to find potential neighbors of a given point

private:
    size_t n; // Resolution of the grid
    std::vector<std::vector<particle_element*>> cells; // n*n*n cells, each cell contains the pointers of the particles contained in the cell
    // Bounding box of all particles
    float xMin; // boundingBox[0]
    float xMax; // boundingBox[1]
    float yMin; // boundingBox[2]
    float yMax; // boundingBox[3]
    float zMin; // boundingBox[4]
    float zMax; // boundingBox[5]

    std::vector<size_t> findCellIndices(const vcl::vec3& p) const; // Find the indices (i, j, k) of the cell containing the point p
    inline size_t findCellIndex(std::vector<size_t> indices, size_t n){ // Find the corresponding index of the cells vector
        return indices[2]+indices[1]*n+indices[0]*n*n;
    }
};

struct scene_exercise : base_scene_exercise
{

    void setup_data(std::map<std::string,GLuint>& shaders, scene_structure& scene, gui_structure& gui);
    void frame_draw(std::map<std::string,GLuint>& shaders, scene_structure& scene, gui_structure& gui);
    void display(std::map<std::string,GLuint>& shaders, scene_structure& scene, gui_structure& gui);

    std::vector<particle_element> particles;
    sph_parameters sph_param;
    uniform_grid sph_grid;
    uniform_grid voxel_grid;


    void update_acceleration(); // Update particles acceleration
    void update_density(); // Update local density, using SPH
    void update_pression(); // Update local pression

    void initialize_sph(); // Choose parameters for SPH simulation
    void define_form(int form, vcl::vec3 base_position); // Define the initial form of the fluid at the begining of the simulation
    void initialize_field_image();
    void set_gui();
    float evaluate_display_field(const vcl::vec3& p);

    gui_parameters gui_param;
    field_display field_image;
    vcl::mesh_drawable sphere;
    vcl::mesh_drawable voxel;
    vcl::segments_drawable borders;
    vcl::mesh liquid;

    vcl::timer_event timer;
};

#endif
