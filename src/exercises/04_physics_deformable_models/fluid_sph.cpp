
#include "fluid_sph.hpp"

#include <random>

#ifdef EXERCISE_FLUID_SPH

using namespace vcl;


// Random value generator
std::default_random_engine generator;
std::uniform_real_distribution<float> distrib(0.0,1.0);
// Counter used to save image on hard drive
int counter_image = 0;

//Kernel
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

void scene_exercise::initialize_sph()
{
    // Influence distance of a particle (size of the kernel)
    const float h = 0.05f;

    // Rest density (consider 1000 Kg/m^3)
    const float rho0 = 1000.0f;

    // Stiffness (consider ~2000 - used in tait equation)
    const float stiffness = 2000.0f;

    // Viscosity parameter
    const float nu = 0.05f;

    // Total mass of a particle (consider rho0 h^2)
    const float m = rho0*h*h;

    // Initial particle spacing (relative to h)
    const float c = 0.95f;


    // Fill a square with particles
    const float epsilon = 1e-3f;
    for(float x=h; x<1.0f-h; x=x+c*h)
    {
        for(float y=-1.0f+h; y<0.0f-h; y=y+c*h)
        {
            particle_element particle;
            particle.p = {x+epsilon*distrib(generator),y,0}; // a zero value in z position will lead to a 2D simulation
            particles.push_back(particle);
        }
    }


    sph_param.h = h;
    sph_param.rho0 = rho0;
    sph_param.nu = nu;
    sph_param.stiffness = stiffness;
    sph_param.m = m;
}


void scene_exercise::setup_data(std::map<std::string,GLuint>& shaders, scene_structure& , gui_structure& gui)
{
    gui.show_frame_camera = false;
    shaders["segment_immediate_mode"] = create_shader_program("shaders/segment_immediate_mode/segment_immediate_mode.vert.glsl","shaders/segment_immediate_mode/segment_immediate_mode.frag.glsl");

    sphere = mesh_drawable( mesh_primitive_sphere(1.0f));
    std::vector<vec3> borders_segments = {{-1,-1,-0.1f},{1,-1,-0.1f}, {1,-1,-0.1f},{1,1,-0.1f}, {1,1,-0.1f},{-1,1,-0.1f}, {-1,1,-0.1f},{-1,-1,-0.1f},
                                          {-1,-1,0.1f} ,{1,-1,0.1f},  {1,-1,0.1f}, {1,1,0.1f},  {1,1,0.1f}, {-1,1,0.1f},  {-1,1,0.1f}, {-1,-1,0.1f},
                                          {-1,-1,-0.1f},{-1,-1,0.1f}, {1,-1,-0.1f},{1,-1,0.1f}, {1,1,-0.1f},{1,1,0.1f},   {-1,1,-0.1f},{-1,1,0.1f}};
    borders = segments_gpu(borders_segments);
    borders.uniform_parameter.color = {0,0,0};


    initialize_sph();
    initialize_field_image();

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
        const std::vector<particle_element*> neighbors = grid.findPotentialNeighbors(part_i);
        std::cout << "Number of particles: " << particles.size() << std::endl;
        std::cout << "Number of neighbors: " << neighbors.size() << std::endl;
        for(particle_element* part_j: neighbors){
            const vec3 v = - sph_param.m * (part_i.pression/pow(part_i.rho, 2) + part_j->pression/pow(part_j->rho, 2)) * gradKernel(part_i.p - part_j->p, sph_param.h);
            f_pression += v;
        }
        part_i.a += f_pression;
    }

    // Add viscosity force
    for(particle_element& part_i: particles){
        vec3 f_viscosity(0., 0., 0.);
        for(particle_element part_j: particles){
            const vec3 v = sph_param.m/part_j.rho * dot(part_i.p-part_j.p, part_i.v-part_j.v)/(dot(part_i.p-part_j.p, part_i.p-part_j.p)+0.01*pow(sph_param.h, 2)) * gradKernel(part_i.p-part_j.p, sph_param.h);
            f_viscosity += v;
        }
        part_i.a += 2*sph_param.nu*f_viscosity;
    }

}


void scene_exercise::frame_draw(std::map<std::string,GLuint>& shaders, scene_structure& scene, gui_structure& gui)
{
    const float dt = timer.update();
    set_gui();

    // Force constant time step
    float h = dt<=1e-6f? 0.0f : timer.scale*0.003f;

    // Create uniform grid
    int n = 10;
    std::vector<float> boundingBox;
    for(int i=0; i<3; i++){
        boundingBox.push_back(-1.0f);
        boundingBox.push_back(1.0f);
    }
    grid = uniform_grid(n, boundingBox, particles);

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

        std::cout << a << std::endl;
//        a = {0,-9.81f,0};

        v = (1-h*damping)*v + h*a;
        p = p + h*v;
    }

    // Collision with border
    for(size_t k=0; k<N; ++k)
    {
        vec3& p = particles[k].p;
        vec3& v = particles[k].v;

        if( p.y<-1 ) {p.y = -1; v.y *= -0.5f;}
        if( p.x<-1 ) {p.x = -1; v.x *= -0.5f;}
        if( p.x>1 )  {p.x = 1;  v.x *= -0.5f;}
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



    // Update field image
    if(gui_param.display_field)
    {
        const size_t im_h = field_image.im.height;
        const size_t im_w = field_image.im.width;
        std::vector<unsigned char>& im_data = field_image.im.data;
        for(size_t ky=0; ky<im_h; ++ky)
        {
            for(size_t kx=0; kx<im_w; ++kx)
            {
                const float x = 2.0f*kx/(im_w-1.0f)-1.0f;
                const float y = 1.0f-2.0f*ky/(im_h-1.0f);

                const float f = evaluate_display_field({x,y,0.0f});
                // adapt this value to set the iso-value of interest for the liquid surface
                const float value = 0.5f*f;

                float r = 1-value;
                float g = 1-value;
                float b = 1;


                im_data[4*(kx+im_w*ky)]   = static_cast<unsigned char>(255*std::max(std::min(r,1.0f),0.0f));
                im_data[4*(kx+im_w*ky)+1] = static_cast<unsigned char>(255*std::max(std::min(g,1.0f),0.0f));
                im_data[4*(kx+im_w*ky)+2] = static_cast<unsigned char>(255*std::max(std::min(b,1.0f),0.0f));
                im_data[4*(kx+im_w*ky)+3] = 255;
            }
        }



        // Display texture
        glBindTexture(GL_TEXTURE_2D, field_image.texture_id);
        glTexSubImage2D(GL_TEXTURE_2D, 0, 0,0, GLsizei(im_w), GLsizei(im_h), GL_RGBA, GL_UNSIGNED_BYTE, &im_data[0]);
        glGenerateMipmap(GL_TEXTURE_2D);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
        field_image.quad.draw(shaders["mesh"],scene.camera);
        glBindTexture(GL_TEXTURE_2D, scene.texture_white);


        // Save texture on hard drive
        if( gui_param.save_field )
        {
            const std::string filename = vcl::zero_fill(std::to_string(counter_image),3);
            image_save_png("output/fluid/file_"+filename+".png",field_image.im);
            ++counter_image;
        }
    }

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
    const size_t N = particles.size();
    for(size_t i=0; i<N; ++i)
    {
        const vec3& pi = particles[i].p;
        const float r  = norm(p-pi);
        const float u = r/d;
        if( u < 4)
            field += std::exp(-u*u);
    }
    return field;
}

// Initialize an image where local density is displayed
void scene_exercise::initialize_field_image()
{
    size_t N = 50; // Image dimension (adapt this value to set the texture precision)

    field_image.quad = mesh_primitive_quad({-1,-1,0},{1,-1,0},{-1,1,0});
    field_image.im.width = N;
    field_image.im.height = N;
    field_image.im.data.resize(4*field_image.im.width*field_image.im.height);
    field_image.texture_id = texture_gpu(field_image.im);

    field_image.quad.uniform_parameter.shading.ambiant = 1.0f;
    field_image.quad.uniform_parameter.shading.diffuse = 0.0f;
    field_image.quad.uniform_parameter.shading.specular = 0.0f;
}

void scene_exercise::update_density(){
    for(particle_element& part_i: particles){
        // Update density of particle i
        float density = 0;
        for(particle_element part_j: particles){
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

std::vector<int> uniform_grid::findCellIndices(particle_element part){
    std::vector<int> indices;
    indices.push_back(floor((part.p.x-xMin)/(xMax-xMin) * n));
    indices.push_back(floor((part.p.y-yMin)/(yMax-yMin) * n));
    indices.push_back(floor((part.p.z-zMin)/(zMax-zMin) * n));
    return indices;
}

uniform_grid::uniform_grid(int n, const std::vector<float>& boundingBox, std::vector<particle_element>& particles){
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
        int index = findCellIndex(findCellIndices(part), n);
        cells[index].push_back(&part);
    }
}

std::vector<particle_element*> uniform_grid::findPotentialNeighbors(particle_element part){
    std::vector<int> indices = findCellIndices(part);
    const int i = indices[0];
    const int j = indices[1];
    const int k = indices[2];

    // TODO clean the code

    std::vector<particle_element*> neighborsPointers = cells[findCellIndex(indices, n)];
    // Corners
    if(i-1>=0 && j-1>=0 && k-1>=0){
        for(particle_element* neighbor: cells[(k-1)+n*(j-1)+n*n*(i-1)]){
            neighborsPointers.push_back(neighbor);
        }
    }
    if(i-1>=0 && j-1>=0 && k+1<n){
        for(particle_element* neighbor: cells[(k+1)+n*(j-1)+n*n*(i-1)]){
            neighborsPointers.push_back(neighbor);
        }
    }
    if(i-1>=0 && j+1<n && k-1>=0){
        for(particle_element* neighbor: cells[(k-1)+n*(j+1)+n*n*(i-1)]){
            neighborsPointers.push_back(neighbor);
        }
    }
    if(i-1>=0 && j+1<n && k+1<n){
        for(particle_element* neighbor: cells[(k+1)+n*(j+1)+n*n*(i-1)]){
            neighborsPointers.push_back(neighbor);
        }
    }
    if(i+1<n && j-1>=0 && k-1>=0){
        for(particle_element* neighbor: cells[(k-1)+n*(j-1)+n*n*(i+1)]){
            neighborsPointers.push_back(neighbor);
        }
    }
    if(i+1<n && j-1>=0 && k+1<n){
        for(particle_element* neighbor: cells[(k+1)+n*(j-1)+n*n*(i+1)]){
            neighborsPointers.push_back(neighbor);
        }
    }
    if(i+1<n && j+1<n && k-1>=0){
        for(particle_element* neighbor: cells[(k-1)+n*(j+1)+n*n*(i+1)]){
            neighborsPointers.push_back(neighbor);
        }
    }
    if(i+1<n && j+1<n && k+1>=0){
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
    if(i-1>=0 && j+1<n){
        for(particle_element* neighbor: cells[k+n*(j+1)+n*n*(i-1)]){
            neighborsPointers.push_back(neighbor);
        }
    }
    if(i+1<n && j-1>=0){
        for(particle_element* neighbor: cells[k+n*(j-1)+n*n*(i+1)]){
            neighborsPointers.push_back(neighbor);
        }
    }
    if(i+1<n && j+1<n){
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
    if(k-1>=0 && j+1<n){
        for(particle_element* neighbor: cells[(k-1)+n*(j+1)+n*n*i]){
            neighborsPointers.push_back(neighbor);
        }
    }
    if(k+1<n && j-1>=0){
        for(particle_element* neighbor: cells[(k+1)+n*(j-1)+n*n*i]){
            neighborsPointers.push_back(neighbor);
        }
    }
    if(k+1<n && j+1<n){
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
    if(k-1>=0 && i+1<n){
        for(particle_element* neighbor: cells[(k-1)+n*j+n*n*(i+1)]){
            neighborsPointers.push_back(neighbor);
        }
    }
    if(k+1<n && i-1>=0){
        for(particle_element* neighbor: cells[(k+1)+n*j+n*n*(i-1)]){
            neighborsPointers.push_back(neighbor);
        }
    }
    if(k+1<n && i+1<n){
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
    if(i+1<n){
        for(particle_element* neighbor: cells[k+n*j+n*n*(i+1)]){
            neighborsPointers.push_back(neighbor);
        }
    }
    if(j-1>=0){
        for(particle_element* neighbor: cells[k+n*(j-1)+n*n*i]){
            neighborsPointers.push_back(neighbor);
        }
    }
    if(j+1<n){
        for(particle_element* neighbor: cells[k+n*(j+1)+n*n*i]){
            neighborsPointers.push_back(neighbor);
        }
    }
    if(k-1>=0){
        for(particle_element* neighbor: cells[(k-1)+n*j+n*n*i]){
            neighborsPointers.push_back(neighbor);
        }
    }
    if(k+1<n){
        for(particle_element* neighbor: cells[(k+1)+n*j+n*n*i]){
            neighborsPointers.push_back(neighbor);
        }
    }

    return neighborsPointers;
}



#endif
