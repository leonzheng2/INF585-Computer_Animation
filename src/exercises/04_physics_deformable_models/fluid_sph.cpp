
#include "fluid_sph.hpp"

#include <random>

#ifdef EXERCISE_FLUID_SPH

using namespace vcl;



// Random value generator
std::default_random_engine generator;
std::uniform_real_distribution<float> distrib(0.0,1.0);
// Counter used to save image on hard drive
int counter_image = 0;


void scene_exercise::initialize_sph()
{
    // Influence distance of a particle (size of the kernel)
    const float h = 0.1f;

    // Rest density (consider 1000 Kg/m^3)
    const float rho0 = 1000.0f;

    // Stiffness (consider ~2000 - used in tait equation)
    const float stiffness = 2000.0f;

    // Viscosity parameter
    const float nu = 2.0f;

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

    std::cout << "coucou" << std::endl;

    // Add pression force
    for(particle_element part_i: particles){
        vec3 f_pression(0., 0., 0.);
        for(particle_element part_j: particles){
            const vec3 v = - sph_param.m * (part_i.pression/pow(part_i.rho, 2) + part_j.pression/pow(part_j.rho, 2)) * gradKernel(part_i.p - part_j.p, sph_param.h);
            std::cout << gradKernel(part_i.p - part_j.p, sph_param.h) << std::endl;
            f_pression += v;
        }
        part_i.a += f_pression;
//        std::cout << f_pression << std::endl;
    }

}



void scene_exercise::frame_draw(std::map<std::string,GLuint>& shaders, scene_structure& scene, gui_structure& gui)
{
    const float dt = timer.update();
    set_gui();

    // Force constant time step
    float h = dt<=1e-6f? 0.0f : timer.scale*0.0003f;

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

        a = {0,-9.81f,0};

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
    for(particle_element part_i: particles){
        // Update density of particle i
        float density = 0;
        for(particle_element part_j: particles){
            density += sph_param.m * smoothKernel(part_j.p, sph_param.h);
        }
        part_i.rho = density;
    }
}

void scene_exercise::update_pression(){
    for(particle_element part_i: particles){
        // Update pression of particle i
        part_i.pression = sph_param.stiffness * pow(part_i.rho/sph_param.rho0, 7);
    }
}

#endif
