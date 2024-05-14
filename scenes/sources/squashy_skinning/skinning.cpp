


#include "skinning.hpp"
#ifdef SCENE_SQUASHY_SKINNING

#include "skinning_loader.hpp"
#include "helper_skeleton.hpp"

#include <fstream>

#ifndef _WIN32
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wfloat-equal"
#pragma GCC diagnostic ignored "-Weffc++"
#endif
#include <glm/glm.hpp>
#include <glm/ext.hpp>
#ifndef _WIN32
#pragma GCC diagnostic pop
#endif

#ifdef _WIN32
#pragma warning(disable : 4267)
#endif


using namespace vcl;

float vcl::generatePoisson(float averageArrivalTime) {
    std::random_device rd;
    std::mt19937 gen(rd());

    std::exponential_distribution<double> exp (1.0/averageArrivalTime);

    return exp(gen);
}


float scene_model::generate_fake_input(float t)  //AICI BANUIESC CA E UNGHIUL PT CAT SE INDOAIE CHESTIA
{
    float T = gui_param.wave_gui.T; // period of the pattern
    float dt = gui_param.wave_gui.dt; // duration of the pattern

    if (t < 0)
        return 0.0f;

    float t_periodic = t - T * int(t / T);
    //std::cout << t<<" " << t_periodic << std::endl;
    float value = 0.0f;
    if (t_periodic < dt) {
        float u = t_periodic / dt;
        float angle_max = 2*3.14159f;
        value = 1 - std::cos(u * angle_max);
    }        

    return value;
}


vec3 deformation_flappy_linear_speed(float w_flappy, vec3 const& speed)
{
    vec3 const deformation = - w_flappy * speed;

    return deformation;
}


vec3 deformation_squashy_linear_speed(float w_squashy, vec3 const& speed, float const speed_norm, vec3 const& p_vertex, vec3 const& center_of_mass) //FOR LINEAR BONE MOTIONS PAGE 6/13
{
    float const squash_factor = w_squashy * speed_norm;  //w_squashy = k_squash in the paper
    float const elongation_scaling = 1+squash_factor;  //SQUASH_FACTOR = s
    float const squeeze_scaling = 1/std::sqrt(1+squash_factor);  //1/sqrt(1 + s)

    mat3 const Id = mat3::identity();
    mat3 const S = {elongation_scaling,0,0, 0,squeeze_scaling,0, 0,0,squeeze_scaling}; // THE SCALING MATRIX
    mat3 const R = rotation_between_vector_mat3({1,0,0}, speed); // Rotation to the correct frame
    mat3 const T = R * S * transpose(R); // Deformation matrix


    vec3 const deformation = (T-Id) * (p_vertex-center_of_mass);  //equality (14) in the paper

    //    T - Id = R S R' - Id     in the paper
    //p_vertex = p^u
    //center_of_mass = c_i

    return deformation;
}

vec3 deformation_flappy_rotation_speed(float w_flappy, vec3 const& p_vertex, vec3 const& p_joint,
                                       vec3 const& un_angular_speed, float vertex_speed_norm, float flappy_max_angle)
{
    vec3 const u_joint_vertex = p_vertex - p_joint;
    mat3 const Id = mat3::identity();

    vec3 const rotation_center = p_joint + dot(u_joint_vertex, un_angular_speed)*un_angular_speed;
    float rotation_angle = vertex_speed_norm * w_flappy ;

    if(rotation_angle>flappy_max_angle)
        rotation_angle = flappy_max_angle + (rotation_angle-flappy_max_angle)/vertex_speed_norm;


    mat3 const R = rotation_from_axis_angle_mat3( un_angular_speed, -rotation_angle );
    vec3 const deformation = (R-Id) * (p_vertex-rotation_center);  //BANUIESC CA E FORMULA (17) DE LA PAG 7/13

    return deformation;
}

vec3 deformation_squashy_rotation_speed(float squashing_power,
                                        vec3 const& p_vertex, vec3 const& p_joint,
                                        vec3 const& un_medial, vec3 const& un_angular_speed,
                                        float vertex_speed_norm, vec3 const& center_of_mass, int squash_type)  //FOR ROTATING BONE MOTIONS
{
    mat3 const Id = mat3::identity();

    // Scaling matrix
    float const squash_factor = squashing_power * vertex_speed_norm;  //squash_factor = s = k_squash * norm
    float const elongation = 1+squash_factor;  //elongation = 1 + s
    float const squeeze = 1/(1+squash_factor);  //squeeze = 1 / (1 + s )
    mat3 const S = { elongation,0,0, 0,squeeze,0, 0,0,squeeze };  //A DOUA MATRICE DE LA PAGINA 6/13

    // Rotate to correct frame
    vec3 const u_elongation = cross(un_medial, un_angular_speed);
    if( norm(u_elongation)< 1e-2f)
        return {0,0,0};

    vec3 const un_elongation = normalize(u_elongation);
    vec3 const un_squeeze = cross(un_medial, un_elongation);
    mat3 const R = mat3(un_elongation, un_squeeze, un_medial);

    // Deformation matrix
    mat3 const T = R*S*transpose(R);  //T = R S R'

    // Center of scaling
    vec3 const p_medial = squash_type==0?  //Pr(p^u) = operator projecting a point into the medial axis
                p_joint + dot(p_vertex-p_joint,un_medial)*un_medial:
                center_of_mass;

    vec3 const deformation = (T-Id) * (p_vertex-p_medial);  //THE FORMULA

    return deformation;
}


vec3 deformation_bending_from_angular_velocity(vec3 const& angular_velocity, vec3 const& p_vertex, vec3 const& p_joint)
{
    vec3 const angular_velocity_direction = normalize(angular_velocity);

    vec3 const u_joint_vertex = p_vertex - p_joint;
    vec3 rotation_center = p_joint + dot(u_joint_vertex, angular_velocity_direction) * angular_velocity_direction;
    float rotation_angle = norm(angular_velocity) * norm(p_vertex - rotation_center);

    mat3 const R = rotation_from_axis_angle_mat3(angular_velocity_direction, -rotation_angle);
    vec3 const deformation = (R - mat3::identity()) * (p_vertex - rotation_center);

    return deformation;
}

// Damien: This is the filter we used initially - it is actually not the correct one (see the paper)
float scene_model::filter(int k)
{
    float dt = 0.017f;
    float frequency = gui_param.oscillation_gui.frequency;
    //float magnitude = gui_param.oscillation_gui.magnitude;
    float attenuation = gui_param.oscillation_gui.attenuation;


    return std::exp(-attenuation * k * dt) * std::cos( frequency * k * dt);

}

void scene_model::convolution( const std::string& filename, int nr_samples ) {
    std::vector<float> convolution_result;
    convolution_result.reserve(nr_samples);
    int i, j;
    float sum;

    for ( i = 0; i < nr_samples; i++ ) {
        float sum = 0.0;

        //for ( j = 0; j <= i; j++ ) {
            //sum = sum + kappa(i - j);
            //sum += (((i-j) >= 0) ? 1.0f : 0.0f) * kappa(j);
        //}

        convolution_result.push_back( kappa(i) );

    }

    std::ofstream output_file(filename );

    if (output_file.is_open()) {
        output_file << "kappa" << std::endl;
        for ( const auto& value : convolution_result ) {
            output_file << value << std::endl;
        }
        output_file.close();
       
    }
}

std::vector<vec3> global_angular_velocities_not_filtered;  //GLOBAL VECTOR TO KEEP THE ANGULAR VELOCITIES NOT FILTERED
std::vector<vec3> global_angular_velocities_filtered;  //GLOBAL VECTOR TO KEEP THE ANGULAR VELOCITIES FILTERED

void scene_model::angularVelocity(const std::string& filename, int nr_samples) {
    //HERE I WILL CODE SUCH THAT WE GET THE ANGULAR VELOCITY THAT WE APPLY ON THE JOINT AND THE ANGULAR VELOCITY OF 
    //THE SAME JOINT BUT AFTER FILTERING

    std::ofstream output_file( filename );

    if (output_file.is_open()) {
        //output_file << "Unfiltered Angular Velocity, Filtered Angular Velocity" << std::endl;
        for (int i = 0; i < global_angular_velocities_not_filtered.size(); ++i) {
            output_file << norm(global_angular_velocities_not_filtered[i]) << "," << norm(global_angular_velocities_filtered[i]) << std::endl;
        }
        output_file.close();
    }
}

// Damien: This is a corrected version of the time filter defined as the response to a unit step
float scene_model::kappa(int k) 
{
    float dt = 0.017; // assume 60 fps (this is hard-coded so it is bad)
    float frequency = gui_param.oscillation_gui.frequency;
    float attenuation = gui_param.oscillation_gui.attenuation;

    float f_slope = gui_param.oscillation_gui.frequency_slope;   //aici e un nr care poate sa creasca ca frcventa sa creasca si ea
    f_slope = 0;
    float df = k*dt*f_slope;
    float u = (frequency+df)*dt*k * 2 * 3.14159f;

    //******************************************************************************************

    //THIS IS FOR A SINUSOIDAL WAVE

    /*float result = std::exp(-attenuation * k * dt) * sin(u);

    return result;   //AICI CRED CA ESTE CE AM DE CAUTAT   */

    //*******************************************************************************************

    //THIS IS FOR THE TRIANGULAR WAVE

    /*float phase = std::fmod(u, 2 * 3.14159f);
    int direction = phase < 3.14159f ? 1 : -1; // 1 for increasing, -1 for decreasing
    float amplitude = std::exp(-attenuation * k * dt);
    float result = (2 / 3.14159f) * std::asin(std::sin(u)) * amplitude * direction;

    return result; */

    //*******************************************************************************************

    //NOW LET S MODIFY FOR THE RECTANGULAR WAVE

    float period = 1.0 / frequency;
    float halfPeriod = period / 2.0;
    float remainder = k  % (int) (period / dt);
    float constant = -attenuation * ( k - remainder) * dt;

    float result = std::exp(constant);
    
    if (sin(u) > 0) {
        return result;
    }
    return -result; 

    //*******************************************************************************************

    //THIS IS FOR A PENDULUM-LIKE WAVE

    /*float period = 1.0 / frequency;
    float halfPeriod = period / 2.0;
    float remainder = k % (int)(period / dt);
    float constant = -attenuation * (k) * dt;

    float result = std::exp(constant) * (sin(u) + sin(3 * u) / 9.0f);

    return result;*/

    //******************************************************************************************

    //BIP BIP COYOTE WAVE
    /*f_slope = 3;
    df = k * dt * f_slope;
    u = (frequency + df) * dt * k * 2 * 3.14159f;
    float result = std::exp(-attenuation * k * dt) * sin(u);

    return result;*/

}
float scene_model::filter_correct(int k)
{
    // Damien: It should be a differentiation of kappa with respect to time
    //  As I am lazy, I do it numerically

    float dt = 0.017; // assume 60 fps (this is hard-coded so it is bad)
    return (kappa(k)-kappa(k-1))/dt;
}



float scene_model::noise_filter(vec3 &curVel, float t) {
        float noisez = powf((perlin((t-0.6)*gui_param.wave_noise_frequency, 8, 0.45, 2.0)*2.0-1.0), 2.0);
        float noisex = powf((perlin((t+0.6)*gui_param.wave_noise_frequency, 8, 0.45, 2.0)*2.0-1.0), 2.0);
        float noise = powf((perlin(t*gui_param.wave_noise_frequency, 8, 0.45, 2.0)*2.0-1.0), 2.0);

        if(noise > 0.8+gui_param.wave_noise_threshold){
            curVel.y = gui_param.wave_noise_max_amplitude*0.8;
        } else if(noise > 0.6+gui_param.wave_noise_threshold) {
            curVel.y = gui_param.wave_noise_max_amplitude*0.6;
        } else if(noise > 0.5+gui_param.wave_noise_threshold) {
            curVel.y = gui_param.wave_noise_max_amplitude*0.5;
        } else if(noise > 0.3+gui_param.wave_noise_threshold) {
            curVel.y = gui_param.wave_noise_max_amplitude*0.3;
        }

        // if(noisez > 0.8+gui_param.wave_noise_threshold){
        //     curVel.z += gui_param.wave_noise_max_amplitude*0.8*0.5;
        // } else if(noisez > 0.6+gui_param.wave_noise_threshold) {
        //     curVel.z += gui_param.wave_noise_max_amplitude*0.6*0.5;
        // } else if(noisez > 0.5+gui_param.wave_noise_threshold) {
        //     curVel.z += gui_param.wave_noise_max_amplitude*0.5*0.5;
        // } else if(noisez > 0.3+gui_param.wave_noise_threshold) {
        //     curVel.z += gui_param.wave_noise_max_amplitude*0.3*0.5;
        // }

        // if(noisex > 0.8+gui_param.wave_noise_threshold){
        //     curVel.x += gui_param.wave_noise_max_amplitude*0.8*0.4;
        // } else if(noisex > 0.6+gui_param.wave_noise_threshold) {
        //     curVel.x += gui_param.wave_noise_max_amplitude*0.6*0.4;
        // } else if(noisex > 0.5+gui_param.wave_noise_threshold) {
        //     curVel.x += gui_param.wave_noise_max_amplitude*0.5*0.4;
        // } else if(noisex > 0.3+gui_param.wave_noise_threshold) {
        //     curVel.x += gui_param.wave_noise_max_amplitude*0.3*0.4;
        // }
        return noise;
}

void scene_model::noise_spike(joint_geometry &curVel, float &t) {

    // Rectangle spike
    // if(t > 0.0 && t < gui_param.wave_noise_threshold){
    //     curVel.p.y = gui_param.wave_noise_max_amplitude*0.8;
    //     // curVel.p.x -= gui_param.wave_noise_max_amplitude*0.8;
    //     // curVel.r.z -= pump_noise_sign*gui_param.wave_noise_max_amplitude*1.0f;
    // }

    // Gaussian noise
    // if(t > 0.0 && t < gui_param.wave_noise_threshold*4){
    //     float val = exp(-((1.0f/gui_param.wave_noise_threshold)*t-2.0f)*((1.0f/gui_param.wave_noise_threshold)*t-2.0f));
    //     curVel.p.y = gui_param.wave_noise_max_amplitude*val;
    //     // curVel.p.x -= gui_param.wave_noise_max_amplitude*0.8;
    //     curVel.r.z -= pump_noise_sign*gui_param.wave_noise_max_amplitude*val;
    // }


    if(t > 0.0 && t < gui_param.wave_noise_threshold){
        float val = sin(t*(3.14/gui_param.wave_noise_threshold));
        curVel.p.y = gui_param.wave_noise_max_amplitude*val;
        // curVel.p.x -= gui_param.wave_noise_max_amplitude*0.8;
        curVel.r.z -= pump_noise_sign*val*gui_param.wave_pump_noise_roll_amplitude;
    }



    // Triangle Spike
    // if(t > 0.0 && t < gui_param.wave_noise_threshold*2){
    //     float val = gui_param.wave_noise_threshold - fabs(gui_param.wave_noise_threshold - t);
    //     curVel.p.y = gui_param.wave_noise_max_amplitude*val;
    //     curVel.p.x -= gui_param.wave_noise_max_amplitude*val;
    //     curVel.r.z -=  pump_noise_sign*gui_param.wave_noise_max_amplitude*val*3.0f;
    // }



    if(t > gui_param.wave_noise_threshold){
        t = 0;
    }
}

void scene_model::velocity_skinning(float magnitude)
{
    float const default_squashy_factor = 0.1f * magnitude;
    float const default_flappy_factor = 0.4f * magnitude;

    float const w_squashy = default_squashy_factor * squashing_power;

    // reset deformation
    deformation_per_vertex.fill({0,0,0});

    assert_vcl_no_msg(weight_flappy.size()==skinning.rest_pose.size());
    size_t const N_joint = skeleton_current.size();

    for(size_t joint=0; joint<N_joint; ++joint)  //ANGULAR VELOCITY
    {
        vec3 const& p_joint = skeleton_speed_per_joint.data[joint].center;
        vec3 linear_speed = (skeleton_speed_per_joint.data[joint].linear_speed / timer_skeleton.scale + skeleton_fake_speed_per_joint.data[joint].linear_speed);

        //vec3 angular_speed = (skeleton_speed_per_joint.data[joint].angular_speed / timer_skeleton.scale + skeleton_fake_speed_per_joint.data[joint].angular_speed);
        auto& saved_rotation = skeleton_joint_rotation_speed[joint].saved_rotation;
        vec3 angular_speed = saved_rotation[saved_rotation.size()-1];  //ANGULAR VELOCITY NOT FILTERED
        if ( joint == picking.selected_joint_true )  //IF IT IS THE SELECTED JOINT
            global_angular_velocities_not_filtered.push_back(angular_speed);  //WE INSERT THE NOT FILTERED ANGULAR VELOCITY IN THE GLOBAL VECTOR

        // This part implements the oscillation skinning
        if (gui_param.oscillation_filter_speedup==false && gui_param.dynamic_type!=6)
        {
            int N_integral = 300; // number of samples used for the numerical convolution (we take 180 samples in the past, roughly 2-3s)

            vec3 v_integral = { 0,0,0 };
            if (saved_rotation.size() > N_integral) {
                float normalization = 0.0f;


                // This step computes the sum of all the filters \int h(t) dt for normalization purpose
                //  This step is actually performed for every vertex and for every joint: this is highly ineffective and it could be optimized in precomputing it only once.
                for (int k = 0; k < N_integral; ++k)
                    normalization += filter(k);

                // The actual convolution loop (v(t) = \int_x v(t-x) h(x) dx), with h(x) the normalized filter such that \int_x h(x) dx = 1
                for (int k = 0; k < N_integral; ++k)
                {
                    // Query the angular velocity at previous time
                    vec3 v = saved_rotation[saved_rotation.size() - 1 - k];

                    // Convolution step
                    v_integral += v * filter(k) / normalization;
                }

                // Final angular velocity is: current_one + oscillation_magnitude * (filtered_angular_velocity - current_one)
                angular_speed = angular_speed + (v_integral-angular_speed)*gui_param.oscillation_gui.magnitude;
            }
        }
        else if (gui_param.oscillation_filter_speedup == false && gui_param.dynamic_type==6)
        {
            int N_integral = 300; // number of samples used for the numerical convolution (we take 180 samples in the past, roughly 2-3s)

            vec3 v_integral = { 0,0,0 };
            if (saved_rotation.size() > N_integral) {

                // The actual convolution loop
                for (int k = 0; k < N_integral; ++k)
                {
                    // Query the angular velocity at previous time
                    vec3 v = saved_rotation[saved_rotation.size() - 1 - k];

                    // Convolution step
                    v_integral += v * filter_correct(k);
                }

                float A = gui_param.oscillation_gui.magnitude/10;  //AICI ESTE A-UL PE CARE IL CAUTAM
                angular_speed =  angular_speed + A * v_integral;  //ANGULAR VELOCITY FILTERED
                if (joint == picking.selected_joint_true)  //IF THE JOINT IS THE SELECTED ONE
                    global_angular_velocities_filtered.push_back(angular_speed); //WE INSERT THE FILTERED ANGULAR VELOCITY IN THE GLOBAL VECTOR
            }
        }
        else if (gui_param.oscillation_filter_speedup) 
        {
            angular_velocity_filtered_storage[joint].update(angular_speed, timer.t);
            if ( gui_param.typeOscillation == 1 )  //sinusoidal
                angular_speed = angular_speed + (angular_velocity_filtered_storage[joint].evaluateSinusoidal(timer.t) - angular_speed) * gui_param.oscillation_gui.magnitude;
            else if ( gui_param.typeOscillation == 2 )  //triangular
                angular_speed = angular_speed + (angular_velocity_filtered_storage[joint].evaluateTriangular(timer.t) - angular_speed) * gui_param.oscillation_gui.magnitude;
            else if( gui_param.typeOscillation == 3 )  //rectangular
                angular_speed = angular_speed + (angular_velocity_filtered_storage[joint].evaluateRectangular(timer.t) - angular_speed) * gui_param.oscillation_gui.magnitude;
            else if ( gui_param.typeOscillation == 4 )  //pendulum
                angular_speed = angular_speed + (angular_velocity_filtered_storage[joint].evaluatePendulum(timer.t) - angular_speed) * gui_param.oscillation_gui.magnitude;
            else if ( gui_param.typeOscillation == 5 ) //bipbip
                angular_speed = angular_speed + (angular_velocity_filtered_storage[joint].evaluateBipBip(timer.t) - angular_speed) * gui_param.oscillation_gui.magnitude;
        }

        // An additional term should be added here to compute the oscillation skinning on linear velocity
        
        // we would need to redo the same as above, but with someting like: skeleton_joint_speed[joint].saved_speed;
        auto& saved_oscillation_lienar = skeleton_joint_speed[joint].saved;
        {
            int N_integral = 300; // number of samples used for the numerical convolution (we take 180 samples in the past, roughly 2-3s)

            vec3 v_integral = { 0,0,0 };
            if (saved_oscillation_lienar.size() > N_integral) {
                float normalization = 0.0f;

                for (int k = 0; k < N_integral; ++k)
                    normalization += filter(k);

                for (int k = 0; k < N_integral; ++k)
                {
                    // Query the angular velocity at previous time
                    vec3 v = saved_oscillation_lienar[saved_oscillation_lienar.size() - 1 - k];

                    // Convolution step
                    v_integral += v * filter(k) / normalization;
                }

                //linear_speed = linear_speed + (v_integral-linear_speed)*gui_param.oscillation_gui.magnitude;
            }
        }


        

        if (joint != 0) {
            quaternion R_parent = skeleton_current[joint].r;
            angular_speed = R_parent.matrix() * angular_speed;
        }



        vec3 un_angular_speed = normalize(angular_speed);

        float const linear_speed_norm = norm(linear_speed);
        float angular_speed_norm = norm(angular_speed);

        // Damien Hack: I remove the squashy effect here
        //vec3 const& center_of_mass = position_center_of_mass(joint);//center_of_mass_per_joint[joint];
        //vec3 const u_medial = center_of_mass - p_joint;
        //if(norm(u_medial)<1e-4f)
        //    std::cout<<"small medial direction"<<std::endl;
        //vec3 const un_medial = normalize(u_medial);

        buffer<int> const& vertices = vertex_depending_on_joint.data[joint];
        buffer<float> const& vertices_weights = vertex_weight_depending_on_joint.data[joint];
        size_t const N_vertex_dependency = vertices.size();


        


        // Linear motion
        if(linear_speed_norm>1e-3f)
        {
            for(size_t k_vertex=0; k_vertex<N_vertex_dependency; ++k_vertex)
            {
                size_t const vertex = vertices.data[k_vertex];
                float const w_skinning = vertices_weights.data[k_vertex];

                vec3 const& p_vertex = save_skinning.data[vertex];
                float const w_flappy = default_flappy_factor * flapping_power * weight_flappy.data[vertex];



                if (gui_param.dynamic_type == 0) // oscillation
                {
                    vec3 const flappy = deformation_flappy_linear_speed(w_flappy, linear_speed);
                    deformation_per_vertex.data[vertex] += w_skinning * (flappy);
                }

                if(gui_param.dynamic_type==1) // wave
                {
                    float distance_vertex_to_joint = distance_to_joint[vertex][joint];
                    float distance_to_time_index_scaling = gui_param.wave_compression_magnitude;

                    float velocity_index_past_float = std::max(std::min(saved_rotation.size() - 2 - distance_vertex_to_joint * distance_to_time_index_scaling + 1, saved_rotation.size() - 2.0f), 0.0f);
                    int velocity_index_past = int(velocity_index_past_float);
                    float relative = velocity_index_past_float - velocity_index_past;


                    vec3 v1 = saved_oscillation_lienar[velocity_index_past + 1];
                    vec3 v0 = saved_oscillation_lienar[velocity_index_past];
                    vec3 v = (1 - relative) * v0 + relative * v1;

                    if (joint != 0) {
                        quaternion R_parent = skeleton_current[joint].r;
                        v = R_parent.matrix() * v;
                    }
                    
                    vec3 const flappy = deformation_flappy_linear_speed(w_flappy, v);
                    deformation_per_vertex.data[vertex] += w_skinning * (flappy);
                }

                if (gui_param.dynamic_type == 3) // oscillation + wave
                {

                    vec3 const flappy_oscillation = deformation_flappy_linear_speed(w_flappy, linear_speed);
                    //deformation_per_vertex.data[vertex] += w_skinning * (flappy_oscillation);

                    float distance_vertex_to_joint = distance_to_joint[vertex][joint];
                    float distance_to_time_index_scaling = gui_param.wave_compression_magnitude;

                    float velocity_index_past_float = std::max(std::min(saved_rotation.size() - 2 - distance_vertex_to_joint * distance_to_time_index_scaling + 1, saved_rotation.size() - 2.0f), 0.0f);
                    int velocity_index_past = int(velocity_index_past_float);
                    float relative = velocity_index_past_float - velocity_index_past;


                    vec3 v1 = saved_oscillation_lienar[velocity_index_past + 1];
                    vec3 v0 = saved_oscillation_lienar[velocity_index_past];
                    vec3 v = (1 - relative) * v0 + relative * v1;

                    if (joint != 0) {
                        quaternion R_parent = skeleton_current[joint].r;
                        v = R_parent.matrix() * v;
                    }

                    vec3 const flappy_wave = deformation_flappy_linear_speed(w_flappy, v);

                    float t = gui_param.propgation_oscillation_weight;
                    deformation_per_vertex.data[vertex] += w_skinning * ((1.0-t) * flappy_wave + t * flappy_oscillation);
                }


                // Damien Hack: I remove the squashy deformation for optimization purpose
                //  If we want to activate it, we also need to reactivate the update of the center of mass
                //vec3 const squashy = deformation_squashy_linear_speed(w_squashy, linear_speed, linear_speed_norm, p_vertex, center_of_mass);
                //deformation_per_vertex.data[vertex] += w_skinning * (flappy  + squashy);
            }
        }


        // Rotation motion
        if(angular_speed_norm>1e-3f)
        {
            for(size_t k_vertex=0; k_vertex<N_vertex_dependency; ++k_vertex)
            {
                size_t const vertex = vertices.data[k_vertex];
                float const w_skinning = vertices_weights.data[k_vertex];

                vec3 const& p_vertex = save_skinning.data[vertex];
                float const w_flappy = default_flappy_factor * flapping_power * weight_flappy.data[vertex];
                float const w_squashy_local = w_squashy * weight_flappy.data[vertex]; // hacky


                if (gui_param.dynamic_type == 0 || gui_param.dynamic_type == 6) // oscillation
                {

                    vec3 const u_joint_vertex = p_vertex - p_joint;
                    vec3 const vertex_speed = cross(angular_speed, u_joint_vertex);
                    float const vertex_speed_norm = norm(vertex_speed);

                    if (angular_speed_norm > 1e-2f) {
                        // Flappy
                        //vec3 const flappy  = deformation_flappy_linear_speed(w_flappy, norm(vertex_speed)*vertex_speed);
                        vec3 const flappy = deformation_flappy_rotation_speed(w_flappy, p_vertex, p_joint, un_angular_speed, vertex_speed_norm, gui_param.flappy_max_angle);

                        deformation_per_vertex.data[vertex] += w_skinning * (flappy);

                        // Damien Hack: I remove the squashy deformation for optimization purpose
                        //  If we want to activate it, we also need to reactivate the update of the center of mass
                        // vec3 const squashy = deformation_squashy_rotation_speed(w_squashy_local, p_vertex, p_joint, un_medial, un_angular_speed, vertex_speed_norm, center_of_mass, gui_param.squash_around);

                        //deformation_per_vertex.data[vertex] += w_skinning * (flappy + squashy);
                    }

                }

                if(gui_param.dynamic_type==1) // wave
                {
                    float distance_vertex_to_joint = distance_to_joint[vertex][joint];
                    float distance_to_time_index_scaling = gui_param.wave_compression_magnitude;

                    float velocity_index_past_float = std::max( std::min(saved_rotation.size() - 2 - distance_vertex_to_joint * distance_to_time_index_scaling+1, saved_rotation.size() - 2.0f), 0.0f);
                    int velocity_index_past = int(velocity_index_past_float);
                    float relative = velocity_index_past_float - velocity_index_past;


                    vec3 v1 = saved_rotation[velocity_index_past + 1];
                    vec3 v0 = saved_rotation[velocity_index_past];
                    vec3 angular_speed = (1 - relative) * v0 + relative * v1;

                    if (joint != 0) {
                        quaternion R_parent = skeleton_current[joint].r;
                        angular_speed = R_parent.matrix() * angular_speed;
                    }

                    float angular_speed_norm = norm(angular_speed);
                    if (angular_speed_norm > 1e-2f) {
                        vec3 un_angular_speed = angular_speed / angular_speed_norm;
                        vec3 u_joint_vertex = p_vertex - p_joint;
                        vec3 vertex_speed = cross(angular_speed, u_joint_vertex);
                        float vertex_speed_norm = norm(vertex_speed);

                        vec3 const flappy = deformation_flappy_rotation_speed(w_flappy, p_vertex, p_joint, un_angular_speed, vertex_speed_norm, gui_param.flappy_max_angle);

                        deformation_per_vertex.data[vertex] += w_skinning * (flappy);

                    }
                }

                if(gui_param.dynamic_type==3){
                    

                    // oscillation
                    vec3 const u_joint_vertex = p_vertex - p_joint;
                    vec3 const vertex_speed = cross(angular_speed, u_joint_vertex);
                    float const vertex_speed_norm = norm(vertex_speed);
                    vec3 flappy_oscillation;
                    vec3 flappy_propgation;

                    if (angular_speed_norm > 1e-2f) {
                        // Flappy
                        //vec3 const flappy  = deformation_flappy_linear_speed(w_flappy, norm(vertex_speed)*vertex_speed);
                        vec3 flappy = deformation_flappy_rotation_speed(w_flappy, p_vertex, p_joint, un_angular_speed, vertex_speed_norm, gui_param.flappy_max_angle);

                        flappy_oscillation = w_skinning * (flappy);

                        // Damien Hack: I remove the squashy deformation for optimization purpose
                        //  If we want to activate it, we also need to reactivate the update of the center of mass
                        // vec3 const squashy = deformation_squashy_rotation_speed(w_squashy_local, p_vertex, p_joint, un_medial, un_angular_speed, vertex_speed_norm, center_of_mass, gui_param.squash_around);

                        //deformation_per_vertex.data[vertex] += w_skinning * (flappy + squashy);
                    }



                    //propgation
                    float distance_vertex_to_joint = distance_to_joint[vertex][joint];
                    float distance_to_time_index_scaling = gui_param.wave_compression_magnitude;

                    float velocity_index_past_float = std::max( std::min(saved_rotation.size() - 2 - distance_vertex_to_joint * distance_to_time_index_scaling+1, saved_rotation.size() - 2.0f), 0.0f);
                    int velocity_index_past = int(velocity_index_past_float);
                    float relative = velocity_index_past_float - velocity_index_past;


                    vec3 v1 = saved_rotation[velocity_index_past + 1];
                    vec3 v0 = saved_rotation[velocity_index_past];
                    vec3 angular_speed = (1 - relative) * v0 + relative * v1;

                    if (joint != 0) {
                        quaternion R_parent = skeleton_current[joint].r;
                        angular_speed = R_parent.matrix() * angular_speed;
                    }

                    float angular_speed_norm = norm(angular_speed);
                    if (angular_speed_norm > 1e-2f) {
                        vec3 un_angular_speed = angular_speed / angular_speed_norm;
                        vec3 u_joint_vertex = p_vertex - p_joint;
                        vec3 vertex_speed = cross(angular_speed, u_joint_vertex);
                        float vertex_speed_norm = norm(vertex_speed);

                        flappy_propgation = deformation_flappy_rotation_speed(w_flappy, p_vertex, p_joint, un_angular_speed, vertex_speed_norm, gui_param.flappy_max_angle);

                    }
                    float t = gui_param.propgation_oscillation_weight;
                    deformation_per_vertex.data[vertex] += w_skinning * ((1.0-t) * flappy_propgation + t * flappy_oscillation);
                }



            }

        }

        //Scaling

    }

    // Damien Hack: This part was implementing the propagating wave: it is currently commented out.
    // 
    //Wave
    //for (size_t joint = 0; joint < N_joint; ++joint)
#ifdef false
    int joint = picking.selected_joint_true;
    if(joint>=0)
    {
        size_t const N_vertex_dependency = vertex_depending_on_joint.data[joint].size();
        for (size_t k_v = 0; k_v < N_vertex_dependency; ++k_v)
        {
            int const vertex = vertex_depending_on_joint.data[joint].data[k_v];

            float time = timer.t;
            //float d = norm(p_vertex-p_joint);
            float d = distance_to_joint[vertex][joint];
            if (d < 0) {
                //std::cout << "Error distance for vertex "<< vertex <<" for join" << joint << std::endl;
            }
            float value = generate_fake_input(time - gui_param.wave_gui.speed * d);// 6 * d - 3 * time);

            float const w_skinning = vertex_weight_depending_on_joint.data[joint].data[k_v];

            vec3 const& p_vertex = save_skinning.data[vertex];
            vec3 const& p_joint = skeleton_speed_per_joint.data[joint].center;

            vec3 deformation;
            {
                float T = gui_param.wave_gui.T; // period of the pattern
                float dt = gui_param.wave_gui.dt; // duration of the pattern

                deformation = { 0,0,0 };
                float t_diff = gui_param.wave_gui.speed * time - d*3;
                float t_periodic = t_diff - T * int(t_diff / T);
                if (t_periodic >0 && t_periodic < dt) {
                    float t_loc = t_periodic / dt;

                    deformation.x = 0.5*sin(2*3.14* t_loc);    
                    deformation.y = 1.0*(1-cos(2*3.14* t_loc));
                }
            }



            //deformation_per_vertex.data[vertex] += w_skinning * gui_param.wave_gui.A * deformation;


            /*
            {
                float value = gui_param.wave_gui.A * norm(velocity);
                vec3 p0 = closest_distance_to_bone.p[vertex];
                mat3 S = mat3::from_scaling({ 1 + value, 1 + value, 1 + value });
                vec3 d = S * (p_vertex - p0) - (p_vertex - p0);
                deformation_per_vertex.data[vertex] += w_skinning * d;
            }*/

        }
    }
#endif


    // Full propagation + oscillation
    if (gui_param.dynamic_type == 4) 
    {
        for (size_t joint = 0; joint < N_joint; ++joint)
        {
            vec3 const& p_joint = skeleton_speed_per_joint.data[joint].center;
            buffer<int> const& vertices = vertex_depending_on_joint.data[joint];
            buffer<float> const& vertices_weights = vertex_weight_depending_on_joint.data[joint];
            size_t const N_vertex_dependency = vertices.size();

            //Update current angular velocity
            vec3 current_angular_velocity = skeleton_joint_rotation_speed[joint].saved_rotation[skeleton_joint_rotation_speed[joint].saved_rotation.size() - 1];
            angular_velocity_filtered_storage[joint].update(current_angular_velocity, timer.t);


            // Change angular velocity to global frame
            if (joint != 0) {
                quaternion R_parent = skeleton_current[joint].r;
                current_angular_velocity = R_parent.matrix() * current_angular_velocity;
            }


            for (size_t k_vertex = 0; k_vertex < N_vertex_dependency; ++k_vertex)
            {
                size_t const vertex = vertices.data[k_vertex];
                float distance_vertex_to_joint = distance_to_joint[vertex][joint];
                float propagation_velocity = gui_param.wave_compression_magnitude;
                float time_delay = distance_vertex_to_joint / propagation_velocity;

                                
                //query delayed angular velocity
                vec3 delayed_angular_velocity;
                {
                    float velocity_index_past_float = std::max(std::min(skeleton_joint_rotation_speed[joint].saved_rotation.size() - 2 - distance_vertex_to_joint * propagation_velocity + 1, skeleton_joint_rotation_speed[joint].saved_rotation.size() - 2.0f), 0.0f);
                    int velocity_index_past = int(velocity_index_past_float);
                    float relative = velocity_index_past_float - velocity_index_past;


                    vec3 v1 = skeleton_joint_rotation_speed[joint].saved_rotation[velocity_index_past + 1];
                    vec3 v0 = skeleton_joint_rotation_speed[joint].saved_rotation[velocity_index_past];
                    delayed_angular_velocity = (1 - relative) * v0 + relative * v1;

                    if (joint != 0) {
                        quaternion R_parent = skeleton_current[joint].r;
                        delayed_angular_velocity = R_parent.matrix() * delayed_angular_velocity;
                    }
                }
                vec3 filtered_angular_velocity = angular_velocity_filtered_storage[joint].evaluateSinusoidal(timer.t);
                vec3 angular_velocity;
                angular_velocity = delayed_angular_velocity +(filtered_angular_velocity - delayed_angular_velocity) * gui_param.oscillation_gui.magnitude;



                if (norm(angular_velocity) > 1e-2f) {



                    float w_flappy = default_flappy_factor * flapping_power * weight_flappy.data[vertex];
                    vec3 const& p_vertex = save_skinning.data[vertex];
                    vec3 un_angular_velocity = normalize(angular_velocity);
                    float vertex_speed_norm = norm(cross(angular_velocity, p_vertex - p_joint));

                    vec3 flappy = deformation_flappy_rotation_speed(w_flappy, p_vertex, p_joint, un_angular_velocity, vertex_speed_norm, gui_param.flappy_max_angle);

                    float const w_skinning = vertices_weights.data[k_vertex];
                    deformation_per_vertex.data[vertex] += w_skinning * (flappy);
                }

            }

            
        }

    }


    // Offset oscillation
    if (gui_param.dynamic_type == 5)
    {
        std::map<int, float> lookup_table_offset_frequency;
        // Cow
        // lookup_table_offset_frequency[11] = 0.5;
        // lookup_table_offset_frequency[12] = -0.5;
        // lookup_table_offset_frequency[10] = 0.5;

        // Girafe
        // lookup_table_offset_frequency[9] = 0.5;
        // lookup_table_offset_frequency[11] = 1.5;
        // lookup_table_offset_frequency[21] = -0.5;
        // lookup_table_offset_frequency[30] = 0.5;
        // lookup_table_offset_frequency[38] = -0.75;
        // lookup_table_offset_frequency[46] = 0.75;

        // Clown
        lookup_table_offset_frequency[19] = -0.5;
        lookup_table_offset_frequency[23] = 0.25;
        lookup_table_offset_frequency[24] = 0.45;



        float default_frequency = max_joint_filtered_value::frequency;
        for (size_t joint = 0; joint < N_joint; ++joint)
        {
            vec3 const& p_joint = skeleton_speed_per_joint.data[joint].center;
            buffer<int> const& vertices = vertex_depending_on_joint.data[joint];
            buffer<float> const& vertices_weights = vertex_weight_depending_on_joint.data[joint];
            size_t const N_vertex_dependency = vertices.size();

            //Update current angular velocity
            vec3 current_angular_velocity = skeleton_joint_rotation_speed[joint].saved_rotation[skeleton_joint_rotation_speed[joint].saved_rotation.size() - 1];
            for (size_t subjoint = 0; subjoint < N_joint; ++subjoint) {
                
                if (lookup_table_offset_frequency.find(subjoint) != lookup_table_offset_frequency.end()) {
                    max_joint_filtered_value::frequency = default_frequency+ lookup_table_offset_frequency[subjoint];
                    angular_velocity_filtered_storage_subjoints[joint][subjoint].update(current_angular_velocity, timer.t);
                    max_joint_filtered_value::frequency = default_frequency;
                }
                else {
                    angular_velocity_filtered_storage_subjoints[joint][subjoint].update(current_angular_velocity, timer.t);
                }
            }
            




            for (size_t k_vertex = 0; k_vertex < N_vertex_dependency; ++k_vertex)
            {
                size_t const vertex = vertices.data[k_vertex];
                int subjoint = joint;
                vec3 filtered_angular_velocity = angular_velocity_filtered_storage_subjoints[joint][subjoint].evaluateSinusoidal(timer.t);
                vec3 angular_velocity = current_angular_velocity + (filtered_angular_velocity - current_angular_velocity) * gui_param.oscillation_gui.magnitude;
                if (joint != 0) {
                    quaternion R_parent = skeleton_current[joint].r;
                    angular_velocity = R_parent.matrix() * angular_velocity;
                }

                auto& joint_dependance = rig_extended_to_ancestor_joint[vertex];
                vec3 angular_subvelocity;
                for (int k_dep = 0; k_dep < joint_dependance.size() - 1; ++k_dep) {
                    int subjoint = joint_dependance[k_dep];
                    if (lookup_table_offset_frequency.find(subjoint) != lookup_table_offset_frequency.end()) {

                        max_joint_filtered_value::frequency = default_frequency + lookup_table_offset_frequency[subjoint];
 
                        vec3 filtered_angular_subvelocity = angular_velocity_filtered_storage_subjoints[joint][subjoint].evaluateSinusoidal(timer.t);
                        max_joint_filtered_value::frequency = default_frequency;

                        vec3 angular_subvelocity = current_angular_velocity + (filtered_angular_subvelocity - current_angular_velocity) * gui_param.oscillation_gui.magnitude;
                        if (joint != 0) {
                            quaternion R_parent = skeleton_current[joint].r;
                            angular_subvelocity = R_parent.matrix() * angular_subvelocity;
                        }


                        float w = rig_extended_to_ancestor_weight[vertex][k_dep];
                        angular_velocity = w* angular_subvelocity+(1-w)* angular_velocity;
                    }
                }


                
                if (norm(angular_velocity) > 1e-3f) {

                    float w_flappy = default_flappy_factor * flapping_power * weight_flappy.data[vertex];
                    vec3 const& p_vertex = save_skinning.data[vertex];
                    vec3 un_angular_velocity = normalize(angular_velocity);
                    float vertex_speed_norm = norm(cross(angular_velocity, p_vertex - p_joint));

                    vec3 flappy = deformation_flappy_rotation_speed(w_flappy, p_vertex, p_joint, un_angular_velocity, vertex_speed_norm, gui_param.flappy_max_angle);

                    float const w_skinning = vertices_weights.data[k_vertex];
                    deformation_per_vertex.data[vertex] += w_skinning * (flappy);
                }

            }


        }

    }


    for(size_t k=0; k<deformation_per_vertex.size(); ++k)
        skinning.deformed.position[k] += deformation_per_vertex[k];




}










void scene_model::frame_draw(std::map<std::string,GLuint>& , scene_structure& _scene, gui_structure& )
{


    timer_measurement["full"].tic();

    scene = _scene;
    timer.update();

    if(gui_param.animation)
        timer_skeleton.update();

    set_gui();

    float const t = timer.t;





    interpolate_skeleton_at_time_with_constraints(skeleton_local_current_before_deformation,
                                                  timer_skeleton.t,
                                                  skeleton.anim,
                                                  gui_param.interpolate,
                                                  record_joint_fixed);

    //    // case of recording -- fix some joints
    //    if(gui_param.record_anim && record_joint_fixed.size()>0)
    //    {
    //        for(int k=0; k<record_joint_fixed.size(); ++k){
    //            int j = record_joint_fixed[k];
    //            skeleton_local_current[j].r = record_joint_fixed_rotation[k];
    //            skeleton_local_current[j].p = record_joint_fixed_position[k];
    //        }
    //    }


    skeleton_local_current = skeleton_local_current_before_deformation;

    // symmetrize
    if(gui_param.symmetry)
    {
        int joint = picking.selected_joint_true;
        if(joint!=-1)
        {
            quaternion const qx = quaternion::axis_angle({1,0,0},3.14f);
            if(symmetrical_joint.find(joint)!=symmetrical_joint.end())
            {
                int const js = symmetrical_joint[joint];
                skeleton_local_interactive_deformation[js].r = qx * skeleton_local_interactive_deformation[joint].r * conjugate(qx);
            }
        }
    }


    // apply deformation to skeleton
    for(size_t k=0; k<skeleton_local_current.size(); ++k) {
        skeleton_local_current[k].p += skeleton_local_interactive_deformation[k].p;
        skeleton_local_current[k].r = skeleton_local_interactive_deformation[k].r * skeleton_local_current[k].r;
    }



   
    skeleton_current = local_to_global(skeleton_local_current, skeleton.connectivity);
    if(is_rotating){
        skeleton_current[0].r = qq;
    }




    // Update per-join velocity
    for(int k=0; k<int(skeleton_current.size()); ++k) {
        skeleton_joint_speed[k].add( skeleton_local_current[k].p, t );
        skeleton_joint_rotation_speed[k].add( skeleton_local_current[k].r, t );

        mat3 R_parent = mat3::identity();
        if(k>0)
            R_parent = skeleton_current[ skeleton.connectivity[k].parent ].r.matrix();

        skeleton_speed_per_joint[k].center = skeleton_current[k].p;
        skeleton_speed_per_joint[k].linear_speed  = R_parent * skeleton_joint_speed[k].avg_speed;
        skeleton_speed_per_joint[k].angular_speed = R_parent * skeleton_joint_rotation_speed[k].avg_rotation_speed;
    }

    vec3 train_scene_velocity = {-3.0f, 0, 0};
    float dt = timer.t - last_t;
    last_t = timer.t;

    start_t_noise += dt;

    if(gui_param.display_type == display_train){
        if(gui_param.wave_noise) {
            for(int joint = 0; joint < skeleton_current.data.size(); joint++){
               noise_filter(skeleton_current.data[joint].p, environment_train["train_track"].uniform.transform.translation.x);
            }
            if(start_t_noise == 0.0){
                start_t_noise = -0.7*generatePoisson(gui_param.wave_noise_frequency);
            }
        } else if(gui_param.pump_noise){
            for(int joint = 0; joint < skeleton_current.data.size(); joint++){
                noise_spike(skeleton_current.data[joint], start_t_noise);
                
            }
            if(start_t_noise == 0.0){
                if(pump_noise_sign == -1.0f){
                    pump_noise_sign = 1.0f;
                } else {
                    pump_noise_sign = -1.0f;
                }
                start_t_noise = -0.7*generatePoisson(gui_param.wave_noise_frequency);
            }
        }
    } else {
        if(gui_param.wave_noise) {
            for(int joint = 0; joint < skeleton_current.data.size(); joint++){
                noise_filter(skeleton_current.data[joint].p, timer.t);
            }
        }
    }
    
    skeleton_previous.add(skeleton_current);


    if(gui_param.display_rest_pose){
        // skeleton_local_interactive_deformation = skeleton_rest_pose;
        // skeleton_local_current_before_deformation = skeleton_rest_pose;
        skeleton_current = skeleton_rest_pose;
    }


    timer_measurement["skinning"].tic();
    if(gui_param.dual_quaternion)
        compute_skinning_dual_quaternion(skinning, skeleton_current, skeleton_rest_pose);
    else {
        compute_skinning(skinning, skeleton_current, skeleton_rest_pose);
        if (gui_param.dynamic_type == 2) {
            compute_skinning_wave_propagation();
        }
    }
    timer_measurement["skinning"].toc();



    character_visual.update_position(skinning.deformed.position);
    character_visual.update_normal(skinning.deformed.normal);


    save_skinning = skinning.deformed.position;

    //    // Update per-vertex velocity
    //    for(int k=0; k<int(skinning.deformed.position.size()); ++k)
    //        vertex_speed[k].add( skinning.deformed.position[k], t );

    /*
    character_visual.uniform.transform.translation = {0,0,-0.8f};
    if(gui_param.display_mesh) {
        glPolygonOffset( 1.0, 1.0 );
        GLuint const texture_id = (gui_param.display_texture? character_visual.texture_id : scene.texture_white);
        draw(character_visual, scene.camera, character_visual.shader, texture_id);
    }
    if(gui_param.display_wireframe) {
        glPolygonOffset( 1.0, 1.0 );
        draw(character_visual, scene.camera, shaders["wireframe_quads"]);
    }
    character_visual.uniform.transform.translation = {0,0,0};
    */


    timer_measurement["center of mass"].tic();
    update_center_of_mass();
    timer_measurement["center of mass"].toc();


    // Apply velocity skinning
    timer_measurement["velocity_skinning"].tic();
    velocity_skinning(1.0f);
    timer_measurement["velocity_skinning"].toc();


    timer_measurement["normals"].tic();
    normal(skinning.deformed.position, skinning.deformed.connectivity, skinning.deformed.normal);
    timer_measurement["normals"].toc();

    if(gui_param.display_deformed_surface)
    {
        character_visual.update_position(skinning.deformed.position);
        character_visual.update_normal(skinning.deformed.normal);
    }




    timer_measurement["drawing"].tic();



    scene.camera.apply_translation_orthogonal_to_screen_plane(-1.0f);
    scene.camera.apply_translation_in_screen_plane(0.5,0.5);
    mat4 lightCam = scene.camera.perspective.matrix()
            *
            scene.camera.view_matrix();
    scene.camera.apply_translation_in_screen_plane(-0.5,-0.5);
    scene.camera.apply_translation_orthogonal_to_screen_plane(1.0f);
    {

        glUseProgram(shaders["depth_map"]);
        glViewport(0, 0, SHADOW_WIDTH, SHADOW_HEIGHT); opengl_debug();
        glBindFramebuffer(GL_FRAMEBUFFER, depthMapFBO); opengl_debug();
        glClear(GL_DEPTH_BUFFER_BIT); opengl_debug();
        //glEnable(GL_DEPTH_TEST); opengl_debug();

        //glClearColor(0.1f, 0.1f, 0.1f, 1.0f); opengl_debug();
        //glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT); // we're not using the stencil buffer now
        //glEnable(GL_DEPTH_TEST); opengl_debug();

        opengl_debug();
        //GLint lightSpaceMatrixLocation = glGetUniformLocation(shaders["depth_map"], "lightSpaceMatrix");
        //glUniformMatrix4fv(lightSpaceMatrixLocation, 1, GL_FALSE, glm::value_ptr(lightSpaceMatrix));

        mat4 Id = mat4::identity();
        //scene.camera.apply_translation_in_screen_plane(0.2,0.2);
        uniform(shaders["depth_map"], "lightSpaceMatrix", lightCam);
        uniform(shaders["depth_map"], "model", Id);


        //glCullFace(GL_FRONT);
        draw(character_visual.data);
        //glCullFace(GL_BACK);

        //draw(character_visual, scene.camera, shaders["depth_map"]); opengl_debug();

        glBindFramebuffer(GL_FRAMEBUFFER, 0); opengl_debug(); // back to default
        glViewport(0, 0, scene.window_width, scene.window_height); opengl_debug();


    }




    //    {


    //        float near_plane = 1.0f, far_plane = 7.5f;
    //        glm::mat4 lightProjection = glm::ortho(-10.0f, 10.0f, -10.0f, 10.0f, near_plane, far_plane);

    //        glm::mat4 lightView = glm::lookAt(glm::vec3(-2.0f, 4.0f, -1.0f),
    //                                          glm::vec3( 0.0f, 0.0f,  0.0f),
    //                                          glm::vec3( 0.0f, 1.0f,  0.0f));


    //        glm::mat4 lightSpaceMatrix = lightProjection * lightView;


    //        glUseProgram(shaders["depth_map"]);
    //        glViewport(0, 0, SHADOW_WIDTH, SHADOW_HEIGHT);
    //        glBindFramebuffer(GL_FRAMEBUFFER, depthMapFBO);
    //        glClear(GL_DEPTH_BUFFER_BIT);
    //        //glActiveTexture(GL_TEXTURE0);


    //        mat4 Id = mat4::identity();

    //        GLint lightSpaceMatrixLocation = glGetUniformLocation(shaders["depth_map"], "lightSpaceMatrix");
    //        glUniformMatrix4fv(lightSpaceMatrixLocation, 1, GL_FALSE, glm::value_ptr(lightSpaceMatrix));
    //        uniform(shaders["depth_map"], "model", Id);

    //        draw(character_visual.data);
    //        //draw(character_visual, scene.camera, shaders["depth_map"]);




    //        glBindFramebuffer(GL_FRAMEBUFFER, 0);

    //        glViewport(0, 0, scene.window_width, scene.window_height);
    //        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    //    }

    /*
    glUseProgram(shaders["mesh_shadow"]);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    quad.uniform.color = {1.0,1.0,1.0};
    quad.texture_id = 0;
    //quad.uniform.shading = {1,0,0};
    //quad.texture_id = depthMap;//character_visual.texture_id;

    //uniform(shaders["mesh_shadow"], "texture_sampler", 0);
    //uniform(shaders["mesh_shadow"], "shadowMap", 1);
    glActiveTexture(GL_TEXTURE0+0);
    glBindTexture(GL_TEXTURE_2D, scene.texture_white);
    glActiveTexture(GL_TEXTURE0+1);
    glBindTexture(GL_TEXTURE_2D, depthMap);
    uniform(shaders["mesh_shadow"],"lightSpaceMatrix",lightCam);

    draw(quad, scene.camera, shaders["mesh_shadow"]);

    glUseProgram(0);
*/



    //recording
    {
        if(gui_param.record_anim){
            if(local_time_record<timer_skeleton.t_max)
                timer_skeleton.t = std::min(local_time_record, timer_skeleton.t_max-0.01f);
        }


        local_time_record += timer_recording.update();
        if(timer_recording.event)
        {

            record_rotation.push_back(skeleton_local_current[recorded_joint].r);
            record_position.push_back(skeleton_local_current[recorded_joint].p);

        }



    }




//character_visual.texture_id=0;
    if(gui_param.display_mesh) {
        //glPolygonOffset( 1.0, 1.0 );

        glUseProgram(shaders["mesh_shadow"]);
        character_visual.shader = shaders["mesh_shadow"];
        //character_visual.texture_id = 0;


        // if(character_visual.texture_id==0 || gui_param.display_texture==false)
        //     character_visual.texture_id = scene.texture_white;

        glActiveTexture(GL_TEXTURE0+0);
        glBindTexture(GL_TEXTURE_2D, character_visual.texture_id==0? scene.texture_white : character_visual.texture_id);

        glActiveTexture(GL_TEXTURE0+1);
        glBindTexture(GL_TEXTURE_2D, depthMap);
        uniform(shaders["mesh_shadow"],"lightSpaceMatrix",lightCam);

        //GLuint id_save = character_visual.texture_id;
        character_visual.texture_id = 0;
        //GLuint const texture_id = (gui_param.display_texture? character_visual.texture_id : scene.texture_white);
        GLuint const id_save = (gui_param.display_texture ? real_texture_id : scene.texture_white);
        draw(character_visual, scene.camera);//, character_visual.shader, texture_id);
        character_visual.texture_id = id_save;



        bool resetScene = false;
        float resetDistance = -50;
        if(environment_train["train_rails"].uniform.transform.translation.x < resetDistance){
            resetScene = true;
        }

        environment_train["train_rails"].uniform.transform.scaling_axis.x = 20.0f;

        if(gui_param.display_type == display_train){
            for (auto& element : environment_train) {
                element.second.uniform.transform.translation += dt*train_scene_velocity;

                // if (element.second.uniform.transform.translation.x < -22)
                //     element.second.uniform.transform.translation += vec3(16, 0, 0);
                
                if((element.first.compare(std::string("train_rails")) == 0 
                || element.first.compare(std::string("train_pillars")) == 0
                || element.first.compare(std::string("train_bars")) == 0) && resetScene) {
                    element.second.uniform.transform.translation = { -20,-0.12f ,0 };
                } 

                if((element.first.compare(std::string("train_rails2")) == 0 
                || element.first.compare(std::string("train_pillars2")) == 0
                || element.first.compare(std::string("train_bars2")) == 0) && resetScene) {
                    element.second.uniform.transform.translation = { 5.0f,-0.12f ,0 };
                }


                glUseProgram(shaders["mesh"]);

                //element.second.texture_id = scene.texture_white;
                //element.second.uniform.shading.diffuse = 0.0;
                //element.second.uniform.shading.specular = 0.0;
                //element.second.uniform.shading.ambiant = 0.8;

                glActiveTexture(GL_TEXTURE0+0);
                glBindTexture(GL_TEXTURE_2D, scene.texture_white);
                //glActiveTexture(GL_TEXTURE0 + 1);
                //glBindTexture(GL_TEXTURE_2D, scene.texture_white);
                //glBindTexture(GL_TEXTURE_2D, element.second.texture_id==0? scene.texture_white : element.second.texture_id);

                //glActiveTexture(GL_TEXTURE0+1);
                //glBindTexture(GL_TEXTURE_2D, depthMap);
                //uniform(shaders["mesh_shadow"],"lightSpaceMatrix", lightCam);

                if(element.first != "train_rails2")
                    draw(element.second, scene.camera, shaders["mesh"]);
                
            }
        }
    }

    glActiveTexture(GL_TEXTURE0+0);
    glBindTexture(GL_TEXTURE_2D, scene.texture_white);
    glActiveTexture(GL_TEXTURE0+1);
    glBindTexture(GL_TEXTURE_2D, scene.texture_white);



    //glUseProgram(shaders["mesh"]);
    //uniform(shaders["mesh"], "texture_sampler", 0);

    if(gui_param.display_wireframe) {
        glPolygonOffset( 1.0, 1.0 );
        draw(character_visual, scene.camera, shaders["wireframe"]);
    }

    glBindTexture(GL_TEXTURE_2D,scene.texture_white);



    // display displacement
    /*
    {
        arrow_line.scaling_width = 0.2f;
        sphere_visual.uniform.transform.scaling = 0.003f;

        for(int k=0; k<int(save_skinning.size()); k+=gui_param.display_deformation_arrow_vertex_offset)
        {
            vec3 const& p = save_skinning[k];
            vec3 const& d = deformation_per_vertex[k];

            if(gui_param.display_deformation_arrows)
                arrow_line.draw(p, p+d, {0.8f,0.8f,1.0f}, scene.camera);

            if(gui_param.display_deformation_target){
                if(norm(d)>1e-5f){
                    sphere_visual.uniform.transform.translation = p+d;
                    draw(sphere_visual, scene.camera);
                }
            }

        }
        arrow_line.scaling_width = 1.0f;
        sphere_visual.uniform.transform.scaling = 1.0f;
    }
    */

    // iterative trajectory
    {
        arrow_line.scaling_width = 0.2f;
        sphere_visual.uniform.transform.scaling = 0.003f;
        sphere_visual.uniform.color = {1,1,1};


        if(gui_param.curved_trajectory && gui_param.display_deformation_arrows)
        {
            buffer<vec3> previous_position = save_skinning;

            int N_alpha = 5; // sampling of the trajectory
            for(int k_alpha=0; k_alpha<N_alpha; ++k_alpha)
            {
                float const alpha = (k_alpha+1.0f)/N_alpha;
                velocity_skinning(alpha);

                for(int k=0; k<int(save_skinning.size()); k+=gui_param.display_deformation_arrow_vertex_offset)
                {
                    vec3 const& p0 = previous_position[k];
                    vec3 const& p1 = save_skinning[k] + deformation_per_vertex[k];

                    arrow_line.draw(p0, p1, {0.8f,0.8f,1.0f}, scene.camera);

                    if(k_alpha==N_alpha-1 && gui_param.display_deformation_target){
                        if(norm(deformation_per_vertex[k])>1e-3f){
                            sphere_visual.uniform.transform.translation = p1;
                            draw(sphere_visual, scene.camera);
                        }
                    }
                }

                previous_position = save_skinning + deformation_per_vertex;
            }
        }

        if(!gui_param.curved_trajectory && gui_param.display_deformation_target)
        {
            for(int k=0; k<int(save_skinning.size()); k+=gui_param.display_deformation_arrow_vertex_offset)
            {
                vec3 const& p = save_skinning[k];
                vec3 const& d = deformation_per_vertex[k];

                arrow_line.draw(p, p+d, {0.8f,0.8f,1.0f}, scene.camera);
            }
        }

        if(gui_param.display_deformation_target)
        {
            for(int k=0; k<int(save_skinning.size()); k+=gui_param.display_deformation_arrow_vertex_offset)
            {
                vec3 const& p = save_skinning[k];
                vec3 const& d = deformation_per_vertex[k];

                if(norm(d)>1e-5f){
                    sphere_visual.uniform.transform.translation = p+d;
                    draw(sphere_visual, scene.camera);
                }
            }
        }

        /*
        for(int k=0; k<int(save_skinning.size()); k+=gui_param.display_deformation_arrow_vertex_offset)
        {
            vec3 const& p = save_skinning[k];
            vec3 const& d = deformation_per_vertex[k];

            if(gui_param.display_deformation_arrows)
                arrow_line.draw(p, p+d, {0.8f,0.8f,1.0f}, scene.camera);

            if(gui_param.display_deformation_target){
                if(norm(d)>1e-5f){
                    sphere_visual.uniform.transform.translation = p+d;
                    draw(sphere_visual, scene.camera);
                }
            }
        }
        */
        arrow_line.scaling_width = 1.0f;
        sphere_visual.uniform.transform.scaling = 1.0f;
    }




    if(gui_param.x_ray)
        glClear(GL_DEPTH_BUFFER_BIT);

    if(gui_param.display_skeleton_bones)
        display_skeleton(skeleton_current, skeleton.connectivity, shaders, scene, segment_drawer);
    if(gui_param.display_skeleton_joints)
        display_joints(skeleton_current, scene, sphere_visual);
    if(gui_param.display_skeleton_frames){
        frame.uniform.transform.scaling = 0.08f*gui_param.frame_scaling;
        display_frames(skeleton_current, scene, frame);
    }

    if(gui_param.display_vertex_to_bone_correspondance) {
        display_vertex_to_bone_correspondance(skinning.deformed.position, skeleton_current, vertex_to_bone_correspondance, segment_drawer, shaders["segment_im"], scene.camera);
    }


    if(gui_param.display_skeleton_pyramid)
    {
        const size_t N = skeleton_current.size();
        for(size_t k=1; k<N; ++k)
        {
            int parent = skeleton.connectivity[k].parent;

            vec3 const& p0 = skeleton_current[parent].p;
            vec3 const& p1 = skeleton_current[k].p;

            vec3 const u = p1-p0;
            const float L = norm(p1-p0);
            const float Lr = std::min(L,0.2f);
            vec3 const un = u/L;

            if(picking.joint_hover==int(k))
                pyramid_skeleton_visual.uniform.color = {0.8f, 0.3f, 0.3f};
            else
                pyramid_skeleton_visual.uniform.color = {0.3f, 0.3f, 0.3f};
            pyramid_skeleton_visual.uniform.transform.translation = p0;
            pyramid_skeleton_visual.uniform.transform.scaling_axis = {Lr,Lr,L};
            pyramid_skeleton_visual.uniform.transform.rotation = rotation_between_vector_mat3({0,0,1}, un);

            draw(pyramid_skeleton_visual, scene.camera);
        }
    }



    {
        segment_drawer.uniform_parameter.color = {1,0,0};
        int N_joint = skeleton_current.size();
        for(int k_joint=0; k_joint<N_joint; ++k_joint)
        {
            vec3 const& p = skeleton_current[k_joint].p;

            if(gui_param.display_joint_linear_speed)
            {
                vec3 const& v = skeleton_speed_per_joint[k_joint].linear_speed + skeleton_fake_speed_per_joint[k_joint].linear_speed;
                arrow.draw(p, p+0.2f*v, {0.3f,0.6f,1.0f}, scene.camera);
            }
            if(gui_param.display_joint_angular_speed)
            {
                vec3 const& v = skeleton_speed_per_joint[k_joint].angular_speed + skeleton_fake_speed_per_joint[k_joint].angular_speed;
                arrow.draw(p, p+0.2f*v, {0.6f,1.0f,0.3f}, scene.camera);
            }
        }


        if(gui_param.display_center_of_mass)
        {
            bool is_translation = (gui_param.type_deformation==0 || picking.selected_joint==0);
            int joint = picking.selected_joint_true;
            if(joint!=-1)
            {
                vec3 const p_com = position_center_of_mass(joint);//center_of_mass_per_joint[joint];
                sphere_visual.uniform.color = {1.0f,0.4f,1.0f};
                sphere_visual.uniform.transform.scaling = 0.02f;
                sphere_visual.uniform.transform.translation = p_com;
                draw(sphere_visual, scene.camera);

                vec3 const& p_joint = skeleton_current[joint].p;
                vec3 const u_joint_com = p_com-p_joint;
                if(!is_translation && gui_param.squash_around==0)
                    arrow_line.draw(p_joint - u_joint_com*0.5f, p_joint + u_joint_com*2.0f, {0.8f,0.2f,0.8f}, scene.camera);

            }
        }
    }



    //  X-RAY display beside this point
    glClear(GL_DEPTH_BUFFER_BIT);

    //    segment_drawer.uniform_parameter.color = {0,0,0};
    //    if(gui_param.display_joint_speed){
    //        segment_drawer.uniform_parameter.color = {0,0,1};
    //        display_joint_speed(skeleton_current, skeleton_joint_speed, segment_drawer, shaders["segment_im"], scene.camera);
    //    }
    //    if(gui_param.display_joint_acceleration){
    //        segment_drawer.uniform_parameter.color = {0,1,1};
    //        display_joint_acceleration(skeleton_current, skeleton_joint_speed, segment_drawer, shaders["segment_im"], scene.camera);
    //    }


    /*
    {
        if(skeleton.connectivity.size()>9)
        {
            //int joint_picked = picking.selected_joint;
            int joint = 9;//skeleton.connectivity[joint_picked].parent;
            if(joint!=-1)
            {
                vec3 const& p_joint = skeleton_current[joint].p;
                vec3 const& bar = center_of_mass_per_joint[joint];
                vec3 const& angular_speed = skeleton_speed_per_joint[joint].angular_speed;

                sphere_visual.uniform.color = {1,1,0};
                sphere_visual.uniform.transform.scaling = 0.03f;
                sphere_visual.uniform.transform.translation = bar;
                draw(sphere_visual, scene.camera);

                vec3 const u_medial = bar - p_joint;
                vec3 const un_medial = normalize(u_medial);

                vec3 const un_angular_speed = normalize(angular_speed);
                vec3 const un_scaling = normalize(cross(un_medial, un_angular_speed));
                vec3 const un_squeeze = cross(un_medial, un_scaling);



                segment_drawer.uniform_parameter.color = {1,0,1};
                segment_drawer.uniform_parameter.p1 = p_joint;
                segment_drawer.uniform_parameter.p2 = p_joint+2*u_medial;
                segment_drawer.draw(shaders["segment_im"], scene.camera);

                segment_drawer.uniform_parameter.color = {0.8,0.8,1};
                segment_drawer.uniform_parameter.p1 = p_joint;
                segment_drawer.uniform_parameter.p2 = p_joint+2*un_scaling;
                segment_drawer.draw(shaders["segment_im"], scene.camera);

                segment_drawer.uniform_parameter.color = {0.2,0.2,0.2};
                segment_drawer.uniform_parameter.p1 = p_joint;
                segment_drawer.uniform_parameter.p2 = p_joint+2*un_squeeze;
                segment_drawer.draw(shaders["segment_im"], scene.camera);


                mat3 S = { 1+0.5f,0,0, 0,1/std::sqrt(1+0.5f),0, 0,0,1.0f };
                mat3 R = mat3(un_scaling, un_squeeze, un_medial);
                mat3 const T = R*S*transpose(R);

                for(int k=0; k<vertex_depending_on_joint[9].size(); ++k)
                {
                    int id = vertex_depending_on_joint[9][k];
                    vec3 const& p_vertex = save_skinning[id];
                    vec3 const p_medial = p_joint + dot(p_vertex-p_joint,un_medial)*un_medial;


                    std::cout<<id<<std::endl;

                    segment_drawer.uniform_parameter.color = {0.0,0.0,0.0};
                    segment_drawer.uniform_parameter.p1 = p_vertex;
                    segment_drawer.uniform_parameter.p2 = p_medial;
                    segment_drawer.draw(shaders["segment_im"], scene.camera);

                    sphere_visual.uniform.transform.translation = T * (p_vertex-p_medial) + p_medial;
                    sphere_visual.uniform.transform.scaling = 0.01f;
                    sphere_visual.uniform.color = {1,0.4,1};
                    draw(sphere_visual, scene.camera);
                }



            }
        }
    }
*/






    segment_drawer.uniform_parameter.color = {0,0,0};

    // Display sphere when hovering over joint ready to be picked
    if(picking.joint_hover!=-1)
        display_sphere_hover(sphere_visual, picking.joint_hover, skeleton_current, scene.camera);

    if(picking.is_selected)
        display_sphere_selected(sphere_visual, picking.p_clicked, picking.p_current, scene.camera);


    if(gui_param.painting.activated)
        if(picking.painting_selected_vertex>-1)
            display_painting_cursor(picking.painting_selected_vertex, skinning, gui_param.painting.radius, gui_param.painting.threshold_percentage, painting_cursor, sphere_visual, scene.camera);



    glClear(GL_DEPTH_BUFFER_BIT);
    if(gui_param.record_anim)
    {
        vec3 e1 = scene.camera.camera_matrix().mat3().row(0);
        vec3 e2 = scene.camera.camera_matrix().mat3().row(1);
        vec3 e3 = scene.camera.camera_matrix().mat3().row(2);

        sphere_visual.uniform.shading = {1.0f, 0.0f, 0.0f};
        sphere_visual.uniform.transform.scaling = 0.02f;
        sphere_visual.uniform.transform.translation = scene.camera.camera_position() - 2.0f*e3 + 0.65f*e2 + 0.75f*e1;

        sphere_visual.uniform.color = {1,0,0};
        draw(sphere_visual, scene.camera);
    }


    timer_measurement["drawing"].toc();

    timer_measurement["full"].toc();




}



void scene_model::diffuse_weight()
{
    if(gui_param.painting.display_weights==1)
        weight_flappy = diffuse_laplacien_weights(weight_flappy, one_ring, 0.1f, 10, this->skinning);
    if(gui_param.painting.display_weights==2)
        weight_squashy = diffuse_laplacien_weights(weight_squashy, one_ring, 0.1f, 10, this->skinning);
    update_painted_color();
}

void scene_model::set_gui()
{
    //HERE I IMPLEMENT THE BUTTON WHICH TRIGGERS THE ANGULAR VELOCITY 
    if (ImGui::CollapsingHeader("Angular Velocity Traker")) {
        bool const startTracker = ImGui::Button("Start");
        if ( startTracker ) {
            angularVelocity("output_angular_velocity.csv", 100);
        }
    }
    if(ImGui::CollapsingHeader("Animation"))
    {
        ImGui::Text("Animation:"); ImGui::SameLine();
        ImGui::Checkbox("Run", &gui_param.animation); ImGui::SameLine();
        bool const stop  = ImGui::Button("Stop"); ImGui::SameLine();
        bool const start = ImGui::Button("Start");
        if(stop)  timer.stop();
        if(start) timer.start();

        ImGui::SliderFloat("Timer",  &timer_skeleton.t, timer_skeleton.t_min, timer_skeleton.t_max, "%.2f s");
        ImGui::SliderFloat("Time scale", &timer_skeleton.scale, 0.005f, 3.0f, "%.3f s");
    }


    if(ImGui::CollapsingHeader("Display Mesh"))
    {
        ImGui::Checkbox("Mesh", &gui_param.display_mesh); ImGui::SameLine();
        ImGui::Checkbox("Wireframe", &gui_param.display_wireframe); ImGui::SameLine();
        ImGui::Checkbox("Texture", &gui_param.display_texture);
    }

    if(ImGui::CollapsingHeader("Sk. display"))
    {
        ImGui::Checkbox("Bones", &gui_param.display_skeleton_bones);  ImGui::SameLine();
        ImGui::Checkbox("Bones pyramid", &gui_param.display_skeleton_pyramid);
        ImGui::Checkbox("Joints", &gui_param.display_skeleton_joints);  ImGui::SameLine();
        ImGui::Checkbox("Frames", &gui_param.display_skeleton_frames);
        ImGui::SliderFloat("Frame Scaling", &gui_param.frame_scaling, 0.1f, 3.0f);
        ImGui::Checkbox("X-Ray", &gui_param.x_ray);

    }

    if(ImGui::CollapsingHeader("Skinning param."))
    {
        ImGui::Checkbox("Rest pose", &gui_param.display_rest_pose);
        ImGui::Checkbox("Interpolate skeleton", &gui_param.interpolate);
        ImGui::Checkbox("Dual quaternion", &gui_param.dual_quaternion);
        if(ImGui::Button("Rest to original pose")){
            skeleton_current = skeleton_rest_pose;
        }
    }





    if(ImGui::CollapsingHeader("Sk. Interaction"))
    {
        ImGui::Checkbox("Symmetry", &gui_param.symmetry);

        ImGui::RadioButton("Translation", &gui_param.type_deformation, 0); ImGui::SameLine();
        ImGui::RadioButton("Rotation", &gui_param.type_deformation, 1);

    }


    if(ImGui::CollapsingHeader("Dynamic Skinning"))
    {
        ImGui::Text("Parameters");

        ImGui::SliderFloat("Drag power", &flapping_power, 0.0f, 2.0f, "%.2f s", 2.0f);

        ImGui::Text("Dynamic type");
        ImGui::RadioButton("Oscillation", &gui_param.dynamic_type, 0); ImGui::SameLine();
        ImGui::RadioButton("Propagation", &gui_param.dynamic_type, 1); ImGui::SameLine();
        ImGui::RadioButton("Propagation-SK", &gui_param.dynamic_type, 2); 
        ImGui::RadioButton("Oscillation & Propagation", &gui_param.dynamic_type, 3);  ImGui::SameLine();
        ImGui::RadioButton("Osc. & Prop. full", &gui_param.dynamic_type, 4); 
        ImGui::RadioButton("Offset oscillation", &gui_param.dynamic_type, 5);
        ImGui::RadioButton("Oscillation Correct and frequency increase", &gui_param.dynamic_type, 6);

        if(ImGui::CollapsingHeader("Noise")){
            ImGui::Checkbox("Perlin Noise", &gui_param.wave_noise);
            ImGui::Checkbox("Pump Noise", &gui_param.pump_noise);
            ImGui::SliderFloat("Noise frequency", &gui_param.wave_noise_frequency, 0.0f, 2.0f);
            ImGui::SliderFloat("Noise threshold", &gui_param.wave_noise_threshold, 0.0f, 1.0f);
            ImGui::SliderFloat("Noise amplitude", &gui_param.wave_noise_max_amplitude, 0.0f, 1.0f);
            if(gui_param.pump_noise){
                ImGui::SliderFloat("Roll Amplitude", &gui_param.wave_pump_noise_roll_amplitude, 0.0f, 2.0f);
            }
        }

        ImGui::Spacing();
        if (gui_param.dynamic_type == 0) {
            ImGui::Checkbox("Speedup filter", &gui_param.oscillation_filter_speedup);
            ImGui::SliderFloat("Oscillation Frequency", &gui_param.oscillation_gui.frequency, 2 * 3.14f * 0.5, 2 * 3.14f * 6);
            ImGui::SliderFloat("Oscillation Damping", &gui_param.oscillation_gui.attenuation, 0.5f, 3.0f);
            ImGui::SliderFloat("Oscillation Magnitude", &gui_param.oscillation_gui.magnitude, 0.0f, 1.5f, "%.4f", 2.0f);

            max_joint_filtered_value::attenuation = gui_param.oscillation_gui.attenuation;
            max_joint_filtered_value::frequency = gui_param.oscillation_gui.frequency;

            ImGui::SliderFloat("Oscillation Frequency Increase Slope (speedup only)", &gui_param.oscillation_gui.frequency_slope, 0.0f, 2*3.14f*1.5f*2);

            max_joint_filtered_value::frequency_slope = gui_param.oscillation_gui.frequency_slope;
        }
        if (gui_param.dynamic_type == 1) {
            ImGui::SliderFloat("Wave compression", &gui_param.wave_compression_magnitude, 0.0f, 20.0f);
        }
        if (gui_param.dynamic_type == 2) {
            ImGui::SliderFloat("Wave compression", &gui_param.wave_compression_magnitude, 0.0f, 20.0f);
        }
        if (gui_param.dynamic_type == 3) {
            ImGui::SliderFloat("Interpolation factor", &gui_param.propgation_oscillation_weight, 0.0f, 1.0f);
            ImGui::SliderFloat("Wave compression", &gui_param.wave_compression_magnitude, 0.0f, 60.0f);
            ImGui::Separator();
            ImGui::Checkbox("Speedup filter", &gui_param.oscillation_filter_speedup);
            ImGui::SliderFloat("Oscillation Frequency", &gui_param.oscillation_gui.frequency, 2 * 3.14f * 0.5, 2 * 3.14f * 6);
            ImGui::SliderFloat("Oscillation Damping", &gui_param.oscillation_gui.attenuation, 0.5f, 3.0f);
            ImGui::SliderFloat("Oscillation Magnitude", &gui_param.oscillation_gui.magnitude, 0.0f, 1.5f, "%.4f", 2.0f);

            max_joint_filtered_value::attenuation = gui_param.oscillation_gui.attenuation;
            max_joint_filtered_value::frequency = gui_param.oscillation_gui.frequency;
        }
        if (gui_param.dynamic_type == 4) {
            ImGui::SliderFloat("Wave compression", &gui_param.wave_compression_magnitude, 0.0f, 60.0f);
            ImGui::Separator();
            ImGui::SliderFloat("Oscillation Frequency", &gui_param.oscillation_gui.frequency, 2 * 3.14f * 0.5, 2 * 3.14f * 6);
            ImGui::SliderFloat("Oscillation Damping", &gui_param.oscillation_gui.attenuation, 0.5f, 3.0f);
            ImGui::SliderFloat("Oscillation Magnitude", &gui_param.oscillation_gui.magnitude, 0.0f, 1.5f, "%.4f", 2.0f);

            max_joint_filtered_value::attenuation = gui_param.oscillation_gui.attenuation;
            max_joint_filtered_value::frequency = gui_param.oscillation_gui.frequency;
        }
        if (gui_param.dynamic_type == 5) {
            ImGui::Checkbox("Speedup filter", &gui_param.oscillation_filter_speedup);
            ImGui::SliderFloat("Oscillation Frequency", &gui_param.oscillation_gui.frequency, 2 * 3.14f * 0.5, 2 * 3.14f * 6);
            ImGui::SliderFloat("Oscillation Damping", &gui_param.oscillation_gui.attenuation, 0.5f, 3.0f);
            ImGui::SliderFloat("Oscillation Magnitude", &gui_param.oscillation_gui.magnitude, 0.0f, 1.5f, "%.4f", 2.0f);

            max_joint_filtered_value::attenuation = gui_param.oscillation_gui.attenuation;
            max_joint_filtered_value::frequency = gui_param.oscillation_gui.frequency;
        }
        if (gui_param.dynamic_type == 6) { // real oscillation + increase frequency over time (use kappa and filter)
            ImGui::Checkbox("Speedup filter", &gui_param.oscillation_filter_speedup);
            //*********************************************************************************************
            ImGui::RadioButton("Speedup filter sinusoidal", &gui_param.typeOscillation, 1 ); 
            ImGui::RadioButton("Speedup filter triangular", &gui_param.typeOscillation, 2);
            ImGui::RadioButton("Speedup filter rectangular", &gui_param.typeOscillation, 3);
            ImGui::RadioButton("Speedup filter pendulum", &gui_param.typeOscillation, 4);
            ImGui::RadioButton("Speedup filter Bip-Bip", &gui_param.typeOscillation, 5);
            //*********************************************************************************************
            ImGui::SliderFloat("Oscillation Frequency Base", &gui_param.oscillation_gui.frequency, 2 * 3.14f * 0.5, 2 * 3.14f * 6);
            ImGui::SliderFloat("Oscillation Frequency Increase Slope", &gui_param.oscillation_gui.frequency_slope, 0.0f, 2*3.14f*1.5f*2);


            ImGui::SliderFloat("Oscillation Damping", &gui_param.oscillation_gui.attenuation, 0.5f, 3.0f);
            ImGui::SliderFloat("Oscillation Magnitude", &gui_param.oscillation_gui.magnitude, 0.0f, 1.5f, "%.4f", 2.0f);


            max_joint_filtered_value::attenuation = gui_param.oscillation_gui.attenuation;
            max_joint_filtered_value::frequency = gui_param.oscillation_gui.frequency;
        }

        
        ImGui::SliderFloat("Flapping max angle", &gui_param.flappy_max_angle, 0, 3.14f);

        if (ImGui::CollapsingHeader("Squash (deactivated)")) {
            ImGui::SliderFloat("Squash power", &squashing_power, 0.0f, 2.0f, "%.2f s", 2.0f);

            ImGui::Text("Squash around:"); ImGui::SameLine();
            ImGui::RadioButton("Axis", &gui_param.squash_around, 0); ImGui::SameLine();
            ImGui::RadioButton("Center", &gui_param.squash_around, 1);
            ImGui::Checkbox("Center of mass", &gui_param.display_center_of_mass);
        }



        ImGui::Spacing();
        ImGui::Spacing();
        ImGui::Text("Debug:");
        //ImGui::Checkbox("Vertex to bone correspondance", &gui_param.display_vertex_to_bone_correspondance);

        ImGui::Checkbox("Fake Speed", &gui_param.fake_speed);
        ImGui::Text("Display deformation:"); ImGui::SameLine();
        ImGui::Checkbox("Surface", &gui_param.display_deformed_surface); ImGui::SameLine();
        ImGui::Checkbox("Arrows",&gui_param.display_deformation_arrows); ImGui::SameLine();
        ImGui::Checkbox("Target",&gui_param.display_deformation_target);
        ImGui::SliderInt("Vertex Offset", &gui_param.display_deformation_arrow_vertex_offset, 1,15);
        ImGui::Checkbox("Curved trajectory", &gui_param.curved_trajectory);

        ImGui::Text("Display joint velocity:"); ImGui::SameLine();
        ImGui::Checkbox("Linear", &gui_param.display_joint_linear_speed); ImGui::SameLine();
        ImGui::Checkbox("Angular", &gui_param.display_joint_angular_speed);



        /*
        ImGui::Checkbox("Flapping", &is_flapping);
        ImGui::Checkbox("Basic Flapping", &basic_flapping);
        ImGui::Checkbox("Cylinder Flapping", &cylinder_flapping);
        ImGui::Checkbox("Display angular speed", &display_angular_speed);
        ImGui::Checkbox("Is rotating", &is_rotating);
        */
    }






    /*
    if(ImGui::Button("Non negative weights"))
    {
        for(int k=0; k<skinning.deformed.position.size(); ++k)
            weight_flappy[k] = std::abs(weight_flappy[k]);
    }

    if( ImGui::Button("Color white") ) {
        for(int k=0; k<skinning.deformed.position.size(); ++k) {
            skinning.deformed.color[k] = {1.0f, 1.0f, 1.0f, 0};
            character_visual.update_color(skinning.deformed.color);
        }
    }
    ImGui::SameLine();
    if( ImGui::Button("Color flappy") ) {
        for(int k=0; k<skinning.deformed.position.size(); ++k) {
            skinning.deformed.color[k] = {1-std::max(weight_flappy[k],0.0f),1+std::min(weight_flappy[k],0.0f),1,0};
            character_visual.update_color(skinning.deformed.color);
        }
    }
    ImGui::SameLine();
    if( ImGui::Button("Color squashy") ) {
        for(int k=0; k<skinning.deformed.position.size(); ++k) {
            skinning.deformed.color[k] = {std::max(weight_squashy[k],0.0f),-std::min(weight_squashy[k],0.0f),0,0};
            character_visual.update_color(skinning.deformed.color);
        }
    }
    */






    if(ImGui::CollapsingHeader("Painting")){
        bool activate = ImGui::Checkbox("Activate Painting Mode", &gui_param.painting.activated);
        if(activate) {
            if(gui_param.painting.activated==true)
                gui_param.painting.display_weights = 1;
            else
                gui_param.painting.display_weights = 0;
            update_painted_color();
        }

        ImGui::Text("Display color");
        bool color_none = ImGui::RadioButton("None", &gui_param.painting.display_weights, 0); ImGui::SameLine();
        bool color_flappy = ImGui::RadioButton("Flappy", &gui_param.painting.display_weights, 1);
        //ImGui::RadioButton("Squashy", &gui_param.painting.display_weights, 2);
        if(color_none || color_flappy)
            update_painted_color();



        ImGui::Text("Value to paint:");
        ImGui::SliderFloat("Value", &gui_param.painting.value, -1.0f, 1.0f); ImGui::SameLine();
        bool zero = ImGui::Button("Zero");
        bool fill = ImGui::Button("Fill");
        if(zero) gui_param.painting.value = 0.0f;
        if(fill) { weight_flappy.fill(gui_param.painting.value); update_painted_color(); }
        ImGui::SliderFloat("Radius", &gui_param.painting.radius, 0.01f, 0.2f);
        ImGui::SliderFloat("Threshold", &gui_param.painting.threshold_percentage, 0.0f, 1.0f);


        bool diffuse = ImGui::Button("Diffuse");
        if(diffuse) diffuse_weight();
    }



    if(ImGui::CollapsingHeader("Model")){

        bool click_sphere    = ImGui::RadioButton("Sphere", &gui_param.display_type, display_sphere); ImGui::SameLine();
        bool click_rondinella = ImGui::RadioButton("Rondinella", &gui_param.display_type, display_rondinella);

        bool click_cylinder_bending   = ImGui::RadioButton("Cylinder B.", &gui_param.display_type, display_cylinder_bending); ImGui::SameLine();
        bool click_cylinder_translate = ImGui::RadioButton("(Tr.)", &gui_param.display_type, display_cylinder_translate); ImGui::SameLine();
        bool click_bar       = ImGui::RadioButton("Bar", &gui_param.display_type, display_bar);
        bool click_long_cylinder = ImGui::RadioButton("Long Cylinder", &gui_param.display_type, display_long_cylinder);

        bool click_character = ImGui::RadioButton("Character", &gui_param.display_type, display_character);
        bool click_girafe = ImGui::RadioButton("Girafe", &gui_param.display_type, display_girafe);
        bool click_spot = ImGui::RadioButton("Spot", &gui_param.display_type, display_spot);
        bool click_dragon = ImGui::RadioButton("Dragon", &gui_param.display_type, display_dragon);
        bool click_snail = ImGui::RadioButton("Snail", &gui_param.display_type, display_snail);
        bool click_custom = ImGui::RadioButton("Custom", &gui_param.display_type, display_custom);
        bool click_rayfish = ImGui::RadioButton("Ray Fish", &gui_param.display_type, display_rayfish);
        bool click_hammerhead = ImGui::RadioButton("Hammerhead", &gui_param.display_type, display_hammerhead);
        bool click_whale = ImGui::RadioButton("Whale", &gui_param.display_type, display_whale);
        bool click_train = ImGui::RadioButton("Train", &gui_param.display_type, display_train);
        bool click_clown = ImGui::RadioButton("Clown", &gui_param.display_type, display_clown);
        bool click_clownSkHat = ImGui::RadioButton("ClownSkHat", &gui_param.display_type, display_clown_skeleton_hat);

        if(click_sphere)    load_sphere_data(skeleton, skinning, weight_flappy, character_visual, timer_skeleton, shader_mesh);
        if(click_cylinder_bending)  load_bending_twisting_cylinder_data(skeleton, skinning, weight_flappy, weight_squashy, character_visual, timer_skeleton, shader_mesh);
        if(click_cylinder_translate)  load_diagonal_translate_cylinder_data(skeleton, skinning, weight_flappy, character_visual, timer_skeleton, shader_mesh);
        if(click_rondinella)  load_rondinella_data(skeleton, skinning, weight_flappy, weight_squashy, character_visual, timer_skeleton, shader_mesh);

        if(click_bar)     load_rectangle_data(skeleton, skinning, character_visual, timer_skeleton, shader_mesh);

        if(click_character) load_character_data(skeleton, skinning, weight_flappy, character_visual, timer_skeleton, shader_mesh);

        data_loaded data_load;
        if(click_girafe)  data_load = load_girafe2_data(shaders["mesh"]);
        if(click_spot)    data_load = load_spot_data(shaders["mesh"]);
        if(click_long_cylinder) data_load = load_long_cylinder_data(shaders["mesh"]);
        if(click_dragon)  data_load = load_dragon_data(shaders["mesh"]);
        if(click_snail)   data_load = load_snail_data(shaders["mesh"]);
        if(click_custom) data_load = load_custom_data(shaders["mesh"]);
        if(click_rayfish) data_load = load_rayfish_data(shaders["mesh"]);
        if(click_whale) data_load = load_whale_data(shaders["mesh"]);
        if(click_hammerhead) data_load = load_hammerhead_data(shaders["mesh"]);
        if(click_train) data_load = load_train_data(shaders["mesh"]);
        if(click_clown) data_load = load_clown_data(shaders["mesh"]);
        if(click_clownSkHat) data_load = load_clown_skeleton_hat_data(shaders["mesh"]);

        if(click_girafe || click_spot || click_long_cylinder || click_dragon || click_snail || click_rayfish || click_whale || click_hammerhead || click_custom || click_train || click_clown || click_clownSkHat)
        {
            skeleton = data_load.skeleton;
            skinning.influence = data_load.skinning_rig;
            skinning.deformed = data_load.shape;
            skinning.rest_pose = data_load.shape.position;
            skinning.rest_pose_normal = data_load.shape.normal;
            symmetrical_joint = data_load.symmetry;
            weight_flappy = data_load.weight_flappy;
            character_visual = data_load.shape;
            character_visual.shader = data_load.shader;
            character_visual.texture_id = data_load.texture_id;
            timer_skeleton.t_max = data_load.anim_time_max;
            real_texture_id = data_load.texture_id;
        }

        if(click_sphere || click_cylinder_bending || click_cylinder_translate || click_rondinella ||
                click_bar || click_character || click_girafe || click_spot || click_long_cylinder ||  click_dragon || click_whale || click_snail || click_hammerhead || click_rayfish || click_custom || click_train || click_clown || click_clownSkHat) {
            resize_structure();
        }
    }

    if(ImGui::CollapsingHeader("Timings")){

#ifndef _WIN32
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wformat-security"
#pragma GCC diagnostic ignored "-Wformat"
#endif

#define S 64

        static float total = 0;
        static timer_event timer_update_timing;
        static std::map<float,std::string> ordered_time; // reversed timing ordered per time
        timer_update_timing.periodic_event_time_step = 2;
        timer_update_timing.update();
        if(timer_update_timing.event)
        {
            total = timer_measurement["full"].t;
            ordered_time.clear();
            for(auto const& timing: timer_measurement)
                ordered_time[timing.second.t] = timing.first;
        }



        char buffer[S];
        snprintf(buffer, S, "Total : %.1fs", double(total));
        ImGui::Text( buffer );
        for(auto it=ordered_time.rbegin(); it!=ordered_time.rend(); ++it)
        {
            std::string const& name = it->second;
            float const t = it->first;
            float const avg = timer_measurement[name].average_timing;
            if(name!="full")
            {
                char buffer_timing[S];
                snprintf(buffer_timing, S, "[%d%] [%.1fms] %s - %.1fs  ", int(t/total*100), double(1000*avg), name.c_str(), double(t));
                ImGui::Text( buffer_timing );
            }
        }

        bool reset = ImGui::Button("Reset timer");
        if(reset)
        {
            timer_measurement.clear();
            ordered_time.clear();
            total = 0;
        }




#undef S
#ifndef _WIN32
#pragma GCC diagnostic pop
#endif

    }



 

    // Currently deactivated
    //if (ImGui::CollapsingHeader("Wave"))
    //{
    //    ImGui::SliderFloat("Period", &gui_param.wave_gui.T, 1.0f, 10.0f);
    //    ImGui::SliderFloat("Duration", &gui_param.wave_gui.dt, 0.5f, 8.0f);
    //    ImGui::SliderFloat("Amplitude", &gui_param.wave_gui.A, 0.001f, 2.0f, "%.4f", 5.0f);
    //    ImGui::SliderFloat("Speed", &gui_param.wave_gui.speed, 1.0f, 5.0f);
    //}

    //    {
    //        static int counter = 0;
    //        if(counter%500==0)
    //        {
    //            std::cout<<"\nTimings : "<<std::endl;
    //            for(auto const& elapsed : timer_measurement)
    //            {
    //                std::cout<<elapsed.first<<" : "<<elapsed.second.t<<" "<<int(elapsed.second.t/timer_measurement["full"].t*100)<<"%" <<std::endl;
    //            }
    //            std::cout<<std::endl;
    //        }
    //        counter ++;
    //    }

}


void scene_model::resize_structure()
{
    picking.selected_joint_true = -1;
    picking.selected_joint = -1;
    picking.is_selected = false;

    skeleton_local_current_before_deformation  = interpolate_skeleton_at_time(0, skeleton.anim, gui_param.interpolate);
    skeleton_local_current = skeleton_local_current_before_deformation;
    int const N_joint = skeleton_local_current.size();
    int const N_vertex = skinning.rest_pose.size();
    skeleton_local_interactive_deformation.resize(N_joint);


    skeleton_current = local_to_global(skeleton_local_current, skeleton.connectivity);
    //skeleton_speed.resize(skeleton_current.size());


    //skeleton_acceleration.resize(skeleton_current.size());
    //skeleton_velocity_tracker.resize(skeleton_current.size());
    skeleton_rest_pose = local_to_global(skeleton.rest_pose, skeleton.connectivity);
    //skeleton_angular_velocity_tracker.resize(skeleton_current.size());


    //previous_angular_velocity.resize(skinning.rest_pose.size());
    //previous_speed.resize(skinning.rest_pose.size());

    if(weight_squashy.size()!=skinning.rest_pose.size()){
        weight_squashy.resize(skinning.rest_pose.size());
        weight_squashy.fill(1.0f);
    }
    if(weight_flappy.size()!=skinning.rest_pose.size()){
        weight_flappy.resize(skinning.rest_pose.size());
        weight_flappy.fill(0.0f);
    }

    //vertex_speed.clear(); vertex_speed.resize(N_vertex);
    skeleton_joint_speed.clear(); skeleton_joint_speed.resize(N_joint);
    skeleton_joint_rotation_speed.resize_clear(skeleton_current.size());
    skeleton_speed_per_joint.resize_clear(N_joint);
    skeleton_fake_speed_per_joint.resize_clear(N_joint);
    skeleton_previous.data.clear();

    for(int k=0; k<N_joint; ++k){
        for(int kt=0; kt<500; ++kt){
            skeleton_joint_speed[k].add({0,0,0},kt*0.1);
            skeleton_joint_rotation_speed[k].add({0,0,0,1},kt*0.1);
        }
    }

    // build one ring
    {
        size_t const N = skinning.deformed.position.size();
        auto const& connectivity = skinning.deformed.connectivity;
        one_ring.clear();
        one_ring.resize(N);

        size_t const N_tri = connectivity.size();
        for(size_t k_tri=0; k_tri<N_tri; ++k_tri)
        {
            unsigned int const u0 = connectivity[k_tri][0];
            unsigned int const u1 = connectivity[k_tri][1];
            unsigned int const u2 = connectivity[k_tri][2];

            one_ring[u0].insert(u1); one_ring[u0].insert(u2);
            one_ring[u1].insert(u0); one_ring[u1].insert(u2);
            one_ring[u2].insert(u0); one_ring[u2].insert(u1);
        }
    }

    // Build triangle around vertex
    {
        size_t const N_triangle = skinning.deformed.connectivity.size();
        triangle_around_vertex.resize_clear(N_vertex);
        for(size_t kt=0; kt<N_triangle; ++kt)
        {
            uint3 const& tri = skinning.deformed.connectivity[kt];
            unsigned int a=tri[0], b=tri[1], c=tri[2];
            triangle_around_vertex[a].push_back(kt);
            triangle_around_vertex[b].push_back(kt);
            triangle_around_vertex[c].push_back(kt);
        }
    }

    update_painted_color();

    skeleton_son_connectivity = compute_joint_sons(skeleton.connectivity);
    vertex_to_bone_correspondance = compute_bone_correspondance(skinning.rest_pose, skinning.influence, skeleton_son_connectivity, skeleton_rest_pose);


    // Rig extended to ancestor
    {
        rig_extended_to_ancestor_joint.resize_clear(N_vertex);
        rig_extended_to_ancestor_weight.resize_clear(N_vertex);

        for(int kv=0; kv<N_vertex; ++kv)
        {
            auto const& influence = skinning.influence[kv];
            int const N_influence = influence.size();

            std::map<int,float> cumulative_weight_per_joint;

            for(int k_influence=0; k_influence<N_influence; ++k_influence)
            {
                int current_joint = influence[k_influence].joint;
                float weight = influence[k_influence].weight;

                while(current_joint != -1)
                {

                    if(cumulative_weight_per_joint.find(current_joint)!=cumulative_weight_per_joint.end())
                        cumulative_weight_per_joint[current_joint] += weight;
                    else
                        cumulative_weight_per_joint[current_joint]  = weight;

                    current_joint = skeleton.connectivity[current_joint].parent;
                }
            }

            for(auto const& it : cumulative_weight_per_joint)
            {
                rig_extended_to_ancestor_joint[kv].push_back(it.first);
                rig_extended_to_ancestor_weight[kv].push_back(it.second);
            }

        }

        //        for(int kv=0; kv<N_vertex; ++kv)
        //        {
        //            for(int k=0; k<rig_extended_to_ancestor_joint[kv].size(); ++k)
        //            {
        //                std::cout<<kv<<" - "<<rig_extended_to_ancestor_joint[kv][k]<<","<<rig_extended_to_ancestor_weight[kv][k]<<std::endl;
        //            }
        //            std::cout<<std::endl;
        //        }
    }

    vertex_depending_on_joint.resize_clear(N_joint);
    vertex_weight_depending_on_joint.resize_clear(N_joint);

    for(int k_vertex=0; k_vertex<N_vertex; ++k_vertex)
    {
        buffer<int> const& joint_dependencies = rig_extended_to_ancestor_joint[k_vertex];
        for(int k_joint=0; k_joint<int(joint_dependencies.size()); ++k_joint)
        {
            int joint = rig_extended_to_ancestor_joint[k_vertex][k_joint];
            float weight = rig_extended_to_ancestor_weight[k_vertex][k_joint];
            vertex_depending_on_joint[joint].push_back(k_vertex);
            vertex_weight_depending_on_joint[joint].push_back(weight);
        }
    }

    skinning_weights_per_joint_per_vertex.resize_clear(N_joint);
    for(int kj=0; kj<N_joint; ++kj)
    {
        skinning_weights_per_joint_per_vertex[kj].resize_clear(N_vertex);
        for(size_t k=0; k<vertex_depending_on_joint[kj].size(); ++k)
        {
            int vertex_id = vertex_depending_on_joint[kj][k];
            float w = vertex_weight_depending_on_joint[kj][k];
            skinning_weights_per_joint_per_vertex[kj][vertex_id] = w;
        }

    }

    angular_velocity_filtered_storage.resize_clear(N_joint);
    angular_velocity_filtered_storage_subjoints.resize_clear(N_joint);
    for (int k = 0; k < N_joint; k++)
        angular_velocity_filtered_storage_subjoints[k].resize_clear(N_joint);



    vcl::buffer<std::set<int> > triangle_set_depending_on_joint(N_joint);
    int const N_triangle = skinning.deformed.connectivity.size();
    for(int k_triangle=0; k_triangle<N_triangle; ++k_triangle)
    {
        uint3 const& tri = skinning.deformed.connectivity[k_triangle];
        unsigned int a = tri[0];
        unsigned int b = tri[1];
        unsigned int c = tri[2];

        buffer<int> const& joints_a = rig_extended_to_ancestor_joint[a];
        buffer<int> const& joints_b = rig_extended_to_ancestor_joint[b];
        buffer<int> const& joints_c = rig_extended_to_ancestor_joint[c];

        for(int j : joints_a)
            triangle_set_depending_on_joint[j].insert(k_triangle);
        for(int j : joints_b)
            triangle_set_depending_on_joint[j].insert(k_triangle);
        for(int j : joints_c)
            triangle_set_depending_on_joint[j].insert(k_triangle);
    }

    triangle_depending_on_joint.resize_clear(N_joint);
    for(int kj=0; kj<N_joint; ++kj)
    {
        for(int t : triangle_set_depending_on_joint[kj])
            triangle_depending_on_joint[kj].push_back(t);
    }


    triangle_area.resize_clear(N_triangle);
    triangle_center.resize_clear(N_triangle);


    center_of_mass_per_joint.resize_clear(N_joint);
    center_of_mass_per_joint_manual_offset.resize_clear(N_joint);

    timer_skeleton.t = timer_skeleton.t_min;

    picking.is_selected = false;
    picking.selected_joint = -1;
    picking.selected_joint_true = -1;
    recorded_joint = -1;

    deformation_per_vertex.resize_clear(N_vertex);


    update_closest_distance_to_bone();
}

vec3 closest_point_to_segment(vec3 const& p, vec3 const& a, vec3 const& b)
{
    vec3 const ap = p - a;
    vec3 const ab = b - a;
    vec3 const u_ab = normalize(ab);

    float proj = dot(ap, u_ab);
    vec3 closest = a + proj * u_ab;

    if (proj <= 0)
        closest = a;
    else if (proj >= norm(ab))
        closest = b;
    
    return closest;
}

void scene_model::update_closest_distance_to_bone()
{
    int const N_joint = skeleton_local_current.size();
    int const N_vertex = skinning.rest_pose.size();

    closest_distance_to_bone.p.resize(N_vertex);
    closest_distance_to_bone.joint.resize(N_vertex);

    for (int k = 0; k < N_vertex; ++k)
    {
        vec3 const& p = skinning.rest_pose[k];

        vec3 closest_point;
        int closest_joint;
        float min_distance = 0.0f;

        int N_dependence = skinning.influence[k].size();
        for (int kj = 0; kj < N_dependence; ++kj)
        {
            int joint = skinning.influence[k][kj].joint;
            int N_son = skeleton_son_connectivity[joint].size();
            for (int k_son = 0; k_son < N_son; ++k_son)
            {
                int joint_son = skeleton_son_connectivity[joint][k_son];

                vec3 const& a = skeleton_current[joint].p;
                vec3 const& b = skeleton_current[joint_son].p;
                vec3 const current_closest = closest_point_to_segment(p, a, b);

                if ((kj == 0 && k_son == 0) || min_distance > norm(current_closest - p)) {
                    closest_point = current_closest;
                    closest_joint = joint;
                    min_distance = norm(current_closest - p);
                }
            }
        }

        closest_distance_to_bone.p[k] = closest_point;
        closest_distance_to_bone.joint[k] = closest_joint;
    }

    distance_to_joint.resize(N_vertex);
    for (int k = 0; k < N_vertex; ++k) {
        distance_to_joint[k].resize(N_joint);

        //set all to cartesian distance first
        for (int kj = 0; kj < skeleton_current.size(); ++kj)
        {
            distance_to_joint[k][kj] = norm(skeleton_current[kj].p- skinning.rest_pose[k]);
        }        
    }


    for (int k = 0; k < N_vertex; ++k)
    {
        int N_dependance = rig_extended_to_ancestor_joint[k].size();
        int joint_closest = 0; //closest_distance_to_bone.joint[k];
        vec3 const& closest_point = closest_distance_to_bone.p[k];


        float L = norm(closest_point - skinning.rest_pose[k]);
        L += norm(closest_point- skeleton_current[joint_closest].p);
        distance_to_joint[k][joint_closest] = L;



        int joint_current = joint_closest;
        while (joint_current >= 0) {
            int joint_parent = skeleton.connectivity[joint_current].parent;
            if (joint_parent >= 0) {

                vec3 const& pj = skeleton_current[joint_current].p;
                vec3 const& ppj = skeleton_current[joint_parent].p;
                L += norm(pj - ppj);
                distance_to_joint[k][joint_parent] = L;
            }
            joint_current = joint_parent;

            if (k == 51) {
                std::cout << " join : " << joint_current << std::endl;
            }
        }
    }

    //exit(0);
}

void scene_model::setup_data(std::map<std::string,GLuint>& _shaders, scene_structure& _scene, gui_structure& gui)
{   //AICI POT SA APELEZ FUNCTIA
    convolution( "output_kappa.csv", 100 );
    //angularVelocity("output_angular_velocity.csv", 100);

    scene = _scene;
    shaders = _shaders;
    shader_mesh = shaders["mesh_bf"];

    shaders["depth_map"] = create_shader_program("scenes/shared_assets/shaders/depth_map/shader.vert.glsl",
                                                 "scenes/shared_assets/shaders/depth_map/shader.frag.glsl");

    shaders["mesh_shadow"] = create_shader_program("scenes/shared_assets/shaders/mesh_shadow/shader.vert.glsl",
                                                   "scenes/shared_assets/shaders/mesh_shadow/shader.frag.glsl");

    //_scene.texture_white = create_texture_gpu(image_load_png("assets/spot/texture.png"));


    glUseProgram(shaders["mesh_shadow"]);
    GLint texLoc = glGetUniformLocation(shaders["mesh_shadow"], "texture_sampler");
    glUniform1i(texLoc, 0);
    texLoc = glGetUniformLocation(shaders["mesh_shadow"], "shadowMap");
    glUniform1i(texLoc, 1);
    glUseProgram(0);


    segment_drawer.init();
    segment_drawer.uniform_parameter.color = {0.0f,0.0f,0.0f};
    glEnable(GL_POLYGON_OFFSET_FILL);


    // Sphere used to display joints
    sphere_visual = mesh_primitive_sphere(1.0f);
    sphere_visual.shader = shader_mesh;

    frame = mesh_primitive_frame();
    frame.uniform.transform.scaling = 0.02f;
    frame.shader = shaders["mesh"];

    // Load initial model
    //load_sphere_data(skeleton, skinning, weight_flappy, character_visual, timer_skeleton, shader_mesh);
    load_bending_cylinder_data(skeleton, skinning, weight_flappy, weight_squashy, character_visual, timer_skeleton, shader_mesh);
    gui_param.display_type = display_cylinder_bending;

    // Setup cursor
    painting_cursor = curve_primitve_circle(40, 1.0f, {0,0,0}, {0,0,1});
    painting_cursor.shader = shaders["curve"];
    painting_cursor.uniform.color = {0.7f, 0.7f, 0.7f};
    painting_cursor.uniform.transform.scaling = gui_param.painting.radius;


    pyramid_skeleton_visual = mesh_primitive_pyramid(0.25f,1.0f);
    pyramid_skeleton_visual.shader = shaders["mesh"];

    gui.show_frame_camera = false;

    quad = mesh_primitive_quad();

    resize_structure();


    {

        glGenFramebuffers(1, &depthMapFBO); opengl_debug();
        glBindFramebuffer(GL_FRAMEBUFFER, depthMapFBO);  opengl_debug();
        glGenTextures(1, &depthMap); opengl_debug();
        glBindTexture(GL_TEXTURE_2D, depthMap); opengl_debug();
        //glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, 1024, 1024, 0, GL_RGB, GL_UNSIGNED_BYTE, NULL); opengl_debug();
        glTexImage2D(GL_TEXTURE_2D, 0, GL_DEPTH_COMPONENT, SHADOW_WIDTH, SHADOW_HEIGHT, 0, GL_DEPTH_COMPONENT, GL_FLOAT, NULL); opengl_debug();
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST); opengl_debug();
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);  opengl_debug();
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);  opengl_debug();
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);  opengl_debug();
        float borderColor[] = { 1.0, 1.0, 1.0, 1.0 };
        glTexParameterfv(GL_TEXTURE_2D, GL_TEXTURE_BORDER_COLOR, borderColor);

        glBindFramebuffer(GL_FRAMEBUFFER, depthMapFBO);  opengl_debug();
        glFramebufferTexture2D(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_TEXTURE_2D, depthMap, 0); opengl_debug();
        //glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, depthMap, 0); opengl_debug();


        //        glGenRenderbuffers(1, &rbo); opengl_debug();
        //        glBindRenderbuffer(GL_RENDERBUFFER, rbo); opengl_debug();
        //        glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH24_STENCIL8, 1024, 1024); opengl_debug();
        //        glBindRenderbuffer(GL_RENDERBUFFER, 0); opengl_debug();
        //        glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_STENCIL_ATTACHMENT, GL_RENDERBUFFER, rbo); opengl_debug();
        //        if(glCheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE) opengl_debug();
        //            std::cout << "ERROR::FRAMEBUFFER:: Framebuffer is not complete!" << std::endl;
        //        glBindFramebuffer(GL_FRAMEBUFFER, 0); opengl_debug();

        //        glBindFramebuffer(GL_FRAMEBUFFER, depthMapFBO);
        //        glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, depthMap, 0);
        glDrawBuffer(GL_NONE);
        glReadBuffer(GL_NONE);
        glBindFramebuffer(GL_FRAMEBUFFER, 0);
    }


    timer_recording.stop();
    record_dt = 0.1f;

    arrow.init(shaders["mesh"]);
    arrow_line.init(shaders["mesh"]);



    // environment_train["cactus1"] = mesh_drawable( mesh_load_file_obj("assets/cactus/cactus.obj") );
    // environment_train["cactus1"].texture_id = scene.texture_green;
    // environment_train["cactus1"].uniform.transform.scaling = 3.0;
    // environment_train["cactus1"].uniform.transform.translation = { -10,0,-1 };
    // environment_train["cactus2"] = mesh_drawable( mesh_load_file_obj("assets/cactus/cactus.obj") );
    // environment_train["cactus2"].texture_id = scene.texture_green;
    // environment_train["cactus2"].uniform.transform.scaling = 3.0;
    // environment_train["cactus2"].uniform.transform.translation = { -20,0,-1 };
    environment_train["train_pillars"] = mesh_drawable(mesh_load_file_obj("assets/train_track/train_pillars.obj"));
    environment_train["train_pillars"].texture_id = scene.texture_pillar;
    environment_train["train_pillars"].uniform.transform.scaling = 0.0025;
    environment_train["train_pillars"].uniform.transform.translation = { -10,-0.12f ,0 };
    environment_train["train_rails"] = mesh_drawable(mesh_load_file_obj("assets/train_track/train_rails.obj"));
    environment_train["train_rails"].texture_id = scene.texture_rails;
    environment_train["train_rails"].uniform.transform.scaling = 0.0025;
    environment_train["train_rails"].uniform.transform.translation = { -10,-0.12f ,0 };
    environment_train["train_rails"].uniform.transform.scaling_axis.x += 3.0f;
    environment_train["train_bars"] = mesh_drawable(mesh_load_file_obj("assets/train_track/train_bars.obj"));
    environment_train["train_bars"].texture_id = scene.texture_bars;
    environment_train["train_bars"].uniform.transform.scaling = 0.0025;
    environment_train["train_bars"].uniform.transform.translation = { -10,-0.12f ,0 };


    environment_train["train_pillars2"] = mesh_drawable(mesh_load_file_obj("assets/train_track/train_pillars.obj"));
    environment_train["train_pillars2"].texture_id = scene.texture_pillar;
    environment_train["train_pillars2"].uniform.transform.scaling = 0.0025;
    environment_train["train_pillars2"].uniform.transform.translation = { 15.0f,-0.12f ,0 };
    // environment_train["train_rails2"] = mesh_drawable(mesh_load_file_obj("assets/train_track/train_rails.obj"));
    // environment_train["train_rails2"].texture_id = scene.texture_rails;
    // environment_train["train_rails2"].uniform.transform.scaling = 0.0025;
    // environment_train["train_rails2"].uniform.transform.translation = { 15.0f,-0.12f ,0 };
    environment_train["train_bars2"] = mesh_drawable(mesh_load_file_obj("assets/train_track/train_bars.obj"));
    environment_train["train_bars2"].texture_id = scene.texture_bars;
    environment_train["train_bars2"].uniform.transform.scaling = 0.0025;
    environment_train["train_bars2"].uniform.transform.translation = { 15.0f,-0.12f ,0 };

}


void scene_model::update_center_of_mass()
{
    // Damien Hack: I add this hack_damien_variable to avoid recomputing the center of mass
    // This destroys the squashy deformation !
    // If we want skashy deformation, we need to re-introduce it and remove this hack.
    static bool hack_damien_variable = false;
    if (hack_damien_variable == false) {
        hack_damien_variable = true;


        size_t const N_triangle = skinning.deformed.connectivity.size();
        for (size_t k_triangle = 0; k_triangle < N_triangle; ++k_triangle)
        {
            uint3 const& tri = skinning.deformed.connectivity[k_triangle];
            unsigned int const u0 = tri[0];
            unsigned int const u1 = tri[1];
            unsigned int const u2 = tri[2];

            vec3 const& p0 = save_skinning[u0];
            vec3 const& p1 = save_skinning[u1];
            vec3 const& p2 = save_skinning[u2];

            triangle_area[k_triangle] = 0.5f * norm(cross(p1 - p0, p2 - p0));
            triangle_center[k_triangle] = (p0 + p1 + p2) / 3.0f;
        }

        size_t const N_joint = skinning_weights_per_joint_per_vertex.size();
        for (size_t k_joint = 0; k_joint < N_joint; ++k_joint)
        {
            size_t N_triangle_dependency = triangle_depending_on_joint.data[k_joint].size();
            vec3 bar = { 0,0,0 };
            float counter = 0.0f;
            for (size_t k_triangle = 0; k_triangle < N_triangle_dependency; ++k_triangle)
            {
                int k_tri = triangle_depending_on_joint.data[k_joint][k_triangle];
                uint3 const& tri = skinning.deformed.connectivity.data[k_tri];

                unsigned int const u0 = tri[0];
                unsigned int const u1 = tri[1];
                unsigned int const u2 = tri[2];

                float const w0 = skinning_weights_per_joint_per_vertex.data[k_joint].data[u0];
                float const w1 = skinning_weights_per_joint_per_vertex.data[k_joint].data[u1];
                float const w2 = skinning_weights_per_joint_per_vertex.data[k_joint].data[u2];

                bar += (w0 + w1 + w2) * triangle_area.data[k_tri] * triangle_center.data[k_tri];
                counter += (w0 + w1 + w2) * triangle_area.data[k_tri];
            }
            if (counter > 1e-6f)
                bar /= counter;
            center_of_mass_per_joint.data[k_joint] = bar;
        }
    }


}

void scene_model::generate_fake_speed()
{
    if(picking.is_selected)
    {

        bool is_translation = (gui_param.type_deformation==0 || picking.selected_joint==0);

        if(is_translation){
            int const joint = picking.selected_joint;
            vec3 const translation = picking.p_current - picking.p_clicked;
            skeleton_fake_speed_per_joint[joint].linear_speed = 2*translation;
        }
        else
        {
            int const joint = skeleton.connectivity[picking.selected_joint].parent;
            vec3 const& p_joint = skeleton_current[joint].p;
            vec3 const& p0 = picking.p_clicked;
            vec3 const& p1 = picking.p_current;

            if(picking.click_button==left) // rotate in plane
            {
                vec3 const u_ref = normalize(p0-p_joint);
                vec3 const u_objective = normalize(p1-p_joint);
                vec3 const axis = normalize(cross(u_ref, u_objective));
                float const angle = std::acos(dot(u_ref,u_objective));
                //quaternion const q = quaternion::axis_angle(axis,angle);
                skeleton_fake_speed_per_joint[joint].angular_speed = 2*axis*angle;
            }
            if(picking.click_button==right) // twist
            {
                vec3 const axis = normalize(p0 - p_joint);
                vec2 const sj = scene.camera.project_position(p_joint);
                vec2 const s0 = scene.camera.project_position(p0);
                vec2 const s1 = scene.camera.project_position(p1);

                vec2 const a = normalize(s0-sj);
                vec2 const b = s1-sj;
                float const angle = -2*3.14f*det(a,b);

                skeleton_fake_speed_per_joint[joint].angular_speed = 2*axis*angle;
            }

        }

    }
}

void scene_model::adapt_skeleton_interactive()
{

    if(picking.is_selected)
    {

        int const parent = skeleton.connectivity[picking.selected_joint].parent;

        quaternion q_parent;
        vec3 p_parent = skeleton_current[0].p;

        bool is_translation = (gui_param.type_deformation==0 || picking.selected_joint==0);

        // Displacement is applied to the current joint
        if(is_translation){
            if(parent>-1) {
                q_parent = skeleton_current[ parent ].r;
            }
            vec3 const translation = picking.p_current - picking.p_clicked;
            skeleton_local_interactive_deformation[picking.selected_joint].p = conjugate(q_parent).apply(translation);
        }

        // Rotation is applied to the parent joint
        if(!is_translation)
        {
            if(parent>-1) {

                p_parent = skeleton_current[ parent ].p;

                int const grand_parent = skeleton.connectivity[parent].parent;
                quaternion q_grand_parent = {0,0,0,1};
                if(grand_parent>-1)
                    q_grand_parent = skeleton_current[grand_parent].r;

                quaternion q;
                if(picking.click_button==left) // rotate in plane
                {
                    vec3 const u_ref = normalize(picking.p_clicked - p_parent);
                    vec3 const u_objective = normalize(picking.p_current - p_parent);

                    const float angle = std::acos(dot(u_ref,u_objective));
                    const vec3 axis = normalize(cross(u_ref, u_objective));
                    q = quaternion::axis_angle(axis,angle);
                }
                if(picking.click_button==right) // twist
                {
                    vec3 const axis = normalize(picking.p_clicked - p_parent);
                    vec2 const p0 = scene.camera.project_position(p_parent);
                    vec2 const p1 = scene.camera.project_position(picking.p_clicked);
                    vec2 const p  = scene.camera.project_position(picking.p_current);

                    vec2 const a = normalize(p1-p0);
                    vec2 const b = p-p0;
                    float const angle = -2*3.14f*det(a,b);
                    q = quaternion::axis_angle(axis,angle);
                }


                skeleton_local_interactive_deformation[parent].r =  conjugate(q_grand_parent) * q * q_grand_parent;;
            }
        }


    }

}



void scene_model::apply_deformation_on_skeleton()
{
    int const parent = skeleton.connectivity[picking.selected_joint].parent;




    // apply deformation on all time steps
    int const N_time = skeleton.anim[picking.selected_joint].size();
    for(int k_time=0; k_time<N_time; ++k_time)
    {
        skeleton.anim[picking.selected_joint][k_time].geometry.p += skeleton_local_interactive_deformation[picking.selected_joint].p;

        if(parent>-1){
            skeleton.anim[parent][k_time].geometry.r = skeleton_local_interactive_deformation[parent].r * skeleton.anim[parent][k_time].geometry.r;

            if(gui_param.symmetry)
            {
                if(symmetrical_joint.find(parent)!=symmetrical_joint.end())
                {
                    int const js = symmetrical_joint[parent];
                    skeleton.anim[js][k_time].geometry.r = skeleton_local_interactive_deformation[js].r * skeleton.anim[js][k_time].geometry.r;
                }
            }
        }
    }



    for(int k=0; k<int(skeleton_local_interactive_deformation.size()); ++k)
    {
        skeleton_local_interactive_deformation[k].p = {0,0,0};
        skeleton_local_interactive_deformation[k].r = {0,0,0,1};
    }

    //    skeleton_local_interactive_deformation[picking.selected_joint].p = {0,0,0};
    //    if(parent>-1)
    //        skeleton_local_interactive_deformation[parent].r = {0,0,0,1};
}

template <typename T, typename F>
int run_picking_selection_generic(scene_structure& scene, GLFWwindow* window, T const& buffer, F const& function_picker)
{
    int index_selected = -1;

    // Create the 3D ray passing by the selected point on the screen
    const vec2 cursor = glfw_cursor_coordinates_window(window);
    const ray r = picking_ray(scene.camera, cursor);

    bool is_selected = false;
    float distance_min = 0.0f;
    const int N = buffer.size();
    for(int k=0; k<N; ++k)
    {
        const vec3& c = function_picker(buffer, k);
        const picking_info info = ray_intersect_sphere(r, c, 0.04f);

        if( info.picking_valid ) // the ray intersects a sphere
        {
            const float distance = norm(info.intersection-r.p); // get the closest intersection
            if( is_selected==false || distance<distance_min ){
                is_selected = true;
                distance_min = distance;

                index_selected = k;
            }
        }
    }

    return index_selected;
}

int run_picking_selection_point_set(scene_structure& scene, GLFWwindow* window, buffer<vec3> const& point_set)
{

    return run_picking_selection_generic(scene, window, point_set, [](buffer<vec3> const& t,int k){return t[k];} );
}

int run_picking_selection_skeleton(scene_structure& scene, GLFWwindow* window, buffer<joint_geometry> const& skeleton_current)
{
    return run_picking_selection_generic(scene, window, skeleton_current, [](buffer<joint_geometry> const& t,int k){return t[k].p;} );
}


void convert_weight_to_color(buffer<float> const& value, buffer<vec4>& color)
{
    assert_vcl(value.size()==color.size(), "Incorrect size");
    size_t const N = value.size();
    for(size_t k=0; k<N; ++k)
        color[k] = { 1-std::max(value[k],0.0f), 1+std::min(value[k],0.0f), 1, 0 };
}

void scene_model::update_painted_color()
{
    switch(gui_param.painting.display_weights)
    {
    case 0:
        skinning.deformed.color.fill({1,1,1,0});
        break;
    case 1:
        convert_weight_to_color(weight_flappy, skinning.deformed.color);
        break;
    case 2:
        convert_weight_to_color(weight_squashy, skinning.deformed.color);
        break;
    }
    character_visual.update_color(skinning.deformed.color);
}

void scene_model::mouse_scroll(scene_structure& , GLFWwindow* window, float , float y_offset)
{
    // Increase/decrease falloff distance when scrolling the mouse
    if(gui_param.painting.activated)
    {
        const bool key_shift = glfw_key_shift_pressed(window);

        if(key_shift)
            gui_param.painting.radius = std::max(gui_param.painting.radius + y_offset*0.01f, 1e-6f);
        else
            gui_param.painting.threshold_percentage = std::min(std::max(gui_param.painting.threshold_percentage + y_offset*0.01f, 0.0f),1.0f);

    }
}

void scene_model::paint_vertices_around_picked_point()
{
    assert_vcl_no_msg(picking.painting_selected_vertex<int(skinning.deformed.position.size()));
    vec3 const& p_picked = skinning.deformed.position[picking.painting_selected_vertex];

    mesh& m = skinning.deformed;
    buffer<vec3>& vertices = m.position;
    size_t const N = vertices.size();
    float const threshold = gui_param.painting.radius * gui_param.painting.threshold_percentage;
    for(size_t k=0; k<N; ++k)
    {
        vec3 const& p = vertices[k];
        float const d = norm(p-p_picked);


        if( d < gui_param.painting.radius ) {
            float const current = weight_flappy[k];
            float const target = gui_param.painting.value;
            float const r      = gui_param.painting.radius;

            if(d<threshold)
            {
                weight_flappy[k] = target;
            }
            else
            {
                float alpha = (d-threshold)/(r-threshold);
                alpha = std::exp(-alpha*alpha*6);
                weight_flappy[k] = (1-alpha) * current + alpha * target;
            }

        }

    }

    update_painted_color();
}

void scene_model::mouse_move(scene_structure& scene_arg, GLFWwindow* window)
{
    const bool key_shift = glfw_key_shift_pressed(window);
    const bool mouse_click_left = glfw_mouse_pressed_left(window);


    if(gui_param.painting.activated && key_shift)
    {
        picking.painting_selected_vertex = run_picking_selection_point_set(scene_arg, window, skinning.deformed.position);
        if(picking.painting_selected_vertex>-1 && mouse_click_left)
            paint_vertices_around_picked_point();
    }


    if(!gui_param.painting.activated)
    {

        // Hovering
        if(!picking.is_selected && key_shift){
            picking.joint_hover = run_picking_selection_skeleton(scene_arg, window, skeleton_current);
        }


        // Displacement
        if(picking.is_selected || picking.picked_center_of_mass)
        {
            const vec2 cursor = glfw_cursor_coordinates_window(window);

            // Get vector orthogonal to camera orientation
            const mat4 M = scene.camera.camera_matrix();
            const vec3 n = {M(0,2),M(1,2),M(2,2)};

            const ray r = picking_ray(scene.camera, cursor);
            const picking_info info = ray_intersect_plane(r, n, picking.p_clicked);

            assert_vcl_no_msg(info.picking_valid);

            picking.p_current = info.intersection;

            if(picking.picked_center_of_mass)
            {
                int joint = picking.selected_joint_true;
                vec3 const& com = center_of_mass_per_joint[joint];
                quaternion q = skeleton_current[joint].r;

                center_of_mass_per_joint_manual_offset[joint] = conjugate(q).apply(picking.p_current-com);
            }
            else
            {
                if(!gui_param.fake_speed)
                    adapt_skeleton_interactive();
                else
                    generate_fake_speed();
            }


        }



    }



}

vcl::vec3 scene_model::position_center_of_mass(int joint)
{
    vec3 const& com_true = center_of_mass_per_joint[joint];
    vec3 const& offset = center_of_mass_per_joint_manual_offset[joint];
    quaternion const& q = skeleton_current[joint].r;

    return com_true + q.apply(offset);
}

void scene_model::mouse_click(scene_structure& scene_arg, GLFWwindow* window, int button, int action, int )
{

    // Check that the mouse is clicked (drag and drop)
    const bool mouse_click_left    = (button==GLFW_MOUSE_BUTTON_LEFT && action==GLFW_PRESS);
    const bool mouse_release_left  = (button==GLFW_MOUSE_BUTTON_LEFT && action==GLFW_RELEASE);
    const bool mouse_click_right    = (button==GLFW_MOUSE_BUTTON_RIGHT && action==GLFW_PRESS);
    const bool mouse_release_right  = (button==GLFW_MOUSE_BUTTON_RIGHT && action==GLFW_RELEASE);
    const bool key_shift = glfw_key_shift_pressed(window);
    const bool mouse_click   = mouse_click_left || mouse_click_right;
    const bool mouse_release = mouse_release_left || mouse_release_right;

    if(mouse_click_left) picking.click_button = left;
    if(mouse_click_right)picking.click_button = right;
    if(mouse_release)    picking.click_button = none;

    if(gui_param.painting.activated && key_shift)
    {
        picking.painting_selected_vertex = run_picking_selection_point_set(scene_arg, window, skinning.deformed.position);
        if(picking.painting_selected_vertex>-1 && mouse_click_left)
            paint_vertices_around_picked_point();
    }

    // Selection
    if(!gui_param.painting.activated){
        if(key_shift && mouse_click)
        {

            if(picking.selected_joint!=-1 && gui_param.display_center_of_mass)
            {
                int joint = picking.selected_joint_true;
                if(joint!=-1)
                {
                    vec3 const& com = position_center_of_mass(joint);
                    int picked_com = run_picking_selection_point_set(scene_arg, window, {com});

                    if(picked_com==0){
                        picking.picked_center_of_mass = true;
                        picking.p_clicked = com;
                        picking.p_current = com;
                    }

                }
            }


            int joint_selected = run_picking_selection_skeleton(scene_arg, window, skeleton_current);
            if(joint_selected!=-1 && picking.picked_center_of_mass==false)
            {
                vec3 const& picked_position = skeleton_current[joint_selected].p;
                picking.is_selected = true;
                picking.selected_joint = joint_selected;
                picking.p_clicked = picked_position;
                picking.p_current = picked_position;
                picking.selected_joint_true = picking.selected_joint;

                if( joint_selected!=0 && gui_param.type_deformation==1 ) // hack to get the actual selected joint which is the parent
                {
                    picking.selected_joint_true = skeleton.connectivity[joint_selected].parent;
                }

                std::cout<<"Picked joint: "<<joint_selected<<std::endl;
            }


        }

        if(mouse_release){
            picking.picked_center_of_mass = false;
            if(picking.is_selected){
                picking.is_selected = false;
                apply_deformation_on_skeleton();
            }
        }
    }

}

void file_write_ascii(std::string const& filename, vcl::buffer<float> values)
{
    std::ofstream fid(filename.c_str());
    
    for (float v : values)
        fid << v << " ";
    fid.close();

}

void scene_model::keyboard_input(scene_structure& , GLFWwindow* , int key, int , int action, int )
{
    //if (key == GLFW_KEY_S) {
    //    std::cout << "Saved data" << std::endl;
    //    file_write_ascii("filtered_speedup.txt", generic_storage["filtered_speedup"]);
    //    file_write_ascii("filtered_real.txt", generic_storage["filtered_real"]);
    //    file_write_ascii("initial_speed.txt", generic_storage["initial_speed"]);
    //}

    if(key==GLFW_KEY_R && action==GLFW_PRESS)
    {
        if(picking.selected_joint!=-1)
        {

            gui_param.record_anim = !gui_param.record_anim;
            std::cout<<"Recording: "<<gui_param.record_anim<<std::endl;

            if(gui_param.record_anim==true)
            {
                record_position.clear();
                record_rotation.clear();

                timer_recording.start();
                timer_recording.update();
                timer_recording.periodic_event_time_step=record_dt;
                local_time_record = 0.0f;
                record_joint_fixed.clear();

                int joint = skeleton.connectivity[picking.selected_joint].parent;
                if( joint!=-1 )
                    recorded_joint = joint;
                else{
                    recorded_joint = 0;
                }



                int j = recorded_joint;
                record_joint_fixed.clear();
                while(j!=-1)
                {
                    record_joint_fixed.insert(j);
                    j = skeleton.connectivity[j].parent;
                }

            }
            else
            {
                timer_recording.stop();
                record_joint_fixed.clear();


                std::cout<<"Recorded joint: "<<recorded_joint<<std::endl;

                if(record_rotation.size()>0)
                {
                    int N_time_record = record_rotation.size();
                    int N_time_prev = skeleton.anim[recorded_joint].size();
                    int N_time_max = std::max(N_time_record, N_time_prev);

                    for(int kj=0; kj<int(skeleton.anim.size()); ++kj)
                    {
                        skeleton.anim[kj].resize(N_time_max);
                        for(int kt=N_time_prev; kt<N_time_record; ++kt)
                            skeleton.anim[kj][kt].geometry = skeleton.anim[kj][N_time_prev-1].geometry;
                    }
                    for(int kt=0; kt<N_time_record; ++kt){
                        skeleton.anim[recorded_joint][kt].geometry.r = record_rotation[kt];
                        skeleton.anim[recorded_joint][kt].geometry.p = record_position[kt];
                    }
                    for(int kt=N_time_record; kt<N_time_prev; ++kt)
                        skeleton.anim[recorded_joint][kt].geometry = skeleton.anim[recorded_joint][N_time_record-1].geometry;



                    for(int kj=0; kj<int(skeleton.anim.size()); ++kj)
                        for(int kt=0; kt<N_time_max; ++kt)
                            skeleton.anim[kj][kt].time = record_dt*kt;

                    timer_skeleton.t_max = record_dt*(N_time_max-1);

                    if(gui_param.symmetry)
                    {
                        if(symmetrical_joint.find(recorded_joint)!=symmetrical_joint.end())
                        {
                            int js = symmetrical_joint[recorded_joint];
                            quaternion const q = quaternion::axis_angle({1,0,0},3.14f);
                            for(int kt=0; kt<N_time_max; ++kt)
                                skeleton.anim[js][kt].geometry.r = q*skeleton.anim[recorded_joint][kt].geometry.r*conjugate(q);
                        }
                    }
                    //skeleton.anim[11] = skeleton.anim[10];
                }



            }
        }
    }


    if(key==GLFW_KEY_S && action==GLFW_PRESS)
    {
        std::cout<<"Export Model ..."<<std::endl;


        {
            std::ofstream fid("scene_web/mesh.json", std::ofstream::out);
            fid << "{\n";
            fid << "\"vertices\": [";
            size_t const N = skinning.rest_pose.size();
            for(size_t k=0; k<N; ++k)
            {
                vec3 const p = skinning.rest_pose[k];
                fid<<p.x<<", "<<p.y<<", "<<p.z;
                if(k<N-1)
                    fid<<", ";
            }
            fid << "],\n";

            fid << "\"uv\": [";
            for(size_t k=0; k<N; ++k)
            {
                vec2 const p = skinning.deformed.texture_uv[k];
                fid<<p.x<<", "<<p.y;
                if(k<N-1)
                    fid<<", ";
            }
            fid << "],\n";

            fid << "\"connectivity\": [";
            size_t const Nt = skinning.deformed.connectivity.size();
            for(size_t kt=0; kt<Nt; ++kt)
            {
                uint3 const f = skinning.deformed.connectivity[kt];
                fid<<f[0]<<", "<<f[1]<<", "<<f[2];
                if(kt<Nt-1)
                    fid<<", ";
            }
            fid << "]\n";

            fid << "}\n";
            fid.close();
        }


        {
            size_t const N = skeleton.connectivity.size();

            std::ofstream fid("scene_web/skeleton.json", std::ofstream::out);
            fid << "{\n";


            fid << "\"names\": [";
            for(size_t k=0; k<N; ++k)
            {
                fid<<"\""<<skeleton.connectivity[k].name<<"\"";
                if(k<N-1) fid<<", ";
            }
            fid << "],\n";


            fid << "\"parent_id\": [";
            for(size_t k=0; k<N; ++k)
            {
                fid<<skeleton.connectivity[k].parent;
                if(k<N-1) fid<<", ";
            }
            fid << "],\n";


            fid << "\"translation\": [";
            for(size_t k=0; k<N; ++k)
            {
                vec3 const p = skeleton.rest_pose[k].p;
                fid<<p.x<<", "<<p.y<<", "<<p.z;
                if(k<N-1) fid<<", ";
            }
            fid << "],\n";

            fid << "\"rotation\": [";
            for(size_t k=0; k<N; ++k)
            {
                quaternion const q = skeleton.rest_pose[k].r;
                fid<<q.x<<", "<<q.y<<", "<<q.z<<", "<<q.w;
                if(k<N-1) fid<<", ";
            }
            fid << "]\n";

            fid << "}\n";
            fid.close();
        }

        {

            std::ofstream fid("scene_web/rig.json", std::ofstream::out);
            fid << "{\n";
            int N_vertex = skinning.influence.size();
            fid << "\"joint\" : ";
            fid<<"[";
            for(int k_vertex=0; k_vertex<N_vertex; ++k_vertex)
            {
                int N_bones = skinning.influence[k_vertex].size();
                fid<<"[";
                for(int k_bone = 0; k_bone<N_bones; ++k_bone)
                {
                    int const j = skinning.influence[k_vertex][k_bone].joint;
                    fid<<j;
                    if(k_bone<N_bones-1) fid<<", ";
                }
                fid<<"]";
                if(k_vertex<N_vertex-1) fid<<", ";
            }
            fid<<"],\n";

            fid << "\"weight\" : ";
            fid<<"[";
            for(int k_vertex=0; k_vertex<N_vertex; ++k_vertex)
            {
                int N_bones = skinning.influence[k_vertex].size();
                fid<<"[";
                for(int k_bone = 0; k_bone<N_bones; ++k_bone)
                {
                    float const w = skinning.influence[k_vertex][k_bone].weight;
                    fid<<w;
                    if(k_bone<N_bones-1) fid<<", ";
                }
                fid<<"]";
                if(k_vertex<N_vertex-1) fid<<", ";
            }
            fid<<"] \n";

            fid << "}\n";
            fid.close();
        }


        // Export anim
        size_t const N_joint = skeleton.rest_pose.size();
        buffer<buffer<vec3>> center_of_mass(N_joint);
        {
            float const t_min = timer_skeleton.t_min;
            float const t_max = timer_skeleton.t_max;
            int N_time = int((t_max-t_min)/record_dt);
            assert_vcl_no_msg(N_time>1);


            buffer<buffer<vec3>> skeleton_position(N_joint);
            buffer<buffer<quaternion>> skeleton_rotation(N_joint);
            for(int k_time=0; k_time<N_time; ++k_time)
            {
                float const dt = (t_max-t_min)/(N_time-1);
                float const t = t_min + k_time * dt;

                interpolate_skeleton_at_time_with_constraints(skeleton_local_current,
                                                              t,
                                                              skeleton.anim,
                                                              gui_param.interpolate,
                                                              record_joint_fixed);
                skeleton_current = local_to_global(skeleton_local_current, skeleton.connectivity);
                for(int k=0; k<int(skeleton_current.size()); ++k) {
                    skeleton_joint_speed[k].add( skeleton_local_current[k].p, t );
                    skeleton_joint_rotation_speed[k].add( skeleton_local_current[k].r, t );



                    mat3 R_parent = mat3::identity();
                    if(k>0)
                        R_parent = skeleton_current[ skeleton.connectivity[k].parent ].r.matrix();

                    skeleton_speed_per_joint[k].center = skeleton_current[k].p;
                    skeleton_speed_per_joint[k].linear_speed  = R_parent * skeleton_joint_speed[k].avg_speed;
                    skeleton_speed_per_joint[k].angular_speed = R_parent * skeleton_joint_rotation_speed[k].avg_rotation_speed;
                }
                update_center_of_mass();



                for(size_t kj=0; kj<N_joint; ++kj)
                {
                    skeleton_position[kj].push_back(skeleton_local_current[kj].p);
                    skeleton_rotation[kj].push_back(skeleton_local_current[kj].r);
                    center_of_mass[kj].push_back(center_of_mass_per_joint[kj]);
                }


            }


            std::ofstream fid("scene_web/anim.json", std::ofstream::out);

            fid << "{\n";

            fid <<"\"time\" : [";
            for(int kt=0; kt<N_time; ++kt)
            {
                float const dt = (t_max-t_min)/(N_time-1);
                float const t = t_min + kt * dt;
                fid << t;
                if(kt<N_time-1) fid << ", ";
            }
            fid<<"], \n";


            fid <<"\"position\" : [";
            for(int kt=0; kt<N_time; ++kt)
            {
                fid<<"[";
                for(size_t kj=0; kj<N_joint; ++kj)
                {
                    vec3 const& p = skeleton_position[kj][kt];
                    fid << p.x<<", "<<p.y<<", "<<p.z;
                    if(kj<N_joint-1) fid << ", ";

                }
                fid<<"]";
                if(kt<N_time-1) fid << ", ";
            }
            fid<<"], \n";

            fid <<"\"rotation\" : [";
            for(int kt=0; kt<N_time; ++kt)
            {
                fid<<"[";
                for(size_t kj=0; kj<N_joint; ++kj)
                {
                    quaternion const& q = skeleton_rotation[kj][kt];
                    fid << q.x<<", "<<q.y<<", "<<q.z<<", "<<q.w;
                    if(kj<N_joint-1) fid << ", ";

                }
                fid<<"]";
                if(kt<N_time-1) fid << ", ";
            }
            fid<<"] \n";

            fid << "}\n";
            fid.close();

        }

        {
            size_t const N_joint = skeleton.connectivity.size();
            size_t const N_vertex = skinning.rest_pose.size();

            std::ofstream fid("scene_web/velocity_skinning.json", std::ofstream::out);
            fid << "{\n";

            fid << "\"flappy_weights\" : [";
            for(size_t k_vertex=0; k_vertex<N_vertex; ++k_vertex)
            {
                float const w = weight_flappy[k_vertex];
                fid << w;
                if(k_vertex<N_vertex-1) fid<<" ,";
            }
            fid << "], \n";
/*
            fid << "\"center of mass\" : [";
            int const N_time = center_of_mass[0].size();
            for(int kt=0; kt<N_time; ++kt)
            {
                fid<<"[";
                for(size_t kj=0; kj<N_joint; ++kj)
                {
                    vec3 const& com = center_of_mass[kj][kt];
                    fid << com.x<<", "<<com.y<<", "<<com.z;
                    if(kj<N_joint-1) fid << ", ";

                }
                fid<<"]";
                if(kt<N_time-1) fid << ", ";
            }
            fid << "], \n";*/


            fid << "\"vertex_depending_on_joint\" : [";
            for(size_t k_joint=0; k_joint<N_joint; ++k_joint)
            {
                size_t const N = vertex_depending_on_joint[k_joint].size();
                fid<<" [";
                for(size_t k_v=0; k_v<N; ++k_v)
                {
                    fid << vertex_depending_on_joint[k_joint][k_v];
                    if(k_v<N-1) fid << ", ";
                }
                fid<<"]";
                if(k_joint<N_joint-1) fid << ", ";
            }
            fid << "],\n";

            fid << "\"vertex_weight_depending_on_joint\" : [";
            for(size_t k_joint=0; k_joint<N_joint; ++k_joint)
            {
                size_t const N = vertex_weight_depending_on_joint[k_joint].size();
                fid<<" [";
                for(size_t k_v=0; k_v<N; ++k_v)
                {
                    fid << vertex_weight_depending_on_joint[k_joint][k_v];
                    if(k_v<N-1) fid << ", ";
                }
                fid<<"]";
                if(k_joint<N_joint-1) fid << ", ";
            }
            ;
            fid << "]\n";

            fid << "}\n";
            fid.close();
        }



        std::cout<<"Export Done"<<std::endl;

    }

    if (key == GLFW_KEY_L && action == GLFW_PRESS)
	{
		size_t const N_joint = skeleton.connectivity.size();
		size_t const N_vertex = skinning.rest_pose.size();
		std::string temp;

		{
			std::cout << "Loading flappy weights..." << std::endl;
            // try{
            std::ifstream fid("scene_web/velocity_skinning.json", std::ofstream::out);
            // } catch (std::invalid_argument) {

            // }
			// throw away '{\n'




			fid >> temp;
			// throw away "\"flappy_weights\" : ["
			fid >> temp;
			fid >> temp;
			//fid >> temp;

			for (size_t k_vertex = 0; k_vertex < N_vertex; ++k_vertex)
			{
				// throw away ' ,'
				fid.get();
				fid.get();
				fid >> temp;
				weight_flappy[k_vertex] = stof(temp);
			}
			fid.close();
			std::cout << "Flappy weights loaded." << std::endl;
		}

		{
			std::cout << "Loading animation..." << std::endl;
			size_t N_key = 0;
			buffer<float> animation_time;
			buffer<joint_geometry_time> animated_joint;
			buffer<buffer<joint_geometry_time> > animation_data_from_file;
			buffer<buffer<joint_geometry_time> > animated_skeleton;

			std::ifstream fid("scene_web/anim.json", std::ofstream::out);
			// throw away '{\n'
			fid >> temp;
			// throw away "\"time\" : ["
			fid >> temp;
			fid >> temp;
			fid.get();
			fid.get();

			// get 1st key frame time
			fid >> temp;

			std::cout << "Reading time data..." << std::endl;
			// read time
			do {
				animation_time.push_back(stof(temp));
				// throw away ',  '
				//fid.get();
				//fid.get();
				fid >> temp;
			} while (temp.find(']') == std::string::npos);
			// store last keyframe
			animation_time.push_back(stof(temp));

			// set number of keys 
			N_key = animation_time.size();
			animated_joint.resize(N_joint);
			animation_data_from_file.resize(N_key);
			animated_skeleton.resize(N_joint);


			// throw away "position :"
			fid >> temp;
			fid >> temp;

			// throw away " "
			fid.get();
			std::cout << "Reading position data..." << std::endl;
			//std::cout << "N_joint: " << N_joint << "; N_key: " << N_key << std::endl;
			// read position
			for (int kt = 0; kt < N_key; ++kt)
			{
				//std::cout << "Reading key " << kt << ": ";
				// throw away " ["
				fid.get();
				fid.get();
				for (size_t kj = 0; kj < N_joint; ++kj)
				{
					//if (kj == 26)
						//std::cout << "stop";
					vec3 p;
					vec4 q;
					//std::cout << kj << ", ";
					/*vec3 const& p = skeleton_position[kj][kt];
					fid << p.x << ", " << p.y << ", " << p.z;
					if (kj < N_joint - 1) fid << ", ";*/

					for (size_t i = 0; i < 3; i++)
					{
						fid >> temp;
						p[i] = stof(temp);
					}
					animated_joint[kj] = { animation_time[kt], {p, q} };
				}
				//std::cout << std::endl;
				animation_data_from_file[kt] = animated_joint;

			}

			// throw away "rotation :"
			fid >> temp;
			fid >> temp;

			// throw away " "
			fid.get();
			std::cout << "Reading rotation data..." << std::endl;
			// read rotation
			for (int kt = 0; kt < N_key; ++kt)
			{
				// throw away " ["
				fid.get();
				fid.get();
				for (size_t kj = 0; kj < N_joint; ++kj)
				{
					vec4 q;

					for (size_t i = 0; i < 4; i++)
					{
						fid >> temp;
						q[i] = stof(temp);
					}
					animation_data_from_file[kt][kj].geometry.r = q;

				}
			}
			fid.close();
			std::cout << "Reconstructing animation data..." << std::endl;
			// resturcture using joints as primary indexing (to match the data structure used in this application)
			for (int kj = 0; kj < N_joint; ++kj)
			{
				animated_skeleton[kj].resize(N_key);
				for (int kt = 0; kt < N_key; ++kt)
				{
					animated_skeleton[kj][kt] = animation_data_from_file[kt][kj];
				}
			}

			skeleton.anim = animated_skeleton;
			timer_skeleton.t_max = find_animation_length(animated_skeleton);
			std::cout << "Animation loaded." << std::endl;

		}
	}
}

void scene_model::compute_skinning_wave_propagation()
{
    const size_t N_vertex = skinning.rest_pose.size();

    std::vector<mat3> T(skeleton_current.data.size());
    std::vector<vec3> tr(skeleton_current.data.size());

    //for (int k = 0; k < T.size(); ++k)
    //{
    //    const quaternion& q = skeleton_current.data[k].r;
    //    const vec3& p = skeleton_current.data[k].p;

    //    const quaternion& q0 = skeleton_rest_pose.data[k].r;
    //    const vec3& p0 = skeleton_rest_pose.data[k].p;

    //    const quaternion qT = q * conjugate(q0);
    //    const mat3 mT = qT.matrix();

    //    tr[k] = p - mT * p0;
    //    T[k] = mT;
    //}
    // noise_filter(skeleton_current.data[0].p, timer.t);








    //vec3 v1 = saved_rotation[velocity_index_past + 1];
    //vec3 v0 = saved_rotation[velocity_index_past];
    //vec3 angular_speed = (1 - relative) * v0 + relative * v1;

    //if (joint != 0) {
    //    quaternion R_parent = skeleton_current[joint].r;
    //    angular_speed = R_parent.matrix() * angular_speed;
    //}
    float max_distance_to_picked = 0.0f;

    for(size_t k = 0; k < N_vertex; ++k){
        if(distance_to_joint[k][picking.selected_joint] > max_distance_to_picked){
            max_distance_to_picked = distance_to_joint[k][picking.selected_joint];
        }
    }


    for (size_t k = 0; k < N_vertex; ++k)
    {
        const buffer<skinning_influence>& influence = skinning.influence.data[k];

        // Transformation matrix for skinning
        mat3 M = mat3::zero();
        vec3 t = { 0,0,0 };

        mat3 M_cur = mat3::zero();
        vec3 t_cur = { 0,0,0 };

        for (size_t kb = 0; kb < influence.size(); ++kb) // Loop over influencing joints
        {
            const int idx = influence.data[kb].joint;
            const float w = influence.data[kb].weight;
            float relative = 1.0;

            int velocity_index_past = skeleton_previous.data.size() - 2;
            
            if (picking.selected_joint != -1)
            {

                float distance_vertex_to_joint = distance_to_joint[k][picking.selected_joint];
                float distance_to_time_index_scaling = gui_param.wave_compression_magnitude;
                

                //float velocity_index_past_float = std::max(std::min(skeleton_previous.data.size() - 2 - distance_vertex_to_joint * distance_to_time_index_scaling * weight_flappy.data[k]  + 1, skeleton_previous.data.size() - 2.0f), 0.0f);
                float velocity_index_past_float = std::max(std::min(skeleton_previous.data.size() - 1 - distance_vertex_to_joint * distance_to_time_index_scaling * weight_flappy.data[k] + 1, skeleton_previous.data.size() - 2.0f), 0.0f);
                velocity_index_past = int(velocity_index_past_float);
                relative = (velocity_index_past_float - velocity_index_past) ;



            }
            joint_geometry const& joint_past = skeleton_previous.data[velocity_index_past + 1][idx];
            const quaternion& r_past = joint_past.r;
            const vec3& p_past = joint_past.p;

            joint_geometry& joint_past_prev = skeleton_previous.data[velocity_index_past][idx]; // the latest one is on the end of the array
            quaternion& r_past_prev = joint_past_prev.r;
            vec3& p_past_prev = joint_past_prev.p;

            
            quaternion r = (1 - relative) * r_past_prev + relative * r_past; //slerp(r_past_prev, r_past, relative); //r_past
            //quaternion r = slerp(r_past_prev, r_past, relative);
            vec3 p = (1-relative)*p_past_prev + relative*p_past;
            // quaternion r = r_past;
            // vec3 p = p_past;

            joint_geometry const& joint_current = skeleton_previous.data[skeleton_previous.data.size() - 1][idx]; // the latest one is on the end of the array
            const quaternion& r_cur = joint_current.r;
            const vec3& p_cur = joint_current.p;

            

            //const quaternion& r = skeleton_current.data[idx].r;
            //const vec3& p = skeleton_current.data[idx].p;

            const quaternion& r0 = skeleton_rest_pose.data[idx].r;
            const quaternion r0_inv = conjugate(r0);
            const vec3& p0 = skeleton_rest_pose.data[idx].p;


            // Convert rotation/translation to matrix
            // mat4 T = mat4::from_mat3_vec3(r.matrix(), p);
            // mat4 T0_inv = mat4::from_mat3_vec3(conjugate(r0).matrix(), conjugate(r0).apply(-p0)); // inverse

            // Skinning
            //M += w*T*T0_inv;

            M += w * (r * r0_inv).matrix();
            t += w * (p - (r * r0_inv).apply(p0));

            M_cur += w * (r_cur * r0_inv).matrix();
            t_cur += w * (p_cur - (r_cur * r0_inv).apply(p0));  
            //M += w * T[idx];
            //t += w * tr[idx];
        }

        // t = t * + norm(t) * distance_to_joint[k][picking.selected_joint];
        // M = M * distance_to_joint[k][picking.selected_joint];


        // Apply skinning transform on vertex
        const vec3& p0 = skinning.rest_pose.data[k];
        vec3 past_deform = M * p0 + t;
        vec3 cur_deform = M_cur * p0 + t_cur;
        // float rel = 1.0f/(5.0f*distance_to_joint[k][picking.selected_joint]+1.0f);
        float rel = distance_to_joint[k][picking.selected_joint]/max_distance_to_picked;
        // rel = rel*rel;
        rel = 1.0f - rel;
        rel = powf(rel+0.077f, 0.6);
        // skinning.deformed.position.data[k] = (1.0f - rel) * cur_deform + rel * past_deform;
        skinning.deformed.position.data[k] = past_deform;

        const vec3& n0 = skinning.rest_pose_normal.data[k];
        vec3 past_deform_n = M * n0;
        // vec3 cur_deform_n = M_cur * n0;
        skinning.deformed.normal.data[k] = past_deform_n;
    }


    // Drag affect
    // --------------------------------------------->    
    // float const default_flappy_factor = 1.0f;
    // size_t const N_joint = skeleton_current.size();
    // for(size_t joint=0; joint<N_joint; ++joint)
    // {
    //     vec3 const& p_joint = skeleton_speed_per_joint.data[joint].center;
    //     vec3 linear_speed = (skeleton_speed_per_joint.data[joint].linear_speed / timer_skeleton.scale + skeleton_fake_speed_per_joint.data[joint].linear_speed);
    //     buffer<int> const& vertices = vertex_depending_on_joint.data[joint];
    //     buffer<float> const& vertices_weights = vertex_weight_depending_on_joint.data[joint];
    //     size_t const N_vertex_dependency = vertices.size();


    //     auto& saved_rotation = skeleton_joint_rotation_speed[joint].saved_rotation;
    //     vec3 angular_speed = saved_rotation[saved_rotation.size()-1];

    //     if (joint != 0) {
    //         quaternion R_parent = skeleton_current[joint].r;
    //         angular_speed = R_parent.matrix() * angular_speed;
    //     }



    //     vec3 un_angular_speed = normalize(angular_speed);

    //     for(size_t k_vertex=0; k_vertex<N_vertex_dependency; ++k_vertex)
    //     {
    //         size_t const vertex = vertices.data[k_vertex];
    //         float const w_skinning = vertices_weights.data[k_vertex];
    //         float const w_flappy = default_flappy_factor * flapping_power * weight_flappy.data[vertex];
    //         vec3 const flappy = deformation_flappy_linear_speed(w_flappy, linear_speed);
    //         skinning.deformed.position[vertex] += w_skinning * (flappy);


    //         // vec3 const& p_vertex = save_skinning.data[vertex];
    //         // vec3 const u_joint_vertex = p_vertex - p_joint;
    //         // vec3 const vertex_speed = cross(angular_speed, u_joint_vertex);
    //         // float const vertex_speed_norm = norm(vertex_speed);

    //         // if (norm(angular_speed) > 1e-2f) {
    //         //     vec3 const flappy_rotation = deformation_flappy_rotation_speed(w_flappy, p_vertex, p_joint, un_angular_speed, vertex_speed_norm, gui_param.flappy_max_angle);

    //         //     skinning.deformed.position[vertex] += w_skinning * flappy_rotation;
    //         // }

    //     }

    // }
    // --------------------------------------------->  

}


float max_joint_filtered_value::attenuation = 1.0f;
float max_joint_filtered_value::frequency = 9.0f;
float max_joint_filtered_value::frequency_slope = 9.0f;

void max_joint_filtered_value::update(vec3 const& new_candidate, float t_current)
{
    float current = norm(new_candidate);
    float damping = evaluate_damping(t_current);
    float previous = damping *norm(value);

    float alpha = 0.6f; // interpolation coefficient between old and new value => limit suddent jumps
    if (current > previous) {
        value = alpha*new_candidate+(1-alpha)*damping*value;
        t0 = alpha*t_current+(1-alpha)*t0;
    }
}
float max_joint_filtered_value::evaluate_damping(float t_current) const
{
    float t = t_current - t0;
    return std::exp(-attenuation * t);
}
vec3 max_joint_filtered_value::evaluateSinusoidal(float t_current) const   //AICI TREBUIE SA MODIFIC PT RECTANGULAR
{
    float t = t_current - t0;

    //SINUSOIDAL WAVE  --DONE

    float filter = evaluate_damping(t_current) * std::cos((frequency + frequency_slope * sqrt(t)) * t);
    return filter * value; 

}

vec3 max_joint_filtered_value::evaluateTriangular(float t_current) const
{
    float t = t_current - t0;

    //THIS IS FOR THE TRIANGULAR WAVE  --DONE

    float period = 1.0f / frequency;
    float halfPeriod = period / 2.0f;
    float remainder = std::fmod(t, period);
    float phase = remainder / halfPeriod;
    float filter = evaluate_damping(t_current) * ((phase < 1.0f) ? (2.0f * phase - 1.0f) : (3.0f - 2.0f * phase));
    return filter * value; 

    /*float phase = std::fmod(u, 2 * 3.14159f);
    int direction = phase < 3.14159f ? 1 : -1; // 1 for increasing, -1 for decreasing
    float amplitude = std::exp(-attenuation * k * dt);
    float result = (2 / 3.14159f) * std::asin(std::sin(u)) * amplitude * direction;

    return result; */
}

vec3 max_joint_filtered_value::evaluateRectangular(float t_current) const
{
    float t = t_current - t0;

    //RECTANGULAR WAVE  --DONE

    float period = 1.0 / frequency;
    float halfPeriod = period / 2.0;
    float remainder = fmod(t, period);
    float filter =  evaluate_damping(t_current) * ((remainder < halfPeriod) ? 1.0f : -1.0f);
    return filter * value; 

    /*float period = 1.0 / frequency;
     float halfPeriod = period / 2.0;
    float remainder = k  % (int) (period / dt);
    float constant = -attenuation * ( k - remainder) * dt;

    float result = std::exp(constant);

    if (sin(u) > 0) {
        return result;
    }
    return -result;*/
}

vec3 max_joint_filtered_value::evaluatePendulum(float t_current) const
{
    float t = t_current - t0;

    //PENDULUM WAVE    --   NU E TERMINAT OFULET

    float period = 1.0f / frequency;
    float dt = 0.017f;
    float k = t / dt;
    float constant = -attenuation * k * dt;
    float result = evaluate_damping(t_current) * (std::cos(frequency * t * 2.0f * 3.14159f) + std::cos(3 * (frequency * t * 2.0f * 3.14159f)) / 9.0f);
    return result * value;

    /*float period = 1.0 / frequency;
    float halfPeriod = period / 2.0;
    float remainder = k % (int)(period / dt);
    float constant = -attenuation * (k) * dt;

    float result = std::exp(constant) * (sin(u) + sin(3 * u) / 9.0f);

    return result;*/
}


vec3 max_joint_filtered_value::evaluateBipBip(float t_current) const
{
    float t = t_current - t0;

    //BIP BIP COYOTE WAVE  --DONE

    float filter = evaluate_damping(t_current) * std::cos((frequency + frequency_slope * sqrt(t)) * t) * 3;
    return filter * value; 

    /*f_slope = 3;
    df = k * dt * f_slope;
    u = (frequency + df) * dt * k * 2 * 3.14159f;
    float result = std::exp(-attenuation * k * dt) * sin(u);

    return result;*/

}

#endif

