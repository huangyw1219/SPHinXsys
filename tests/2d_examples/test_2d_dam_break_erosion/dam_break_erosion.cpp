/**
 * @file 	dam_break_erosion.cpp
 * @brief 	2D dam break with soil erosion/deposition.
 */
#include "dam_break_erosion.h"

#include <algorithm>
#include <iomanip>

namespace
{
class WaterVelocityFreeze : public LocalDynamics
{
  public:
    explicit WaterVelocityFreeze(SPHBody &water_body)
        : LocalDynamics(water_body),
          vel_(particles_->getVariableDataByName<Vecd>("Velocity")) {};

    void update(size_t index_i, Real dt = 0.0)
    {
        vel_[index_i] = Vecd::Zero();
    };

  protected:
    Vecd *vel_;
};

void erodeSoilParticles(RealBody &soil_body, FluidBody &eroded_body)
{
    auto &soil_particles = soil_body.getBaseParticles();
    auto &eroded_particles = eroded_body.getBaseParticles();
    int *erosion_flag = soil_particles.getVariableDataByName<int>("ErosionFlag");
    Vecd *soil_pos = soil_particles.getVariableDataByName<Vecd>("Position");
    Vecd *soil_vel = soil_particles.getVariableDataByName<Vecd>("Velocity");
    Real *soil_rho = soil_particles.getVariableDataByName<Real>("Density");

    Vecd *eroded_pos = eroded_particles.getVariableDataByName<Vecd>("Position");
    Vecd *eroded_vel = eroded_particles.getVariableDataByName<Vecd>("Velocity");
    Real *eroded_rho = eroded_particles.getVariableDataByName<Real>("Density");
    Real *eroded_mass = eroded_particles.getVariableDataByName<Real>("Mass");
    Real *eroded_vol = eroded_particles.getVariableDataByName<Real>("VolumetricMeasure");
    int *eroded_state = eroded_particles.getVariableDataByName<int>("ErosionState");

    StdVec<size_t> to_erode;
    for (size_t i = 0; i < soil_particles.TotalRealParticles(); ++i)
    {
        if (erosion_flag[i] == 1)
        {
            to_erode.push_back(i);
        }
    }
    if (to_erode.empty())
        return;

    std::sort(to_erode.begin(), to_erode.end(), std::greater<size_t>());
    for (size_t index_i : to_erode)
    {
        eroded_particles.checkEnoughReserve();
        size_t buffer_index = eroded_particles.TotalRealParticles();
        size_t new_index = eroded_particles.createRealParticleFrom(buffer_index);

        eroded_pos[new_index] = soil_pos[index_i];
        eroded_vel[new_index] = soil_vel[index_i];
        eroded_rho[new_index] = soil_rho[index_i];
        eroded_vol[new_index] = particle_volume;
        eroded_mass[new_index] = eroded_rho[new_index] * eroded_vol[new_index];
        eroded_state[new_index] = 1;

        soil_particles.switchToBufferParticle(index_i);
    }
}

void depositSoilParticles(RealBody &soil_body, FluidBody &eroded_body)
{
    auto &soil_particles = soil_body.getBaseParticles();
    auto &eroded_particles = eroded_body.getBaseParticles();
    int *deposition_flag = eroded_particles.getVariableDataByName<int>("DepositionFlag");
    int *erosion_state = soil_particles.getVariableDataByName<int>("ErosionState");
    Vecd *soil_pos = soil_particles.getVariableDataByName<Vecd>("Position");
    Vecd *soil_vel = soil_particles.getVariableDataByName<Vecd>("Velocity");
    Real *soil_rho = soil_particles.getVariableDataByName<Real>("Density");
    Real *soil_mass = soil_particles.getVariableDataByName<Real>("Mass");
    Real *soil_vol = soil_particles.getVariableDataByName<Real>("VolumetricMeasure");
    Mat3d *soil_stress = soil_particles.getVariableDataByName<Mat3d>("StressTensor3D");
    Vecd *soil_pos0 = soil_particles.getVariableDataByName<Vecd>("InitialPosition");
    Vecd *soil_disp = soil_particles.getVariableDataByName<Vecd>("Displacement");

    Vecd *eroded_pos = eroded_particles.getVariableDataByName<Vecd>("Position");
    Vecd *eroded_vel = eroded_particles.getVariableDataByName<Vecd>("Velocity");
    Real *eroded_rho = eroded_particles.getVariableDataByName<Real>("Density");

    StdVec<size_t> to_deposit;
    for (size_t i = 0; i < eroded_particles.TotalRealParticles(); ++i)
    {
        if (deposition_flag[i] == 1)
        {
            to_deposit.push_back(i);
        }
    }
    if (to_deposit.empty())
        return;

    std::sort(to_deposit.begin(), to_deposit.end(), std::greater<size_t>());
    for (size_t index_i : to_deposit)
    {
        soil_particles.checkEnoughReserve();
        size_t buffer_index = soil_particles.TotalRealParticles();
        size_t new_index = soil_particles.createRealParticleFrom(buffer_index);

        soil_pos[new_index] = eroded_pos[index_i];
        soil_pos0[new_index] = eroded_pos[index_i];
        soil_disp[new_index] = Vecd::Zero();
        soil_vel[new_index] = eroded_vel[index_i];
        soil_rho[new_index] = eroded_rho[index_i];
        soil_vol[new_index] = particle_volume;
        soil_mass[new_index] = soil_rho[new_index] * soil_vol[new_index];
        soil_stress[new_index] = Mat3d::Zero();
        erosion_state[new_index] = 2;

        eroded_particles.switchToBufferParticle(index_i);
    }
}
} // namespace

int main(int ac, char *av[])
{
    SPHSystem sph_system(system_domain_bounds, particle_spacing_ref);
    sph_system.handleCommandlineOptions(ac, av);

    RealBody soil_block(sph_system, makeShared<Soil>("SoilBody"));
    soil_block.defineMaterial<PlasticContinuum>(rho0_s, c_s, Youngs_modulus, poisson, friction_angle, cohesion);
    ParticleBuffer<ReserveSizeFactor> soil_buffer(0.2);
    soil_block.generateParticlesWithReserve<BaseParticles, Lattice>(soil_buffer);

    FluidBody water_block(sph_system, makeShared<WaterBlock>("WaterBody"));
    water_block.defineMaterial<WeaklyCompressibleFluid>(rho0_f, c_f);
    water_block.generateParticles<BaseParticles, Lattice>();

    FluidBody eroded_soil(sph_system, makeShared<Soil>("ErodedSoil"));
    eroded_soil.defineClosure<WeaklyCompressibleFluid, HerschelBulkleyViscosity>(
        ConstructArgs(rho0_s, c_s), ConstructArgs(hb_min_shear_rate, hb_max_shear_rate, hb_consistency, hb_power_index, hb_yield_stress));
    ParticleBuffer<ReserveSizeFactor> eroded_buffer(0.2);
    eroded_soil.generateParticlesWithReserve<BaseParticles, Lattice>(eroded_buffer);

    SolidBody wall_boundary(sph_system, makeShared<WallBoundary>("WallBoundary"));
    wall_boundary.defineMaterial<Solid>();
    wall_boundary.generateParticles<BaseParticles, Lattice>();

    auto &soil_particles = soil_block.getBaseParticles();
    int *soil_erosion_state = soil_particles.registerStateVariableData<int>("ErosionState");
    soil_particles.addEvolvingVariable<int>("ErosionState");
    size_t initial_soil_particle_count = soil_particles.TotalRealParticles();
    for (size_t i = 0; i < soil_particles.TotalRealParticles(); ++i)
        soil_erosion_state[i] = 0;

    auto &eroded_particles = eroded_soil.getBaseParticles();
    int *eroded_state = eroded_particles.registerStateVariableData<int>("ErosionState");
    eroded_particles.addEvolvingVariable<int>("ErosionState");
    for (size_t i = 0; i < eroded_particles.TotalRealParticles(); ++i)
        eroded_state[i] = 1;

    // Remove initially generated eroded particles; keep as buffer.
    if (eroded_particles.TotalRealParticles() > 0)
    {
        for (size_t i = eroded_particles.TotalRealParticles(); i-- > 0;)
        {
            eroded_particles.switchToBufferParticle(i);
        }
    }

    InnerRelation soil_block_inner(soil_block);
    ContactRelation soil_block_contact(soil_block, {&wall_boundary});
    ContactRelation soil_water_contact(soil_block, {&water_block});
    ComplexRelation soil_block_complex(soil_block_inner, soil_block_contact);

    InnerRelation water_block_inner(water_block);
    ContactRelation water_wall_contact(water_block, {&wall_boundary});
    ContactRelation water_fluid_contact(water_block, {&eroded_soil, &soil_block});
    ComplexRelation water_block_complex(water_block_inner, {&water_fluid_contact, &water_wall_contact});

    InnerRelation eroded_inner(eroded_soil);
    ContactRelation eroded_wall_contact(eroded_soil, {&wall_boundary});
    ContactRelation eroded_fluid_contact(eroded_soil, {&water_block, &soil_block});
    ComplexRelation eroded_complex(eroded_inner, {&eroded_fluid_contact, &eroded_wall_contact});

    Gravity gravity(Vecd(0.0, -gravity_g));
    SimpleDynamics<GravityForce<Gravity>> soil_gravity(soil_block, gravity);
    SimpleDynamics<GravityForce<Gravity>> water_gravity(water_block, gravity);
    SimpleDynamics<GravityForce<Gravity>> eroded_gravity(eroded_soil, gravity);

    SimpleDynamics<NormalDirectionFromBodyShape> wall_boundary_normal_direction(wall_boundary);
    SimpleDynamics<NormalDirectionFromBodyShape> soil_boundary_normal_direction(soil_block);
    SimpleDynamics<SoilInitialCondition> soil_initial_condition(soil_block);
    SimpleDynamics<WaterInitialCondition> water_initial_condition(water_block, 0.4);
    InteractionWithUpdate<LinearGradientCorrectionMatrixComplex> soil_correction_matrix(soil_block_inner, soil_block_contact);
    Dynamics1Level<continuum_dynamics::PlasticIntegration1stHalfWithWallRiemann> soil_stress_relaxation(soil_block_inner, soil_block_contact);
    Dynamics1Level<continuum_dynamics::PlasticIntegration2ndHalfWithWallRiemann> soil_density_relaxation(soil_block_inner, soil_block_contact);
    InteractionWithUpdate<fluid_dynamics::DensitySummationComplexFreeSurface> soil_density_by_summation(soil_block_inner, soil_block_contact);
    InteractionDynamics<continuum_dynamics::StressDiffusion> soil_stress_diffusion(soil_block_inner);
    InteractionWithUpdate<FreeSurfaceIndicationComplex> soil_surface_indicator(soil_block_inner, soil_block_contact);
    InteractionWithUpdate<TransportVelocityCorrectionComplex<AllParticles>> soil_transport_velocity_correction(soil_block_inner, soil_block_contact);
    InteractionWithUpdate<FreeSurfaceNormalComplex> soil_free_surface_normal(soil_block_inner, soil_block_contact);
    ReduceDynamics<fluid_dynamics::AcousticTimeStep> soil_acoustic_time_step(soil_block, 0.4);

    InteractionWithUpdate<LinearGradientCorrectionMatrixComplex> water_correction_matrix(water_block_inner, water_wall_contact);
    Dynamics1Level<fluid_dynamics::MultiPhaseIntegration1stHalfWithWallRiemann> water_pressure_relaxation(
        water_block_inner, water_fluid_contact, water_wall_contact);
    Dynamics1Level<fluid_dynamics::MultiPhaseIntegration2ndHalfWithWallRiemann> water_density_relaxation(
        water_block_inner, water_fluid_contact, water_wall_contact);
    InteractionWithUpdate<fluid_dynamics::BaseDensitySummationComplex<Inner<>, Contact<>, Contact<>>>
        water_density_by_summation(water_block_inner, water_fluid_contact, water_wall_contact);
    DampingWithRandomChoice<InteractionSplit<DampingPairwiseWithWall<Vec2d, FixedDampingRate>>>
        water_damping(0.2, DynamicsArgs(water_block_inner, "Velocity", mu_f), DynamicsArgs(water_wall_contact, "Velocity", mu_f));
    InteractionDynamics<fluid_dynamics::BoundingFromWall> water_near_wall_bounding(water_wall_contact);
    ReduceDynamics<fluid_dynamics::AcousticTimeStep> water_acoustic_time_step(water_block, 0.4);

    InteractionWithUpdate<LinearGradientCorrectionMatrixComplex> eroded_correction_matrix(eroded_inner, eroded_wall_contact);
    Dynamics1Level<fluid_dynamics::MultiPhaseIntegration1stHalfWithWallRiemann> eroded_pressure_relaxation(
        eroded_inner, eroded_fluid_contact, eroded_wall_contact);
    Dynamics1Level<fluid_dynamics::MultiPhaseIntegration2ndHalfWithWallRiemann> eroded_density_relaxation(
        eroded_inner, eroded_fluid_contact, eroded_wall_contact);
    InteractionWithUpdate<fluid_dynamics::BaseDensitySummationComplex<Inner<>, Contact<>, Contact<>>>
        eroded_density_by_summation(eroded_inner, eroded_fluid_contact, eroded_wall_contact);
    InteractionDynamics<fluid_dynamics::DistanceFromWall> eroded_distance_to_wall(eroded_wall_contact);
    InteractionDynamics<fluid_dynamics::BoundingFromWall> eroded_near_wall_bounding(eroded_wall_contact);
    InteractionWithUpdate<fluid_dynamics::VelocityGradientWithWall<LinearGradientCorrection>> eroded_vel_grad(
        eroded_inner, eroded_wall_contact);
    SimpleDynamics<fluid_dynamics::ShearRateDependentViscosity> eroded_shear_viscosity(eroded_soil);
    InteractionWithUpdate<fluid_dynamics::NonNewtonianViscousForceWithWall<AngularConservative>> eroded_viscous_force(
        eroded_inner, eroded_wall_contact);
    ReduceDynamics<fluid_dynamics::AcousticTimeStep> eroded_acoustic_time_step(eroded_soil, 0.4);

    InteractionDynamics<ErosionIdentification> erosion_identification(soil_water_contact, erosion_velocity_threshold);
    SimpleDynamics<WaterVelocityFreeze> water_velocity_freeze(water_block);
    SimpleDynamics<DepositionIdentification> deposition_identification(eroded_soil, deposition_velocity_threshold);
    SimpleDynamics<UpdateDisplacement> update_displacement(soil_block);

    SimpleDynamics<continuum_dynamics::VerticalStress> vertical_stress(soil_block);
    SimpleDynamics<continuum_dynamics::AccDeviatoricPlasticStrain> accumulated_deviatoric_plastic_strain(soil_block);

    BodyStatesRecordingToVtp body_states_recording(sph_system);
    body_states_recording.addToWrite<Real>(soil_block, "Pressure");
    body_states_recording.addToWrite<Real>(soil_block, "Density");
    body_states_recording.addToWrite<Real>(soil_block, "VerticalStress");
    body_states_recording.addToWrite<Real>(soil_block, "AccDeviatoricPlasticStrain");
    body_states_recording.addToWrite<int>(soil_block, "Indicator");
    body_states_recording.addToWrite<int>(soil_block, "ErosionState");
    body_states_recording.addToWrite<Vecd>(soil_block, "Velocity");
    body_states_recording.addToWrite<Vecd>(soil_block, "Displacement");
    body_states_recording.addToWrite<Vecd>(water_block, "Velocity");
    body_states_recording.addToWrite<Vecd>(eroded_soil, "Velocity");

    sph_system.initializeSystemCellLinkedLists();
    sph_system.initializeSystemConfigurations();
    wall_boundary_normal_direction.exec();
    soil_boundary_normal_direction.exec();
    soil_gravity.exec();
    water_gravity.exec();
    eroded_gravity.exec();
    soil_initial_condition.exec();
    water_initial_condition.exec();
    soil_correction_matrix.exec();
    water_correction_matrix.exec();
    eroded_correction_matrix.exec();
    update_displacement.exec();

    Real &physical_time = *sph_system.getSystemVariableDataByName<Real>("PhysicalTime");
    size_t number_of_iterations = 0;
    int initial_velocity_freeze_steps = 200;
    int screen_output_interval = 500;
    Real End_Time = 2.0;
    Real D_Time = End_Time / 50;

    body_states_recording.writeToFile();

    while (physical_time < End_Time)
    {
        Real integration_time = 0.0;
        while (integration_time < D_Time)
        {
            soil_density_by_summation.exec();
            soil_surface_indicator.exec();
            soil_free_surface_normal.exec();
            soil_transport_velocity_correction.exec();

            water_density_by_summation.exec();
            water_damping.exec();
            water_near_wall_bounding.exec();
            eroded_density_by_summation.exec();

            if (number_of_iterations < static_cast<size_t>(initial_velocity_freeze_steps))
            {
                water_velocity_freeze.exec();
            }

            Real dt = SMIN(soil_acoustic_time_step.exec(), water_acoustic_time_step.exec());
            bool has_eroded_particles = eroded_soil.getBaseParticles().TotalRealParticles() > 0;
            if (has_eroded_particles)
            {
                dt = SMIN(dt, eroded_acoustic_time_step.exec());
            }

            soil_stress_diffusion.exec();
            soil_stress_relaxation.exec(dt);
            soil_density_relaxation.exec(dt);

            water_pressure_relaxation.exec(dt);
            water_density_relaxation.exec(dt);

            if (has_eroded_particles)
            {
                eroded_distance_to_wall.exec();
                eroded_near_wall_bounding.exec();
                eroded_vel_grad.exec();
                eroded_shear_viscosity.exec();
                eroded_viscous_force.exec();
                eroded_pressure_relaxation.exec(dt);
                eroded_density_relaxation.exec(dt);
            }

            erosion_identification.exec();
            deposition_identification.exec();
            erodeSoilParticles(soil_block, eroded_soil);
            depositSoilParticles(soil_block, eroded_soil);

            integration_time += dt;
            physical_time += dt;

            if (number_of_iterations % screen_output_interval == 0)
            {
                std::cout << std::fixed << std::setprecision(9) << "N=" << number_of_iterations << std::setprecision(4)
                          << "\tTime = " << physical_time << std::scientific << "\tdt = " << dt << "\n";
                size_t current_soil_particles = soil_particles.TotalRealParticles();
                size_t current_eroded_particles = eroded_particles.TotalRealParticles();
                std::cout << std::fixed << std::setprecision(0)
                          << "SoilInit=" << initial_soil_particle_count
                          << "\tSoilRemain=" << current_soil_particles
                          << "\tEroded=" << current_eroded_particles
                          << "\tSum=" << current_soil_particles + current_eroded_particles << "\n";
            }
            number_of_iterations++;

            soil_block.updateCellLinkedList();
            water_block.updateCellLinkedList();
            eroded_soil.updateCellLinkedList();
            soil_block_complex.updateConfiguration();
            water_block_complex.updateConfiguration();
            eroded_complex.updateConfiguration();
            soil_water_contact.updateConfiguration();
            soil_correction_matrix.exec();
            water_correction_matrix.exec();
            eroded_correction_matrix.exec();
            update_displacement.exec();
        }
        vertical_stress.exec();
        accumulated_deviatoric_plastic_strain.exec();
        body_states_recording.writeToFile();
    }

    return 0;
}
