/**
 * @file 	water_soil_erosion.cpp
 * @brief 	2D water-soil coupling with erosion and deposition.
 */
#include "sphinxsys.h"

#include <fstream>
#include <iomanip>
#include <map>
#include <sstream>
using namespace SPH;
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real DL = 4.0;                        /**< Tank length. */
Real DH = 0.5;                        /**< Tank height. */
Real soil_height = 0.3;
Real soil_tail_height = 0.1;
Real soil_flat_length = 2.0;
Real soil_slope_end = 3.0;
Real water_length = 1.0;
Real water_height = 0.1;
Real particle_spacing_ref = 0.02; /**< Initial reference particle spacing. */
Real BW = particle_spacing_ref * 4.0;
BoundingBoxd system_domain_bounds(Vec2d(-BW, -BW), Vec2d(DL + BW, DH + BW));
//----------------------------------------------------------------------
//	Material properties.
//----------------------------------------------------------------------
Real gravity_g = 9.8;
// water
Real rho0_f = 1000.0;
Real c_f = 50.0;
Real U_f = 1.0;
// soil (DP)
Real rho0_s = 1850.0;
Real Youngs_modulus = 1.8e6;
Real poisson = 0.3;
Real c_s = sqrt(Youngs_modulus / (rho0_s * 3.0 * (1.0 - 2.0 * poisson)));
Real cohesion = 5.0e3;
Real friction_angle = 25.0 * Pi / 180.0;
// eroded soil (HBP)
Real hb_min_shear_rate = 1.0e-4;
Real hb_max_shear_rate = 1.0e3;
Real hb_consistency_index = 15.0;
Real hb_power_index = 0.45;
Real hb_yield_stress = 50.0;
// erosion / deposition criteria
Real erosion_tau_crit = 15.0;
Real erosion_coeff = 1.0;
Real deposition_velocity = 0.05;
//----------------------------------------------------------------------
//	Geometric shapes.
//----------------------------------------------------------------------
class WaterBlock : public MultiPolygonShape
{
  public:
    explicit WaterBlock(const std::string &shape_name) : MultiPolygonShape(shape_name)
    {
        std::vector<Vecd> water_shape{
            Vecd(0.0, soil_height),
            Vecd(0.0, soil_height + water_height),
            Vecd(water_length, soil_height + water_height),
            Vecd(water_length, soil_height),
            Vecd(0.0, soil_height)};
        multi_polygon_.addAPolygon(water_shape, ShapeBooleanOps::add);
    }
};

class SoilBlock : public MultiPolygonShape
{
  public:
    explicit SoilBlock(const std::string &shape_name) : MultiPolygonShape(shape_name)
    {
        std::vector<Vecd> soil_shape{
            Vecd(0.0, 0.0),
            Vecd(0.0, soil_height),
            Vecd(soil_flat_length, soil_height),
            Vecd(soil_slope_end, soil_tail_height),
            Vecd(DL, soil_tail_height),
            Vecd(DL, 0.0),
            Vecd(0.0, 0.0)};
        multi_polygon_.addAPolygon(soil_shape, ShapeBooleanOps::add);
    }
};

class WallBoundary : public ComplexShape
{
  public:
    explicit WallBoundary(const std::string &shape_name) : ComplexShape(shape_name)
    {
        add<GeometricShapeBox>(Transform(Vec2d(0.5 * DL, 0.5 * DH)), Vec2d(0.5 * DL + BW, 0.5 * DH + BW));
        subtract<GeometricShapeBox>(Transform(Vec2d(0.5 * DL, 0.5 * DH)), Vec2d(0.5 * DL, 0.5 * DH));
    }
};
//----------------------------------------------------------------------
//	Soil initial condition with in-situ stress.
//----------------------------------------------------------------------
class SoilInitialCondition : public continuum_dynamics::ContinuumInitialCondition
{
  public:
    explicit SoilInitialCondition(RealBody &soil_body)
        : continuum_dynamics::ContinuumInitialCondition(soil_body) {};

  protected:
    void update(size_t index_i, Real dt)
    {
        Real y = pos_[index_i][1];
        Real gama = 1.0 - sin(friction_angle);
        Real stress_yy = -rho0_s * gravity_g * y;
        stress_tensor_3D_[index_i](1, 1) = stress_yy;
        stress_tensor_3D_[index_i](0, 0) = stress_yy * gama;
        stress_tensor_3D_[index_i](2, 2) = stress_yy * gama;
    };
};
//----------------------------------------------------------------------
//	Displacement tracking for soil.
//----------------------------------------------------------------------
class SoilDisplacementUpdate : public LocalDynamics
{
  public:
    explicit SoilDisplacementUpdate(RealBody &soil_body)
        : LocalDynamics(soil_body),
          pos0_(particles_->registerStateVariableDataFrom<Vecd>("InitialPosition", "Position")),
          pos_(particles_->getVariableDataByName<Vecd>("Position")),
          displacement_(particles_->registerStateVariableData<Vecd>("Displacement"))
    {
        particles_->addEvolvingVariable<Vecd>("Displacement");
    }

    void update(size_t index_i, Real dt)
    {
        displacement_[index_i] = pos_[index_i] - pos0_[index_i];
    }

  protected:
    Vecd *pos0_;
    Vecd *pos_;
    Vecd *displacement_;
};
//----------------------------------------------------------------------
//	PVD writer for ParaView time series.
//----------------------------------------------------------------------
class PvdWriter
{
  public:
    explicit PvdWriter(const std::string &output_folder, const StdVec<std::string> &body_names)
        : output_folder_(output_folder), body_names_(body_names) {}

    void record(Real physical_time, const std::string &sequence)
    {
        for (const auto &body_name : body_names_)
        {
            records_[body_name].push_back({physical_time, body_name + "_" + sequence + ".vtp"});
            writeBodyPvd(body_name);
        }
    }

    void record(Real physical_time)
    {
        record(physical_time, sequenceFromTime(physical_time));
    }

  private:
    struct RecordEntry
    {
        Real time;
        std::string file;
    };

    std::string output_folder_;
    StdVec<std::string> body_names_;
    std::map<std::string, StdVec<RecordEntry>> records_;

    std::string padValueWithZeros(size_t value, size_t max_string_width = 10)
    {
        std::ostringstream s_time;
        s_time << std::setw(max_string_width) << std::setfill('0') << value;
        return s_time.str();
    }

    std::string sequenceFromTime(Real physical_time)
    {
        size_t i_time = size_t(physical_time * 1.0e6);
        return padValueWithZeros(i_time);
    }

    void writeBodyPvd(const std::string &body_name)
    {
        std::string filefullpath = output_folder_ + "/" + body_name + ".pvd";
        std::ofstream out_file(filefullpath.c_str(), std::ios::trunc);
        out_file << "<?xml version=\"1.0\"?>\n";
        out_file << "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
        out_file << " <Collection>\n";
        for (const auto &entry : records_[body_name])
        {
            out_file << "  <DataSet timestep=\"" << std::fixed << std::setprecision(9) << entry.time
                     << "\" group=\"\" part=\"0\" file=\"" << entry.file << "\"/>\n";
        }
        out_file << " </Collection>\n";
        out_file << "</VTKFile>\n";
        out_file.close();
    }
};
//----------------------------------------------------------------------
//	Unified transport velocity correction (copied from cohesive soil case).
//----------------------------------------------------------------------
template <typename... T>
class TransportVelocityCorrection;

template <class AdaptationType, class LimiterType, typename... CommonControlTypes>
class TransportVelocityCorrection<Inner<AdaptationType, LimiterType>, CommonControlTypes...>
    : public fluid_dynamics::TransportVelocityCorrection<Base, DataDelegateInner, CommonControlTypes...>
{
    using SmoothingLengthRatioType = typename AdaptationType::SmoothingLengthRatioType;

  public:
    explicit TransportVelocityCorrection(BaseInnerRelation &inner_relation, Real coefficient = 0.2)
        : fluid_dynamics::TransportVelocityCorrection<Base, DataDelegateInner, CommonControlTypes...>(inner_relation),
          h_ref_(this->getSPHAdaptation().ReferenceSmoothingLength()),
          correction_scaling_(coefficient * h_ref_ * h_ref_),
          Vol_(this->particles_->template getVariableDataByName<Real>("VolumetricMeasure")),
          pos_div_(this->particles_->template registerStateVariableData<Real>("PositionDivergence")),
          pos_(this->particles_->template getVariableDataByName<Vecd>("Position")),
          h_ratio_(DynamicCast<AdaptationType>(this, this->getSPHAdaptation())), limiter_(h_ref_ * h_ref_),
          indicator_(this->particles_->template registerStateVariableData<int>("Indicator")),
          corner_indicator_(this->particles_->template registerStateVariableData<int>("CornerIndicator")),
          surface_normal_(this->particles_->template registerStateVariableData<Vecd>("SurfaceNormal"))
    {
        static_assert(std::is_base_of<Limiter, LimiterType>::value,
                      "Limiter is not the base of LimiterType!");
    }
    virtual ~TransportVelocityCorrection() {};
    void interaction(size_t index_i, Real dt = 0.0)
    {
        if (this->within_scope_(index_i))
        {
            Vecd inconsistency = Vecd::Zero();
            const Neighborhood &inner_neighborhood = this->inner_configuration_[index_i];
            for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
            {
                size_t index_j = inner_neighborhood.j_[n];
                inconsistency -= (this->kernel_correction_(index_i) + this->kernel_correction_(index_j)) *
                                 inner_neighborhood.dW_ij_[n] * this->Vol_[index_j] * inner_neighborhood.e_ij_[n];
            }
            this->kernel_gradient_integral_[index_i] = inconsistency;
        }
    };
    void update(size_t index_i, Real dt = 0.0)
    {
        if (this->within_scope_(index_i))
        {
            Real inv_h_ratio = 1.0 / h_ratio_(index_i);
            Real squared_norm = this->kernel_gradient_integral_[index_i].squaredNorm();
            Vecd pos_transport = correction_scaling_ * limiter_(squared_norm) *
                                 this->kernel_gradient_integral_[index_i] * inv_h_ratio * inv_h_ratio;
            if (this->indicator_[index_i])
            {
                pos_transport = pos_transport - pos_transport.dot(this->surface_normal_[index_i]) * this->surface_normal_[index_i];
                if (this->pos_div_[index_i] < 0.6 * Dimensions)
                    pos_transport = Vecd::Zero();
            }
            pos_[index_i] += pos_transport;
        }
    };

  protected:
    const Real h_ref_, correction_scaling_;
    Real *Vol_, *pos_div_;
    Vecd *pos_;
    SmoothingLengthRatioType h_ratio_;
    LimiterType limiter_;
    int *indicator_, *corner_indicator_;
    Vecd *surface_normal_;
};

template <class LimiterType, class ParticleScope>
using TransportVelocityCorrectionInner =
    TransportVelocityCorrection<Inner<SPHAdaptation, LimiterType>, NoKernelCorrection, ParticleScope>;

template <typename... CommonControlTypes>
class TransportVelocityCorrection<Contact<Boundary>, CommonControlTypes...>
    : public fluid_dynamics::TransportVelocityCorrection<Base, DataDelegateContact, CommonControlTypes...>
{
  public:
    explicit TransportVelocityCorrection(BaseContactRelation &contact_relation)
        : fluid_dynamics::TransportVelocityCorrection<Base, DataDelegateContact, CommonControlTypes...>(contact_relation)
    {
        for (size_t k = 0; k != this->contact_particles_.size(); ++k)
        {
            wall_Vol_.push_back(this->contact_particles_[k]->template getVariableDataByName<Real>("VolumetricMeasure"));
        }
    };
    virtual ~TransportVelocityCorrection() {};
    void interaction(size_t index_i, Real dt = 0.0)
    {
        if (this->within_scope_(index_i))
        {
            Vecd inconsistency = Vecd::Zero();
            for (size_t k = 0; k < this->contact_configuration_.size(); ++k)
            {
                Real *wall_Vol_k = wall_Vol_[k];
                Neighborhood &contact_neighborhood = (*this->contact_configuration_[k])[index_i];
                for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
                {
                    size_t index_j = contact_neighborhood.j_[n];
                    inconsistency -= 2.0 * this->kernel_correction_(index_i) * contact_neighborhood.dW_ij_[n] *
                                     wall_Vol_k[index_j] * contact_neighborhood.e_ij_[n];
                }
            }
            this->kernel_gradient_integral_[index_i] += inconsistency;
        }
    };

  protected:
    StdVec<Real *> wall_Vol_;
};

template <class AdaptationType, class LimiterType, typename... CommonControlTypes>
using BaseTransportVelocityCorrectionComplex =
    ComplexInteraction<TransportVelocityCorrection<Inner<AdaptationType, LimiterType>, Contact<Boundary>>, CommonControlTypes...>;

template <class ParticleScope>
using TransportVelocityCorrectionComplex =
    BaseTransportVelocityCorrectionComplex<SPHAdaptation, NoLimiter, NoKernelCorrection, ParticleScope>;
//----------------------------------------------------------------------
//	Free surface normal direction (copied from cohesive soil case).
//----------------------------------------------------------------------
template <typename... InteractionTypes>
class FreeSurfaceNormal;

template <class DataDelegationType>
class FreeSurfaceNormal<DataDelegationType>
    : public LocalDynamics, public DataDelegationType
{
  public:
    template <class BaseRelationType>
    explicit FreeSurfaceNormal(BaseRelationType &base_relation)
        : LocalDynamics(base_relation.getSPHBody()), DataDelegationType(base_relation),
          surface_normal_(particles_->registerStateVariableData<Vecd>("SurfaceNormal")),
          color_gradient_(particles_->registerStateVariableData<Vecd>("ColorGradient")),
          B_(particles_->getVariableDataByName<Matd>("LinearGradientCorrectionMatrix")),
          indicator_(particles_->registerStateVariableData<int>("Indicator")),
          Vol_(particles_->getVariableDataByName<Real>("VolumetricMeasure"))
    {
        particles_->addEvolvingVariable<Vecd>("SurfaceNormal");
        particles_->addEvolvingVariable<Vecd>("ColorGradient");
    };
    virtual ~FreeSurfaceNormal() {};

  protected:
    Vecd *surface_normal_, *color_gradient_;
    Matd *B_;
    int *indicator_;
    Real *Vol_;
};

template <>
class FreeSurfaceNormal<Inner<>>
    : public FreeSurfaceNormal<DataDelegateInner>
{
  public:
    explicit FreeSurfaceNormal(BaseInnerRelation &inner_relation) : FreeSurfaceNormal<DataDelegateInner>(inner_relation) {};
    virtual ~FreeSurfaceNormal() {};
    void interaction(size_t index_i, Real dt = 0.0)
    {
        if (indicator_[index_i])
        {
            Vecd color_gradient = ZeroData<Vecd>::value;
            const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
            for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
            {
                size_t index_j = inner_neighborhood.j_[n];
                color_gradient -= inner_neighborhood.dW_ij_[n] * Vol_[index_j] * inner_neighborhood.e_ij_[n];
            }
            color_gradient_[index_i] = color_gradient;
        }
    };
    void update(size_t index_i, Real dt = 0.0)
    {
        if (indicator_[index_i])
        {
            surface_normal_[index_i] = B_[index_i] * color_gradient_[index_i] / (B_[index_i] * color_gradient_[index_i]).norm();
        }
        else
        {
            surface_normal_[index_i] = ZeroData<Vecd>::value;
        }
    };
};
using FreeSurfaceNormalInner = FreeSurfaceNormal<Inner<>>;

template <>
class FreeSurfaceNormal<Contact<>>
    : public FreeSurfaceNormal<DataDelegateContact>
{
  public:
    explicit FreeSurfaceNormal(BaseContactRelation &contact_relation)
        : FreeSurfaceNormal<DataDelegateContact>(contact_relation)
    {
        for (size_t k = 0; k != this->contact_particles_.size(); ++k)
        {
            contact_Vol_.push_back(contact_particles_[k]->getVariableDataByName<Real>("VolumetricMeasure"));
        }
    };
    virtual ~FreeSurfaceNormal() {};
    void interaction(size_t index_i, Real dt = 0.0)
    {
        if (indicator_[index_i])
        {
            for (size_t k = 0; k < contact_configuration_.size(); ++k)
            {
                Vecd color_gradient = ZeroData<Vecd>::value;
                Real *Vol_k = contact_Vol_[k];
                Neighborhood &contact_neighborhood = (*contact_configuration_[k])[index_i];
                for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
                {
                    size_t index_j = contact_neighborhood.j_[n];
                    color_gradient -= contact_neighborhood.dW_ij_[n] * Vol_k[index_j] * contact_neighborhood.e_ij_[n];
                }
                color_gradient_[index_i] += color_gradient;
            }
        }
    };

  protected:
    StdVec<Real *> contact_Vol_;
};
using FreeSurfaceNormalComplex =
    ComplexInteraction<FreeSurfaceNormal<Inner<>, Contact<>>>;
//----------------------------------------------------------------------
//	Erosion and deposition transfer utilities.
//----------------------------------------------------------------------
class SoilToErodedTransfer : public LocalDynamics, public DataDelegateContact
{
  public:
    SoilToErodedTransfer(BaseContactRelation &soil_contact, BaseParticles &eroded_particles)
        : LocalDynamics(soil_contact.getSPHBody()), DataDelegateContact(soil_contact),
          eroded_particles_(eroded_particles),
          indicator_(particles_->getVariableDataByName<int>("Indicator")),
          erosion_state_(particles_->getVariableDataByName<int>("ErosionState")),
          pos_(particles_->getVariableDataByName<Vecd>("Position")),
          vel_(particles_->getVariableDataByName<Vecd>("Velocity")),
          Vol_(particles_->getVariableDataByName<Real>("VolumetricMeasure"))
    {
        for (size_t k = 0; k < contact_particles_.size(); ++k)
        {
            contact_vel_.push_back(contact_particles_[k]->getVariableDataByName<Vecd>("Velocity"));
            contact_rho_.push_back(contact_particles_[k]->getVariableDataByName<Real>("Density"));
        }
    }

    void exec(Real dt)
    {
        Vecd *eroded_pos = eroded_particles_.getVariableDataByName<Vecd>("Position");
        Vecd *eroded_vel = eroded_particles_.getVariableDataByName<Vecd>("Velocity");
        Real *eroded_rho = eroded_particles_.getVariableDataByName<Real>("Density");
        Real *eroded_mass = eroded_particles_.getVariableDataByName<Real>("Mass");
        Real *eroded_Vol = eroded_particles_.getVariableDataByName<Real>("VolumetricMeasure");
        int *eroded_state = eroded_particles_.getVariableDataByName<int>("ErosionState");

        size_t index_i = 0;
        while (index_i < particles_->TotalRealParticles())
        {
            if (indicator_[index_i] && erosion_state_[index_i] == 0)
            {
                Real tau_max = 0.0;
                for (size_t k = 0; k < contact_configuration_.size(); ++k)
                {
                    Vecd *vel_k = contact_vel_[k];
                    Real *rho_k = contact_rho_[k];
                    Neighborhood &contact_neighborhood = (*contact_configuration_[k])[index_i];
                    for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
                    {
                        size_t index_j = contact_neighborhood.j_[n];
                        Vecd v_rel = vel_k[index_j] - vel_[index_i];
                        Real tau = rho_k[index_j] * erosion_coeff * v_rel.squaredNorm();
                        tau_max = SMAX(tau_max, tau);
                    }
                }
                if (tau_max > erosion_tau_crit)
                {
                    UnsignedInt new_index = eroded_particles_.createRealParticleFrom(0);
                    eroded_pos[new_index] = pos_[index_i];
                    eroded_vel[new_index] = vel_[index_i];
                    eroded_rho[new_index] = rho0_s;
                    eroded_Vol[new_index] = Vol_[index_i];
                    eroded_mass[new_index] = rho0_s * Vol_[index_i];
                    eroded_state[new_index] = 1;

                    particles_->switchToBufferParticle(index_i);
                    continue;
                }
            }
            index_i++;
        }
    }

  protected:
    BaseParticles &eroded_particles_;
    int *indicator_;
    int *erosion_state_;
    Vecd *pos_;
    Vecd *vel_;
    Real *Vol_;
    StdVec<Vecd *> contact_vel_;
    StdVec<Real *> contact_rho_;
};

class ErodedToSoilTransfer : public LocalDynamics, public DataDelegateContact
{
  public:
    ErodedToSoilTransfer(BaseContactRelation &eroded_soil_contact, BaseParticles &soil_particles)
        : LocalDynamics(eroded_soil_contact.getSPHBody()), DataDelegateContact(eroded_soil_contact),
          soil_particles_(soil_particles),
          eroded_state_(particles_->getVariableDataByName<int>("ErosionState")),
          pos_(particles_->getVariableDataByName<Vecd>("Position")),
          vel_(particles_->getVariableDataByName<Vecd>("Velocity")),
          Vol_(particles_->getVariableDataByName<Real>("VolumetricMeasure"))
    {}

    void exec(Real dt)
    {
        Vecd *soil_pos = soil_particles_.getVariableDataByName<Vecd>("Position");
        Vecd *soil_vel = soil_particles_.getVariableDataByName<Vecd>("Velocity");
        Real *soil_rho = soil_particles_.getVariableDataByName<Real>("Density");
        Real *soil_mass = soil_particles_.getVariableDataByName<Real>("Mass");
        Real *soil_Vol = soil_particles_.getVariableDataByName<Real>("VolumetricMeasure");
        Mat3d *soil_stress = soil_particles_.getVariableDataByName<Mat3d>("StressTensor3D");
        Vecd *soil_pos0 = soil_particles_.getVariableDataByName<Vecd>("InitialPosition");
        int *soil_state = soil_particles_.getVariableDataByName<int>("ErosionState");

        size_t index_i = 0;
        while (index_i < particles_->TotalRealParticles())
        {
            bool has_soil_contact = false;
            Neighborhood &contact_neighborhood = (*contact_configuration_[0])[index_i];
            if (contact_neighborhood.current_size_ > 0)
            {
                has_soil_contact = true;
            }

            if (has_soil_contact && eroded_state_[index_i] == 1 && vel_[index_i].norm() < deposition_velocity)
            {
                UnsignedInt new_index = soil_particles_.createRealParticleFrom(0);
                soil_pos[new_index] = pos_[index_i];
                soil_vel[new_index] = vel_[index_i];
                soil_rho[new_index] = rho0_s;
                soil_Vol[new_index] = Vol_[index_i];
                soil_mass[new_index] = rho0_s * Vol_[index_i];
                soil_pos0[new_index] = pos_[index_i];
                soil_state[new_index] = 2;
                soil_stress[new_index] = Mat3d::Zero();

                particles_->switchToBufferParticle(index_i);
                continue;
            }
            index_i++;
        }
    }

  protected:
    BaseParticles &soil_particles_;
    int *eroded_state_;
    Vecd *pos_;
    Vecd *vel_;
    Real *Vol_;
};
//----------------------------------------------------------------------
//	Main program starts here.
//----------------------------------------------------------------------
int main(int ac, char *av[])
{
    SPHSystem sph_system(system_domain_bounds, particle_spacing_ref);
    sph_system.handleCommandlineOptions(ac, av);
    //------------------------------------------------------------------
    //	Creating bodies with corresponding materials and particles.
    //------------------------------------------------------------------
    FluidBody water_block(sph_system, makeShared<WaterBlock>("WaterBody"));
    water_block.defineMaterial<WeaklyCompressibleFluid>(rho0_f, c_f);
    water_block.generateParticles<BaseParticles, Lattice>();

    SolidBody wall_boundary(sph_system, makeShared<WallBoundary>("WallBoundary"));
    wall_boundary.defineMaterial<Solid>();
    wall_boundary.generateParticles<BaseParticles, Lattice>();

    ParticleBuffer<ReserveSizeFactor> soil_buffer(0.5);
    RealBody soil_block(sph_system, makeShared<SoilBlock>("SoilBody"));
    soil_block.defineMaterial<PlasticContinuum>(rho0_s, c_s, Youngs_modulus, poisson, friction_angle, cohesion);
    soil_block.generateParticlesWithReserve<BaseParticles, Lattice>(soil_buffer);

    ParticleBuffer<ReserveSizeFactor> eroded_buffer(1.0);
    FluidBody eroded_soil(sph_system, makeShared<SoilBlock>("ErodedSoil"));
    eroded_soil.defineClosure<WeaklyCompressibleFluid, HerschelBulkleyViscosity>(
        ConstructArgs(rho0_s, c_f), ConstructArgs(hb_min_shear_rate, hb_max_shear_rate, hb_consistency_index, hb_power_index, hb_yield_stress));
    eroded_soil.generateParticlesWithReserve<BaseParticles, Lattice>(eroded_buffer);

    BaseParticles &soil_particles = soil_block.getBaseParticles();
    BaseParticles &eroded_particles = eroded_soil.getBaseParticles();
    soil_particles.registerStateVariableData<int>("ErosionState");
    eroded_particles.registerStateVariableData<int>("ErosionState");

    // switch all eroded soil particles to buffer, keep index 0 as template
    while (eroded_particles.TotalRealParticles() > 0)
    {
        eroded_particles.switchToBufferParticle(0);
    }
    int *soil_state = soil_particles.getVariableDataByName<int>("ErosionState");
    for (size_t i = 0; i < soil_particles.TotalRealParticles(); ++i)
    {
        soil_state[i] = 0;
    }
    //------------------------------------------------------------------
    //	Define body relation map.
    //------------------------------------------------------------------
    InnerRelation water_inner(water_block);
    ContactRelation water_eroded_contact(water_block, {&eroded_soil});
    ContactRelation water_wall_contact(water_block, {&wall_boundary, &soil_block});

    InnerRelation eroded_inner(eroded_soil);
    ContactRelation eroded_water_contact(eroded_soil, {&water_block});
    ContactRelation eroded_wall_contact(eroded_soil, {&wall_boundary, &soil_block});
    ContactRelation eroded_soil_contact(eroded_soil, {&soil_block});

    InnerRelation soil_inner(soil_block);
    ContactRelation soil_wall_contact(soil_block, {&wall_boundary});
    ContactRelation soil_fluid_contact(soil_block, {&water_block, &eroded_soil});
    ContactRelation soil_water_contact(soil_block, {&water_block});
    ContactRelation soil_eroded_contact(soil_block, {&eroded_soil});

    //------------------------------------------------------------------
    //	Define the main numerical methods used in the simulation.
    //------------------------------------------------------------------
    Gravity gravity(Vecd(0.0, -gravity_g));
    SimpleDynamics<GravityForce<Gravity>> gravity_to_water(water_block, gravity);
    SimpleDynamics<GravityForce<Gravity>> gravity_to_eroded(eroded_soil, gravity);
    SimpleDynamics<GravityForce<Gravity>> gravity_to_soil(soil_block, gravity);

    SimpleDynamics<NormalDirectionFromBodyShape> wall_boundary_normal_direction(wall_boundary);
    SimpleDynamics<SoilInitialCondition> soil_initial_condition(soil_block);
    SimpleDynamics<SoilDisplacementUpdate> soil_displacement(soil_block);

    InteractionWithUpdate<LinearGradientCorrectionMatrixComplex> correction_matrix(soil_inner, soil_wall_contact);
    Dynamics1Level<continuum_dynamics::PlasticIntegration1stHalfWithWallRiemann> soil_stress_relaxation(soil_inner, soil_wall_contact);
    Dynamics1Level<continuum_dynamics::PlasticIntegration2ndHalfWithWallRiemann> soil_density_relaxation(soil_inner, soil_wall_contact);
    InteractionWithUpdate<fluid_dynamics::DensitySummationComplexFreeSurface> soil_density_by_summation(soil_inner, soil_wall_contact);
    InteractionDynamics<continuum_dynamics::StressDiffusion> soil_stress_diffusion(soil_inner);
    InteractionWithUpdate<FreeSurfaceIndicationComplex> soil_surface_indicator(soil_inner, soil_wall_contact);
    InteractionWithUpdate<TransportVelocityCorrectionComplex<AllParticles>> soil_transport_correction(soil_inner, soil_wall_contact);
    InteractionWithUpdate<FreeSurfaceNormalComplex> soil_free_surface_normal(soil_inner, soil_wall_contact);
    ReduceDynamics<fluid_dynamics::AcousticTimeStep> soil_acoustic_time_step(soil_block, 0.4);

    Dynamics1Level<fluid_dynamics::MultiPhaseIntegration1stHalfWithWallRiemann>
        water_pressure_relaxation(water_inner, water_eroded_contact, water_wall_contact);
    Dynamics1Level<fluid_dynamics::MultiPhaseIntegration2ndHalfWithWallRiemann>
        water_density_relaxation(water_inner, water_eroded_contact, water_wall_contact);

    Dynamics1Level<fluid_dynamics::MultiPhaseIntegration1stHalfWithWallRiemann>
        eroded_pressure_relaxation(eroded_inner, eroded_water_contact, eroded_wall_contact);
    Dynamics1Level<fluid_dynamics::MultiPhaseIntegration2ndHalfWithWallRiemann>
        eroded_density_relaxation(eroded_inner, eroded_water_contact, eroded_wall_contact);

    InteractionWithUpdate<fluid_dynamics::DensitySummationComplexFreeSurface>
        water_density_by_summation(water_inner, water_wall_contact);
    InteractionWithUpdate<fluid_dynamics::BaseDensitySummationComplex<Inner<>, Contact<>, Contact<>>>
        eroded_density_by_summation(eroded_inner, eroded_water_contact, eroded_wall_contact);

    InteractionWithUpdate<fluid_dynamics::VelocityGradientWithWall<NoKernelCorrection>> eroded_velocity_gradient(eroded_inner, eroded_wall_contact);
    SimpleDynamics<fluid_dynamics::ShearRateDependentViscosity> eroded_shear_rate_viscosity(eroded_soil);
    InteractionWithUpdate<fluid_dynamics::NonNewtonianViscousForceWithWall<AngularConservative>> eroded_viscous_acceleration(eroded_inner, eroded_wall_contact);

    InteractionWithUpdate<solid_dynamics::PressureForceFromFluid<decltype(water_density_relaxation)>> water_pressure_on_soil(soil_water_contact);
    InteractionWithUpdate<solid_dynamics::PressureForceFromFluid<decltype(eroded_density_relaxation)>> eroded_pressure_on_soil(soil_eroded_contact);
    solid_dynamics::AverageVelocityAndAcceleration soil_average_velocity(soil_block);

    ReduceDynamics<fluid_dynamics::AdvectionTimeStep> water_advection_time_step(water_block, U_f);
    ReduceDynamics<fluid_dynamics::AdvectionTimeStep> eroded_advection_time_step(eroded_soil, U_f);
    ReduceDynamics<fluid_dynamics::AcousticTimeStep> water_acoustic_time_step(water_block);
    ReduceDynamics<fluid_dynamics::AcousticTimeStep> eroded_acoustic_time_step(eroded_soil);

    SoilToErodedTransfer soil_to_eroded(soil_fluid_contact, eroded_particles);
    ErodedToSoilTransfer eroded_to_soil(eroded_soil_contact, soil_particles);

    //------------------------------------------------------------------
    //	Define the methods for I/O operations.
    //------------------------------------------------------------------
    BodyStatesRecordingToVtp body_states_recording(sph_system);
    body_states_recording.addToWrite<Vecd>(water_block, "Velocity");
    body_states_recording.addToWrite<Vecd>(soil_block, "Velocity");
    body_states_recording.addToWrite<Vecd>(soil_block, "Displacement");
    body_states_recording.addToWrite<Vecd>(eroded_soil, "Velocity");
    body_states_recording.addToWrite<int>(soil_block, "ErosionState");
    body_states_recording.addToWrite<int>(eroded_soil, "ErosionState");
    body_states_recording.addToWrite<Real>(soil_block, "Pressure");
    PvdWriter pvd_writer(sph_system.getIOEnvironment().OutputFolder(),
                         {"WaterBody", "SoilBody", "ErodedSoil"});

    //------------------------------------------------------------------
    //	Prepare the simulation with cell linked list and configuration.
    //------------------------------------------------------------------
    sph_system.initializeSystemCellLinkedLists();
    sph_system.initializeSystemConfigurations();
    wall_boundary_normal_direction.exec();
    gravity_to_water.exec();
    gravity_to_eroded.exec();
    gravity_to_soil.exec();
    soil_initial_condition.exec();
    soil_displacement.exec();
    correction_matrix.exec();

    //------------------------------------------------------------------
    //	Setup for time-stepping control.
    //------------------------------------------------------------------
    Real &physical_time = *sph_system.getSystemVariableDataByName<Real>("PhysicalTime");
    size_t number_of_iterations = 0;
    int screen_output_interval = 200;
    Real End_Time = 2.0;
    Real D_Time = End_Time / 50.0;

    //------------------------------------------------------------------
    //	First output.
    //------------------------------------------------------------------
    body_states_recording.writeToFile();
    pvd_writer.record(0.0);

    //------------------------------------------------------------------
    //	Main loop starts here.
    //------------------------------------------------------------------
    while (physical_time < End_Time)
    {
        Real integration_time = 0.0;
        while (integration_time < D_Time)
        {
            Real Dt = SMIN(water_advection_time_step.exec(), eroded_advection_time_step.exec());

            water_density_by_summation.exec();
            eroded_density_by_summation.exec();
            soil_density_by_summation.exec();
            soil_surface_indicator.exec();
            soil_free_surface_normal.exec();
            soil_transport_correction.exec();

            Real relaxation_time = 0.0;
            while (relaxation_time < Dt)
            {
                Real dt = SMIN(SMIN(water_acoustic_time_step.exec(), eroded_acoustic_time_step.exec()), soil_acoustic_time_step.exec());

                water_pressure_relaxation.exec(dt);
                eroded_pressure_relaxation.exec(dt);

                water_density_relaxation.exec(dt);
                eroded_density_relaxation.exec(dt);

                eroded_velocity_gradient.exec();
                eroded_shear_rate_viscosity.exec();
                eroded_viscous_acceleration.exec();

                soil_stress_diffusion.exec();
                water_pressure_on_soil.exec();
                eroded_pressure_on_soil.exec();
                soil_stress_relaxation.exec(dt);
                soil_density_relaxation.exec(dt);
                soil_displacement.exec(dt);

                soil_to_eroded.exec(dt);
                eroded_to_soil.exec(dt);

                relaxation_time += dt;
                integration_time += dt;
                physical_time += dt;
            }

            if (number_of_iterations % screen_output_interval == 0)
            {
                std::cout << std::fixed << std::setprecision(9) << "N=" << number_of_iterations << "\tTime = "
                          << physical_time << "\tDt = " << Dt << "\n";
            }
            number_of_iterations++;

            water_block.updateCellLinkedList();
            eroded_soil.updateCellLinkedList();
            soil_block.updateCellLinkedList();
            water_eroded_contact.updateConfiguration();
            water_wall_contact.updateConfiguration();
            eroded_water_contact.updateConfiguration();
            eroded_wall_contact.updateConfiguration();
            eroded_soil_contact.updateConfiguration();
            soil_fluid_contact.updateConfiguration();
            soil_wall_contact.updateConfiguration();
            soil_eroded_contact.updateConfiguration();
            correction_matrix.exec();
        }
        body_states_recording.writeToFile();
        pvd_writer.record(physical_time);
    }

    return 0;
}
