/**
 * @file 	cohesive_soil_failure.h
 * @brief 	2D cohesive soil failure.
 * @details This case employs the SPH to simulate cohesive granular materials.
 *          The transport velocity formulation is adopted to address the tensile instability.
 * @author Shuaihao Zhang
 */
#include "sphinxsys.h" //SPHinXsys Library.
using namespace SPH;   // Namespace cite here.
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup (对应论文 4.8 节泥沙侵蚀案例).
//----------------------------------------------------------------------
Real Ls = 2.0;   /**< 砂土床长度 Ls (m). */
Real Hs = 0.06;  /**< 砂土床高度 Hs (m). */
Real Lw = 1.0;   /**< 水柱长度 Lw (m). */
Real Hw = 0.10;  /**< 水柱高度 Hw (m). */
Real LL = Ls;    /**< 仍沿用原变量名，表示土体长度。 */
Real LH = Hs;    /**< 仍沿用原变量名，表示土体高度。 */
Real DL = Ls;    /**< 计算域长度与土床长度一致。 */
Real DH = Hw + Hs + 0.10; /**< 计算域高度，预留自由表面空间。 */
Real particle_spacing_ref = 0.0035; /**< 初始粒子间距 dp (m). */
Real BW = particle_spacing_ref * 4; /**< 边界扩展宽度。 */
BoundingBoxd system_domain_bounds(Vec2d(-BW, -BW), Vec2d(DL + BW, DH + BW));
//----------------------------------------------------------------------
//	Material properties of the soil (按表 4.1 设置).
//----------------------------------------------------------------------
Real rho0_s = 1540.0;                                                   // 土体密度 rho_s (kg/m^3)
Real gravity_g = 9.8;                                                   // 重力加速度 (m/s^2)
Real Youngs_modulus = 1.0e6;                                            // 弹性模量 (文献未给出，保持合理数值)
Real poisson = 0.2;                                                     // 泊松比 e
Real c_s = sqrt(Youngs_modulus / (rho0_s * 3.0 * (1.0 - 2.0 * poisson))); // 声速
Real cohesion = 0.0;                                                    // 内聚力 C (kPa -> Pa，这里取 0)
Real friction_angle = 31.0 * Pi / 180.0;                                // 内摩擦角 φ (°)
// Herschel–Bulkley–Papanastasiou (HBP) 模型参数（泥浆/泥流，见表 4.1）
Real hbp_yield_stress = 0.0;  // 屈服应力 τ0 (Pa)，泥浆稀薄时可取 0
Real hbp_consistency = 0.001; // 黏度系数 Vs (Pa·s^n)
Real hbp_flow_index = 1.0;    // 幂律指数 n
Real hbp_regularization = 0.0; // Papanastasiou 正则化参数 m (文献给出 0)
Real d50 = 0.0035;             // 中值粒径 d50 (m)
Real erosion_velocity_eps = 1.0e-4; // 沉积判据中的速度阈值 (m/s)
// 侵蚀/淤积状态阈值（用于 DP/HBP 切换，需根据实验标定）
Real erosion_shear_rate = 2.0;    // 侵蚀触发剪切率 (1/s)
Real deposition_shear_rate = 0.5; // 回沉/淤积剪切率 (1/s)
//----------------------------------------------------------------------
//	Material properties of water (按图示参数设置).
//----------------------------------------------------------------------
Real rho0_f = 1000.0;                                     // 水密度 rho_w (kg/m^3)
Real mu_f = 0.001;                                        // 动力黏度 μ (Pa·s)
Real c_f = 10.0 * sqrt(gravity_g * (Hw + Hs));            // 人工声速，保证弱可压缩条件
Real U_f = 1.0;                                           // 参考速度，用于时间步估计
Real water_soil_force_scale = 1.0;                        // 水-土力耦合系数
Real erosion_drag_coeff = 5.0 * rho0_f;                   // 侵蚀粒子随流拖曳系数
Real water_soil_repulsion_strength = 5.0e3;               // 非侵蚀土体水体惩罚反力系数
Real water_soil_damping = 5.0;                            // 非侵蚀土体水体阻尼系数
Real erosion_velocity_relaxation = 1.0;                   // 侵蚀粒子速度松弛系数(0-1)
//----------------------------------------------------------------------
//	Geometric shapes used in this case.
//----------------------------------------------------------------------
Vec2d soil_block_halfsize = Vec2d(0.5 * LL, 0.5 * LH); // 土体几何尺寸
Vec2d soil_block_translation = soil_block_halfsize;   // 土体左下角位于 (0,0)
Vec2d outer_wall_halfsize = Vec2d(0.5 * DL + BW, 0.5 * DH + BW);
Vec2d outer_wall_translation = Vec2d(-BW, -BW) + outer_wall_halfsize;
Vec2d inner_wall_halfsize = Vec2d(0.5 * DL, 0.5 * DH);
Vec2d inner_wall_translation = inner_wall_halfsize;
//----------------------------------------------------------------------
//	Complex for wall boundary
//----------------------------------------------------------------------
class WallBoundary : public ComplexShape
{
  public:
    explicit WallBoundary(const std::string &shape_name) : ComplexShape(shape_name)
    {
        add<GeometricShapeBox>(Transform(outer_wall_translation), outer_wall_halfsize);
        subtract<GeometricShapeBox>(Transform(inner_wall_translation), inner_wall_halfsize);
    }
};
std::vector<Vecd> soil_shape{
    Vecd(0, 0), Vecd(0, LH), Vecd(LL, LH), Vecd(LL, 0), Vecd(0, 0)};
std::vector<Vecd> water_shape{
    Vecd(0, LH), Vecd(0, LH + Hw), Vecd(Lw, LH + Hw), Vecd(Lw, LH), Vecd(0, LH)};

class Soil : public MultiPolygonShape
{
  public:
    explicit Soil(const std::string &shape_name) : MultiPolygonShape(shape_name)
    {
        multi_polygon_.addAPolygon(soil_shape, ShapeBooleanOps::add);
    }
};

class WaterBlock : public MultiPolygonShape
{
  public:
    explicit WaterBlock(const std::string &shape_name) : MultiPolygonShape(shape_name)
    {
        multi_polygon_.addAPolygon(water_shape, ShapeBooleanOps::add);
    }
};
//----------------------------------------------------------------------
//	未侵蚀土体粒子集合：作为水体动态壁面
//----------------------------------------------------------------------
class NonErodedSoilPart : public BodyPartByParticle
{
  public:
    explicit NonErodedSoilPart(RealBody &soil_body)
        : BodyPartByParticle(soil_body),
          erosion_state_(soil_body.getBaseParticles().registerStateVariableData<int>("ErosionState"))
    {
        TaggingParticleMethod tagging_particle_method =
            std::bind(&NonErodedSoilPart::tagByErosionState, this, std::placeholders::_1);
        tagParticles(tagging_particle_method);
    }

    bool tagByErosionState(size_t particle_index)
    {
        return erosion_state_[particle_index] == 0;
    }

    void updateTags()
    {
        TaggingParticleMethod tagging_particle_method =
            std::bind(&NonErodedSoilPart::tagByErosionState, this, std::placeholders::_1);
        tagParticles(tagging_particle_method);
    }

  protected:
    int *erosion_state_;
};

//---------------------------------------------------------------------- 
//	侵蚀土体粒子集合：作为水体可穿透区域
//----------------------------------------------------------------------
class ErodedSoilPart : public BodyPartByParticle
{
  public:
    explicit ErodedSoilPart(RealBody &soil_body)
        : BodyPartByParticle(soil_body),
          erosion_state_(soil_body.getBaseParticles().registerStateVariableData<int>("ErosionState"))
    {
        TaggingParticleMethod tagging_particle_method =
            std::bind(&ErodedSoilPart::tagByErosionState, this, std::placeholders::_1);
        tagParticles(tagging_particle_method);
    }

    bool tagByErosionState(size_t particle_index)
    {
        return erosion_state_[particle_index] == 1;
    }

    void updateTags()
    {
        TaggingParticleMethod tagging_particle_method =
            std::bind(&ErodedSoilPart::tagByErosionState, this, std::placeholders::_1);
        tagParticles(tagging_particle_method);
    }

  protected:
    int *erosion_state_;
};
//----------------------------------------------------------------------
//	application dependent initial condition
//----------------------------------------------------------------------
class SoilInitialCondition : public continuum_dynamics::ContinuumInitialCondition
{
  public:
    explicit SoilInitialCondition(RealBody &granular_column)
        : continuum_dynamics::ContinuumInitialCondition(granular_column) {};

  protected:
    void update(size_t index_i, Real dt)
    {
        /** initial stress */
        Real y = pos_[index_i][1];
        Real gama = 1 - sin(friction_angle);
        Real stress_yy = -rho0_s * gravity_g * y;
        stress_tensor_3D_[index_i](1, 1) = stress_yy;
        stress_tensor_3D_[index_i](0, 0) = stress_yy * gama;
        stress_tensor_3D_[index_i](2, 2) = stress_yy * gama;
    };
};
//----------------------------------------------------------------------
//	Water-soil interaction forces (scheme A: pressure + viscous coupling)
//	公式：F_p = -∑ p_j ∇W_ij V_j, F_v = 2μ (v_j - v_i)/(r_ij + 0.01h) ∇W_ij V_j
//	      其中 ∇W_ij = dW_ij * e_ij, 最终力乘以粒子体积作为 ForcePrior.
//----------------------------------------------------------------------
class SoilForceFromWater : public ForcePrior, public DataDelegateContact
{
  public:
    explicit SoilForceFromWater(BaseContactRelation &contact_relation)
        : ForcePrior(contact_relation.getSPHBody(), "SoilWaterForce"), DataDelegateContact(contact_relation),
          vel_(particles_->registerStateVariableData<Vecd>("Velocity")),
          erosion_state_(particles_->registerStateVariableData<int>("ErosionState")),
          Vol_(particles_->getVariableDataByName<Real>("VolumetricMeasure"))
    {
        for (size_t k = 0; k != contact_particles_.size(); ++k)
        {
            contact_vel_.push_back(contact_particles_[k]->registerStateVariableData<Vecd>("Velocity"));
            contact_p_.push_back(contact_particles_[k]->registerStateVariableData<Real>("Pressure"));
            contact_Vol_.push_back(contact_particles_[k]->getVariableDataByName<Real>("VolumetricMeasure"));
            smoothing_length_.push_back(contact_bodies_[k]->getSPHAdaptation().ReferenceSmoothingLength());
        }
    }

    void update(size_t index_i, Real dt)
    {
        Vecd force = Vecd::Zero();
        for (size_t k = 0; k < contact_configuration_.size(); ++k)
        {
            Vecd *vel_k = contact_vel_[k];
            Real *p_k = contact_p_[k];
            Real *Vol_k = contact_Vol_[k];
            Real smoothing_length_k = smoothing_length_[k];
            Neighborhood &contact_neighborhood = (*contact_configuration_[k])[index_i];
            for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
            {
                size_t index_j = contact_neighborhood.j_[n];
                Vecd gradW = contact_neighborhood.dW_ij_[n] * contact_neighborhood.e_ij_[n];
                Real r_ij = contact_neighborhood.r_ij_[n];
                force -= p_k[index_j] * gradW * Vol_k[index_j];
                Vecd vel_derivative = 2.0 * (vel_k[index_j] - vel_[index_i]) / (r_ij + 0.01 * smoothing_length_k);
                force += 2.0 * mu_f * vel_derivative * contact_neighborhood.dW_ij_[n] * Vol_k[index_j];
                if (erosion_state_[index_i] == 1)
                {
                    force += erosion_drag_coeff * (vel_k[index_j] - vel_[index_i]) *
                             contact_neighborhood.W_ij_[n] * Vol_k[index_j];
                }
            }
        }
        current_force_[index_i] = water_soil_force_scale * force * Vol_[index_i];
        ForcePrior::update(index_i, dt);
    }

  protected:
    Vecd *vel_;
    int *erosion_state_;
    Real *Vol_;
    StdVec<Vecd *> contact_vel_;
    StdVec<Real *> contact_p_;
    StdVec<Real *> contact_Vol_;
    StdVec<Real> smoothing_length_;
};

class WaterForceFromSoil : public ForcePrior, public DataDelegateContact
{
  public:
    explicit WaterForceFromSoil(BaseContactRelation &contact_relation)
        : ForcePrior(contact_relation.getSPHBody(), "WaterSoilForce"), DataDelegateContact(contact_relation),
          vel_(particles_->registerStateVariableData<Vecd>("Velocity")),
          Vol_(particles_->getVariableDataByName<Real>("VolumetricMeasure"))
    {
        for (size_t k = 0; k != contact_particles_.size(); ++k)
        {
            contact_vel_.push_back(contact_particles_[k]->registerStateVariableData<Vecd>("Velocity"));
            contact_p_.push_back(contact_particles_[k]->registerStateVariableData<Real>("Pressure"));
            contact_Vol_.push_back(contact_particles_[k]->getVariableDataByName<Real>("VolumetricMeasure"));
            contact_erosion_state_.push_back(contact_particles_[k]->getVariableDataByName<int>("ErosionState"));
            smoothing_length_.push_back(contact_bodies_[k]->getSPHAdaptation().ReferenceSmoothingLength());
        }
    }

    void update(size_t index_i, Real dt)
    {
        Vecd force = Vecd::Zero();
        for (size_t k = 0; k < contact_configuration_.size(); ++k)
        {
            Vecd *vel_k = contact_vel_[k];
            Real *p_k = contact_p_[k];
            Real *Vol_k = contact_Vol_[k];
            int *erosion_state_k = contact_erosion_state_[k];
            Real smoothing_length_k = smoothing_length_[k];
            Neighborhood &contact_neighborhood = (*contact_configuration_[k])[index_i];
            for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
            {
                size_t index_j = contact_neighborhood.j_[n];
                Vecd gradW = contact_neighborhood.dW_ij_[n] * contact_neighborhood.e_ij_[n];
                Real r_ij = contact_neighborhood.r_ij_[n];
                force -= p_k[index_j] * gradW * Vol_k[index_j];
                Vecd vel_derivative = 2.0 * (vel_k[index_j] - vel_[index_i]) / (r_ij + 0.01 * smoothing_length_k);
                force += 2.0 * mu_f * vel_derivative * contact_neighborhood.dW_ij_[n] * Vol_k[index_j];
            }
        }
        current_force_[index_i] = water_soil_force_scale * force * Vol_[index_i];
        ForcePrior::update(index_i, dt);
    }

  protected:
    Vecd *vel_;
    Real *Vol_;
    StdVec<Vecd *> contact_vel_;
    StdVec<Real *> contact_p_;
    StdVec<Real *> contact_Vol_;
    StdVec<int *> contact_erosion_state_;
    StdVec<Real> smoothing_length_;
};
//----------------------------------------------------------------------
//	水体-未侵蚀土体的动态壁面处理（近似 Contact<Wall>）
//	仅当土体粒子为未侵蚀状态（ErosionState==0）时施加边界压力贡献
//----------------------------------------------------------------------
class WaterWallBoundaryFromSoil : public LocalDynamics, public DataDelegateContact
{
  public:
    explicit WaterWallBoundaryFromSoil(BaseContactRelation &contact_relation)
        : LocalDynamics(contact_relation.getSPHBody()), DataDelegateContact(contact_relation),
          rho_(particles_->registerStateVariableData<Real>("Density")),
          p_(particles_->registerStateVariableData<Real>("Pressure")),
          mass_(particles_->registerStateVariableData<Real>("Mass")),
          force_(particles_->registerStateVariableData<Vecd>("Force")),
          force_prior_(particles_->registerStateVariableData<Vecd>("ForcePrior")),
          drho_dt_(particles_->registerStateVariableData<Real>("DensityChangeRate"))
    {
        for (size_t k = 0; k != contact_particles_.size(); ++k)
        {
            wall_Vol_.push_back(contact_particles_[k]->getVariableDataByName<Real>("VolumetricMeasure"));
        }
    }

    void interaction(size_t index_i, Real dt = 0.0)
    {
        Vecd force = Vecd::Zero();
        Real rho_dissipation(0);
        for (size_t k = 0; k < contact_configuration_.size(); ++k)
        {
            Real *wall_Vol_k = wall_Vol_[k];
            Neighborhood &wall_neighborhood = (*contact_configuration_[k])[index_i];
            for (size_t n = 0; n != wall_neighborhood.current_size_; ++n)
            {
                size_t index_j = wall_neighborhood.j_[n];
                Vecd &e_ij = wall_neighborhood.e_ij_[n];
                Real dW_ijV_j = wall_neighborhood.dW_ij_[n] * wall_Vol_k[index_j];
                Real r_ij = wall_neighborhood.r_ij_[n];
                Real face_wall_external_acceleration = (force_prior_[index_i] / mass_[index_i]).dot(-e_ij);
                Real p_j_in_wall = p_[index_i] + rho_[index_i] * r_ij * SMAX(Real(0), face_wall_external_acceleration);
                force -= (p_[index_i] + p_j_in_wall) * dW_ijV_j * e_ij;
                rho_dissipation += 0.0;
            }
        }
        force_[index_i] += force * particles_->getVariableDataByName<Real>("VolumetricMeasure")[index_i];
        drho_dt_[index_i] += rho_dissipation * rho_[index_i];
    }

  protected:
    Real *rho_, *p_, *mass_, *drho_dt_;
    Vecd *force_, *force_prior_;
    StdVec<Real *> wall_Vol_;
};
//----------------------------------------------------------------------
//	水体-未侵蚀土体的不可穿透处理（惩罚反力）
//	用于阻止水体穿透非侵蚀土体（ErosionState==0）
//----------------------------------------------------------------------
class WaterRepulsionFromSoil : public LocalDynamics, public DataDelegateContact
{
  public:
    explicit WaterRepulsionFromSoil(BaseContactRelation &contact_relation)
        : LocalDynamics(contact_relation.getSPHBody()), DataDelegateContact(contact_relation),
          vel_(particles_->registerStateVariableData<Vecd>("Velocity")),
          force_(particles_->registerStateVariableData<Vecd>("Force"))
    {
        for (size_t k = 0; k != contact_particles_.size(); ++k)
        {
            contact_Vol_.push_back(contact_particles_[k]->getVariableDataByName<Real>("VolumetricMeasure"));
            contact_vel_.push_back(contact_particles_[k]->registerStateVariableData<Vecd>("Velocity"));
        }
    }

    void interaction(size_t index_i, Real dt = 0.0)
    {
        Vecd force = Vecd::Zero();
        for (size_t k = 0; k < contact_configuration_.size(); ++k)
        {
            Real *soil_Vol_k = contact_Vol_[k];
            Vecd *soil_vel_k = contact_vel_[k];
            Neighborhood &soil_neighborhood = (*contact_configuration_[k])[index_i];
            for (size_t n = 0; n != soil_neighborhood.current_size_; ++n)
            {
                size_t index_j = soil_neighborhood.j_[n];
                Vecd &e_ij = soil_neighborhood.e_ij_[n];
                Real r_ij = soil_neighborhood.r_ij_[n];
                Real penetration = particle_spacing_ref - r_ij;
                if (penetration > 0.0)
                {
                    Real normal_speed = (vel_[index_i] - soil_vel_k[index_j]).dot(e_ij);
                    Real penalty = water_soil_repulsion_strength * penetration;
                    force += (penalty - water_soil_damping * SMIN(Real(0), normal_speed)) *
                             e_ij * soil_Vol_k[index_j];
                }
            }
        }
        force_[index_i] += force;
    }

  protected:
    Vecd *vel_;
    Vecd *force_;
    StdVec<Real *> contact_Vol_;
    StdVec<Vecd *> contact_vel_;
};
//----------------------------------------------------------------------
//	侵蚀判据与沉积判据（对应论文 4.4.2 / 4.4.3）
//	交界面判定：2h 支持域内是否存在水粒子
//	流速估计：v_i = sum(v_j W_ij V_j) / sum(W_ij V_j)
//	侵蚀判据：E = v_i - v_ic, 其中 v_ic 按 d50 分段（式 4.3/4.4）
//	沉积判据：D > dp 且 E <= 0 且 v_is = 0（式 4.10-4.12）
//----------------------------------------------------------------------
class ErosionStateByVelocity : public LocalDynamics, public DataDelegateContact
{
  public:
    explicit ErosionStateByVelocity(BaseContactRelation &contact_relation)
        : LocalDynamics(contact_relation.getSPHBody()), DataDelegateContact(contact_relation),
          pos_(particles_->getVariableDataByName<Vecd>("Position")),
          vel_(particles_->registerStateVariableData<Vecd>("Velocity")),
          erosion_state_(particles_->registerStateVariableData<int>("ErosionState")),
          erosion_start_pos_(particles_->registerStateVariableData<Vecd>("ErosionStartPosition")),
          interface_indicator_(particles_->registerStateVariableData<int>("InterfaceIndicator"))
    {
        particles_->addEvolvingVariable<Vecd>("ErosionStartPosition");
        particles_->addEvolvingVariable<int>("InterfaceIndicator");
        for (size_t k = 0; k != contact_particles_.size(); ++k)
        {
            contact_vel_.push_back(contact_particles_[k]->registerStateVariableData<Vecd>("Velocity"));
            contact_Vol_.push_back(contact_particles_[k]->getVariableDataByName<Real>("VolumetricMeasure"));
        }
    }

    void update(size_t index_i, Real dt)
    {
        Real weight_sum = 0.0;
        Vecd velocity_sum = Vecd::Zero();
        interface_indicator_[index_i] = 0;
        for (size_t k = 0; k < contact_configuration_.size(); ++k)
        {
            Vecd *vel_k = contact_vel_[k];
            Real *Vol_k = contact_Vol_[k];
            Neighborhood &contact_neighborhood = (*contact_configuration_[k])[index_i];
            if (contact_neighborhood.current_size_ > 0)
                interface_indicator_[index_i] = 1;
            for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
            {
                size_t index_j = contact_neighborhood.j_[n];
                Real W_ij = contact_neighborhood.W_ij_[n];
                Real weight = W_ij * Vol_k[index_j];
                weight_sum += weight;
                velocity_sum += vel_k[index_j] * weight;
            }
        }

        if (!interface_indicator_[index_i])
        {
            erosion_state_[index_i] = 0;
            return;
        }

        Real v_i = (weight_sum > TinyReal) ? velocity_sum.norm() / weight_sum : 0.0;
        Real d50_mm = d50 * 1000.0;
        Real v_ic = (d50_mm < 0.1) ? 0.1 * pow(d50_mm, -0.2) : 0.35 * pow(d50_mm, 0.45);
        Real E = v_i - v_ic;

        if (E > 0.0)
        {
            if (erosion_state_[index_i] == 0)
                erosion_start_pos_[index_i] = pos_[index_i];
            erosion_state_[index_i] = 1;
        }
        else if (erosion_state_[index_i] == 1)
        {
            Real D = (pos_[index_i] - erosion_start_pos_[index_i]).norm();
            if (D > particle_spacing_ref && vel_[index_i].norm() <= erosion_velocity_eps)
            {
                erosion_state_[index_i] = 0;
            }
        }
    }

  protected:
    Vecd *pos_, *vel_;
    int *erosion_state_, *interface_indicator_;
    Vecd *erosion_start_pos_;
    StdVec<Vecd *> contact_vel_;
    StdVec<Real *> contact_Vol_;
};
//----------------------------------------------------------------------
//	侵蚀土体随水流运动（速度松弛）
//	侵蚀粒子速度向水体局部平均速度收敛
//----------------------------------------------------------------------
class ErodedSoilVelocityRelaxation : public LocalDynamics, public DataDelegateContact
{
  public:
    explicit ErodedSoilVelocityRelaxation(BaseContactRelation &contact_relation)
        : LocalDynamics(contact_relation.getSPHBody()), DataDelegateContact(contact_relation),
          vel_(particles_->registerStateVariableData<Vecd>("Velocity")),
          erosion_state_(particles_->registerStateVariableData<int>("ErosionState"))
    {
        for (size_t k = 0; k != contact_particles_.size(); ++k)
        {
            contact_vel_.push_back(contact_particles_[k]->registerStateVariableData<Vecd>("Velocity"));
            contact_Vol_.push_back(contact_particles_[k]->getVariableDataByName<Real>("VolumetricMeasure"));
        }
    }

    void update(size_t index_i, Real dt)
    {
        if (erosion_state_[index_i] == 0)
            return;

        Real weight_sum = 0.0;
        Vecd velocity_sum = Vecd::Zero();
        for (size_t k = 0; k < contact_configuration_.size(); ++k)
        {
            Vecd *vel_k = contact_vel_[k];
            Real *Vol_k = contact_Vol_[k];
            Neighborhood &contact_neighborhood = (*contact_configuration_[k])[index_i];
            for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
            {
                size_t index_j = contact_neighborhood.j_[n];
                Real weight = contact_neighborhood.W_ij_[n] * Vol_k[index_j];
                weight_sum += weight;
                velocity_sum += vel_k[index_j] * weight;
            }
        }

        if (weight_sum > TinyReal)
        {
            Vecd target_vel = velocity_sum / weight_sum;
            vel_[index_i] = (1.0 - erosion_velocity_relaxation) * vel_[index_i] +
                            erosion_velocity_relaxation * target_vel;
        }
    }

  protected:
    Vecd *vel_;
    int *erosion_state_;
    StdVec<Vecd *> contact_vel_;
    StdVec<Real *> contact_Vol_;
};
//----------------------------------------------------------------------
//	Unified transport velocity correction
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
                // acceleration for transport velocity
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
                    // acceleration for transport velocity
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
//	Free surface normal direction
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
