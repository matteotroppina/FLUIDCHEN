#pragma once

#include <vector>

#include "Cell.hpp"
#include "Fields.hpp"
#include <iostream>
#include <map>

/**
 * @brief Abstract of boundary conditions.
 *
 * This class patches the physical values to the given field.
 */
class Boundary {
  public:
    /**
     * @brief Method to patch the velocity boundary conditions to the given field.
     *
     * @param[in] Field to be applied
     */
    virtual void applyVelocity(Fields &field) = 0;

    /**
     * @brief Method to patch the pressure boundary conditions to the given field.
     *
     * @param[in] Field to be applied
     */
    virtual void applyPressure(Fields &field) = 0;

    /**
     * @brief Method to patch the flux (F & G) boundary conditions to the given field.
     *
     * @param[in] Field to be applied
     */
    virtual void applyFlux(Fields &field);

    virtual void applyTemperature(Fields &field);

    /**
     * @brief Method to patch the K boundary conditions to the given field (K-epsilon turbulence model).
     *
     * @param[in] Field to be applied
     */
    virtual void applyK(Fields &field);

    /**
     * @brief Method to patch the Epsilon boundary conditions to the given field (K-eosilon turbulence model).
     *
     * @param[in] Field to be applied
     */
    virtual void applyEpsilon(Fields &field);

    /**
     * @brief Method to patch the NuT boundary conditions to the given field (K-eosilon turbulence model).
     *
     * @param[in] Field to be applied
     */
    virtual void applyTurbulence(Fields &field);

    virtual ~Boundary() = default;

  protected:
    Boundary(std::vector<Cell *> cells);
    std::vector<Cell *> _cells;
    /// C_nu
    double C0{0.09};
    double l0{0.1};
};

/**
 * @brief Fixed wall boundary condition for the outer boundaries of the domain.
 * Dirichlet for velocities, which is zero, Neumann for pressure
 */
class FixedWallBoundary : public Boundary {
  public:
    FixedWallBoundary(std::vector<Cell *> cells);
    FixedWallBoundary(std::vector<Cell *> cells, std::map<int, double> wall_temperature);
    FixedWallBoundary(std::vector<Cell *> cells, double wall_temperature);
    virtual ~FixedWallBoundary() = default;
    virtual void applyVelocity(Fields &field);
    virtual void applyPressure(Fields &field);
    void applyTemperature(Fields &field);

    void applyK(Fields &field);
    void applyEpsilon(Fields &field);

  private:
        double _wall_temperature;

};

class InnerObstacle : public Boundary {
  public:
    InnerObstacle(std::vector<Cell *> cells);
    virtual ~InnerObstacle() = default;
    virtual void applyVelocity(Fields &field); // do nothing
    virtual void applyPressure(Fields &field); // do nothing
    virtual void applyFlux(Fields &field); // do nothing

    // virtual void applyK(Fields &field);
    // virtual void applyEpsilon(Fields &field);
};

/**
 * @brief Moving wall boundary condition for the outer boundaries of the domain.
 * Dirichlet for velocities for the given velocity parallel to the fluid,
 * Neumann for pressure
 */
class MovingWallBoundary : public Boundary {
  public:
    MovingWallBoundary(std::vector<Cell *> cells, double wall_velocity);
    MovingWallBoundary(std::vector<Cell *> cells, std::map<int, double> wall_velocity,
                       std::map<int, double> wall_temperature);
    virtual ~MovingWallBoundary() = default;
    virtual void applyVelocity(Fields &field);
    virtual void applyPressure(Fields &field);

    // virtual void applyK(Fields &field);
    // virtual void applyEpsilon(Fields &field);
    
  private:
    std::map<int, double> _wall_velocity;
    std::map<int, double> _wall_temperature;
};

class FixedVelocityBoundary : public Boundary {
  public:
    FixedVelocityBoundary(std::vector<Cell *> cells, double inflow_u_velocity, double inflow_v_velocity);
    FixedVelocityBoundary(std::vector<Cell *> cells, std::map<int, double> inflow_u_velocity, std::map<int, double> inflow_v_velocity,
                       std::map<int, double> wall_temperature);
    virtual ~FixedVelocityBoundary() = default;
    virtual void applyVelocity(Fields &field);
    virtual void applyPressure(Fields &field);

    virtual void applyK(Fields &field);
    // virtual void applyEpsilon(Fields &field);

  private:
    std::map<int, double> _inflow_u_velocity;
    std::map<int, double> _inflow_v_velocity;
    std::map<int, double> _wall_temperature;
};

class ZeroGradientBoundary : public Boundary {
  public:
    ZeroGradientBoundary(std::vector<Cell *> cells);
    ZeroGradientBoundary(std::vector<Cell *> cells, std::map<int, double> wall_temperature);
    virtual ~ZeroGradientBoundary() = default;
    virtual void applyVelocity(Fields &field);
    virtual void applyPressure(Fields &field);

    virtual void applyK(Fields &field);
    // virtual void applyEpsilon(Fields &field);

  private:
    std::map<int, double> _wall_temperature;
};

