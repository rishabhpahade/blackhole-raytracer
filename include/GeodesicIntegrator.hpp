#ifndef BLACKHOLE_GEODESICINTEGRATOR_HPP
#define BLACKHOLE_GEODESICINTEGRATOR_HPP

#include "Vector3.hpp"
#include "Ray.hpp"
#include "BlackHole.hpp"
#include "AccretionDisk.hpp"

/**
 * @brief Evaluates 6D state derivatives and integrates null geodesics using
 * an explicit 4th-Order Runge-Kutta (RK4) scheme with curvature-adaptive stepping.
 */
class GeodesicIntegrator {
public:
    BlackHole blackHole;
    AccretionDisk disk;

    double initialStepSize{0.06};
    double minStepSize{1e-5};
    double maxStepSize{2.0};
    double escapeRadius{60.0};
    int maxSteps{15000};
    bool testAccretionDisk{true};

    GeodesicIntegrator() noexcept = default;
    GeodesicIntegrator(const BlackHole& bh, const AccretionDisk& d,
                       double step = 0.06, double escape = 60.0, int maxIter = 15000) noexcept
        : blackHole(bh), disk(d), initialStepSize(step), escapeRadius(escape), maxSteps(maxIter) {}

    /**
     * @brief Structure representing the 6D phase space state: position and spatial velocity.
     */
    struct State6D {
        Vector3 pos;
        Vector3 vel;

        constexpr State6D operator+(const State6D& rhs) const noexcept {
            return {pos + rhs.pos, vel + rhs.vel};
        }

        constexpr State6D operator*(double scalar) const noexcept {
            return {pos * scalar, vel * scalar};
        }
    };

    /**
     * @brief Computes phase-space derivative d/dλ [x, v] = [v, a(x, v)].
     * 
     * In Schwarzschild spacetime, the Cartesian null-geodesic acceleration is:
     * a = - 3 * M * |x × v|^2 / |x|^5 * x
     */
    [[nodiscard]] State6D evaluateDerivative(const State6D& state) const noexcept;

    /**
     * @brief Advances a 6D state by step size h using standard Runge-Kutta 4th Order.
     */
    [[nodiscard]] State6D rk4Step(const State6D& state, double h) const noexcept;

    /**
     * @brief Computes adaptive step size based on radial distance to the photon sphere and horizon.
     */
    [[nodiscard]] double computeAdaptiveStepSize(double r) const noexcept;

    /**
     * @brief Integrates a single light ray until termination condition is met.
     * @param ray Ray object to evolve and populate with final state and metrics.
     */
    void integrateRay(Ray& ray) const noexcept;
};

#endif // BLACKHOLE_GEODESICINTEGRATOR_HPP
