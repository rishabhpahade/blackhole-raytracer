#ifndef BLACKHOLE_BLACKHOLE_HPP
#define BLACKHOLE_BLACKHOLE_HPP

#include "Vector3.hpp"
#include <cmath>

/**
 * @brief Represents a stationary, uncharged Schwarzschild Black Hole.
 * 
 * In geometrized units where G = c = 1:
 * - Event Horizon radius:  r_h  = 2M
 * - Photon Sphere radius:  r_ph = 3M
 * - Critical impact param: b_c  = 3 * sqrt(3) * M
 */
class BlackHole {
public:
    double mass{1.0};
    Vector3 position{0.0, 0.0, 0.0};
    double horizonTolerance{1e-3};

    // Constructors
    constexpr BlackHole() noexcept = default;
    constexpr explicit BlackHole(double mass_, const Vector3& pos = {0.0, 0.0, 0.0}, double tolerance = 1e-3) noexcept
        : mass(mass_), position(pos), horizonTolerance(tolerance) {}

    // Schwarzschild Event Horizon radius: r_h = 2M
    [[nodiscard]] constexpr double eventHorizonRadius() const noexcept {
        return 2.0 * mass;
    }

    // Photon Sphere radius: r_ph = 3M
    [[nodiscard]] constexpr double photonSphereRadius() const noexcept {
        return 3.0 * mass;
    }

    // Critical impact parameter for capture: b_c = 3 * sqrt(3) * M
    [[nodiscard]] double criticalImpactParameter() const noexcept {
        return 3.0 * std::sqrt(3.0) * mass;
    }

    // ISCO (Innermost Stable Circular Orbit for test particles): r_ISCO = 6M
    [[nodiscard]] constexpr double iscoRadius() const noexcept {
        return 6.0 * mass;
    }

    /**
     * @brief Tests whether a given position vector is inside or at the event horizon boundary.
     * @param pos 3D spatial coordinate
     * @return true if distance to center <= r_h * (1.0 + horizonTolerance)
     */
    [[nodiscard]] bool isInsideHorizon(const Vector3& pos) const noexcept {
        const double dist = (pos - position).magnitude();
        return dist <= eventHorizonRadius() * (1.0 + horizonTolerance);
    }

    /**
     * @brief Computes distance from a point to the black hole center.
     */
    [[nodiscard]] double distanceTo(const Vector3& pos) const noexcept {
        return (pos - position).magnitude();
    }
};

#endif // BLACKHOLE_BLACKHOLE_HPP
