#ifndef BLACKHOLE_RAY_HPP
#define BLACKHOLE_RAY_HPP

#include "Vector3.hpp"
#include <string>
#include <limits>

/**
 * @brief Categorical outcome status of a numerically traced null geodesic.
 */
enum class RayStatus {
    Captured,            // Plunged across the event horizon (r <= 2M)
    Escaped,             // Escaped to asymptotic infinity (r >= r_escape)
    HitAccretionDisk,    // Crossed the equatorial accretion disk plane (r_in <= r <= r_out)
    MaximumStepsReached, // Exceeded iteration limit before resolving
    NumericalFailure     // Encountered NaN, infinity, or step underflow
};

/**
 * @brief Converts RayStatus enum to a readable string.
 */
inline std::string rayStatusToString(RayStatus status) {
    switch (status) {
        case RayStatus::Captured:            return "Captured";
        case RayStatus::Escaped:             return "Escaped";
        case RayStatus::HitAccretionDisk:    return "HitAccretionDisk";
        case RayStatus::MaximumStepsReached: return "MaximumStepsReached";
        case RayStatus::NumericalFailure:    return "NumericalFailure";
        default:                             return "Unknown";
    }
}

/**
 * @brief Represents a photon (null geodesic) undergoing relativistic propagation.
 */
class Ray {
public:
    Vector3 initialPosition{0.0, 0.0, 0.0};
    Vector3 initialDirection{0.0, 0.0, -1.0}; // Normalized initial wavevector
    Vector3 position{0.0, 0.0, 0.0};          // Current spatial position
    Vector3 velocity{0.0, 0.0, -1.0};         // Current spatial wavevector dx/dλ

    RayStatus status{RayStatus::MaximumStepsReached};
    int steps{0};
    double closestApproach{std::numeric_limits<double>::infinity()};
    double diskHitRadius{0.0};
    Vector3 diskHitPoint{0.0, 0.0, 0.0};
    Vector3 color{0.0, 0.0, 0.0};             // Linear RGB [0..1]

    constexpr Ray() noexcept = default;

    Ray(const Vector3& origin, const Vector3& direction) noexcept
        : initialPosition(origin),
          initialDirection(direction.normalized()),
          position(origin),
          velocity(direction.normalized()),
          closestApproach(origin.magnitude()) {}

    /**
     * @brief Resets the ray with new initial coordinates.
     */
    void reset(const Vector3& origin, const Vector3& direction) noexcept {
        initialPosition = origin;
        initialDirection = direction.normalized();
        position = origin;
        velocity = initialDirection;
        status = RayStatus::MaximumStepsReached;
        steps = 0;
        closestApproach = origin.magnitude();
        diskHitRadius = 0.0;
        diskHitPoint = {0.0, 0.0, 0.0};
        color = {0.0, 0.0, 0.0};
    }

    /**
     * @brief Computes total scattering / deflection angle between initial and current directions.
     * @return Angle in radians in [0, pi]
     */
    [[nodiscard]] double deflectionAngle() const noexcept {
        const double d = initialDirection.dot(velocity.normalized());
        const double clamped = std::max(-1.0, std::min(1.0, d));
        return std::acos(clamped);
    }
};

#endif // BLACKHOLE_RAY_HPP
