#ifndef BLACKHOLE_ACCRETIONDISK_HPP
#define BLACKHOLE_ACCRETIONDISK_HPP

#include "Vector3.hpp"
#include <cmath>
#include <algorithm>

/**
 * @brief Represents an geometrically thin, optically glowing equatorial accretion disk.
 * 
 * Located at z = 0 with inner radius r_in (typically ISCO = 6M) and outer radius r_out.
 * Radiates with effective temperature profile T(r) ∝ r^(-3/4) * (1 - sqrt(r_in/r))^(1/4)
 * or power-law T(r) ∝ r^(-3/4).
 */
class AccretionDisk {
public:
    double innerRadius{6.0};  // ISCO for Schwarzschild M=1 is 6M
    double outerRadius{20.0}; // Outer disk extent
    double peakIntensity{2.5};
    bool enableDopplerBeaming{false}; // For non-rotating Schwarzschild, disk orbital velocity can be Keplerian

    constexpr AccretionDisk() noexcept = default;
    constexpr AccretionDisk(double rInner, double rOuter, double intensity = 2.5) noexcept
        : innerRadius(rInner), outerRadius(rOuter), peakIntensity(intensity) {}

    /**
     * @brief Checks if a radial coordinate r falls within the accretion disk boundaries.
     */
    [[nodiscard]] constexpr bool containsRadius(double r) const noexcept {
        return r >= innerRadius && r <= outerRadius;
    }

    /**
     * @brief Determines whether a segment between two points crosses the z=0 equatorial disk plane.
     * @param p0 Previous position
     * @param p1 Current position
     * @param outIntersectionPoint Resulting 3D intersection point on the disk
     * @param outRadius Resulting radial distance r = sqrt(x^2 + y^2)
     * @return true if the segment intersected the disk between r_in and r_out
     */
    [[nodiscard]] bool checkIntersection(const Vector3& p0, const Vector3& p1,
                                         Vector3& outIntersectionPoint, double& outRadius) const noexcept;

    /**
     * @brief Computes radiated linear RGB color based on radial temperature model T(r) ∝ r^(-3/4).
     * @param r Radial distance from black hole center
     * @param phi Azimuthal angle around disk
     * @return Linear RGB color (can exceed 1.0 before tone mapping)
     */
    [[nodiscard]] Vector3 computeColor(double r, double phi = 0.0) const noexcept;

    /**
     * @brief Converts blackbody / thermal temperature to RGB color representation.
     */
    [[nodiscard]] static Vector3 temperatureToRGB(double normalizedTemp) noexcept;
};

#endif // BLACKHOLE_ACCRETIONDISK_HPP
