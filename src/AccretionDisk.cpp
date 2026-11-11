#include "AccretionDisk.hpp"
#include <cmath>
#include <algorithm>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

bool AccretionDisk::checkIntersection(const Vector3& p0, const Vector3& p1,
                                     Vector3& outIntersectionPoint, double& outRadius) const noexcept {
    // Check if the segment spans across the equatorial z = 0 plane
    if ((p0.z > 0.0 && p1.z > 0.0) || (p0.z < 0.0 && p1.z < 0.0)) {
        return false;
    }

    const double dz = p1.z - p0.z;
    if (std::abs(dz) < 1e-12) {
        // Parallel or residing directly on disk plane
        const double r0 = std::sqrt(p0.x * p0.x + p0.y * p0.y);
        if (containsRadius(r0)) {
            outIntersectionPoint = p0;
            outRadius = r0;
            return true;
        }
        return false;
    }

    // Linear interpolation parameter t ∈ [0, 1] where z(t) = 0
    const double t = -p0.z / dz;
    if (t < 0.0 || t > 1.0) {
        return false;
    }

    outIntersectionPoint = p0 + (p1 - p0) * t;
    outRadius = std::sqrt(outIntersectionPoint.x * outIntersectionPoint.x +
                          outIntersectionPoint.y * outIntersectionPoint.y);

    return containsRadius(outRadius);
}

Vector3 AccretionDisk::temperatureToRGB(double normalizedTemp) noexcept {
    // High-fidelity Planckian blackbody gradient approximation
    // Range: [0.0 (cool dark red) -> 0.5 (warm gold) -> 1.0 (blazing white-blue)]
    const double t = std::clamp(normalizedTemp, 0.0, 1.0);

    double r = 0.0, g = 0.0, b = 0.0;

    if (t < 0.25) {
        // Deep infrared/red to crimson
        const double s = t / 0.25;
        r = 0.3 * s + 0.6 * (s * s);
        g = 0.05 * s * s;
        b = 0.01 * s;
    } else if (t < 0.6) {
        // Crimson to fiery amber/orange
        const double s = (t - 0.25) / 0.35;
        r = 0.9 + 0.1 * s;
        g = 0.05 + 0.55 * s;
        b = 0.01 + 0.15 * (s * s);
    } else if (t < 0.85) {
        // Amber/gold to brilliant warm white
        const double s = (t - 0.6) / 0.25;
        r = 1.0;
        g = 0.6 + 0.35 * s;
        b = 0.16 + 0.64 * s;
    } else {
        // Intense bluish white
        const double s = (t - 0.85) / 0.15;
        r = 1.0 - 0.05 * s;
        g = 0.95 + 0.05 * s;
        b = 0.8 + 0.2 * s;
    }

    return {r, g, b};
}

Vector3 AccretionDisk::computeColor(double r, double phi, double g) const noexcept {
    if (!containsRadius(r)) {
        return {0.0, 0.0, 0.0};
    }

    // Page-Thorne / Shakura-Sunyaev standard thin disk temperature profile:
    // T(r) ∝ r^(-3/4) * [1 - sqrt(r_in / r)]^(1/4)
    const double rRatio = innerRadius / r;
    const double torqueCutoff = std::max(0.0, 1.0 - std::sqrt(rRatio));
    const double tempProfile = std::pow(rRatio, 0.75) * std::pow(torqueCutoff, 0.25);

    // Subtle spiral density waves and azimuthal texture modulation
    const double spiralPattern = 0.85 + 0.15 * std::sin(4.0 * phi - 1.2 * std::log(r / innerRadius + 1.0));
    const double radialFade = 1.0 - std::pow((r - innerRadius) / (outerRadius - innerRadius), 2.0);

    // Blueshifted gas radiates as if hotter, redshifted gas as if cooler.
    const double shift = (g > 0.0 && std::isfinite(g)) ? g : 1.0;
    const double effectiveTemp = std::clamp(tempProfile * 1.5 * shift, 0.0, 1.0);
    const Vector3 baseColor = temperatureToRGB(effectiveTemp);

    // Relativistic beaming: the bolometric intensity scales as g^4, so the
    // approaching side of the disk is dramatically brighter than the receding
    // side.
    const double beaming = shift * shift * shift * shift;
    const double emissionIntensity =
        peakIntensity * (tempProfile * 2.8 + 0.2) * spiralPattern * radialFade * beaming;

    return baseColor * emissionIntensity;
}
