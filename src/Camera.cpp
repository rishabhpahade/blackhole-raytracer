#include "Camera.hpp"
#include <cmath>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

Camera::Camera() noexcept {
    updateBasis();
}

Camera::Camera(const Vector3& pos, const Vector3& lookAt, const Vector3& worldUp,
               double fovDeg, int width, int height, int spp) noexcept
    : position(pos),
      target(lookAt),
      upVector(worldUp),
      fovDegrees(fovDeg),
      imageWidth(width),
      imageHeight(height),
      samplesPerPixel(spp) {
    updateBasis();
}

void Camera::updateBasis() noexcept {
    forward = (target - position).normalized();
    if (forward.magnitudeSquared() < 1e-12) {
        forward = {0.0, 0.0, 1.0};
    }

    // A degenerate reference up (zero length, or (anti)parallel to the view
    // direction) would collapse the basis to all-zero vectors and make every
    // pixel produce the same ray. Fall back to a different world axis.
    Vector3 upRef = upVector.normalized();
    if (upRef.magnitudeSquared() < 1e-12) {
        upRef = {0.0, 1.0, 0.0};
    }
    if (std::abs(forward.dot(upRef)) > 0.9999) {
        upRef = (std::abs(forward.z) < 0.9) ? Vector3{0.0, 0.0, 1.0} : Vector3{1.0, 0.0, 0.0};
    }

    right = forward.cross(upRef).normalized();
    up = right.cross(forward).normalized();

    const double fovRad = fovDegrees * (M_PI / 180.0);
    const double halfHeight = std::tan(fovRad * 0.5);
    const double halfWidth = aspectRatio() * halfHeight;

    viewPlaneWidth = 2.0 * halfWidth;
    viewPlaneHeight = 2.0 * halfHeight;
}

Ray Camera::generateRay(double px, double py) const noexcept {
    // Map raster pixel coordinate to [-halfWidth, +halfWidth] and [-halfHeight, +halfHeight]
    const double halfWidth = viewPlaneWidth * 0.5;
    const double halfHeight = viewPlaneHeight * 0.5;

    const double u = ((px + 0.5) / static_cast<double>(imageWidth) - 0.5) * 2.0 * halfWidth;
    const double v = (0.5 - (py + 0.5) / static_cast<double>(imageHeight)) * 2.0 * halfHeight;

    // Ray direction vector in 3D camera coordinate frame
    const Vector3 dir = (forward + right * u + up * v).normalized();
    return Ray(position, dir);
}
