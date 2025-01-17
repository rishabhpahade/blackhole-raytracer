#ifndef BLACKHOLE_CAMERA_HPP
#define BLACKHOLE_CAMERA_HPP

#include "Vector3.hpp"
#include "Ray.hpp"
#include <cmath>

/**
 * @brief Pinhole virtual camera generating initial null geodesics for each image pixel.
 */
class Camera {
public:
    Vector3 position{0.0, 3.5, -24.0}; // Default viewing position with slight inclination
    Vector3 target{0.0, 0.0, 0.0};     // Looking at black hole center
    Vector3 upVector{0.0, 1.0, 0.0};
    
    double fovDegrees{42.0};
    int imageWidth{800};
    int imageHeight{600};
    int samplesPerPixel{1}; // 1 = standard, >1 = anti-aliasing sub-pixel jitter

    // Orthonormal basis vectors
    Vector3 forward{0.0, 0.0, 1.0};
    Vector3 right{1.0, 0.0, 0.0};
    Vector3 up{0.0, 1.0, 0.0};

    double viewPlaneWidth{0.0};
    double viewPlaneHeight{0.0};

    Camera() noexcept;
    Camera(const Vector3& pos, const Vector3& lookAt, const Vector3& worldUp,
           double fovDeg, int width, int height, int spp = 1) noexcept;

    /**
     * @brief Recalculates orthonormal basis vectors and viewing plane bounds.
     */
    void updateBasis() noexcept;

    /**
     * @brief Generates an initial light ray for a specific (sub)pixel coordinate.
     * @param px Normalized horizontal coordinate in [0, width)
     * @param py Normalized vertical coordinate in [0, height)
     * @return Initialized Ray originating at camera position pointed into the scene
     */
    [[nodiscard]] Ray generateRay(double px, double py) const noexcept;

    [[nodiscard]] double aspectRatio() const noexcept {
        return static_cast<double>(imageWidth) / static_cast<double>(imageHeight);
    }
};

#endif // BLACKHOLE_CAMERA_HPP
