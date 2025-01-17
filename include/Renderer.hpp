#ifndef BLACKHOLE_RENDERER_HPP
#define BLACKHOLE_RENDERER_HPP

#include "Camera.hpp"
#include "GeodesicIntegrator.hpp"
#include "BlackHole.hpp"
#include "AccretionDisk.hpp"
#include "Vector3.hpp"
#include <vector>
#include <string>
#include <chrono>

/**
 * @brief Rendering mode options for scene generation.
 */
enum class BackgroundMode {
    ProceduralStarfield, // Deep space starry sky with nebula clouds
    CoordinateGrid,      // Lensed celestial grid for lensing visualization
    RayClassification    // Diagnostic false-color map of RayStatus
};

/**
 * @brief Aggregated simulation and performance metrics.
 */
struct SimulationMetrics {
    int imageWidth{0};
    int imageHeight{0};
    long long totalRays{0};
    long long capturedCount{0};
    long long escapedCount{0};
    long long diskHitCount{0};
    long long maxStepsCount{0};
    long long failureCount{0};
    
    double capturedPercentage{0.0};
    double escapedPercentage{0.0};
    double diskHitPercentage{0.0};

    double averageSteps{0.0};
    int maxStepsObserved{0};
    double renderTimeSeconds{0.0};

    [[nodiscard]] std::string formatSummary() const;
};

/**
 * @brief Coordinates pixel ray generation, multi-threaded geodesic integration,
 * shading, color mapping, metrics gathering, and image export.
 */
class Renderer {
public:
    Camera camera;
    GeodesicIntegrator integrator;
    BackgroundMode bgMode{BackgroundMode::ProceduralStarfield};
    bool enableToneMapping{true};
    double exposure{1.0};

    SimulationMetrics lastMetrics;

    Renderer(const Camera& cam, const GeodesicIntegrator& gi,
             BackgroundMode mode = BackgroundMode::ProceduralStarfield) noexcept
        : camera(cam), integrator(gi), bgMode(mode) {}

    /**
     * @brief Renders the scene and stores the resulting image buffer.
     * @param outLinearBuffer Floating-point RGB pixel buffer (width * height * 3)
     */
    void render(std::vector<float>& outLinearBuffer);

    /**
     * @brief Computes procedural background color for escaped rays.
     * @param dir Asymptotic normalized ray direction
     */
    [[nodiscard]] Vector3 sampleBackground(const Vector3& dir) const noexcept;

    /**
     * @brief Converts linear HDR RGB to 8-bit sRGB with optional Reinhard/ACES tone mapping.
     */
    [[nodiscard]] static std::vector<unsigned char> convertTo8BitSRGB(
        const std::vector<float>& linearBuffer, int width, int height, bool applyToneMap, double exp);

    /**
     * @brief Writes image buffer to PNG format using stb_image_write.
     */
    static bool savePNG(const std::string& filename, const std::vector<unsigned char>& byteBuffer,
                        int width, int height, int channels = 3);

    /**
     * @brief Runs complete batch generation for all required deliverable artifacts.
     */
    void generateDeliverables(const std::string& outputDir);

    /**
     * @brief Generates deflection angle versus impact parameter dataset (CSV).
     */
    static void generateDeflectionData(const GeodesicIntegrator& integrator,
                                      const std::string& outputPath,
                                      double bMin = 3.0, double bMax = 50.0, int samples = 500);
};

#endif // BLACKHOLE_RENDERER_HPP
