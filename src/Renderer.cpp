#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "external/stb_image_write.h"

#include "Renderer.hpp"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <algorithm>
#include <random>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

std::string SimulationMetrics::formatSummary() const {
    std::ostringstream oss;
    oss << "=====================================================\n"
        << "      SCHWARZSCHILD BLACK-HOLE RAY TRACER METRICS    \n"
        << "=====================================================\n"
        << "Image Resolution:          " << imageWidth << " x " << imageHeight << " (" << totalRays << " rays)\n"
        << "Render Execution Time:     " << std::fixed << std::setprecision(3) << renderTimeSeconds << " seconds\n"
        << "Ray Throughput:            " << static_cast<long long>(totalRays / std::max(0.001, renderTimeSeconds)) << " rays/sec\n"
        << "-----------------------------------------------------\n"
        << "RAY CLASSIFICATION BREAKDOWN:\n"
        << "  * Captured (Shadow):     " << capturedCount << " (" << std::fixed << std::setprecision(2) << capturedPercentage << "%)\n"
        << "  * Escaped (Lensed):      " << escapedCount << " (" << std::fixed << std::setprecision(2) << escapedPercentage << "%)\n"
        << "  * Hit Accretion Disk:    " << diskHitCount << " (" << std::fixed << std::setprecision(2) << diskHitPercentage << "%)\n"
        << "  * Max Steps Reached:     " << maxStepsCount << "\n"
        << "  * Numerical Failures:    " << failureCount << "\n"
        << "-----------------------------------------------------\n"
        << "INTEGRATION DEPTH STATISTICS:\n"
        << "  * Average Steps per Ray: " << std::fixed << std::setprecision(1) << averageSteps << "\n"
        << "  * Maximum Steps Observed:" << maxStepsObserved << "\n"
        << "=====================================================\n";
    return oss.str();
}

Vector3 Renderer::sampleBackground(const Vector3& dir) const noexcept {
    const Vector3 n = dir.normalized();

    // Map asymptotic direction to celestial spherical coordinates
    const double theta = std::acos(std::clamp(n.y, -1.0, 1.0)); // [0, pi]
    const double phi = std::atan2(n.z, n.x) + M_PI;             // [0, 2pi]

    if (bgMode == BackgroundMode::CoordinateGrid) {
        // Celestial coordinate grid lines every 15 degrees
        const double gridSpacingRad = 15.0 * (M_PI / 180.0);
        const double thetaMod = std::fmod(theta, gridSpacingRad);
        const double phiMod = std::fmod(phi, gridSpacingRad);

        const double lineWidth = 0.02;
        if (thetaMod < lineWidth || phiMod < lineWidth) {
            return {0.3, 0.7, 1.0}; // Light cyan grid lines
        }

        // Soft background checkered hue
        const int checker = (static_cast<int>(theta / gridSpacingRad) + static_cast<int>(phi / gridSpacingRad)) % 2;
        return checker ? Vector3{0.02, 0.03, 0.06} : Vector3{0.05, 0.06, 0.10};
    }

    // Procedural Starfield with multi-spectral stars and subtle nebulae
    // 1. Soft deep-space nebula background
    const double nebula1 = std::sin(3.0 * phi) * std::cos(4.0 * theta);
    const double nebula2 = std::sin(5.0 * phi + 1.2) * std::sin(3.0 * theta);
    Vector3 bg = {
        0.015 + 0.025 * std::max(0.0, nebula1),
        0.010 + 0.020 * std::max(0.0, nebula2),
        0.030 + 0.050 * std::max(0.0, (nebula1 + nebula2) * 0.5)
    };

    // 2. High-frequency deterministic star hash
    const double starSeed = std::sin(n.x * 127.1 + n.y * 311.7 + n.z * 74.7) * 43758.5453123;
    const double starFract = starSeed - std::floor(starSeed);

    if (starFract > 0.985) {
        // Star intensity with varied color temperatures
        const double intensity = std::pow((starFract - 0.985) / 0.015, 6.0) * 8.0;
        
        // Spectral distribution based on hash sub-fraction
        const double colorSeed = (starSeed * 10.0) - std::floor(starSeed * 10.0);
        if (colorSeed < 0.25) {
            // Blue Supergiant
            bg += Vector3{0.6, 0.8, 1.0} * intensity;
        } else if (colorSeed < 0.7) {
            // Solar white / yellow
            bg += Vector3{1.0, 0.95, 0.85} * intensity;
        } else {
            // Red / Orange Dwarf
            bg += Vector3{1.0, 0.5, 0.3} * intensity;
        }
    }

    return bg;
}

void Renderer::render(std::vector<float>& outLinearBuffer) {
    const int w = camera.imageWidth;
    const int h = camera.imageHeight;
    outLinearBuffer.assign(w * h * 3, 0.0f);

    SimulationMetrics metrics;
    metrics.imageWidth = w;
    metrics.imageHeight = h;
    metrics.totalRays = static_cast<long long>(w) * h;

    long long totalStepsAccum = 0;
    int maxStepsFound = 0;

    const auto startTime = std::chrono::high_resolution_clock::now();

    #pragma omp parallel for schedule(dynamic, 16) reduction(+:totalStepsAccum) reduction(max:maxStepsFound)
    for (int y = 0; y < h; ++y) {
        for (int x = 0; x < w; ++x) {
            Vector3 pixelColor{0.0, 0.0, 0.0};
            RayStatus lastStatus = RayStatus::MaximumStepsReached;
            int lastSteps = 0;

            const int spp = std::max(1, camera.samplesPerPixel);
            for (int s = 0; s < spp; ++s) {
                // Sub-pixel jitter for anti-aliasing if spp > 1
                double offsetX = 0.0;
                double offsetY = 0.0;
                if (spp > 1) {
                    offsetX = ((s % 2) - 0.5) * 0.5;
                    offsetY = ((s / 2) - 0.5) * 0.5;
                }

                Ray ray = camera.generateRay(x + offsetX, y + offsetY);
                integrator.integrateRay(ray);

                lastStatus = ray.status;
                lastSteps = ray.steps;

                if (bgMode == BackgroundMode::RayClassification) {
                    // False-color classification map
                    switch (ray.status) {
                        case RayStatus::Captured:
                            pixelColor += Vector3{0.0, 0.0, 0.0}; // Black
                            break;
                        case RayStatus::HitAccretionDisk:
                            pixelColor += Vector3{1.0, 0.6, 0.0}; // Amber Gold
                            break;
                        case RayStatus::Escaped:
                            pixelColor += Vector3{0.1, 0.45, 1.0}; // Electric Blue
                            break;
                        case RayStatus::MaximumStepsReached:
                            pixelColor += Vector3{0.3, 0.3, 0.3}; // Gray
                            break;
                        case RayStatus::NumericalFailure:
                            pixelColor += Vector3{1.0, 0.0, 1.0}; // Magenta
                            break;
                    }
                } else {
                    // Physical visual shading
                    if (ray.status == RayStatus::Captured) {
                        pixelColor += Vector3{0.0, 0.0, 0.0};
                    } else if (ray.status == RayStatus::HitAccretionDisk) {
                        pixelColor += ray.color;
                    } else if (ray.status == RayStatus::Escaped) {
                        pixelColor += sampleBackground(ray.velocity);
                    } else {
                        pixelColor += ray.color;
                    }
                }
            }

            pixelColor = pixelColor / static_cast<double>(spp);

            const int idx = (y * w + x) * 3;
            outLinearBuffer[idx + 0] = static_cast<float>(pixelColor.x);
            outLinearBuffer[idx + 1] = static_cast<float>(pixelColor.y);
            outLinearBuffer[idx + 2] = static_cast<float>(pixelColor.z);

            totalStepsAccum += lastSteps;
            maxStepsFound = std::max(maxStepsFound, lastSteps);

            #pragma omp critical
            {
                switch (lastStatus) {
                    case RayStatus::Captured:            metrics.capturedCount++; break;
                    case RayStatus::Escaped:             metrics.escapedCount++; break;
                    case RayStatus::HitAccretionDisk:    metrics.diskHitCount++; break;
                    case RayStatus::MaximumStepsReached: metrics.maxStepsCount++; break;
                    case RayStatus::NumericalFailure:    metrics.failureCount++; break;
                }
            }
        }
    }

    const auto endTime = std::chrono::high_resolution_clock::now();
    const std::chrono::duration<double> elapsed = endTime - startTime;

    metrics.renderTimeSeconds = elapsed.count();
    metrics.capturedPercentage = (100.0 * metrics.capturedCount) / metrics.totalRays;
    metrics.escapedPercentage = (100.0 * metrics.escapedCount) / metrics.totalRays;
    metrics.diskHitPercentage = (100.0 * metrics.diskHitCount) / metrics.totalRays;
    metrics.averageSteps = static_cast<double>(totalStepsAccum) / metrics.totalRays;
    metrics.maxStepsObserved = maxStepsFound;

    lastMetrics = metrics;
}

std::vector<unsigned char> Renderer::convertTo8BitSRGB(
    const std::vector<float>& linearBuffer, int width, int height, bool applyToneMap, double exp) {
    
    std::vector<unsigned char> byteBuffer(width * height * 3);

    for (size_t i = 0; i < linearBuffer.size(); ++i) {
        float val = linearBuffer[i] * static_cast<float>(exp);

        if (applyToneMap) {
            // Extended Reinhard Tone Mapping: L_d = (L * (1 + L / L_white^2)) / (1 + L)
            const float lWhiteSq = 16.0f;
            val = (val * (1.0f + val / lWhiteSq)) / (1.0f + val);
        }

        // Gamma correction to sRGB (gamma ≈ 2.2)
        val = std::pow(std::max(0.0f, val), 1.0f / 2.2f);
        val = std::clamp(val * 255.0f + 0.5f, 0.0f, 255.0f);

        byteBuffer[i] = static_cast<unsigned char>(val);
    }

    return byteBuffer;
}

bool Renderer::savePNG(const std::string& filename, const std::vector<unsigned char>& byteBuffer,
                       int width, int height, int channels) {
    const int result = stbi_write_png(filename.c_str(), width, height, channels, byteBuffer.data(), width * channels);
    return result != 0;
}

void Renderer::generateDeflectionData(const GeodesicIntegrator& integrator,
                                      const std::string& outputPath,
                                      double bMin, double bMax, int samples) {
    std::ofstream csv(outputPath);
    if (!csv.is_open()) {
        std::cerr << "Failed to open CSV output path: " << outputPath << "\n";
        return;
    }

    csv << "ImpactParameter_b,NumericalDeflection_rad,EinsteinWeakField_rad,ClosestApproach_r,Status\n";

    const double M = integrator.blackHole.mass;
    const double zStart = -60.0 * M;
    const Vector3 initialDir{0.0, 0.0, 1.0}; // Parallel beam traveling in +z direction

    for (int i = 0; i < samples; ++i) {
        const double t = static_cast<double>(i) / (samples - 1);
        const double b = bMin + (bMax - bMin) * t;

        // Position ray at x = -b, y = 0, z = zStart
        Ray ray(Vector3{-b, 0.0, zStart}, initialDir);
        integrator.integrateRay(ray);

        const double alphaNumerical = ray.deflectionAngle();
        const double alphaEinstein = (4.0 * M) / b;

        csv << std::fixed << std::setprecision(6)
            << b << ","
            << alphaNumerical << ","
            << alphaEinstein << ","
            << ray.closestApproach << ","
            << rayStatusToString(ray.status) << "\n";
    }

    csv.close();
}

void Renderer::generateDeliverables(const std::string& outputDir) {
    std::cout << "[1/4] Rendering primary visual simulation (black_hole.png)...\n";
    bgMode = BackgroundMode::ProceduralStarfield;
    std::vector<float> heroLinearBuffer;
    render(heroLinearBuffer);

    auto heroBytes = convertTo8BitSRGB(heroLinearBuffer, camera.imageWidth, camera.imageHeight, enableToneMapping, exposure);
    const std::string heroPath = outputDir + "/black_hole.png";
    if (savePNG(heroPath, heroBytes, camera.imageWidth, camera.imageHeight, 3)) {
        std::cout << "  -> Saved: " << heroPath << "\n";
    }

    // Save metrics
    const std::string metricsPath = outputDir + "/simulation_metrics.txt";
    std::ofstream metricsFile(metricsPath);
    if (metricsFile.is_open()) {
        metricsFile << lastMetrics.formatSummary();
        metricsFile.close();
        std::cout << "[2/4] Saved simulation metrics to " << metricsPath << "\n";
    }

    std::cout << "[3/4] Rendering ray classification diagnostic map (ray_classification.png)...\n";
    bgMode = BackgroundMode::RayClassification;
    std::vector<float> classLinearBuffer;
    render(classLinearBuffer);

    auto classBytes = convertTo8BitSRGB(classLinearBuffer, camera.imageWidth, camera.imageHeight, false, 1.0);
    const std::string classPath = outputDir + "/ray_classification.png";
    if (savePNG(classPath, classBytes, camera.imageWidth, camera.imageHeight, 3)) {
        std::cout << "  -> Saved: " << classPath << "\n";
    }

    std::cout << "[4/4] Computing deflection vs impact parameter curve...\n";
    GeodesicIntegrator cleanIntegrator = integrator;
    cleanIntegrator.testAccretionDisk = false; // Pure light deflection in empty vacuum Schwarzschild
    const std::string csvPath = outputDir + "/deflection_vs_impact_parameter.csv";
    generateDeflectionData(cleanIntegrator, csvPath, 3.0, 50.0, 600);
    std::cout << "  -> Saved: " << csvPath << "\n";

    std::cout << "\n" << lastMetrics.formatSummary() << "\n";
}
