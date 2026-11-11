#include "Vector3.hpp"
#include "BlackHole.hpp"
#include "AccretionDisk.hpp"
#include "Camera.hpp"
#include "GeodesicIntegrator.hpp"
#include "Renderer.hpp"
#include <iostream>
#include <string>
#include <filesystem>
#include <cmath>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

int main(int argc, char* argv[]) {
    std::cout << "==========================================================\n";
    std::cout << "     SCHWARZSCHILD BLACK-HOLE NUMERICAL RAY TRACER        \n";
    std::cout << "  Simulating Null Geodesics & Relativistic Light Bending  \n";
    std::cout << "==========================================================\n\n";

    // Configurable rendering resolution
    int width = 1280;
    int height = 720;
    int samplesPerPixel = 2;

    if (argc >= 3) {
        width = std::stoi(argv[1]);
        height = std::stoi(argv[2]);
    }
    if (argc >= 4) {
        samplesPerPixel = std::stoi(argv[3]);
    }

    std::cout << "Target Configuration:\n";
    std::cout << "  * Resolution:        " << width << " x " << height << "\n";
    std::cout << "  * Samples Per Pixel: " << samplesPerPixel << "\n";

    // 1. Spacetime & Black Hole parameters (G = c = 1)
    const double M = 1.0;
    BlackHole blackHole(M, Vector3{0.0, 0.0, 0.0}, 1e-3);
    std::cout << "  * Black Hole Mass M: " << blackHole.mass << "\n";
    std::cout << "  * Event Horizon r_h: " << blackHole.eventHorizonRadius() << " M\n";
    std::cout << "  * Photon Sphere r_ph:" << blackHole.photonSphereRadius() << " M\n";
    std::cout << "  * Critical Impact b_c:" << blackHole.criticalImpactParameter() << " M\n";

    // 2. Accretion Disk in equatorial plane (z = 0)
    AccretionDisk disk(6.0 * M, 22.0 * M, 3.0);
    std::cout << "  * Accretion Disk:    r_in = " << disk.innerRadius 
              << " M (ISCO), r_out = " << disk.outerRadius << " M\n";

    // 3. Virtual Observer Camera
    // The disk lies in the z = 0 plane, so the camera's |z| is its height above
    // the disk and its in-plane radius is |x,y|. A small height relative to the
    // in-plane distance gives the classic near-edge-on view where the far side
    // of the disk is lensed up over the shadow. The previous default sat 26M
    // above the plane with only 3.2M of in-plane offset — a nearly face-on
    // view in which the disk filled ~77% of the frame.
    const double inclinationDegrees = 15.0;
    const double orbitRadius = 56.0;
    const double inclinationRad = inclinationDegrees * (M_PI / 180.0);
    const Vector3 camPos{
        0.0,
        orbitRadius * std::cos(inclinationRad),
        orbitRadius * std::sin(inclinationRad)
    };
    const Vector3 camTarget{0.0, 0.0, 0.0};
    // The disk normal is +z, so +z is "up" in the image at a shallow
    // inclination. Camera::updateBasis falls back to another axis if a caller
    // passes an up vector parallel to the view direction.
    const Vector3 camUp{0.0, 0.0, 1.0};
    const double fov = 32.0;

    Camera camera(camPos, camTarget, camUp, fov, width, height, samplesPerPixel);

    // 4. Numerical Geodesic Integrator
    const double initialStep = 0.06;
    // Keep the escape sphere well outside the camera so the asymptotic
    // direction is sampled far from the observer.
    const double escapeRadius = 120.0 * M;
    const int maxSteps = 15000;
    GeodesicIntegrator integrator(blackHole, disk, initialStep, escapeRadius, maxSteps);

    // 5. Output directory setup
    const std::string outputDir = "output";
    std::filesystem::create_directories(outputDir);

    // 6. Execute rendering pipeline
    Renderer renderer(camera, integrator, BackgroundMode::ProceduralStarfield);
    renderer.exposure = 0.18;
    renderer.enableToneMapping = true;

    std::cout << "\nStarting numerical simulation and deliverable generation...\n";
    renderer.generateDeliverables(outputDir);

    std::cout << "Simulation completed successfully. Deliverables generated in ./" << outputDir << "/\n";
    return 0;
}
