#include "Vector3.hpp"
#include "BlackHole.hpp"
#include "AccretionDisk.hpp"
#include "Camera.hpp"
#include "GeodesicIntegrator.hpp"
#include "Renderer.hpp"
#include <iostream>
#include <string>
#include <filesystem>

int main(int argc, char* argv[]) {
    std::cout << "==========================================================\n";
    std::cout << "     SCHWARZSCHILD BLACK-HOLE NUMERICAL RAY TRACER        \n";
    std::cout << "  Simulating Null Geodesics & Relativistic Light Bending  \n";
    std::cout << "==========================================================\n\n";

    // Configurable rendering resolution
    int width = 1280;
    int height = 720;
    int samplesPerPixel = 1;

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
    // Position observer at (0, 3.2, -26) to see dramatic gravitational lensing of the disk
    const Vector3 camPos{0.0, 3.2, -26.0};
    const Vector3 camTarget{0.0, 0.0, 0.0};
    const Vector3 camUp{0.0, 1.0, 0.0};
    const double fov = 42.0;

    Camera camera(camPos, camTarget, camUp, fov, width, height, samplesPerPixel);

    // 4. Numerical Geodesic Integrator
    const double initialStep = 0.06;
    const double escapeRadius = 60.0 * M;
    const int maxSteps = 15000;
    GeodesicIntegrator integrator(blackHole, disk, initialStep, escapeRadius, maxSteps);

    // 5. Output directory setup
    const std::string outputDir = "output";
    std::filesystem::create_directories(outputDir);

    // 6. Execute rendering pipeline
    Renderer renderer(camera, integrator, BackgroundMode::ProceduralStarfield);
    renderer.exposure = 1.1;
    renderer.enableToneMapping = true;

    std::cout << "\nStarting numerical simulation and deliverable generation...\n";
    renderer.generateDeliverables(outputDir);

    std::cout << "Simulation completed successfully. Deliverables generated in ./" << outputDir << "/\n";
    return 0;
}
