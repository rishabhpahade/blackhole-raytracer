#include "Vector3.hpp"
#include "BlackHole.hpp"
#include "Ray.hpp"
#include "GeodesicIntegrator.hpp"
#include "AccretionDisk.hpp"
#include "Camera.hpp"
#include "Renderer.hpp"
#include <iostream>
#include <cmath>
#include <vector>
#include <iomanip>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define TEST_ASSERT(cond, msg) \
    do { \
        if (!(cond)) { \
            std::cerr << "  [FAIL] " << msg << " (Line " << __LINE__ << ")\n"; \
            testsFailed++; \
        } else { \
            std::cout << "  [PASS] " << msg << "\n"; \
            testsPassed++; \
        } \
    } while(0)

int main() {
    std::cout << "==========================================================\n";
    std::cout << "   RUNNING SCHWARZSCHILD RAY TRACER SCIENTIFIC TESTS     \n";
    std::cout << "==========================================================\n\n";

    int testsPassed = 0;
    int testsFailed = 0;

    const double M = 1.0;
    BlackHole bh(M);
    AccretionDisk disk(6.0 * M, 20.0 * M);

    // ---------------------------------------------------------
    // TEST 1: Event Horizon & Radial Infall (r <= 2M)
    // ---------------------------------------------------------
    std::cout << "--- Test 1: Event Horizon Capture (r_h = 2M) ---\n";
    {
        GeodesicIntegrator integrator(bh, disk, 0.04, 80.0, 30000);
        integrator.testAccretionDisk = false;

        // Ray launched directly towards center
        Ray ray(Vector3{0.0, 0.0, -30.0}, Vector3{0.0, 0.0, 1.0});
        integrator.integrateRay(ray);

        TEST_ASSERT(ray.status == RayStatus::Captured, "Direct radial ray is captured by horizon");
        TEST_ASSERT(ray.closestApproach <= 2.005 * M, "Capture occurs at/inside r <= 2M");
        TEST_ASSERT(ray.position.isFinite(), "Captured ray coordinates remain finite");
    }

    // ---------------------------------------------------------
    // TEST 2: Critical Impact Parameter (b_c = 3*sqrt(3)*M ≈ 5.19615M)
    // ---------------------------------------------------------
    std::cout << "\n--- Test 2: Critical Impact Parameter (b_c = 3*sqrt(3)*M) ---\n";
    {
        GeodesicIntegrator integrator(bh, disk, 0.03, 80.0, 40000);
        integrator.testAccretionDisk = false;
        const double bc = bh.criticalImpactParameter();
        const double zStart = -60.0;

        // Sub-critical ray (b < bc): must plunge and be captured
        const double bSub = bc - 0.08;
        Ray raySub(Vector3{-bSub, 0.0, zStart}, Vector3{0.0, 0.0, 1.0});
        integrator.integrateRay(raySub);
        TEST_ASSERT(raySub.status == RayStatus::Captured,
                    "Sub-critical ray (b = " + std::to_string(bSub) + ") is captured");

        // Super-critical ray (b > bc): must orbit/scatter and escape
        const double bSuper = bc + 0.08;
        Ray raySuper(Vector3{-bSuper, 0.0, zStart}, Vector3{0.0, 0.0, 1.0});
        integrator.integrateRay(raySuper);
        TEST_ASSERT(raySuper.status == RayStatus::Escaped,
                    "Super-critical ray (b = " + std::to_string(bSuper) + ") escapes");

        // Verify closest approach of escaped near-critical ray is close to photon sphere r ≈ 3M
        TEST_ASSERT(raySuper.closestApproach < 3.8 * M && raySuper.closestApproach > 2.5 * M,
                    "Near-critical ray grazes near photon sphere r_min = " + std::to_string(raySuper.closestApproach));
    }

    // ---------------------------------------------------------
    // TEST 3: Relativistic Deflection (Einstein 4M/b + 2nd Order Post-Newtonian)
    // ---------------------------------------------------------
    std::cout << "\n--- Test 3: Relativistic Light Deflection & Asymptotic 4M/b ---\n";
    {
        const double zStart = -600.0;
        GeodesicIntegrator asymptoticIntegrator(bh, disk, 0.2, 600.0, 80000);
        asymptoticIntegrator.testAccretionDisk = false;

        // Test across weak-field regime where b << escape radius (600M)
        const std::vector<double> testB = {50.0, 100.0, 150.0};

        for (double b : testB) {
            Ray ray(Vector3{-b, 0.0, zStart}, Vector3{0.0, 0.0, 1.0});
            asymptoticIntegrator.integrateRay(ray);

            const double alphaNumerical = ray.deflectionAngle();
            const double alphaEinstein = (4.0 * M) / b;
            const double secondOrderTerm = (15.0 * M_PI / 4.0) * (M * M) / (b * b);
            const double alphaPostNewtonian = alphaEinstein + secondOrderTerm;

            const double relError = std::abs(alphaNumerical - alphaPostNewtonian) / alphaPostNewtonian;

            TEST_ASSERT(relError < 0.02,
                        "Deflection at b = " + std::to_string(b) + " matches General Relativity theory within "
                        + std::to_string(relError * 100.0) + "% (num: "
                        + std::to_string(alphaNumerical) + " rad, PN theory: " + std::to_string(alphaPostNewtonian) + " rad)");
        }
    }

    // ---------------------------------------------------------
    // TEST 4: Spacetime Axisymmetry & Reflection Symmetry
    // ---------------------------------------------------------
    std::cout << "\n--- Test 4: Axisymmetric Spacetime Invariance ---\n";
    {
        GeodesicIntegrator integrator(bh, disk, 0.04, 80.0, 30000);
        integrator.testAccretionDisk = false;
        const double b = 10.0;
        const double zStart = -50.0;

        Ray rayLeft(Vector3{-b, 0.0, zStart}, Vector3{0.0, 0.0, 1.0});
        Ray rayRight(Vector3{b, 0.0, zStart}, Vector3{0.0, 0.0, 1.0});
        Ray rayTop(Vector3{0.0, b, zStart}, Vector3{0.0, 0.0, 1.0});

        integrator.integrateRay(rayLeft);
        integrator.integrateRay(rayRight);
        integrator.integrateRay(rayTop);

        const double diffLR = std::abs(rayLeft.deflectionAngle() - rayRight.deflectionAngle());
        const double diffLT = std::abs(rayLeft.deflectionAngle() - rayTop.deflectionAngle());

        TEST_ASSERT(diffLR < 1e-10, "Left/Right deflection symmetry error < 1e-10");
        TEST_ASSERT(diffLT < 1e-10, "Planar/Vertical deflection symmetry error < 1e-10");
    }

    // ---------------------------------------------------------
    // TEST 5: RK4 Numerical Convergence
    // ---------------------------------------------------------
    std::cout << "\n--- Test 5: RK4 Numerical Convergence ---\n";
    {
        GeodesicIntegrator integrator(bh, disk);
        GeodesicIntegrator::State6D s0{Vector3{-8.0, 0.0, -20.0}, Vector3{0.0, 0.0, 1.0}};

        // Integrate for fixed total affine parameter distance Δλ = 10.0
        const double totalLambda = 10.0;

        auto advanceState = [&](double h) {
            GeodesicIntegrator::State6D s = s0;
            const int steps = static_cast<int>(totalLambda / h);
            for (int i = 0; i < steps; ++i) {
                s = integrator.rk4Step(s, h);
            }
            return s;
        };

        const auto sRef = advanceState(0.0025);
        const auto sCoarse = advanceState(0.08);
        const auto sMedium = advanceState(0.04);
        const auto sFine = advanceState(0.02);

        const double errCoarse = (sCoarse.pos - sRef.pos).magnitude();
        const double errMedium = (sMedium.pos - sRef.pos).magnitude();
        const double errFine = (sFine.pos - sRef.pos).magnitude();

        TEST_ASSERT(errMedium < errCoarse * 0.15,
                    "Medium error (" + std::to_string(errMedium) + ") < Coarse error (" + std::to_string(errCoarse) + ") * 0.15");
        TEST_ASSERT(errFine < errMedium * 0.15,
                    "Fine error (" + std::to_string(errFine) + ") < Medium error (" + std::to_string(errMedium) + ") * 0.15");
    }

    // ---------------------------------------------------------
    // TEST 6: Accretion Disk Crossing & Gravitational Inward Lensing
    // ---------------------------------------------------------
    std::cout << "\n--- Test 6: Accretion Disk Detection & Inward Bending ---\n";
    {
        GeodesicIntegrator diskIntegrator(bh, disk, 0.04, 60.0, 20000);
        diskIntegrator.testAccretionDisk = true;

        // Ray launched from z = +5 above the disk downward towards z = -5 at x=10, y=0
        Ray rayDisk(Vector3{10.0, 0.0, 5.0}, Vector3{0.0, 0.0, -1.0});
        diskIntegrator.integrateRay(rayDisk);

        TEST_ASSERT(rayDisk.status == RayStatus::HitAccretionDisk, "Ray crossing z=0 hits accretion disk");
        // Gravitational attraction pulls ray from x=10.0 to x ≈ 9.72M at disk crossing
        TEST_ASSERT(rayDisk.diskHitRadius < 10.0 && rayDisk.diskHitRadius > 9.0,
                    "Gravitational bending pulls intersection inward to r = " + std::to_string(rayDisk.diskHitRadius) + " M");
        TEST_ASSERT(rayDisk.color.isFinite() && rayDisk.color.x > 0.0, "Disk emission color is finite and positive");
    }

    // ---------------------------------------------------------
    // Summary
    // ---------------------------------------------------------
    std::cout << "\n==========================================================\n";
    std::cout << "TEST RESULTS: " << testsPassed << " Passed, " << testsFailed << " Failed\n";
    std::cout << "==========================================================\n";

    return (testsFailed == 0) ? 0 : 1;
}
