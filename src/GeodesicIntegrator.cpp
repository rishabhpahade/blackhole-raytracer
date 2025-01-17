#include "GeodesicIntegrator.hpp"
#include <cmath>
#include <algorithm>

GeodesicIntegrator::State6D GeodesicIntegrator::evaluateDerivative(const State6D& state) const noexcept {
    const Vector3 relPos = state.pos - blackHole.position;
    const double rSq = relPos.magnitudeSquared();
    const double r = std::sqrt(rSq);

    // Safeguard against division by zero or extreme proximity singularity
    if (r < 1e-8) {
        return {state.vel, {0.0, 0.0, 0.0}};
    }

    // Conserved orbital angular momentum vector L = r × v
    const Vector3 angMom = relPos.cross(state.vel);
    const double L2 = angMom.magnitudeSquared();

    // Relativistic acceleration: a = - (3 * M * L^2 / r^5) * r
    const double r5 = rSq * rSq * r;
    const double factor = -3.0 * blackHole.mass * L2 / r5;
    const Vector3 accel = relPos * factor;

    return {state.vel, accel};
}

GeodesicIntegrator::State6D GeodesicIntegrator::rk4Step(const State6D& state, double h) const noexcept {
    const State6D k1 = evaluateDerivative(state);
    const State6D k2 = evaluateDerivative(state + k1 * (0.5 * h));
    const State6D k3 = evaluateDerivative(state + k2 * (0.5 * h));
    const State6D k4 = evaluateDerivative(state + k3 * h);

    const State6D combined = k1 + k2 * 2.0 + k3 * 2.0 + k4;
    return state + combined * (h / 6.0);
}

double GeodesicIntegrator::computeAdaptiveStepSize(double r) const noexcept {
    const double rh = blackHole.eventHorizonRadius();
    const double rph = blackHole.photonSphereRadius();

    // As ray approaches the photon sphere (r ≈ 3M) or horizon (r ≈ 2M),
    // decrease step size smoothly to resolve sharp orbital bending.
    double curvatureFactor = 1.0;

    if (r < rph * 2.0) {
        // Distance to horizon relative to characteristic radius 2M
        const double deltaR = std::max(0.01, (r - rh) / (2.0 * blackHole.mass));
        curvatureFactor = std::clamp(deltaR, 0.05, 1.0);
    } else if (r > escapeRadius * 0.5) {
        // In the asymptotic weak field, take larger steps
        curvatureFactor = std::min(2.5, r / (escapeRadius * 0.5));
    }

    const double h = initialStepSize * curvatureFactor;
    return std::clamp(h, minStepSize, maxStepSize);
}

void GeodesicIntegrator::integrateRay(Ray& ray) const noexcept {
    State6D state{ray.position, ray.velocity};
    ray.closestApproach = (state.pos - blackHole.position).magnitude();

    for (int step = 0; step < maxSteps; ++step) {
        const double r = (state.pos - blackHole.position).magnitude();
        ray.closestApproach = std::min(ray.closestApproach, r);

        // 1. Horizon capture condition: r <= 2M * (1 + tolerance)
        if (blackHole.isInsideHorizon(state.pos)) {
            ray.status = RayStatus::Captured;
            ray.position = state.pos;
            ray.velocity = state.vel;
            ray.steps = step;
            ray.color = {0.0, 0.0, 0.0}; // Shadow is completely black
            return;
        }

        // 2. Numerical validity check
        if (!state.pos.isFinite() || !state.vel.isFinite()) {
            ray.status = RayStatus::NumericalFailure;
            ray.steps = step;
            ray.color = {1.0, 0.0, 1.0}; // Magenta warning color for diagnostics
            return;
        }

        // 3. Adaptive step size calculation
        const double h = computeAdaptiveStepSize(r);

        // 4. Perform Runge-Kutta 4th-Order step
        const State6D nextState = rk4Step(state, h);

        // 5. Equatorial Accretion Disk intersection detection
        if (testAccretionDisk) {
            Vector3 hitPoint;
            double hitRadius;
            if (disk.checkIntersection(state.pos, nextState.pos, hitPoint, hitRadius)) {
                ray.status = RayStatus::HitAccretionDisk;
                ray.diskHitPoint = hitPoint;
                ray.diskHitRadius = hitRadius;
                const double phi = std::atan2(hitPoint.y, hitPoint.x);
                ray.color = disk.computeColor(hitRadius, phi);
                ray.position = hitPoint;
                ray.velocity = nextState.vel;
                ray.steps = step;
                return;
            }
        }

        // 6. Escape condition: r >= escapeRadius moving outward (r · v > 0)
        if (r >= escapeRadius && (state.pos - blackHole.position).dot(state.vel) > 0.0) {
            ray.status = RayStatus::Escaped;
            ray.position = nextState.pos;
            ray.velocity = nextState.vel;
            ray.steps = step;
            return;
        }

        // Advance state
        state = nextState;
    }

    // Maximum iterations reached without resolving
    ray.status = RayStatus::MaximumStepsReached;
    ray.position = state.pos;
    ray.velocity = state.vel;
    ray.steps = maxSteps;
    ray.color = {0.05, 0.05, 0.05};
}
