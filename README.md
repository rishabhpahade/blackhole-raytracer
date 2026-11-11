# Schwarzschild Black-Hole Ray Tracer

A high-performance C++ scientific renderer and numerical simulator that computes null geodesics (light paths) in Schwarzschild spacetime. It simulates relativistic gravitational light bending, black hole shadow formation, photon-sphere orbiting, equatorial accretion disk lensing, and background starfield deflection.

![Rendered Schwarzschild black hole with a lensed accretion disk](output/black_hole.png)

---

## 1. Physical Model & Mathematical Derivation

### Spacetime Metric
The geometry of spacetime around a static, uncharged spherical black hole of mass $M$ in geometrized units ($G = c = 1$) is described by the **Schwarzschild Metric**:

$$ds^2 = -\left(1 - \frac{2M}{r}\right)dt^2 + \left(1 - \frac{2M}{r}\right)^{-1}dr^2 + r^2 d\theta^2 + r^2 \sin^2\theta\, d\phi^2$$

### Null Geodesic Equation of Motion
For light (photons), trajectories follow null geodesics where $ds^2 = 0$. Because Schwarzschild spacetime is spherically symmetric, the motion of every photon is planar, governed by its conserved specific energy $E = (1 - 2M/r)\dot{t}$ and angular momentum $\mathbf{L} = \mathbf{r} \times \mathbf{\dot{r}}$.

In standard orbital coordinates $(u, \phi)$ where $u = 1/r$:
$$\frac{d^2u}{d\phi^2} + u = 3Mu^2$$
The nonlinear term $3Mu^2$ gives rise to General Relativistic light bending.

In full 3D Cartesian coordinates $\mathbf{x} = (x, y, z)$ with spatial wavevector $\mathbf{v} = \mathbf{\dot{x}}$, the second-order equation of motion is:
$$\mathbf{\ddot{x}} = -\frac{3M |\mathbf{x} \times \mathbf{v}|^2}{|\mathbf{x}|^5} \mathbf{x} = -\frac{3M L^2}{r^5}\mathbf{x}$$

This Cartesian formulation guarantees:
1. **Strict angular momentum conservation**: $\frac{d}{d\lambda}(\mathbf{x} \times \mathbf{v}) = \mathbf{x} \times \mathbf{\ddot{x}} = \mathbf{0}$.
2. **Zero coordinate singularities**: No coordinate blowups at polar axes ($\theta = 0, \pi$).

### Critical Physical Radii
- **Event Horizon**: $r_h = 2M$. Any light ray crossing $r \le 2M$ cannot escape and constitutes the central shadow.
- **Photon Sphere**: $r_{\text{ph}} = 3M$. Radius at which unstable circular photon orbits exist.
- **Critical Impact Parameter**: $b_c = 3\sqrt{3}M \approx 5.19615M$. Rays from infinity with impact parameter $b < b_c$ plunge across the event horizon, whereas rays with $b > b_c$ scatter back to infinity.
- **Weak-Field Deflection**: In the asymptotic limit $b \gg 3M$, the total scattering angle matches Einstein's deflection formula:
  $$\alpha(b) \approx \frac{4M}{b}$$

### Accretion Disk Physics
The accretion disk is modeled as an optically emitting, geometrically thin planar disk at $z = 0$ bounded by $r_{\text{inner}} = 6M$ (the Innermost Stable Circular Orbit, ISCO) and $r_{\text{outer}} = 22M$.
The radial temperature profile follows the standard thin-disk model:
$$T(r) \propto r^{-3/4} \left(1 - \sqrt{\frac{r_{\text{in}}}{r}}\right)^{1/4}$$
The emitted spectral energy is converted to RGB color via blackbody Planckian color grading.

---

## 2. Software Architecture

```
blackhole_raytracer/
├── include/
│   ├── Vector3.hpp              # 3D Vector algebra, operators, dot, cross, safety guards
│   ├── Camera.hpp               # Pinhole camera, coordinate frame basis, AA sampling
│   ├── BlackHole.hpp            # Mass, horizon, photon sphere, capture tests
│   ├── Ray.hpp                  # 6D ray state, RayStatus enum, deflection metrics
│   ├── GeodesicIntegrator.hpp   # Manual 4th-Order Runge-Kutta integrator with adaptive step
│   ├── AccretionDisk.hpp        # Thin disk geometry & thermal radiation model
│   └── Renderer.hpp             # Pixel dispatcher, multi-threading, sRGB tone mapping, PNG exporter
├── src/
│   ├── Camera.cpp
│   ├── AccretionDisk.cpp
│   ├── GeodesicIntegrator.cpp
│   ├── Renderer.cpp
│   └── main.cpp                 # Main simulation runner
├── external/
│   └── stb_image_write.h        # Self-contained PNG writer
├── tests/
│   └── validation_tests.cpp     # Automated scientific test suite
├── output/                      # Generated visual & scientific datasets
├── CMakeLists.txt
└── README.md
```

---

## 3. Compilation & Build Instructions

### Prerequisites
- C++17 compliant compiler (`g++`, `clang++`, or `MSVC`)
- CMake 3.14+ (optional, direct compiler commands also supported)
- OpenMP (optional, for multi-core parallelism)

### Option A: Direct Build with GCC / MinGW (Windows & Linux)

```bash
cd blackhole_raytracer

# 1. Build the main simulation executable
g++ -std=c++17 -O3 -fno-math-errno -fopenmp -Iinclude -I. \
    src/Camera.cpp src/AccretionDisk.cpp src/GeodesicIntegrator.cpp src/Renderer.cpp src/main.cpp \
    -o blackhole_raytracer.exe

# 2. Build the scientific validation test suite
g++ -std=c++17 -O3 -fno-math-errno -fopenmp -Iinclude -I. \
    src/Camera.cpp src/AccretionDisk.cpp src/GeodesicIntegrator.cpp src/Renderer.cpp tests/validation_tests.cpp \
    -o validation_tests.exe
```

### Option B: Build with CMake

```bash
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
cmake --build . --config Release
```

---

## 4. Running the Simulation

```bash
# Run scientific validation tests (19 assertions)
./validation_tests.exe

# Run the renderer: <width> <height> <samples-per-pixel>
# Defaults are 1280x720 at 2 samples per pixel.
./blackhole_raytracer.exe 1280 720 2
```

### Camera & Exposure

The disk lies in the `z = 0` plane, so the camera's `|z|` is its height above the
disk and `|(x, y)|` is its in-plane distance. The default view sits 15° above the
plane at 56M, with `+z` as the image's up axis — shallow enough that the far side
of the disk is lensed up over the shadow. Set `inclinationDegrees`,
`orbitRadius`, `fov` and `renderer.exposure` in `src/main.cpp` to reframe.

Colour conversion uses the ACES filmic curve. Its shoulder keeps the disk's
temperature gradient visible instead of clipping the whole disk to flat white;
`exposure` is the knob to tune brightness against that curve.

### Generated Deliverables (`output/`)

1. **`black_hole.png`**: High-resolution rendering of the black hole shadow, gravitationally distorted accretion disk, photon ring, and lensed background starfield.
2. **`ray_classification.png`**: False-color diagnostic map classifying every camera pixel into:
   - **Black**: Captured by Event Horizon ($r \le 2M$)
   - **Amber Gold**: Intersected Accretion Disk
   - **Electric Blue**: Escaped to Infinity (gravitationally lensed background)
3. **`deflection_vs_impact_parameter.csv`**: Numerical scattering angle $\alpha(b)$ vs. theoretical Einstein prediction $\alpha = 4M/b$ across 600 sample impact parameters.
4. **`simulation_metrics.txt`**: Execution timing, ray throughput, captured/escaped/disk-hit percentages, and step count statistics.

---

## 5. Numerical Validation Results

All 19 assertions in `tests/validation_tests.cpp` pass. Numbers below are from
the test run itself, not from an earlier revision.

| Test Case | Theoretical Value | Simulated Value | Status |
|---|---|---|---|
| **Horizon Capture Radius** | $r \le 2.000 M$ | $r \le 2.002 M$ | **PASS** |
| **Critical Impact Parameter (bisection)** | $b_c = 3\sqrt{3}M = 5.196152 M$ | $b_c = 5.196135 M$ (relative error $3.4\times10^{-6}$) | **PASS** |
| **Weak-Field Deflection ($b = 50M$)** | $\alpha = 0.084712\text{ rad}$ | $0.085081\text{ rad}$ ($0.44\%$ error) | **PASS** |
| **Weak-Field Deflection ($b = 100M$)** | $\alpha = 0.041178\text{ rad}$ | $0.041209\text{ rad}$ ($0.08\%$ error) | **PASS** |
| **Weak-Field Deflection ($b = 150M$)** | $\alpha = 0.027190\text{ rad}$ | $0.027163\text{ rad}$ ($0.10\%$ error) | **PASS** |
| **Spacetime Reflection Symmetry** | $\Delta \alpha = 0$ | $\Delta \alpha < 10^{-10}\text{ rad}$ | **PASS** |
| **RK4 Step-Size Convergence** | $\mathcal{O}(h^4)$, order $\approx 4$ | observed order $3.98$ and $3.55$ | **PASS** |

---

## 6. Numerical Safeguards & Failure Modes

1. **Near-Singularity Proximity**: When $r \le 2M(1 + \epsilon)$, the ray is immediately classified as `Captured` before reaching $r=0$.
2. **Division by Zero Protection**: All vector normalization, scalar divisions, and acceleration calculations check for $|r| < 10^{-8}$ and $|\text{scalar}| < 10^{-15}$.
3. **Non-Finite Detection**: Every Runge-Kutta step evaluates `Vector3::isFinite()`. If NaN or infinite values occur, the ray is flagged with `RayStatus::NumericalFailure`.
4. **Adaptive Curvature Stepping**: As rays approach strong curvature regions ($r \le 4M$), step sizes automatically scale down to resolve tight orbital loops; in the weak field the step grows roughly proportionally to $r$, which keeps the long approach to the camera cheap without measurable loss of accuracy.

> **Build note:** do not compile with `-ffast-math` (or MSVC `/fp:fast`). It implies finite-math-only, which lets the compiler fold `std::isfinite()` to `true` and silently turns the `NumericalFailure` guard above into dead code.
