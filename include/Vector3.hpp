#ifndef BLACKHOLE_VECTOR3_HPP
#define BLACKHOLE_VECTOR3_HPP

#include <cmath>
#include <iostream>
#include <sstream>

/**
 * @brief Represents a 3D Euclidean vector with double precision.
 * 
 * Provides vector arithmetic, dot/cross products, normalization,
 * and safety checks against division-by-zero and non-finite values.
 */
class Vector3 {
public:
    double x{0.0};
    double y{0.0};
    double z{0.0};

    // Constructors
    constexpr Vector3() noexcept = default;
    constexpr Vector3(double x_, double y_, double z_) noexcept : x(x_), y(y_), z(z_) {}
    constexpr explicit Vector3(double scalar) noexcept : x(scalar), y(scalar), z(scalar) {}

    // Vector addition and subtraction
    constexpr Vector3 operator+(const Vector3& rhs) const noexcept {
        return {x + rhs.x, y + rhs.y, z + rhs.z};
    }

    constexpr Vector3 operator-(const Vector3& rhs) const noexcept {
        return {x - rhs.x, y - rhs.y, z - rhs.z};
    }

    constexpr Vector3 operator-() const noexcept {
        return {-x, -y, -z};
    }

    // Component-wise multiplication
    constexpr Vector3 operator*(const Vector3& rhs) const noexcept {
        return {x * rhs.x, y * rhs.y, z * rhs.z};
    }

    // Scalar multiplication and division
    constexpr Vector3 operator*(double scalar) const noexcept {
        return {x * scalar, y * scalar, z * scalar};
    }

    Vector3 operator/(double scalar) const {
        if (std::abs(scalar) < 1e-15) {
            return {0.0, 0.0, 0.0};
        }
        const double inv = 1.0 / scalar;
        return {x * inv, y * inv, z * inv};
    }

    // Compound assignment operators
    Vector3& operator+=(const Vector3& rhs) noexcept {
        x += rhs.x;
        y += rhs.y;
        z += rhs.z;
        return *this;
    }

    Vector3& operator-=(const Vector3& rhs) noexcept {
        x -= rhs.x;
        y -= rhs.y;
        z -= rhs.z;
        return *this;
    }

    Vector3& operator*=(double scalar) noexcept {
        x *= scalar;
        y *= scalar;
        z *= scalar;
        return *this;
    }

    Vector3& operator/=(double scalar) {
        if (std::abs(scalar) < 1e-15) {
            x = y = z = 0.0;
            return *this;
        }
        const double inv = 1.0 / scalar;
        x *= inv;
        y *= inv;
        z *= inv;
        return *this;
    }

    // Equality checks
    constexpr bool operator==(const Vector3& rhs) const noexcept {
        return x == rhs.x && y == rhs.y && z == rhs.z;
    }

    constexpr bool operator!=(const Vector3& rhs) const noexcept {
        return !(*this == rhs);
    }

    // Dot product
    [[nodiscard]] constexpr double dot(const Vector3& rhs) const noexcept {
        return x * rhs.x + y * rhs.y + z * rhs.z;
    }

    // Cross product
    [[nodiscard]] constexpr Vector3 cross(const Vector3& rhs) const noexcept {
        return {
            y * rhs.z - z * rhs.y,
            z * rhs.x - x * rhs.z,
            x * rhs.y - y * rhs.x
        };
    }

    // Length calculations
    [[nodiscard]] double magnitudeSquared() const noexcept {
        return dot(*this);
    }

    [[nodiscard]] double magnitude() const noexcept {
        return std::sqrt(magnitudeSquared());
    }

    [[nodiscard]] double norm() const noexcept {
        return magnitude();
    }

    // Normalization with safe handling of near-zero magnitudes
    [[nodiscard]] Vector3 normalized() const noexcept {
        const double magSq = magnitudeSquared();
        if (magSq < 1e-30) {
            return {0.0, 0.0, 0.0};
        }
        const double invMag = 1.0 / std::sqrt(magSq);
        return {x * invMag, y * invMag, z * invMag};
    }

    void normalize() noexcept {
        *this = normalized();
    }

    // Vector reflection
    [[nodiscard]] Vector3 reflect(const Vector3& normal) const noexcept {
        return *this - normal * (2.0 * dot(normal));
    }

    // Numerical validity check
    [[nodiscard]] bool isFinite() const noexcept {
        return std::isfinite(x) && std::isfinite(y) && std::isfinite(z);
    }

    // String formatting
    [[nodiscard]] std::string toString() const {
        std::ostringstream oss;
        oss << "(" << x << ", " << y << ", " << z << ")";
        return oss.str();
    }
};

// Global helper operators
inline constexpr Vector3 operator*(double scalar, const Vector3& v) noexcept {
    return v * scalar;
}

inline constexpr double dot(const Vector3& a, const Vector3& b) noexcept {
    return a.dot(b);
}

inline constexpr Vector3 cross(const Vector3& a, const Vector3& b) noexcept {
    return a.cross(b);
}

inline std::ostream& operator<<(std::ostream& os, const Vector3& v) {
    os << v.toString();
    return os;
}

#endif // BLACKHOLE_VECTOR3_HPP
