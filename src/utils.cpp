#include "utils.hpp"

std::mt19937 gen(42);

double gen_angle() {
    std::uniform_real_distribution<> dis(0.0, 2.0 * PI);

    return dis(gen);
}

double gen_radius(double r, double h) {
    std::uniform_real_distribution<> dis(0.0, 1.0);
    double u = dis(gen);

    return std::sqrt(r * r + u * ((r + h) * (r + h) - r * r));
}

double cross_prod(const glm::dvec2 &a, const glm::dvec2 &b) {
    return a.x * b.y - b.x * a.y;
}

double point_to_segment_dist(const glm::dvec2 &a, const glm::dvec2 &b, const glm::dvec2 &p) {
    glm::dvec2 ap = p - a;
    glm::dvec2 bp = p - b;
    glm::dvec2 ab = b - a;
    glm::dvec2 ba = a - b;

    if (glm::dot(ap, ab) >= 0 && glm::dot(bp, ba) >= 0) {
        return 0.5 * glm::abs(cross_prod(ap, ab)) / glm::length(ab);
    } else {
        return glm::min(glm::distance(a, p), glm::distance(b, p));
    }
}

// A = U * D * transpose(V)
void svd(const glm::dmat2& A, glm::dmat2& U, glm::dmat2& D, glm::dmat2& V) {
    double theta = std::atan2(
        A[0][1] - A[1][0],
        A[0][0] + A[1][1]
    );

    double ct = std::cos(theta), st = std::sin(theta);
    glm::dmat2 R(
        glm::dvec2(ct,  st),
        glm::dvec2(-st, ct)
    );

    glm::dmat2 S = glm::transpose(R) * A;

    double phi = 0.5 * std::atan2(
        S[1][0] + S[0][1], 
        S[0][0] - S[1][1]
    );

    double cp = std::cos(phi), sp = std::sin(phi);
    V = glm::dmat2(
        glm::dvec2(cp, sp), 
        glm::dvec2(-sp, cp)
    );

    U = R * V;

    double tr = trace(S);
    double diff  = std::sqrt(
        0.25 * (S[0][0] - S[1][1]) * (S[0][0] - S[1][1]) + 
        S[1][0] * S[0][1]
    );
    double sv1 = 0.5 * tr + diff;
    double sv2 = 0.5 * tr - diff;

    if (sv1 < 0.0) { sv1 = -sv1; U[0] = -U[0]; }
    if (sv2 < 0.0) { sv2 = -sv2; U[1] = -U[1]; }

    D = glm::dmat2(sv1, 0.0, 0.0, sv2);
}

glm::dmat2 extract_rotation(const glm::dmat2& A) {
    glm::dmat2 U, D, V;
    svd(A, U, D, V);

    return U * glm::transpose(V);
}

glm::dmat2 moore_penrose_inverse(const glm::dmat2& A, double eps) {
    glm::dmat2 U, D, V;
    svd(A, U, D, V);

    double d1 = D[0][0];
    double d2 = D[1][1];

    double inv1 = (d1 > eps) ? 1.0 / d1 : 0.0;
    double inv2 = (d2 > eps) ? 1.0 / d2 : 0.0;

    glm::dmat2 Dinv(inv1, 0.0, 0.0, inv2);

    return V * Dinv * glm::transpose(U);
}

double trace(const glm::dmat2 &matrix) {
    return matrix[0][0] + matrix[1][1];
}

glm::dmat2 tensor_product(const glm::dvec2 &a, const glm::dvec2 &b) {
    return glm::dmat2(
        a.x * b.x, a.y * b.x, 
        a.x * b.y, a.y * b.y
    );
}

double calc_orbital_velocity(double M, double R) {
    return sqrt(G * M / R);
}