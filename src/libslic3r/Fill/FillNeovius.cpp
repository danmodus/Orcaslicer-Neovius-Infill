#include "../ClipperUtils.hpp"
#include "../ShortestPath.hpp"
#include "../Surface.hpp"
#include <cmath>
#include <algorithm>
#include <iostream>
#include "FillBase.hpp"
#include "FillNeovius.hpp"

namespace Slic3r {

// Neovius Function: 3*(cos(x)+cos(y)+cos(z)) + 4*cos(x)*cos(y)*cos(z) = 0
// We solve for y given x, z.
// let A = cos(x), C = cos(z), B = cos(y)
// 3(A + B + C) + 4ABC = 0
// 3A + 3B + 3C + 4ABC = 0
// B(3 + 4AC) = -(3A + 3C)
// B = -3(A + C) / (3 + 4AC)
// y = acos(B)

static inline double f(double x, double z_cos, bool vertical, bool flip)
{
    // Neovius is symmetric in x, y, z so functional form is similar.
    // If vertical, we swap variables conceptually, but sine/cosine symmetry holds for Neovius?
    // Actually, Neovius uses Cosines everywhere.
    // standard: 3(cx + cy + cz) + 4(cx*cy*cz) = 0

    // We pass z_cos directly.
    double cx = cos(x);
    double cz = z_cos;

    // Calculate RHS: B = -3(A+C) / (3+4AC)
    double numerator   = -3.0 * (cx + cz);
    double denominator = 3.0 + 4.0 * cx * cz;

    // Handle singularity or numerical issues?
    // Denominator is 0 if cx*cz = -0.75.

    double val = 0.0;
    if (std::abs(denominator) < 1e-9) {
        // Singularity, slope implies asymptote.
        // Return something that fits? Or clamp?
        // If denom -> 0, B -> infinity, so no solution for real y.
        return 0; // Return 0 or NAN?
    } else {
        val = numerator / denominator;
    }

    // Clamp to valid range for acos
    if (val > 1.0)
        val = 1.0;
    if (val < -1.0)
        val = -1.0;

    double y = acos(val);

    // Return y. Note: acos returns [0, pi].
    // Surface symmetry: cos(y) = val. y = +/- acos(val).
    // The previous implementation added phase shifts.
    // Gyroid: asin(a/r) + ...

    // For Neovius, we need to trace the continuous line.
    // If we just return acos(val), we get the principal value.
    // We might need to handle the flip to alternate between the positive and negative branch
    // or phase shift to make it continuous across the period.

    // Simple approach:
    if (flip)
        return -y;
    return y;
}

static inline Polyline make_wave(const std::vector<Vec2d>& one_period,
                                 double                    width,
                                 double                    height,
                                 double                    offset,
                                 double                    scaleFactor,
                                 double                    z_cos,
                                 bool                      vertical,
                                 bool                      flip)
{
    std::vector<Vec2d> points = one_period;
    double             period = 2 * PI; // Standard 2pi period

    // The input 'one_period' is usually generated for [0, 2pi].
    // But points.back()(0) gives exact period.
    if (!points.empty())
        period = points.back()(0);

    if (width != period && !points.empty()) // do not extend if already truncated
    {
        points.reserve(one_period.size() * size_t(floor(width / period)));
        // Logic to extend points...
        // For Neovius, periodicity is 2pi.
        // Copy paste loop from Gyroid:

        points.pop_back();

        size_t n = points.size();
        if (n > 0) {
            do {
                points.emplace_back(points[points.size() - n].x() + period, points[points.size() - n].y());
            } while (points.back()(0) < width - EPSILON);
        }

        // Add final point
        // points.emplace_back(Vec2d(width, f(width, z_cos, vertical, flip)));
        // We can't easily call f(width) without exact context, but assuming periodicity:
        // Actually, we should probably re-eval f if we want exact end.
        points.emplace_back(points[points.size() - n].x() + period, points[points.size() - n].y());
    }

    // and construct the final polyline to return:
    Polyline polyline;
    polyline.points.reserve(points.size());
    for (auto& point : points) {
        // Apply offset (Y-shift) ?
        // For Gyroid, y is shifted by M_PI steps.
        // For Neovius, basic period is 2pi.

        double y_val = point.y() + offset;

        // Handling flip/mirroring if needed?

        point(1) = y_val;

        // Clamp logic
        // point(1) = std::clamp(double(point.y()), 0., height); // Slicer clamping

        // Swap for vertical
        Vec2d p_out = point;
        if (vertical)
            std::swap(p_out(0), p_out(1));

        polyline.points.emplace_back((p_out * scaleFactor).cast<coord_t>());
    }

    return polyline;
}

static std::vector<Vec2d> make_one_period(double width, double scaleFactor, double z_cos, bool vertical, bool flip, double tolerance)
{
    std::vector<Vec2d> points;
    double             dx    = PI / 18.0; // 10 degree steps
    double             limit = std::min(2 * PI, width);
    points.reserve(coord_t(ceil(limit / dx)));

    for (double x = 0.; x < limit + EPSILON; x += dx) {
        points.emplace_back(Vec2d(x, f(x, z_cos, vertical, flip)));
    }

    // Refinement steps (simpler version than Gyroid for now)

    return points;
}

static Polylines make_neovius_waves(double gridZ, double density_adjusted, double line_spacing, double width, double height)
{
    const double scaleFactor = scale_(line_spacing) / density_adjusted;

    // Neovius tolerance
    const double tolerance = std::min(line_spacing / 2, FillNeovius::PatternTolerance) / unscale<double>(scaleFactor);

    // Z in domain
    const double z     = gridZ / scaleFactor;
    const double z_cos = cos(z);

    // Vertical / Horizontal preference logic
    // Neovius is symmetric, but we need to choose an axis for functional parameterization.
    // If |cos(z)| is large, B denominator (3+4AC) might cross zero?
    // Denom = 0 when A*C = -0.75.
    // We want to avoid denom approx 0.

    bool vertical = false;
    // Optimization: check which projection is better?
    // Only matters if we hit the singularity.

    // Setup period vectors
    // Neovius repeats every 2pi

    // We generate "Upper" and "Lower" branches of acos?
    // acos returns [0, pi]. -acos returns [-pi, 0].
    // The surface exists as bubbles?
    // Actually, network Neovius is continuous.

    // Let's generate one period for [0, 2pi].
    bool               flip           = false;
    std::vector<Vec2d> one_period_pos = make_one_period(width, scaleFactor, z_cos, vertical, false, tolerance);
    std::vector<Vec2d> one_period_neg = make_one_period(width, scaleFactor, z_cos, vertical, true, tolerance);

    Polylines result;

    // Periodicity in Y is also 2Pi?
    // We tile these waves.

    for (double y0 = -2 * PI; y0 < height + 2 * PI; y0 += 2 * PI) {
        // Add positive branch
        result.emplace_back(make_wave(one_period_pos, width, height, y0, scaleFactor, z_cos, vertical, false));
        // Add negative branch at same y0 (since acos is +/-)
        result.emplace_back(make_wave(one_period_neg, width, height, y0, scaleFactor, z_cos, vertical, true));
    }

    return result;
}

constexpr double FillNeovius::PatternTolerance;

void FillNeovius::_fill_surface_single(const FillParams&              params,
                                       unsigned int                   thickness_layers,
                                       const std::pair<float, Point>& direction,
                                       ExPolygon                      expolygon,
                                       Polylines&                     polylines_out)
{
    auto infill_angle = float(this->angle + (CorrectionAngle * 2 * PI) / 360.);
    if (std::abs(infill_angle) >= EPSILON)
        expolygon.rotate(-infill_angle);

    BoundingBox bb               = expolygon.contour.bounding_box();
    double      density_adjusted = std::max(0., params.density * DensityAdjust / params.multiline);
    coord_t     distance         = coord_t(scale_(this->spacing) / density_adjusted);

    bb.merge(align_to_grid(bb.min, Point(2 * PI * distance, 2 * PI * distance)));
    coord_t expand = 10 * (scale_(this->spacing));
    bb.offset(expand);

    Polylines polylines = make_neovius_waves(scale_(this->z), density_adjusted, this->spacing, ceil(bb.size()(0) / distance) + 1.,
                                             ceil(bb.size()(1) / distance) + 1.);

    for (Polyline& pl : polylines)
        pl.translate(bb.min);

    multiline_fill(polylines, params, spacing);
    polylines = intersection_pl(polylines, expolygon);

    if (!polylines.empty()) {
        const double minlength = scale_(0.8 * this->spacing);
        polylines.erase(std::remove_if(polylines.begin(), polylines.end(),
                                       [minlength](const Polyline& pl) { return pl.length() < minlength; }),
                        polylines.end());
    }

    if (!polylines.empty()) {
        size_t polylines_out_first_idx = polylines_out.size();
        chain_or_connect_infill(std::move(polylines), expolygon, polylines_out, this->spacing, params);

        if (std::abs(infill_angle) >= EPSILON) {
            for (auto it = polylines_out.begin() + polylines_out_first_idx; it != polylines_out.end(); ++it)
                it->rotate(infill_angle);
        }
    }
}

} // namespace Slic3r
