#include "../ClipperUtils.hpp"
#include "../ShortestPath.hpp"
#include "../Surface.hpp"
#include <cmath>
#include <algorithm>
#include <iostream>
#include "FillBase.hpp"
#include "FillNeovius.hpp"

namespace Slic3r {

// Neovius Surface Equation:
// 3(cos(x) + cos(y) + cos(z)) + 4 cos(x) cos(y) cos(z) = 0
//
// We solve for one variable given the other two.
// e.g. solve for y:
// 3 cos(y) + 4 cos(x) cos(z) cos(y) = -3(cos(x) + cos(z))
// cos(y) * (3 + 4 cos(x) cos(z)) = -3(cos(x) + cos(z))
// cos(y) = -3(cos(x) + cos(z)) / (3 + 4 cos(x) cos(z))
//
// Singularity:
// Denominator D = 3 + 4 cos(x) cos(z)
// D = 0 when cos(x) cos(z) = -0.75
// This occurs when cos(x) and cos(z) have opposite signs and magnitude product is 0.75.
// Near this singularity, cos(y) -> infinity, which means no real solution for y, OR the surface is vertical.
// However, since cos(y) must be in [-1, 1], the valid domain is limited.
// If |RHS| > 1, there is no surface at that (x, z).
//
// We want to trace continuous isolines.
// We can switch axes to avoid the "vertical" slopes.
// The normal vector is Gradient F.
// F = 3(cx + cy + cz) + 4 cx cy cz
// dF/dx = -sin(x) * (3 + 4 cy cz)
// dF/dy = -sin(y) * (3 + 4 cx cz)
// dF/dz = -sin(z) * (3 + 4 cx cy)
//
// If we are tracing in the XY plane (fixed z), we are following curvature.
// We want to express y(x) or x(y).
// We choose the variable to solve for such that the denominator is furthest from 0?
// Actually simpler: we just detect if calculating y(x) is stable.
// Stability breakdown happens when Slope dy/dx -> infinity.
// Implicit differentiation:
// dy/dx = - (dF/dx) / (dF/dy)
//        = - (sin(x)(3+4cy cz)) / (sin(y)(3+4cx cz))
//
// If (3+4cx cz) is close to 0, slope is infinite -> vertical tangent -> switch to x(y).

static inline double solve_neovius(double v1, double v2_cos, double v2_sin, bool& valid)
{
    // Solve for v3 (return angle) given v1 (angle) and v2 info.
    // cos(v3) = -3(cos(v1) + cos(v2)) / (3 + 4 cos(v1) cos(v2))

    double cv1 = cos(v1);
    double cv2 = v2_cos;

    double num = -3.0 * (cv1 + cv2);
    double den = 3.0 + 4.0 * cv1 * cv2;

    if (std::abs(den) < 1e-6) {
        // Singularity.
        valid = false;
        return 0.0;
    }

    double val = num / den;
    if (val > 1.0 + EPSILON || val < -1.0 - EPSILON) {
        valid = false; // no solution
        return 0.0;
    }

    // Clamp
    if (val > 1.0)
        val = 1.0;
    if (val < -1.0)
        val = -1.0;

    valid = true;
    return acos(val);
}

// Function pointer logic for subdivision?
// We need to pass context.
struct NeoviusContext
{
    double z_cos;
    double z_sin;
    // We are generating periodic waves.
    // Basic equation: 3(cx+cy+cz) + 4 cx cy cz = 0
};

static inline Vec2d evaluate_neovius(double t, const NeoviusContext& ctx, bool upper_branch)
{
    // We parameterize the curve. Ideally by arc length, but here 't' is likely x or y.
    // Let's assume we primarily drive with x, but handle steep slopes via subdivision.
    // However, for strict x-driving, we can't handle vertical segments.
    // Gyroid implementation uses a clever analytical approximation or just `asin` composition.
    // For Neovius, we don't have a clean `asin(sin...)` form.
    //
    // We will just return (x, y(x)).

    bool   valid;
    double y = solve_neovius(t, ctx.z_cos, ctx.z_sin, valid);
    // If invalid, we are in a void or singularity.
    // But Neovius is a connected minimal surface. Holes appear if we slice it.
    // Yes, isolines can terminate or form loops.

    if (!valid) {
        // Return NaN or special marks?
        // In the valid domain of the TPMS, existence is guaranteed?
        // Actually for TPMS P, D, G, holes exist. Neovius (Cellular) also has voids boundaries.
        // It separates two sub-volumes.
        // So hitting the boundary |val|=1 is normal. That's a turning point (vertical tangent).
        // If val > 1, we are outside the surface.
        return Vec2d(t, std::nan(""));
    }

    if (!upper_branch)
        y = -y;
    // y is in [-pi, pi] technically? acos gives [0, pi].
    // symmetry gives -y.

    return Vec2d(t, y);
}

// Check for slope/curvature to subdivide
static void subdivide_segment(std::vector<Vec2d>& points, Vec2d p1, Vec2d p2, const NeoviusContext& ctx, bool upper_branch, double tolerance)
{
    // Midpoint in parameter space (x)
    double t_mid       = (p1.x() + p2.x()) * 0.5;
    Vec2d  p_mid_curve = evaluate_neovius(t_mid, ctx, upper_branch);

    if (std::isnan(p_mid_curve.y()))
        return;

    // Linear midpoint
    Vec2d  p_mid_segment = (p1 + p2) * 0.5;
    double distinct_sq   = (p_mid_curve - p_mid_segment).squaredNorm();

    if (distinct_sq > tolerance * tolerance) {
        if (std::abs(p1.x() - p2.x()) > 1e-5) { // limit recursion
            subdivide_segment(points, p1, p_mid_curve, ctx, upper_branch, tolerance);
            points.push_back(p_mid_curve);
            subdivide_segment(points, p_mid_curve, p2, ctx, upper_branch, tolerance);
        }
    }
}

static std::vector<std::vector<Vec2d>> make_one_period(double width, double scaleFactor, double z_cos, double z_sin, double tolerance)
{
    std::vector<std::vector<Vec2d>> segments;
    std::vector<Vec2d>              current_segment;
    double                          dx    = PI / 60.0; // Finer base resolution
    double                          limit = std::min(2 * PI, width);

    NeoviusContext ctx{z_cos, z_sin};

    auto add_point = [&](double x) {
        bool   valid;
        double y = solve_neovius(x, z_cos, z_sin, valid);
        if (valid) {
            Vec2d p_curr(x, y);
            if (current_segment.empty()) {
                current_segment.push_back(p_curr);
            } else {
                // Recursive subdivision for smoothness
                subdivide_segment(current_segment, current_segment.back(), p_curr, ctx, true, tolerance);
                current_segment.push_back(p_curr);
            }
        } else if (!current_segment.empty()) {
            segments.push_back(std::move(current_segment));
            current_segment.clear();
        }
    };

    for (double x = 0; x <= limit + EPSILON; x += dx) {
        add_point(std::min(x, limit));
    }
    if (!current_segment.empty())
        segments.push_back(std::move(current_segment));

    return segments;
}

static inline void make_wave(Polylines&                             result,
                             const std::vector<std::vector<Vec2d>>& segments,
                             double                                 width,
                             double                                 height,
                             double                                 offset_y,
                             double                                 scaleFactor,
                             bool                                   vertical,
                             bool                                   flip_y)
{
    // Tiling logic
    if (segments.empty())
        return;

    double period      = 2 * PI;
    int    num_periods = int(ceil(width / period + 0.1));

    for (const auto& segment : segments) {
        Polyline polyline;
        for (int i = 0; i < num_periods; ++i) {
            double shift = i * period;
            for (const auto& p : segment) {
                double x = p.x() + shift;
                if (x > width + EPSILON)
                    break;

                double y_raw = p.y();
                if (flip_y)
                    y_raw = -y_raw;

                double y_final = y_raw + offset_y;
                Vec2d  out     = vertical ? Vec2d(y_final, x) : Vec2d(x, y_final);
                polyline.points.emplace_back((out * scaleFactor).cast<coord_t>());
            }
            // If the segment ends exactly at PI or 2PI, we could connect to the next tile.
            // But breaking is safer to avoid gaps/spikes when solve_neovius hits limits.
            if (!polyline.points.empty() && segment.back().x() < period - EPSILON) {
                result.push_back(std::move(polyline));
                polyline = Polyline();
            }
        }
        if (!polyline.points.empty())
            result.push_back(std::move(polyline));
    }
}

static Polylines make_neovius_waves(double gridZ, double density_adjusted, double line_spacing, double width, double height)
{
    const double scaleFactor = scale_(line_spacing) / density_adjusted;
    const double tolerance   = std::min(line_spacing / 2, FillNeovius::PatternTolerance) / unscale<double>(scaleFactor);

    // gridZ should be the absolute Z.
    const double z_rad = gridZ / scaleFactor;
    const double z_cos = cos(z_rad);
    const double z_sin = sin(z_rad);

    // Detect if we are in a predominantly "Vertical" configuration?
    // Neovius is cubic symmetric.
    // If we map X->X, Y->Y, Z->Z, we just generate standard orientation.
    // We should not swap axes arbitrarily unless it helps solve the equation.
    // But since we solve for y(x), we require |dy/dx| to be manageable.
    // If tangent is vertical, y(x) fails.
    // But we handle this via gaps or we could swap x/y if needed.
    // For specific layers (Z fixed), the isolines are level sets.

    // For Neovius, we switch axes when |cos(z)| is large to stay away from the singularities
    // of the y(x,z) solver which occur when 3 + 4 cx cz = 0.
    bool vertical = (std::abs(z_cos) > 0.5);

    auto segments = make_one_period(width, scaleFactor, z_cos, z_sin, tolerance);

    Polylines result;
    double    start_y = -2 * PI;
    double    end_y   = height + 2 * PI;

    for (double y_base = start_y; y_base < end_y; y_base += 2 * PI) {
        make_wave(result, segments, width, height, y_base, scaleFactor, vertical, false);
        make_wave(result, segments, width, height, y_base, scaleFactor, vertical, true);
    }

    return result;
}

constexpr double FillNeovius::PatternTolerance;
constexpr double FillNeovius::DensityAdjust;

void FillNeovius::_fill_surface_single(const FillParams&              params,
                                       unsigned int                   thickness_layers,
                                       const std::pair<float, Point>& direction,
                                       ExPolygon                      expolygon,
                                       Polylines&                     polylines_out)
{
    auto infill_angle = float(this->angle + (CorrectionAngle * 2 * PI) / 360.);
    if (std::abs(infill_angle) >= EPSILON)
        expolygon.rotate(-infill_angle);

    BoundingBox bb = expolygon.contour.bounding_box();

    // params.layer_height is not always populated correctly for infill?
    // params usually carries generic info.
    // However, FillBase has 'z' member.
    // We should use 'this->z' which is the print Z coordinate.

    double density_adjusted = std::max(0., params.density * DensityAdjust / params.multiline);

    // Avoid division by zero
    if (density_adjusted <= 0)
        return;

    coord_t distance = coord_t(scale_(this->spacing) / density_adjusted);
    if (distance <= 0)
        return;

    // Align grid
    bb.merge(align_to_grid(bb.min, Point(coord_t(2 * PI * distance), coord_t(2 * PI * distance))));

    coord_t expand = 10 * (scale_(this->spacing));
    bb.offset(expand);

    // Grid Z should be unscaled Z?
    // The gyroid implementation uses scale_(this->z).
    // Let's verify: Gyroid uses input z as `z` member (unscaled).
    // `scale_(this->z)` converts mm to internal units.
    // Inside make_neovius_waves, we rely on `scaleFactor = scale_(spacing)/density`.
    // And `z_rad = gridZ / scaleFactor`.
    // So `z_rad = scale_(z) / (scale_(spacing)/density) = z * density / spacing`.
    // This makes Z dimensionless relative to the cell size. Correct.

    Polylines polylines = make_neovius_waves(scale_(this->z), density_adjusted, this->spacing,
                                             ceil(bb.size()(0) / (double) distance) + 1., // Width in dimensionless? No, wait.
                                             ceil(bb.size()(1) / (double) distance) + 1.  // Height
    );

    // Note: The make_neovius_waves Width/Height arguments above logic seems to assume 'distance' is the unit scale?
    // In Gyroid:
    //   distance = scale(spacing)/density
    //   width passed = ceil(bb.size / distance) ...
    // Inside make_gyroid:
    //   scaleFactor = scale(spacing)/density = distance (in value)
    //   points calculated in dimensionless units [0, width]
    //   then multiplied by scaleFactor.
    // So YES, we pass dimensionless width/height.

    for (Polyline& pl : polylines)
        pl.translate(bb.min);

    multiline_fill(polylines, params, spacing);
    {
        ExPolygons expolygons_off;
        if (this->overlap > 0)
            expolygons_off = offset_ex(expolygon, float(scale_(this->overlap)));
        else
            expolygons_off.push_back(expolygon);
        polylines = intersection_pl(polylines, expolygons_off);
    }

    if (!polylines.empty()) {
        const double minlength = scale_(0.2 * this->spacing); // More lenient pruning
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
