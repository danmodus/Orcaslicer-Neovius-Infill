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
// This is recursive
static void subdivide_segment(
    std::vector<Vec2d>& points, Vec2d p1, Vec2d p2, const NeoviusContext& ctx, bool upper_branch, double tolerance_sq)
{
    // Midpoint in parameter space (x)
    double t_mid       = (p1.x() + p2.x()) * 0.5;
    Vec2d  p_mid_curve = evaluate_neovius(t_mid, ctx, upper_branch);

    if (std::isnan(p_mid_curve.y())) {
        // If midpoint is invalid, we might be crossing a boundary.
        // We can't easily refine this with simple recursion assuming connectedness.
        // But assuming we are inside a valid period, let's just stop or assume linear.
        return;
    }

    // Check distance from linear segment midpoint to curve midpoint
    Vec2d  p_mid_segment = (p1 + p2) * 0.5;
    double distinct_sq   = (p_mid_curve - p_mid_segment).squaredNorm();

    if (distinct_sq > tolerance_sq) {
        // Recurse
        if (std::abs(t_mid - p1.x()) > 1e-4) { // limit recursion depth
            subdivide_segment(points, p1, p_mid_curve, ctx, upper_branch, tolerance_sq);
            points.push_back(p_mid_curve);
            subdivide_segment(points, p_mid_curve, p2, ctx, upper_branch, tolerance_sq);
        }
    }
}

static std::vector<Vec2d> make_one_period(double width, double scaleFactor, double z_cos, double z_sin, double tolerance)
{
    // Generate Neovius profile for x in [0, 2pi]
    // Due to symmetry, we can generate just one branch and flip?
    // Equation is symmetric Ox, Oy.

    std::vector<Vec2d> points;
    double             dx    = PI / 9.0; // 20 degrees coarse steps
    double             limit = std::min(2 * PI, width);

    // We only generate if the slice of surface exists here.
    // Check z condition?
    // 3(cx+cy+cz) + 4cx cy cz = 0
    // At x=0: 3(1+cy+cz) + 4 cy cz = 0 => 3+3cz + cy(3+4cz) = 0 => cy = -3(1+cz)/(3+4cz).

    NeoviusContext ctx{z_cos, z_sin};

    // We need to trace the "Upper" branch (y > 0) and "Lower" branch (y < 0)?
    // Or just one and handle the other by symmetry in make_wave.
    // Let's generate the positive branch (acos output).
    // Note: this may be discontinuous if the isoline breaks.

    bool  last_valid = false;
    Vec2d last_p;

    // Initial point
    {
        bool   valid;
        double y = solve_neovius(0, z_cos, z_sin, valid);
        if (valid) {
            last_p = Vec2d(0, y);
            points.push_back(last_p);
            last_valid = true;
        }
    }

    for (double x = dx; x <= limit + EPSILON; x += dx) {
        bool   valid;
        double y = solve_neovius(x, z_cos, z_sin, valid);

        if (valid) {
            Vec2d p_curr(x, y);

            if (last_valid) {
                // Determine if we should subdivide
                // We use the tolerance check similar to Gyroid

                // But wait, Gyroid generates 'points' first then refines passes.
                // We'll use the recursive approach inline or post-pass.
                // Let's use the post-pass approach from Gyroid for consistency and speed.
                points.push_back(p_curr);
            } else {
                // Gap in validity? Start new segment or jump?
                // For infill, we usually want continuous lines.
                // If there's a gap, it means the surface doesn't exist there.
                // effectively we reset.
                points.push_back(p_curr);
            }
            last_p     = p_curr;
            last_valid = true;
        } else {
            last_valid = false;
            // What if valid region ends?
            // We just skip.
        }
    }

    // Ensure end point at exactly 2pi if needed for periodicity
    if (width >= 2 * PI - EPSILON && points.size() > 1) {
        // Force 2pi point to match 0 point y (periodicity)
        bool   valid;
        double y0 = solve_neovius(0, z_cos, z_sin, valid);
        if (valid) {
            Vec2d endP(2 * PI, y0);
            if ((points.back() - endP).norm() > 1e-4) {
                points.push_back(endP);
            }
        }
    }

    // Refinement pass
    // We iterate the points vector and insert points where curvature is high.
    // However, vector insertion is slow. Better rebuild or use list.
    // Gyroid.cpp uses a loop with insertions.

    for (int pass = 0; pass < 2; ++pass) // 2 passes of refinement
    {
        size_t size = points.size();
        for (size_t i = 1; i < size; ++i) {
            Vec2d& p_left  = points[i - 1];
            Vec2d& p_right = points[i];

            // Check middle
            double x_mid       = (p_left.x() + p_right.x()) * 0.5;
            Vec2d  p_mid_curve = evaluate_neovius(x_mid, ctx, true); // true = upper branch logic

            if (std::isnan(p_mid_curve.y()))
                continue; // In hole

            // Linear midpoint
            Vec2d p_mid_seg = (p_left + p_right) * 0.5;

            if ((p_mid_curve - p_mid_seg).squaredNorm() > scale_(tolerance) * scale_(tolerance)) {
                // Insert
                points.emplace_back(p_mid_curve);
            }
        }

        // Sort back into place
        std::sort(points.begin(), points.end(), [](const Vec2d& a, const Vec2d& b) { return a.x() < b.x(); });
    }

    return points;
}

static inline Polyline make_wave(
    const std::vector<Vec2d>& one_period, double width, double height, double offset_y, double scaleFactor, bool vertical, bool flip_y)
{
    // Tiling logic
    Polyline polyline;
    if (one_period.empty())
        return polyline;

    double period = 2 * PI;

    // We assemble the wave by repeating 'one_period'
    // For Neovius, y(x) + 2*pi is also a solution?
    // The cos(y) is periodic. So y + 2*k*pi is valid.

    std::vector<Vec2d> points;
    double             current_x_base = 0;

    // How many periods to cover width?
    int num_periods = int(ceil(width / period + 0.1));

    points.reserve(one_period.size() * num_periods);

    for (int i = 0; i < num_periods; ++i) {
        double shift = i * period;
        for (const auto& p : one_period) {
            double x = p.x() + shift;
            if (x > width + EPSILON)
                break;

            double y_raw = p.y();
            if (flip_y)
                y_raw = -y_raw; // The other branch of acos

            // Apply Offset (Y-position of the wave center)
            double y_final = y_raw + offset_y;

            // Clamp to bounding box height?
            // Actually slicing usually clips later, but clamping helps robustness
            // y_final = std::clamp(y_final, 0.0, height);

            // Transform to output coords
            Vec2d out = vertical ? Vec2d(y_final, x) : Vec2d(x, y_final);

            points.emplace_back(out);
        }

        // Add a break between periods? No, it should be continuous.
        // p.back() of period i should match p.front() of period i+1 (mod 2pi)
    }

    // Convert to scaled points
    polyline.points.reserve(points.size());
    for (const auto& pt : points) {
        polyline.points.emplace_back((pt * scaleFactor).cast<coord_t>());
    }

    return polyline;
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

    bool vertical = false; // Standard orientation

    // Generate one period of the "positive" branch y = +acos(...)
    std::vector<Vec2d> period_pos = make_one_period(width, scaleFactor, z_cos, z_sin, tolerance);

    Polylines result;

    // Neovius cells repeat every 2*PI in space.
    // We need to cover the 'height' with these waves.
    // For each Y-period (2PI), we have a positive branch centered at 0,
    // and potentially other branches?
    // cos(y) = val. Solutions are y = +/- acos(val) + 2*k*PI.
    // So we need:
    // 1. y = acos(val) + 2kPI
    // 2. y = -acos(val) + 2kPI

    double start_y = -2 * PI; // Start a bit outside
    double end_y   = height + 2 * PI;

    for (double y_base = start_y; y_base < end_y; y_base += 2 * PI) {
        // Positive branch
        Polyline p1 = make_wave(period_pos, width, height, y_base, scaleFactor, vertical, false);
        if (!p1.points.empty())
            result.emplace_back(std::move(p1));

        // Negative branch
        Polyline p2 = make_wave(period_pos, width, height, y_base, scaleFactor, vertical, true);
        if (!p2.points.empty())
            result.emplace_back(std::move(p2));
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
