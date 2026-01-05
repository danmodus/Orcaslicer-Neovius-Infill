#ifndef slic3r_FillNeovius_hpp_
#define slic3r_FillNeovius_hpp_

#include "../libslic3r.h"

#include "FillBase.hpp"

namespace Slic3r {

class FillNeovius : public Fill
{
public:
    FillNeovius() {}
    Fill* clone() const override { return new FillNeovius(*this); }

    // require bridge flow since most of this pattern hangs in air
    bool use_bridge_flow() const override { return false; }
    bool is_self_crossing() override { return false; }

    // Correction applied to regular infill angle to maximize printing
    // speed in default configuration (degrees)
    static constexpr float CorrectionAngle = -45.;

    // Density adjustment to have a good %of weight.
    // Neovius has higher surface area density than Gyroid, so we decrease this factor
    // (relative to Gyroid's 2.44) to increase spacing and maintain correct solid volume fraction.
    // Approx factor: 2.44 (Gyroid) * (Area_Gyroid / Area_Neovius) ~= 2.44 * (2.4 / 3.2) ~= 1.8. 
    // We choose 2.0 to be slightly robust.
    static constexpr double DensityAdjust = 2.0;

    // Neovius upper resolution tolerance (mm^-2)
    static constexpr double PatternTolerance = 0.5;

protected:
    void _fill_surface_single(const FillParams&              params,
                              unsigned int                   thickness_layers,
                              const std::pair<float, Point>& direction,
                              ExPolygon                      expolygon,
                              Polylines&                     polylines_out) override;
};

} // namespace Slic3r

#endif // slic3r_FillNeovius_hpp_
