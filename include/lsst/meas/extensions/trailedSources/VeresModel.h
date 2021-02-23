#ifndef LSST_MEAS_EXTENSIONS_TRAILEDSOURCES_VERES_H
#define LSST_MEAS_EXTENSIONS_TRAILEDSOURCES_VERES_H

#include "lsst/pex/config.h"
#include "lsst/afw/image.h"
#include "lsst/afw/geom.h"
#include "lsst/geom.h"

namespace lsst {
namespace meas {
namespace extensions {
namespace trailedSources {

class VeresModel {
public:
    typedef afw::image::Image<float> Image;
    typedef afw::image::Exposure<float> Exposure;
    typedef afw::image::Image<float>::Array Array;

    explicit VeresModel(Exposure const& data, std::vector<double> const& params);

    // Update parameters, and compute model image and chi-squared (for passing to optimizer)
    double operator()(std::vector<double> const& params) const;

    /// Return the current model image
    std::shared_ptr<Image> getModelImage() const { return _image; }

    /// Return the PSF sigma
    double getSigma() const { return _sigma; }

private:
    double _computeModel(double x, double y, double xc, double yc,
                         double F, double L, double theta) const;

    double _sigma;
    std::vector<double> _params;
    lsst::geom::Box2I _bbox;
    lsst::geom::Extent2I _dims;
    std::shared_ptr<Image> _image;
    Array _model;
    Array _data;
    Array _variance;
};

}}}} // lsst::meas::extensions::trailedSources

#endif // LSST_MEAS_EXTENSIONS_TRAILEDSOURCES_VERES_H