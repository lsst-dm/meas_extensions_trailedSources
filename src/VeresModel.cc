#include "lsst/meas/extensions/trailedSources/VeresModel.h"
#include "lsst/geom.h"
#include "lsst/afw/image.h"
#include "lsst/afw/detection.h"

namespace lsst {
namespace meas {
namespace extensions {
namespace trailedSources {

typedef afw::image::Image<float> Image;
typedef afw::image::Exposure<float> Exposure;
typedef afw::image::Image<float>::Array Array;

VeresModel::VeresModel(
    Exposure const& data
) : _sigma(data.getPsf()->computeShape().getTraceRadius()),
    _params({0.,0.,0.,0.,0.}),
    _bbox(data.getBBox()),
    _dims(data.getDimensions()),
    _image(new Image(_bbox)),
    _model(_image->getArray()),
    _data(data.getMaskedImage().getImage()->getArray()),
    _variance(data.getMaskedImage().getVariance()->getArray()) {}

double VeresModel::operator()(std::vector<double> const& params) const {

    // Unpack params
    double xc = params[0];    // Centroid x
    double yc = params[1];    // Centroid y
    double F = params[2];     // Flux
    double L = params[3];     // Trail length
    double theta = params[4]; // Angle from +x-axis

    // reset params
    // _setParams(params);

    // From computeKernelImage()
    // Compute model image and chi-squared
    double chiSq = 0.0;
    double tmp = 0.0;
    for (int yIndex = 0, yp = _image->getY0(); yIndex < _dims.getY(); ++yIndex, ++yp) {
        Array::Reference modelRow = _model[yIndex];
        Array::Reference dataRow = _data[yIndex];
        Array::Reference varRow = _variance[yIndex];
        for (int xIndex = 0, xp = _image->getX0(); xIndex < _dims.getX(); ++xIndex, ++xp) {
            modelRow[xIndex] = _computeModel(xp,yp,xc,yc,F,L,theta);
            tmp = dataRow[xIndex] - modelRow[xIndex];
            chiSq += tmp*tmp/varRow[xIndex];
        }
    }

    return chiSq;
}

double VeresModel::_computeModel(double x, double y, double xc, double yc,
                                 double F, double L, double theta) const {
    // Computes the Veres et al model at a given position (pixel)
    double xp = (x-xc)*std::cos(theta) + (y-yc)*std::sin(theta);
    double yp = (x-xc)*std::sin(theta) - (y-yc)*std::cos(theta);
    double A = std::exp(-0.5 * yp*yp / (_sigma*_sigma));
    double B = std::erf((xp+L/2) / (std::sqrt(2.0) * _sigma));
    double C = std::erf((xp-L/2) / (std::sqrt(2.0) * _sigma));
    return F * A * (B - C) / (L * 2 * std::sqrt(2.0 * geom::PI) * _sigma);
}

}}}} // lsst::meas::extensions::trailedSources