#include "lsst/meas/extensions/trailedSources/VeresModel.h"
#include "lsst/geom.h"
#include "lsst/afw/image.h"
#include "lsst/afw/detection.h"

namespace lsst {
namespace meas {
namespace extensions {
namespace trailedSources {

typedef afw::image::Image<float> Image;
typedef afw::image::MaskedImage<float> MaskedImage;

VeresModel::VeresModel(
    Exposure const& data,
    std::vector<double> const& params
) : _sigma(data.getPsf()->computeShape().getTraceRadius()),
    _params(params),
    _bbox(data.getBBox()),
    _dims(data.getDimensions()),
    _image(new Image(_bbox)),
    _data(std::make_shared<MaskedImage>(data.getMaskedImage())) {}

double VeresModel::operator()(std::vector<double> const& params) const {
    typedef afw::image::Image<float>::Array Array;

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
    Array model = _image->getArray();
    Array data = _data->getImage()->getArray();
    Array var = _data->getVariance()->getArray();
    double chiSq = 0.0;
    double tmp = 0.0;
    for (int yIndex = 0, yp = _image->getY0(); yIndex < _dims.getY(); ++yIndex, ++yp) {
        Array::Reference modelRow = model[yIndex];
        Array::Reference dataRow = data[yIndex];
        Array::Reference varRow = var[yIndex];
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
    double xp = (x-xc)*std::cos(theta) - (y-yc)*std::sin(theta);
    double yp = (x-xc)*std::sin(theta) + (y-yc)*std::cos(theta);
    double A = std::exp(-0.5 * yp*yp / (_sigma*_sigma));
    double B = std::erf((xp+L/2) / (std::sqrt(2.0) * _sigma));
    double C = std::erf((xp-L/2) / (std::sqrt(2.0) * _sigma));
    return F * A * (B - C) / (L * 2);
}

}}}} // lsst::meas::extensions::trailedSources