// -*- LSST-C++ -*-
/*
 * This file is part of meas_extensions_trailedSources.
 *
 * Developed for the LSST Data Management System.
 * This product includes software developed by the LSST Project
 * (https://www.lsst.org).
 * See the COPYRIGHT file at the top-level directory of this distribution
 * for details of code ownership.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

#include "lsst/pex/exceptions.h"
#include "lsst/geom/Box.h"
#include "lsst/afw/geom/ellipses/Quadrupole.h"
#include "lsst/afw/image.h"
#include "lsst/afw/detection/Psf.h"
#include "lsst/afw/table/Source.h"
#include "lsst/afw/math/Function.h"
#include "lsst/afw/math/minimize.h"
#include "lsst/afw/math/Statistics.h"

#include "lsst/meas/extensions/trailedSources/ConvolvedTrailedSource.h"

namespace lsst {
namespace meas {
namespace extensions {
namespace trailedSources {

namespace {
lsst::meas::base::FlagDefinitionList flagDefinitions;
}

base::FlagDefinition const ConvolvedTrailedSourceAlgorithm::FAILURE = flagDefinitions.addFailureFlag();

base::FlagDefinitionList const& ConvolvedTrailedSourceAlgorithm::getFlagDefinitions() { return flagDefinitions; }

ConvolvedTrailedSourceAlgorithm::ConvolvedTrailedSourceAlgorithm(
    Control const& ctrl,
    std::string const& name,
    afw::table::Schema& schema
) : _ctrl(ctrl),
    _doc("Convolved trailed source (Veres et al 2012)"),
    _schema(schema),
    _xHeadKey(schema.addField<double>(name + "_x0", _doc)),
    _yHeadKey(schema.addField<double>(name + "_y0", _doc)),
    _xTailKey(schema.addField<double>(name + "_x1", _doc)),
    _yTailKey(schema.addField<double>(name + "_y1", _doc)),
    _totalFluxKey(schema.addField<double>(name + "_totalFlux", _doc)),
    _sourceFluxKey(schema.addField<double>(name + "_sourceFlux", _doc)),
    _chiSqKey(schema.addField<double>(name + "_chiSq", _doc)),
    _flagHandler(base::FlagHandler::addFields(schema, name, getFlagDefinitions())),
    _centroidExtractor(schema, name) {}

template <typename ReturnT>
double ConvolvedTrailedSourceFunction2<ReturnT>::_computeModel(
                                                  std::vector<double> const& params,
                                                  double sigma,
                                                  double x,
                                                  double y) const {
    // Unpack params
    double x0 = params[0];
    double y0 = params[1];
    double F = params[2];
    double L = params[3];
    double theta = params[4];

    // Computes the Veres et al model at a given position (pixel)
    double xp = (x-x0)*std::cos(theta) - (y-y0)*std::sin(theta);
    double yp = (x-x0)*std::sin(theta) + (y-y0)*std::cos(theta);
    double A = std::exp(-0.5 * yp*yp / (sigma*sigma));
    double B = std::erf((xp+L/2) / (std::sqrt(2.0) * sigma));
    double C = std::erf((xp-L/2) / (std::sqrt(2.0) * sigma));
    return 2.0 * F * A * (B - C) / L;
}

double ConvolvedTrailedSourceAlgorithm::computeModel(double x, double y, double x0, double y0, double F,
                                                      double L, double theta, double sigma) const {
    // Computes the Veres et al model at a given position (pixel)
    double xp = (x-x0)*std::cos(theta) - (y-y0)*std::sin(theta);
    double yp = (x-x0)*std::sin(theta) + (y-y0)*std::cos(theta);
    double A = std::exp(-0.5 * yp*yp / (sigma*sigma));
    double B = std::erf((xp+L/2) / (std::sqrt(2.0) * sigma));
    double C = std::erf((xp-L/2) / (std::sqrt(2.0) * sigma));
    return 2.0 * F * A * (B - C) / L;
}

std::shared_ptr<afw::image::Image<double>> ConvolvedTrailedSourceAlgorithm::computeModelImage(
    afw::table::SourceRecord& measRecord, afw::image::Exposure<float> const& exposure) const {
    // Get model parameters (ASSUMES NAIVE HAS BEEN COMPUTED)
    geom::Point2D center = _centroidExtractor(measRecord, _flagHandler);
    double xc = center.getX();
    double yc = center.getY();
    std::string CON = "ext_trailedSources_Convolved";
    double F = measRecord[_schema.find<double>(CON + "_sourceFlux").key];
    double x0 = measRecord[_schema.find<double>(CON + "_x0").key];
    double y0 = measRecord[_schema.find<double>(CON + "_y0").key];
    double x1 = measRecord[_schema.find<double>(CON + "_x1").key];
    double y1 = measRecord[_schema.find<double>(CON + "_y1").key];
    double x_m_x = x1 - x0;
    double y_m_y = y1 - y0;
    double L = std::sqrt(x_m_x*x_m_x + y_m_y*y_m_y); // Length of trail
    double theta = -std::atan2(y_m_y, x_m_x); // Angle measured from the image frame +x-axis (why flipped?)
    double sigma = exposure.getPsf()->computeShape(center).getTraceRadius(); // Psf sigma.

    // Generate blank image, size of footprint
    // Add templating later
    geom::Box2I sourceBB = measRecord.getFootprint()->getBBox();
    geom::Extent2I dims = sourceBB.getDimensions();
    std::shared_ptr<afw::image::Image<double>> im(new afw::image::Image<double>(sourceBB));
    afw::image::Image<double>::Array array = im->getArray();

    // Loop over pixels and evaluate computeModel()
    // Based on GuassianPsf::doComputeKernelImage()
    double sum = 0.0;
    for (int yIndex = 0, yp = im->getY0(); yIndex < dims.getY(); ++yIndex, ++yp) {
        afw::image::Image<double>::Array::Reference row = array[yIndex];
        for (int xIndex = 0, xp = im->getX0(); xIndex < dims.getX(); ++xIndex, ++xp) {
            sum += row[xIndex] = computeModel(xp,yp,xc,yc,F,L,theta,sigma);
        }
    }
    //ndarray::asEigenMatrix(array) /= sum;
    //ndarray::asEigenMatrix(array) *= F;
    return im;
}

void ConvolvedTrailedSourceAlgorithm::measure(afw::table::SourceRecord& measRecord,
                                              afw::image::Exposure<float> const& exposure) const {

    typedef afw::image::MaskedImage<float> MaskedImage;
    // Get model parameters (ASSUMES NAIVE HAS BEEN COMPUTED)
    geom::Point2D center = _centroidExtractor(measRecord, _flagHandler);
    double xc = center.getX();
    double yc = center.getY();
    std::string NAIVE = "ext_trailedSources_Naive";
    double F = measRecord[_schema.find<double>(NAIVE + "_flux").key];
    double x0 = measRecord[_schema.find<double>(NAIVE + "_x0").key];
    double y0 = measRecord[_schema.find<double>(NAIVE + "_y0").key];
    double x1 = measRecord[_schema.find<double>(NAIVE + "_x1").key];
    double y1 = measRecord[_schema.find<double>(NAIVE + "_y1").key];
    double x_m_x = x1 - x0;
    double y_m_y = y1 - y0;
    double L = std::sqrt(x_m_x*x_m_x + y_m_y*y_m_y); // Length of trail
    double theta = -std::atan2(y_m_y, x_m_x); // Angle measured from the image frame +x-axis (why flipped?)
    double sigma = exposure.getPsf()->computeShape(center).getTraceRadius(); // Psf sigma.

    // Get the exposures masked image and variance plane
    MaskedImage im = exposure.getMaskedImage();

    // Get data and errors
    const int height = im.getHeight();
    const int width = im.getHeight();
    std::vector<double> values(height * width);
    std::vector<double> errors(height * width);
    std::vector<double> x_position(height * width);
    std::vector<double> y_position(height * width);

    // From afw masked image examples
    for (int y = 0; y != height; ++y) {
        int x = 0;
        for (MaskedImage::x_iterator ptr = im.row_begin(y), end = im.row_end(y); ptr != end; ++ptr, ++x) {
            values[x + (width * y)] = ptr.image();
            errors[x + (width * y)] = ptr.variance();
            x_position[x + (width * y)] = x;
            y_position[x + (width * y)] = y;
        }
    }

    // Construct Function2
    std::vector<double> params = {xc, yc, F/(L * sigma * std::sqrt(2*3.1415)), L, theta};
    ConvolvedTrailedSourceFunction2<double> modelFunction(params, sigma);

    // Minimize
    // double ssize = _ctrl.stepsize; // Why doesn't this work??
    double ssize = 0.1;
    std::vector<double> stepsize(modelFunction.getNParameters(), ssize);
    //for (int i = 0; i < modelFunction.getNParameters(); ++i) {
    //    stepsize[i] = ssize * params[i];
    //}

    double sigma_min = 1.0; // I don't know what this does...

    afw::math::FitResults fitResults = afw::math::minimize(
        modelFunction, params, stepsize, values, errors, x_position, y_position, sigma_min);

    // Make an assert for a valid solution here

    // Save optimized parameters
    // modelFunction->setParameters(fitResults.parameterList);
    std::vector<double> fitParams(fitResults.parameterList);
    double xc_fit = fitParams[0];
    double yc_fit = fitParams[1];
    double F_fit = fitParams[2];
    double L_fit = fitParams[3];
    double theta_fit = fitParams[4];

    // Transform to endpoints
    double x0_fit = xc_fit - L_fit/2 * std::cos(theta_fit);
    double y0_fit = yc_fit + L_fit/2 * std::sin(theta_fit);
    double x1_fit = xc_fit + L_fit/2 * std::cos(theta_fit);
    double y1_fit = yc_fit - L_fit/2 * std::sin(theta_fit);

    // set keys
    measRecord.set(_xHeadKey, x0_fit);
    measRecord.set(_yHeadKey, y0_fit);
    measRecord.set(_xTailKey, x1_fit);
    measRecord.set(_yTailKey, y1_fit);
    measRecord.set(_totalFluxKey, F);
    measRecord.set(_sourceFluxKey, F_fit);
    measRecord.set(_chiSqKey, fitResults.chiSq / (values.size() - modelFunction.getNParameters() - 1));
    _flagHandler.setValue(measRecord, FAILURE.number, false);

}

void ConvolvedTrailedSourceAlgorithm::fail(afw::table::SourceRecord& measRecord,
                                           base::MeasurementError* error) const {
    _flagHandler.handleFailure(measRecord, error);
}

}}}} // lsst::meas::extensions::trailedSources
