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

#include "lsst/meas/extensions/trailedSources/trailedSources.h"

namespace lsst {
namespace meas {
namespace extensions {
namespace trailedSources {

namespace {
lsst::meas::base::FlagDefinitionList flagDefinitions;
}

base::FlagDefinition const TrailedSourceAlgorithm::FAILURE = flagDefinitions.addFailureFlag();

base::FlagDefinitionList const& TrailedSourceAlgorithm::getFlagDefinitions() { return flagDefinitions; }

TrailedSourceAlgorithm::TrailedSourceAlgorithm(
    Control const& ctrl,
    std::string const& name,
    std::string const& doc,
    afw::table::Schema& schema
) : _ctrl(ctrl),
    _doc(doc),
    _schema(schema),
    _xHeadKey(schema.addField<double>(name + "_x0", _doc)),
    _yHeadKey(schema.addField<double>(name + "_y0", _doc)),
    _xTailKey(schema.addField<double>(name + "_x1", _doc)),
    _yTailKey(schema.addField<double>(name + "_y1", _doc)),
    _fluxKey(schema.addField<double>(name + "_flux", _doc)),
    _flagHandler(base::FlagHandler::addFields(schema, name, getFlagDefinitions())),
    _centroidExtractor(schema, name) {}

void TrailedSourceAlgorithm::fail(afw::table::SourceRecord& measRecord,
                                  base::MeasurementError* error) const {
    _flagHandler.handleFailure(measRecord, error);
}

void NaiveTrailedSourceAlgorithm::measure(afw::table::SourceRecord& measRecord,
                                          afw::image::Exposure<float> const& exposure) const {
    // Make error flag for if no psf [ ]
    // get centroid
    geom::Point2D center = _centroidExtractor(measRecord, _flagHandler);
    double x = center.getX();
    double y = center.getY();
    // get quadrupole
    // afw::geom::ellipses::Quadrupole quad = exposure.getPsf()->computeShape(center);
    double Ixx = measRecord.getIxx();
    double Iyy = measRecord.getIyy();
    double Ixy = measRecord.getIxy();
    // calculate ellipse parameters (from afw::geom::ellipses::BaseCore)
    double xx_p_yy = Ixx + Iyy;
    double xx_m_yy = Ixx - Iyy;
    double t = std::sqrt(xx_m_yy * xx_m_yy + 4 * Ixy * Ixy);
    double a = std::sqrt(0.5 * (xx_p_yy + t));
    double theta = 0.5 * std::atan2(2.0 * Ixy, xx_m_yy);
    // Calculate end points
    double x0 = -a*std::cos(theta);
    double y0 = -a*std::sin(theta);
    double x1 = a*std::cos(theta);
    double y1 = a*std::sin(theta);
    // calculate flux
    double F = measRecord.getApInstFlux(); // Change this later
    // set keys
    measRecord.set(_xHeadKey, x + x0);
    measRecord.set(_yHeadKey, y + y0);
    measRecord.set(_xTailKey, x + x1);
    measRecord.set(_yTailKey, y + y1);
    measRecord.set(_fluxKey, F);
    _flagHandler.setValue(measRecord, FAILURE.number, false);
}

double ConvolvedTrailedSourceAlgorithm::_computeModel(double x, double y, double x0, double y0,
                                                      double L, double theta, double sigma) const {
    // Computes the Veres et al model at a given position (pixel)
    double xp = (x-x0)*std::cos(theta) - (y-y0)*std::sin(theta);
    double yp = (x-x0)*std::sin(theta) + (y-y0)*std::cos(theta);
    double A = std::exp(-0.5 * yp*yp / (sigma*sigma));
    double B = std::erf((xp+L/2) / (std::sqrt(2.0) * sigma));
    double C = std::erf((xp-L/2) / (std::sqrt(2.0) * sigma));
    return 2.0 * A * (B - C) / L;
}

std::shared_ptr<afw::image::Image<double>> ConvolvedTrailedSourceAlgorithm::computeModelImage(
    afw::table::SourceRecord& measRecord, afw::image::Exposure<float> const& exposure) const {
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
    double theta = -std::atan2(y_m_y, x_m_x); // Angle measured from the image frame +x-axis
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
            sum += row[xIndex] = _computeModel(xp,yp,xc,yc,L,theta,sigma);
        }
    }
    ndarray::asEigenMatrix(array) /= sum;
    ndarray::asEigenMatrix(array) *= F;
    return im;
}

void ConvolvedTrailedSourceAlgorithm::measure(afw::table::SourceRecord& measRecord,
                                              afw::image::Exposure<float> const& exposure) const {
    /*/ set keys
    measRecord.set(_xHeadKey, x + x0);
    measRecord.set(_yHeadKey, y + y0);
    measRecord.set(_xTailKey, x + x1);
    measRecord.set(_yTailKey, y + y1);
    measRecord.set(_fluxKey, F);
    */
    _flagHandler.setValue(measRecord, FAILURE.number, false);

}

}}}} // lsst::meas::extensions::trailedSources
