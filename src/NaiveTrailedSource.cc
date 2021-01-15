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

#include "lsst/meas/extensions/trailedSources/NaiveTrailedSource.h"

namespace lsst {
namespace meas {
namespace extensions {
namespace trailedSources {

namespace {
lsst::meas::base::FlagDefinitionList flagDefinitions;
}

base::FlagDefinition const NaiveTrailedSourceAlgorithm::FAILURE = flagDefinitions.addFailureFlag();

base::FlagDefinitionList const& NaiveTrailedSourceAlgorithm::getFlagDefinitions() { return flagDefinitions; }

NaiveTrailedSourceAlgorithm::NaiveTrailedSourceAlgorithm(
    Control const& ctrl,
    std::string const& name,
    afw::table::Schema& schema
) : _ctrl(ctrl),
    _doc("Naive trailed source"),
    _schema(schema),
    _xHeadKey(schema.addField<double>(name + "_x0", _doc)),
    _yHeadKey(schema.addField<double>(name + "_y0", _doc)),
    _xTailKey(schema.addField<double>(name + "_x1", _doc)),
    _yTailKey(schema.addField<double>(name + "_y1", _doc)),
    _fluxKey(schema.addField<double>(name + "_flux", _doc)),
    _flagHandler(base::FlagHandler::addFields(schema, name, getFlagDefinitions())),
    _centroidExtractor(schema, name) {}

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

void NaiveTrailedSourceAlgorithm::fail(afw::table::SourceRecord& measRecord,
                                       base::MeasurementError* error) const {
    _flagHandler.handleFailure(measRecord, error);
}

}}}} // lsst::meas::extensions::trailedSources
