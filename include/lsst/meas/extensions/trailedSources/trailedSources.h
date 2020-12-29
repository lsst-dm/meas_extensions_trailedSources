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

#ifndef LSST_MEAS_EXTENSIONS_TRAILEDSOURCES_H
#define LSST_MEAS_EXTENSIONS_TRAILEDSOURCES_H

#include "lsst/pex/config.h"
#include "lsst/meas/base/Algorithm.h"
#include "lsst/meas/base/CentroidUtilities.h"
#include "lsst/meas/base/FlagHandler.h"
#include "lsst/meas/base/InputUtilities.h"

namespace lsst {
namespace meas {
namespace extensions {
namespace trailedSources {

/**
 * Control class to handle TrailedSourceAlgorithm's configuration
 */
class TrailedSourceControl {
public:
    /**
     * Default constructor
     */
    TrailedSourceControl() : _name("") {}

    TrailedSourceControl(std::string const& name) : _name(name) {}

private:
    std::string _name;
};

class NaiveTrailedSourceControl : public TrailedSourceControl {
public:
    NaiveTrailedSourceControl() : TrailedSourceControl("ext_trailedSources_Naive") {}
};

class ConvolvedTrailedSourceControl : public TrailedSourceControl {
public:
    ConvolvedTrailedSourceControl() : TrailedSourceControl("ext_trailedSources_Convolved") {}
};

/*
 * Main Trailed-source measurement algorithm
 */
class TrailedSourceAlgorithm : public base::SimpleAlgorithm {
public:

    static base::FlagDefinitionList const& getFlagDefinitions();
    static base::FlagDefinition const FAILURE;

    typedef TrailedSourceControl Control;

    TrailedSourceAlgorithm(Control const& ctrl, std::string const& name,
                           std::string const& doc, afw::table::Schema& schema);

    virtual void measure(afw::table::SourceRecord& measRecord,
                         afw::image::Exposure<float> const& exposure) const {}

    virtual void fail(afw::table::SourceRecord& measRecord,
                      meas::base::MeasurementError* error = nullptr) const;

protected:
    Control _ctrl;
    std::string _doc;
    afw::table::Schema _schema;
    afw::table::Key<double> _xHeadKey;
    afw::table::Key<double> _yHeadKey;
    afw::table::Key<double> _xTailKey;
    afw::table::Key<double> _yTailKey;
    afw::table::Key<double> _fluxKey;
    base::FlagHandler _flagHandler;
    base::SafeCentroidExtractor _centroidExtractor;
};

class NaiveTrailedSourceAlgorithm : public TrailedSourceAlgorithm {
public:
    typedef NaiveTrailedSourceControl Control;
    NaiveTrailedSourceAlgorithm(Control const& ctrl, std::string const& name, afw::table::Schema& schema) :
        TrailedSourceAlgorithm(ctrl, name, "Naive trailed source", schema) {}

    void measure(afw::table::SourceRecord& measRecord,
                 afw::image::Exposure<float> const& exposure) const;
};

class ConvolvedTrailedSourceAlgorithm : public TrailedSourceAlgorithm {
public:
    typedef ConvolvedTrailedSourceControl Control;
    ConvolvedTrailedSourceAlgorithm(Control const& ctrl, std::string const& name, afw::table::Schema& schema):
        TrailedSourceAlgorithm(ctrl, name, "Convolved trailed source (Veres et al 2012)", schema) {}

    // Make this private??
    double _computeModel(double x, double y, double x0, double y0,
                         double L, double theta, double sigma) const;

    std::shared_ptr<afw::image::Image<double>> computeModelImage(
        afw::table::SourceRecord& measRecord, afw::image::Exposure<float> const& exposure) const;

    void measure(afw::table::SourceRecord& measRecord,
                 afw::image::Exposure<float> const& exposure) const;

};

}}}} // namespace lsst::meas::extensions::trailedSources

#endif // LSST_MEAS_EXTENSIONS_TRAILEDSOURCES_H
