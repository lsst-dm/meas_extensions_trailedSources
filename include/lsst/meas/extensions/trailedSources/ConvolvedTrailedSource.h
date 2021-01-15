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

#ifndef LSST_MEAS_EXTENSIONS_TRAILEDSOURCES_CONVOLVED_H
#define LSST_MEAS_EXTENSIONS_TRAILEDSOURCES_CONVOVLED_H

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
class ConvolvedTrailedSourceControl {
public:
    /**
     * Default constructor
     */
    explicit ConvolvedTrailedSourceControl() : _name("ext_trailedSources_Convolved") {}

private:
    std::string _name;

};

/*
 * Main Trailed-source measurement algorithm
 */
class ConvolvedTrailedSourceAlgorithm : public base::SimpleAlgorithm {
public:

    static base::FlagDefinitionList const& getFlagDefinitions();
    static base::FlagDefinition const FAILURE;

    typedef ConvolvedTrailedSourceControl Control;

    explicit ConvolvedTrailedSourceAlgorithm(Control const& ctrl, std::string const& name,
                                            afw::table::Schema& schema);

    double computeModel(double x, double y, double x0, double y0, double F,
                        double L, double theta, double sigma) const;

    std::shared_ptr<afw::image::Image<double>> computeModelImage(
        afw::table::SourceRecord& measRecord, afw::image::Exposure<float> const& exposure) const;

    void measure(afw::table::SourceRecord& measRecord,
                 afw::image::Exposure<float> const& exposure) const;

    void fail(afw::table::SourceRecord& measRecord,
              meas::base::MeasurementError* error = nullptr) const;

private:
    Control _ctrl;
    std::string _doc;
    afw::table::Schema _schema;
    afw::table::Key<double> _xHeadKey;
    afw::table::Key<double> _yHeadKey;
    afw::table::Key<double> _xTailKey;
    afw::table::Key<double> _yTailKey;
    afw::table::Key<double> _totalFluxKey;
    afw::table::Key<double> _sourceFluxKey;
    afw::table::Key<double> _chiSqKey;
    base::FlagHandler _flagHandler;
    base::SafeCentroidExtractor _centroidExtractor;
};

template <typename ReturnT>
class ConvolvedTrailedSourceFunction2 : public afw::math::Function2<ReturnT> {
public:
    explicit ConvolvedTrailedSourceFunction2(std::vector<double> const& params, double sigma)
    : afw::math::Function2<ReturnT>(params), _sigma(sigma) {};

    // Needed for some reason...
    std::shared_ptr<afw::math::Function2<ReturnT>> clone() const override {
        return std::shared_ptr<afw::math::Function2<ReturnT>>(
            new ConvolvedTrailedSourceFunction2(this->_params, _sigma)
        );
    }

    ReturnT operator()(double x, double y) const noexcept(afw::math::IS_NOTHROW_INIT<ReturnT>) override {
        return _computeModel(this->_params, _sigma, x, y);
    }

private:
    double _computeModel(std::vector<double> const& params, double sigma, double x, double y) const;

    double _sigma;
};

}}}} // namespace lsst::meas::extensions::trailedSources

#endif // LSST_MEAS_EXTENSIONS_TRAILEDSOURCES_H
