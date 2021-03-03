// -*- lsst-c++ -*-
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

/**
 * Implementation of an axisymmetric 2D Gaussian convolved with a line -- a model
 * for a fast-moving, trailed-source (Veres et al. 2012).
 *
 * VeresModel is designed to compute the chi-squared of the model given an
 * `lsst::afw::image::Exposure`, and to be passed to scipy.optimize for minimization.
 */
class VeresModel {
public:
    typedef afw::image::Image<float> Image;
    typedef afw::image::Exposure<float> Exposure;
    typedef afw::image::Image<float>::Array Array;

    /**
     * Constructor for VeresModel.
     *
     * @param data Exposure passed from the measurement task.
     */
    explicit VeresModel(Exposure const& data);

    /**
     * Compute the model and chi-squared given the data.
     *
     * @param params Model parameters. [centroid x, centroid y, trail flux,
     *     trail length, trail angle].
     *
     * @return The chi-squared of the model, given the data.
     */
    double operator()(std::vector<double> const& params) const;

    /// Return the current model image.
    std::shared_ptr<Image> getModelImage() const { return _image; }

    /// Return the PSF sigma.
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