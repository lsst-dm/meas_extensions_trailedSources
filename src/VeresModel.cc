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

VeresModel::VeresModel(
    Exposure const& data
) : _sigma(data.getPsf()->computeShape().getTraceRadius()),
    _bbox(data.getBBox()),
    _model(Image::Array(_bbox.getHeight(),_bbox.getWidth())),
    _data(data.getMaskedImage().getImage()->getArray()),
    _variance(data.getMaskedImage().getVariance()->getArray()),
    _gradModel({0.0,0.0,0.0,0.0,0.0}) {}

double VeresModel::operator()(std::vector<double> const& params) const {

    double xc = params[0];    // Centroid x
    double yc = params[1];    // Centroid y
    double F = params[2];     // Flux
    double L = params[3];     // Trail length
    double theta = params[4]; // Angle from +x-axis

    // From computeKernelImage()
    // Compute model image and chi-squared
    double chiSq = 0.0;
    double tmp = 0.0;
    for (int yIndex = 0, yp = _bbox.getBeginY(); yIndex < _bbox.getHeight(); ++yIndex, ++yp) {
        Image::Array::Reference modelRow = _model[yIndex];
        Image::Array::Reference dataRow = _data[yIndex];
        Image::Array::Reference varRow = _variance[yIndex];
        for (int xIndex = 0, xp = _bbox.getBeginX(); xIndex < _bbox.getWidth(); ++xIndex, ++xp) {
            modelRow[xIndex] = _computeModel(xp,yp,xc,yc,F,L,theta);
            tmp = dataRow[xIndex] - modelRow[xIndex];
            chiSq += tmp*tmp/varRow[xIndex];
        }
    }

    return chiSq;
}

std::vector<double> VeresModel::gradient(std::vector<double> const& params) {

    double xc = params[0];    // Centroid x
    double yc = params[1];    // Centroid y
    double F = params[2];     // Flux
    double L = params[3];     // Trail length
    double theta = params[4]; // Angle from +x-axis

    // Compute gradients of the model and of chi-squared
    std::vector<double> gradChiSq = {0.0,0.0,0.0,0.0,0.0};
    double tmp = 0.0;
    for (int yIndex = 0, yp = _bbox.getBeginY(); yIndex < _bbox.getHeight(); ++yIndex, ++yp) {
        Image::Array::Reference modelRow = _model[yIndex];
        Image::Array::Reference dataRow = _data[yIndex];
        Image::Array::Reference varRow = _variance[yIndex];
        for (int xIndex = 0, xp = _bbox.getBeginX(); xIndex < _bbox.getWidth(); ++xIndex, ++xp) {
            modelRow[xIndex] = _computeModel(xp,yp,xc,yc,F,L,theta);
            _computeGradient(xp,yp,xc,yc,F,L,theta);
            tmp = -2.0 * (dataRow[xIndex] - modelRow[xIndex]) / varRow[xIndex];
            for (int k=0; k<5; ++k) {
                gradChiSq[k] += _gradModel[k] * tmp;
            }
        }
    }
    return gradChiSq;
}

double VeresModel::_computeModel(double x, double y, double xc, double yc,
                                 double F, double L, double theta) const {
    double xp = (x-xc)*std::cos(theta) + (y-yc)*std::sin(theta);
    double yp = (x-xc)*std::sin(theta) - (y-yc)*std::cos(theta);
    double A = std::exp(-0.5 * yp*yp / (_sigma*_sigma));
    double B = std::erf((xp+L/2) / (std::sqrt(2.0) * _sigma));
    double C = std::erf((xp-L/2) / (std::sqrt(2.0) * _sigma));
    return F * A * (B - C) / (L * 2 * std::sqrt(2.0 * geom::PI) * _sigma);
}

void VeresModel::_computeGradient(double x, double y, double xc, double yc,
                                  double F, double L, double theta) {
    double xp = (x-xc)*std::cos(theta) + (y-yc)*std::sin(theta);
    double yp = (x-xc)*std::sin(theta) - (y-yc)*std::cos(theta);

    // Duplicated quantities
    double F2L = F/(2.0*L);
    double ypSq = yp*yp;
    double sqrt2 = std::sqrt(2.0);
    double sqrt2Pi = std::sqrt(2.0*geom::PI);
    double sigmaSq = _sigma*_sigma;
    double sigmaSq8 = sigmaSq * 8.0;
    double eypSq =  std::exp(-ypSq/(2.0*sigmaSq));
    double LPlus = L+2.0*xp;
    double LMinus= L-2.0*xp;
    double erfPlus = std::erf(LPlus/(2.0*sqrt2*_sigma));
    double erfMinus = std::erf(LMinus/(2.0*sqrt2*_sigma));
    double expPlus = std::exp(-LPlus*LPlus/sigmaSq8);

    // Compute partials wrt the transformed coordinates
    double dfdxp = F2L/(geom::PI*sigmaSq)*std::exp(-4.0*ypSq/sigmaSq8)*expPlus*
        (1.0 - std::exp(L*xp/sigmaSq));
    double dfdyp = -F2L*yp/(sqrt2Pi*_sigma*sigmaSq)*eypSq*(erfMinus+erfPlus);

    // Use the chain rule to get partials wrt the centroid and rotation angle
    double dxpdxc = -std::cos(theta);
    double dxpdyc = -std::sin(theta);
    double dxpdTheta = -yp;
    double dypdxc = -std::sin(theta);
    double dypdyc = std::cos(theta);
    double dypdTheta = xp;
    double dfdxc = dfdxp*dxpdxc + dfdyp*dypdxc;
    double dfdyc = dfdxp*dxpdyc + dfdyp*dypdyc;
    double dfdTheta = dfdxp*dxpdTheta + dfdyp*dypdTheta;

    double dfdF = _computeModel(x,y,xc,yc,1.0,L,theta); // dfdF = f / F

    double dfdL = F2L/(L*sqrt2Pi*_sigma)*eypSq*(L/(sqrt2Pi*_sigma)*
        (std::exp(-LMinus*LMinus/sigmaSq8)+expPlus) - erfMinus - erfPlus);

    _gradModel[0] = dfdxc;
    _gradModel[1] = dfdyc;
    _gradModel[2] = dfdF;
    _gradModel[3] = dfdL;
    _gradModel[4] = dfdTheta;
}

}}}} // lsst::meas::extensions::trailedSources