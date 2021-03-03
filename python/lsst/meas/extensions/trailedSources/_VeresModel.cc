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

#include "pybind11/pybind11.h"
#include "lsst/utils/python.h"
#include "pybind11/stl.h"

#include "lsst/meas/extensions/trailedSources/VeresModel.h"
#include "lsst/geom.h"
#include "lsst/afw/image.h"
#include "lsst/afw/detection.h"

namespace py = pybind11;
using namespace pybind11::literals;

namespace lsst {
namespace meas {
namespace extensions {
namespace trailedSources {

void wrapVeresModel(utils::python::WrapperCollection& wrappers) {

    wrappers.addSignatureDependency("lsst.geom");
    wrappers.addSignatureDependency("lsst.afw.image");
    wrappers.addSignatureDependency("lsst.afw.detection");

    wrappers.wrapType(
        py::class_<VeresModel, std::shared_ptr<VeresModel>>(wrappers.module, "VeresModel"),
        [](auto & mod, auto & cls) {
            // cls.def(py::init<lsst::geom::Box2I const&, double>(), "bbox"_a, "sigma"_a);
            cls.def(py::init<afw::image::Exposure<float> const&>(), "data"_a);
            cls.def("__call__",
                py::overload_cast<std::vector<double> const&>(&VeresModel::operator(), py::const_));
            cls.def("getSigma", &VeresModel::getSigma);
            cls.def("getModelImage", &VeresModel::getModelImage);
    });
}
}}}} // lsst::meas::extensions::trailedSources