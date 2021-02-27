// _VeresModel.cc

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