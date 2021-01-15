#include "pybind11/pybind11.h"

#include "lsst/meas/extensions/trailedSources/NaiveTrailedSource.h"
#include "lsst/meas/extensions/trailedSources/ConvolvedTrailedSource.h"
#include "lsst/pex/config/python.h"

namespace py = pybind11;
using namespace pybind11::literals;

namespace lsst {
namespace meas {
namespace extensions {
namespace trailedSources {

PYBIND11_MODULE(trailedSources, mod) {
    py::module::import("lsst.meas.base");
    py::module::import("lsst.afw.table");
    py::module::import("lsst.afw.image");

    /* Module level */
    py::class_<NaiveTrailedSourceAlgorithm, std::shared_ptr<NaiveTrailedSourceAlgorithm>, base::SimpleAlgorithm>
            clsNaiveTrailedSourceAlgorithm(mod, "NaiveTrailedSourceAlgorithm");
    py::class_<NaiveTrailedSourceControl> clsNaiveTrailedSourceControl(mod, "NaiveTrailedSourceControl");

    py::class_<ConvolvedTrailedSourceAlgorithm, std::shared_ptr<ConvolvedTrailedSourceAlgorithm>, base::SimpleAlgorithm>
            clsConvolvedTrailedSourceAlgorithm(mod, "ConvolvedTrailedSourceAlgorithm");
    py::class_<ConvolvedTrailedSourceControl> clsConvolvedTrailedSourceControl(mod, "ConvolvedTrailedSourceControl");

    /* Constructors */
    clsNaiveTrailedSourceAlgorithm.def(py::init<NaiveTrailedSourceAlgorithm::Control const&,
                                  std::string const&, afw::table::Schema&>(),
                                  "ctrl"_a, "name"_a, "schema"_a);
    clsNaiveTrailedSourceControl.def(py::init<>());

    clsConvolvedTrailedSourceAlgorithm.def(py::init<ConvolvedTrailedSourceAlgorithm::Control const&,
                                  std::string const&, afw::table::Schema&>(),
                                  "ctrl"_a, "name"_a, "schema"_a);
    clsConvolvedTrailedSourceControl.def(py::init<>());

    /* Members */
    clsNaiveTrailedSourceAlgorithm.def("measure", &NaiveTrailedSourceAlgorithm::measure);
    clsNaiveTrailedSourceAlgorithm.def("fail", &NaiveTrailedSourceAlgorithm::fail);

    clsConvolvedTrailedSourceAlgorithm.def("measure", &ConvolvedTrailedSourceAlgorithm::measure);
    clsConvolvedTrailedSourceAlgorithm.def("computeModelImage",
                                        &ConvolvedTrailedSourceAlgorithm::computeModelImage);
    clsConvolvedTrailedSourceAlgorithm.def("fail", &ConvolvedTrailedSourceAlgorithm::fail);
}
}}}} // lsst::meas::extensions::trailedSources