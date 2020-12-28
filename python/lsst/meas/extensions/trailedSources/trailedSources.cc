#include "pybind11/pybind11.h"

#include "lsst/meas/extensions/trailedSources/trailedSources.h"
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

    /* Module level */
    py::class_<TrailedSourceAlgorithm, std::shared_ptr<TrailedSourceAlgorithm>, base::SimpleAlgorithm>
            clsTrailedSourceAlgorithm(mod, "TrailedSourceAlgorithm");
    py::class_<TrailedSourceControl> clsTrailedSourceControl(mod, "TrailedSourceControl");

    py::class_<NaiveTrailedSourceAlgorithm, std::shared_ptr<NaiveTrailedSourceAlgorithm>, TrailedSourceAlgorithm>
            clsNaiveTrailedSourceAlgorithm(mod, "NaiveTrailedSourceAlgorithm");
    py::class_<NaiveTrailedSourceControl, TrailedSourceControl>
            clsNaiveTrailedSourceControl(mod, "NaiveTrailedSourceControl");

    /* Constructors */
    clsNaiveTrailedSourceAlgorithm.def(py::init<NaiveTrailedSourceAlgorithm::Control const&,
                                  std::string const&, afw::table::Schema&>(),
                                  "ctrl"_a, "name"_a, "schema"_a);
    clsNaiveTrailedSourceControl.def(py::init<>());

    /* Members */
    clsTrailedSourceAlgorithm.def("fail", &TrailedSourceAlgorithm::fail);
    clsNaiveTrailedSourceAlgorithm.def("measure", &NaiveTrailedSourceAlgorithm::measure);
}
}}}} // lsst::meas::extensions::trailedSources