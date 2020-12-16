#include "pybind11/pybind11.h"

#include "lsst/pex/config/python.h"
#include "lsst/meas/base/python.h"

#include "lsst/meas/extensions/trailedSources/trailedSources.h"

namespace py = pybind11;
using namespace pybind11::literals;

namespace lsst {
namespace meas {
namespace extensions {
namespace trailedSources {

PYBIND11_MODULE(TrailedSourceControl, mod) {
    py::module::import("lsst.meas.base");
    py::module::import("lsst.afw.table");
    
    /* Module level */
    py::class_<TrailedSourceAlgorithm, std::shared_ptr<TrailedSourceAlgorithm>, base::SimpleAlgorithm>
            clsTrailedSourceAlgorithm(mod, "TrailedSourceAlgorithm");
    py::class_<TrailedSourceControl> clsTrailedSourceControl(mod, "TrailedSourceControl");
    
    py::class_<NaiveTrailedSourceAlgorithm, std::shared_ptr<NaiveTrailedSourceAlgorithm, base::SimpleAlgorithm>>
            clsNaiveTrailedSourceAlgorithm(mod, "NaiveTrailedSourceAlgorithm");
    py::class_<NaiveTrailedSourceControl, TrailedSourceControl> 
            clsNaiveTrailedSourceControl(mod, "NaiveTrailedSourceControl");
    
    /* Constructors */
    clsNaiveTrailedSourceAlgorithm.def(
            py::init<NaiveTrailedSourceAlgorithm::Control const&, std::string const&, afw::table::Schema&>(),
            "ctrl"_a, "name"_a, "schema"_a);
    clsNaiveTrailedSourceControl.def(py::init<>());

    /* Members */
    clsTrailedSourceAlgorithm.def("measure", &TrailedSourceAlgorithm::measure);
    clsTrailedSourceAlgorithm.def("fail", &TrailedSourceAlgorithm::fail);
}
}}}} // lsst::meas::extensions::trailedSources