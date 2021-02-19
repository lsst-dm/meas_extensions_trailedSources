#include "pybind11/pybind11.h"

#include "lsst/utils/python.h"

namespace lsst {
namespace meas {
namespace extensions {
namespace trailedSources {

void wrapVeresModel(lsst::utils::python::WrapperCollection & wrappers);

PYBIND11_MODULE(_trailedSources, mod) {
    lsst::utils::python::WrapperCollection wrappers(mod, "lsst.meas.extensions.trailedSources");
    wrapVeresModel(wrappers);
    wrappers.finish();
}

}}}} // lsst::meas::extensions::trailedSources