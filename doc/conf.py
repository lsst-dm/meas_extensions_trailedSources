"""Sphinx configuration file for an LSST stack package.

This configuration only affects single-package Sphinx documentation builds.
"""

from documenteer.sphinxconfig.stackconf import build_package_configs
import lsst.meas.extensions.trailedSources


_g = globals()
_g.update(build_package_configs(
    project_name='meas_extensions_trailedSources',
    version=lsst.meas.extensions.trailedSources.version.__version__))
