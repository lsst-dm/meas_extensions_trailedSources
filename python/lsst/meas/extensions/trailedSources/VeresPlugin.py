#
# This file is part of meas_extensions_trailedSources.
#
# Developed for the LSST Data Management System.
# This product includes software developed by the LSST Project
# (http://www.lsst.org).
# See the COPYRIGHT file at the top-level directory of this distribution
# for details of code ownership.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
""" Veres trailed source measurement plugin.

Computes the length, angle from +x-axis, end points, flux, and centroid of an
extended source. Base on the model from Veres el al. 2012.
"""

import numpy as np
from scipy.optimize import minimize

from lsst.meas.base.pluginRegistry import register
from lsst.meas.base import SingleFramePlugin, SingleFramePluginConfig
from lsst.meas.base import FlagHandler, FlagDefinitionList, SafeCentroidExtractor

from ._trailedSources import VeresModel

__all__ = ("SingleFrameVeresTrailConfig", "SingleFrameVeresTrailPlugin")


class SingleFrameVeresTrailConfig(SingleFramePluginConfig):
    """Config class for SingleFrameVeresTrailPlugin
    """

    def setDefaults(self):
        super().setDefaults()


@register("ext_trailedSources_Veres")
class SingleFrameVeresTrailPlugin(SingleFramePlugin):
    """Veres trailed source characterization plugin.

    Parameters
    ----------
    config: `SingleFrameNaiveTrailConfig`
        Plugin configuration.
    name: `str`
        Plugin name.
    schema: `lsst.afw.table.Schema`
        Schema for the output catalog.
    metadata: `lsst.daf.base.PropertySet`
        Metadata to be attached to output catalog.

    See also
    --------
    lsst.meas.base.SingleFramePlugin
    """

    ConfigClass = SingleFrameVeresTrailConfig

    @classmethod
    def getExecutionOrder(cls):
        return cls.APCORR_ORDER + 0.2  # Needs to happen after Naive

    def __init__(self, config, name, schema, metadata):
        super().__init__(config, name, schema, metadata)

        self.keyX0 = schema.addField(name + "_x0", type="D", doc="Trail head X coordinate.", units="pixel")
        self.keyY0 = schema.addField(name + "_y0", type="D", doc="Trail head Y coordinate.", units="pixel")
        self.keyX1 = schema.addField(name + "_x1", type="D", doc="Trail tail X coordinate.", units="pixel")
        self.keyY1 = schema.addField(name + "_y1", type="D", doc="Trail tail Y coordinate.", units="pixel")
        self.keyL = schema.addField(name + "_length", type="D", doc="Length of trail.", units="pixel")
        self.keyTheta = schema.addField(name + "_angle", type="D", doc="Angle of trail from +x-axis.")
        self.keyFlux = schema.addField(name + "_flux", type="D", doc="Trailed source flux.", units="count")
        self.keyRChiSq = schema.addField(name + "_rChiSq", type="D", doc="Reduced chi-squared of fit")

        flagDefs = FlagDefinitionList()
        flagDefs.addFailureFlag("No trailed-sources measured")
        self.flagHandler = FlagHandler.addFields(schema, name, flagDefs)

        self.centroidExtractor = SafeCentroidExtractor(schema, name)

    def measure(self, measRecord, exposure):
        """Run the Veres trailed source measurement plugin.

        Parameters
        ----------
        measRecord : `lsst.afw.table.SourceRecord`
            Record describing the object being measured.
        exposure : `lsst.afw.image.Exposure`
            Pixel data to be measured.

        See also
        --------
        lsst.meas.base.SingleFramePlugin.measure
        """
        xc, yc = self.centroidExtractor(measRecord, self.flagHandler)
        # Look at measRecord for Naive end points ##
        # ASSUMES NAIVE ALREADY RAN
        # Should figure out an assert for this
        x0 = measRecord.get("ext_trailedSources_Naive_x0")
        y0 = measRecord.get("ext_trailedSources_Naive_y0")
        x1 = measRecord.get("ext_trailedSources_Naive_x1")
        y1 = measRecord.get("ext_trailedSources_Naive_y1")
        F = measRecord.get("ext_trailedSources_Naive_flux")
        L = np.sqrt((x1-x0)**2 + (y1-y0)**2)
        theta = np.arctan2(y1-y0, x1-x0)

        # Make VeresModel
        model = VeresModel(exposure)

        # Do optimization with scipy
        params = np.array([xc, yc, F, L, theta])
        results = minimize(model, params, method='Powell')  # Should allow user to choose alg

        # Calculate end points and reduced chi-squared
        xc_fit, yc_fit, F_fit, L_fit, theta_fit = results.x
        a = L_fit/2
        x0_fit = xc_fit - a * np.cos(theta_fit)
        y0_fit = yc_fit - a * np.sin(theta_fit)
        x1_fit = xc_fit + a * np.cos(theta_fit)
        y1_fit = yc_fit + a * np.sin(theta_fit)
        rChiSq = results.fun / (exposure.image.array.size - 6)

        # Set keys
        measRecord.set(self.keyX0, x0_fit)
        measRecord.set(self.keyY0, y0_fit)
        measRecord.set(self.keyX1, x1_fit)
        measRecord.set(self.keyY1, y1_fit)
        measRecord.set(self.keyFlux, F_fit)
        measRecord.set(self.keyL, L_fit)
        measRecord.set(self.keyTheta, theta_fit)
        measRecord.set(self.keyRChiSq, rChiSq)

    def fail(self, measRecord, error=None):
        """Record failure

        See also
        --------
        lsst.meas.base.SingleFramePlugin.fail
        """
        self.flagHandler.handleFailure(measRecord)
