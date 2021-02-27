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
    """

    Trailed source characterization plugin using the Veres et al. 2012 model.
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
        """
        Measure trailed source end points and flux,
        and record reduced chi-squared for fit.
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
        """
        self.flagHandler.handleFailure(measRecord)
