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
        x0 = measRecord["ext_trailedSources_Naive_x0"]
        y0 = measRecord["ext_trailedSources_Naive_y0"]
        x1 = measRecord["ext_trailedSources_Naive_x1"]
        y1 = measRecord["ext_trailedSources_Naive_y1"]
        F = measRecord["ext_trailedSources_Naive_flux"]
        L = np.sqrt((x1-x0)**2 + (y1-y0)**2)
        theta = np.arctan2(y1-y0, x1-x0)

        # Make VeresModel
        params = np.array([xc, yc, F, L, theta])
        model = VeresModel(exposure, params)

        # Bounds
        bbox = exposure.getBBox()
        xmin = bbox.getY0()
        xmax = bbox.getX()
        ymin = bbox.getY0()
        ymax = bbox.getY()
        bounds = ((xmin, xmax), (ymin, ymax), (1.0, None),
                  (0.0, xmax), (-2*np.pi, 2*np.pi))
        # Do optimization with scipy
        results = minimize(model, params, bound=bounds,
                           method="Nelder-Mead")  # Allow user to choose alg

        # Calculate end points and reduced chi-squared
        xc_fit, yc_fit, L_fit, F_fit, theta_fit = results.x
        x0_fit = xc_fit - L_fit/2 * np.cos(theta_fit)
        y0_fit = yc_fit - L_fit/2 * np.sin(theta_fit)
        x1_fit = xc_fit + L_fit/2 * np.cos(theta_fit)
        y1_fit = yc_fit + L_fit/2 * np.sin(theta_fit)
        rChiSq = results.fun / (exposure.image.array.size - 6)

        # Set keys
        measRecord.set(self.keyX0, x0_fit)
        measRecord.set(self.keyY0, y0_fit)
        measRecord.set(self.keyX1, x1_fit)
        measRecord.set(self.keyY1, y1_fit)
        measRecord.set(self.keyFlux, F_fit)
        measRecord.set(self.keyRChiSq, rChiSq)

    def fail(self, measRecord, error=None):
        """Record failure
        """
        self.flagHandler.handleFailure(measRecord)
