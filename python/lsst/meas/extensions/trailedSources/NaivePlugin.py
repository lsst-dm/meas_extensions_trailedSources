import numpy as np

from lsst.meas.base.pluginRegistry import register
from lsst.meas.base import SingleFramePlugin, SingleFramePluginConfig
from lsst.meas.base import FlagHandler, FlagDefinitionList, SafeCentroidExtractor

__all__ = ("SingleFrameNaiveTrailConfig", "SingleFrameNaiveTrailPlugin")


class SingleFrameNaiveTrailConfig(SingleFramePluginConfig):
    """Config class for SingleFrameNaiveTrailPlugin
    """

    def setDefaults(self):
        super().setDefaults()


@register("ext_trailedSources_Naive")
class SingleFrameNaiveTrailPlugin(SingleFramePlugin):
    """Naive trailed source characterization plugin
    """

    ConfigClass = SingleFrameNaiveTrailConfig

    @classmethod
    def getExecutionOrder(cls):
        return cls.APCORR_ORDER + 0.1

    def __init__(self, config, name, schema, metadata):
        super().__init__(config, name, schema, metadata)

        # Measurement Keys
        self.keyX0 = schema.addField(name + "_x0", type="D", doc="Trail head X coordinate.", units="pixel")
        self.keyY0 = schema.addField(name + "_y0", type="D", doc="Trail head Y coordinate.", units="pixel")
        self.keyX1 = schema.addField(name + "_x1", type="D", doc="Trail tail X coordinate.", units="pixel")
        self.keyY1 = schema.addField(name + "_y1", type="D", doc="Trail tail Y coordinate.", units="pixel")
        self.keyFlux = schema.addField(name + "_flux", type="D", doc="Trailed source flux.", units="count")
        self.keyL = schema.addField(name + "_length", type="D", doc="Trail length.", units="pixel")
        self.keyAngle = schema.addField(name + "_angle", type="D", doc="Angle measured from +x-axis.")

        # Measurement Error Keys
        self.keyX0Err = schema.addField(name + "_x0Err", type="D",
                                        doc="Trail head X coordinate error.", units="pixel")
        self.keyY0Err = schema.addField(name + "_y0Err", type="D",
                                        doc="Trail head Y coordinate error.", units="pixel")
        self.keyX1Err = schema.addField(name + "_x1Err", type="D",
                                        doc="Trail tail X coordinate error.", units="pixel")
        self.keyY1Err = schema.addField(name + "_y1Err", type="D",
                                        doc="Trail tail Y coordinate error.", units="pixel")
        # self.keyFluxErr = schema.addField(name + "_fluxErr", type="D",
        #                      doc="Trailed source flux error.", units="count")

        flagDefs = FlagDefinitionList()
        flagDefs.addFailureFlag("No trailed-source measured")
        self.flagHandler = FlagHandler.addFields(schema, name, flagDefs)

        self.centriodExtractor = SafeCentroidExtractor(schema, name)

    def measure(self, measRecord, exposure):
        """
        Measure trailed source end points and flux
        """
        xc, yc = self.centriodExtractor(measRecord, self.flagHandler)
        Ixx, Iyy, Ixy = measRecord.getShape().getParameterVector()
        xmy = Ixx - Iyy
        xpy = Ixx + Iyy
        xmy2 = xmy*xmy
        # xpy2 = xpy*xpy
        xy2 = Ixy*Ixy
        a = np.sqrt(0.5 * (xpy + np.sqrt((xmy)**2 + 4.0*xy2)))
        theta = 0.5 * np.arctan2(2.0 * Ixy, xmy)

        x0 = xc-a*np.cos(theta)
        y0 = yc-a*np.sin(theta)
        x1 = xc+a*np.cos(theta)
        y1 = yc+a*np.sin(theta)
        F = measRecord.get("base_SdssShape_instFlux")  # Should this even be done?

        # Calculate errors (Clean this up later...)
        xcErr2, ycErr2 = np.diag(measRecord.getCentroidErr())  # Is there a "safe" error extractor?
        IxxErr2, IyyErr2, IxyErr2 = np.diag(measRecord.getShapeErr())
        desc = np.sqrt(xmy2 + 4.0*xy2)  # Descriminant^1/2 of EV equation
        denom = 2*np.sqrt(2.0*(Ixx + np.sqrt(4.0*xy2 + xmy2 + Iyy)))  # Denominator for dadIxx and dadIyy
        dadIxx = (1.0 + (xmy/desc)) / denom
        dadIyy = (1.0 - (xmy/desc)) / denom
        dadIxy = (4.0*Ixy) / (desc * denom)
        aErr2 = IxxErr2*dadIxx*dadIxx + IyyErr2*dadIyy*dadIyy + IxyErr2*dadIxy*dadIxy
        thetaErr2 = ((IxxErr2 + IyyErr2)*xy2 + xmy2*IxyErr2) / (desc*desc*desc*desc)

        dxda = np.cos(theta)
        dxdt = a*np.sin(theta)
        dyda = np.sin(theta)
        dydt = a*np.cos(theta)
        xErr2 = aErr2*dxda*dxda + thetaErr2*dxdt*dxdt
        yErr2 = aErr2*dyda*dyda + thetaErr2*dydt*dydt
        x0Err = np.sqrt(xErr2 + xcErr2)  # Same for x1
        y0Err = np.sqrt(yErr2 + ycErr2)  # Same for y1

        # Set flags
        measRecord.set(self.keyX0, x0)
        measRecord.set(self.keyY0, y0)
        measRecord.set(self.keyX1, x1)
        measRecord.set(self.keyY1, y1)
        measRecord.set(self.keyFlux, F)
        measRecord.set(self.keyL, 2*a)
        measRecord.set(self.keyAngle, theta)
        measRecord.set(self.keyX0Err, x0Err)
        measRecord.set(self.keyY0Err, y0Err)
        measRecord.set(self.keyX1Err, x0Err)
        measRecord.set(self.keyY1Err, y0Err)
        # measRecord.set(self.keyFluxErr, fluxErr)

    def fail(self, measRecord, error=None):
        """Record failure
        """
        self.flagHandler.handleFailure(measRecord)
