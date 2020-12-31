# Trailed source plugins
import numpy as np
from scipy.optimize import minimize

import lsst.meas.base
from lsst.meas.base import SingleFramePlugin, SingleFramePluginConfig
from lsst.meas.base.pluginRegistry import register
from lsst.pex.config import Field

__all__ = ("SingleFrameTrailedSourceConfig", "SingleFrameTrailedSourcePlugin")


class SingleFrameTrailedSourceConfig(SingleFramePluginConfig):
    """Config class for SingleFrameTrailedSourcePlugin
    """
    optimizerMethod = Field(
        doc="Optimization method",
        dtype=str,
        default="Nelder-Mead",
    )
    optimizerMaxIter = Field(
        doc="Maximum number of optimization iterations",
        dtype=int,
        default=1000,
    )
    psfSigma = Field(
        doc="PSF sigma for trailed-source model",
        dtype=float,
        default=1.5,
    )
    modelNumPoints = Field(
        doc="Number of sumation steps in the model",
        dtype=int,
        default=10,
    )

    def setDefaults(self):
        super().setDefaults()


@register("ext_trailedSources_TrailedSource")
class SingleFrameTrailedSourcePlugin(SingleFramePlugin):
    """
    Measurement plugin to calculate the endpoints
    and total flux of a trailed source.
    """

    ConfigClass = SingleFrameTrailedSourceConfig

    @classmethod
    def getExecutionOrder(cls):
        return cls.APCORR_ORDER + 0.1

    def __init__(self, config, name, schema, metadata):
        super().__init__(config, name, schema, metadata)
        self.keyNHeadX = schema.addField(name + "_head_x_naive", type="D", doc="X coordinate of trail head.",
                                         units="pixel")
        self.keyNHeadY = schema.addField(name + "_head_y_naive", type="D", doc="Y coordinate of trail head.",
                                         units="pixel")
        self.keyNTailX = schema.addField(name + "_tail_x_naive", type="D", doc="X coordinate of trail tail.",
                                         units="pixel")
        self.keyNTailY = schema.addField(name + "_tail_y_naive", type="D", doc="Y coordinate of trail tail.",
                                         units="pixel")
        self.keyNFlux = schema.addField(name + "_flux_naive", type="D", doc="Flux of trailed source (n)",
                                        units="count")
        self.keyHeadX = schema.addField(name + "_head_x", type="D", doc="X coordinate of trail head",
                                        units="pixel")
        self.keyHeadY = schema.addField(name + "_head_y", type="D", doc="Y coordinate of trail head",
                                        units="pixel")
        self.keyTailX = schema.addField(name + "_tail_x", type="D", doc="X coordinate of trail tail",
                                        units="pixel")
        self.keyTailY = schema.addField(name + "_tail_y", type="D", doc="Y coordinate of trail tail",
                                        units="pixel")
        self.keyFlux = schema.addField(name + "_flux", type="D", doc="Flux of trailed source",
                                       units="count")
        self.keyRChi2 = schema.addField(name + "_reduced_chi_squared", type="D",
                                        doc="Reduced chi-squared of best fit")

        flagDefs = lsst.meas.base.FlagDefinitionList()
        flagDefs.addFailureFlag("No trailed-sources measured")
        self.flagHandler = lsst.meas.base.FlagHandler.addFields(schema, name, flagDefs)

        self.centroidExtractor = lsst.meas.base.SafeCentroidExtractor(schema, name)

    def measure(self, measRecord, exposure):
        """
        Measure trailed source end points and flux.
        Also record chi-squared and reduced chi-squared.
        """
        x0, y0 = self.centroidExtractor(measRecord, self.flagHandler)
        x_h, y_h, x_t, y_t = self.getEndPoints(measRecord)
        F = measRecord.getApInstFlux()
        measRecord.set(self.keyNHeadX, x0 + x_h)
        measRecord.set(self.keyNHeadY, y0 + y_h)
        measRecord.set(self.keyNTailX, x0 + x_t)
        measRecord.set(self.keyNTailY, y0 + y_t)
        measRecord.set(self.keyNFlux, F)

        params = np.array([F, x0 + x_h, y0 + y_h, x0 + x_t, y0 + y_t])
        res = self.optimize(measRecord, exposure, params)
        F, x_h, y_h, x_t, y_t = res.x
        reduced_chi2 = res.fun / (exposure.image.array.size - 5 - 1)
        measRecord.set(self.keyHeadX, x_h)
        measRecord.set(self.keyHeadY, y_h)
        measRecord.set(self.keyTailX, x_t)
        measRecord.set(self.keyTailY, y_t)
        measRecord.set(self.keyFlux, F)
        measRecord.set(self.keyRChi2, reduced_chi2)

    def fail(self, measRecord, error=None):
        """Record failure
        """
        self.flagHandler.handleFailure(measRecord)

    def getEllipseParams(self, measRecord):
        """
        Calculate semi-major axis and angle from horizonal of ellipse
        using quadrupole moments.
        """
        Ixx, Iyy, Ixy = measRecord.getShape().getParameterVector()
        a = np.sqrt(0.5 * (Ixx + Iyy + np.sqrt((Ixx - Iyy)**2.0 + 4.0*(Ixy**2.0))))
        theta = 0.5 * np.arctan2(2.0 * Ixy, Ixx - Iyy)
        # if theta < 0.0: theta += np.pi
        return a, theta

    def getEndPoints(self, measRecord):
        """Calculate trailed source end points.
        """
        a, theta = self.getEllipseParams(measRecord)
        x_h = -a*np.cos(theta)
        y_h = -a*np.sin(theta)
        x_t = a*np.cos(theta)
        y_t = a*np.sin(theta)
        return x_h, y_h, x_t, y_t

    def psf(self, exposure, x0, y0, sigma):
        """
        Point-spread function convolution kernel.
        Double guassian, centered at x0,y0.
        """
        bb = exposure.getBBox()
        x = np.arange(bb.beginX, bb.endX) - x0
        y = np.arange(bb.beginY, bb.endY) - y0

        r = np.sqrt((x[:, None])**2 + (y[None, :])**2)
        psf = np.exp(-r**2./2./sigma**2.) / (2*np.pi*sigma**2.)

        return psf

    def makeTrailedSourceImage(self, exposure, params, npts=10, sigma=1.5):
        """
        Generate trailed source image.
        Integrate a psf over path.
        """
        F, xh, yh, xt, yt = params
        image = np.zeros(exposure.getDimensions())

        xs = np.linspace(xh, xt, npts)
        ys = np.linspace(yh, yt, npts)

        for i in range(xs.size):
            image += self.psf(exposure, xs[i], ys[i], sigma)

        image = F * (image / image.sum())
        return image

    def chiSquared(self, parameters, exposure, sigma, npts):
        data = exposure.image.array
        error = exposure.getVariance().array
        model = self.makeTrailedSourceImage(exposure, parameters, sigma=sigma, npts=npts)
        return np.sum(((data - model)**2.0 / error))

    def optimize(self, measRecord, exposure, params):
        sigma = self.config.psfSigma
        npts = self.config.modelNumPoints
        method = self.config.optimizerMethod

        res = minimize(self.chiSquared, params, args=(exposure, sigma, npts),
                       method=method, options={"maxiter": self.config.optimizerMaxIter})
        return res