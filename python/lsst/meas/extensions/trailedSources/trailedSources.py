# Trailed source detection plugins
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
        return cls.APCORR_ORDER

    def __init__(self, config, name, schema, metadata):
        super().__init__(config, name, schema, metadata)
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
        res = self.optimize(measRecord, exposure)
        F, x_h, y_h, x_t, y_t = res.x
        x0, y0 = self.centroidExtractor(measRecord, self.flagHandler)
        reduced_chi2 = res.fun / (exposure.image.array.size - 5 - 1)
        measRecord.set(self.keyHeadX, x0 + x_h)
        measRecord.set(self.keyHeadY, y0 + y_h)
        measRecord.set(self.keyTailX, x0 + x_t)
        measRecord.set(self.keyTailY, y0 + y_t)
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
        if theta < 0.0:
            theta += np.pi
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

    def psf(self, x0, y0, xdim, ydim, sigma):
        """
        Point-spread function convolution kernel.
        Double guassian, centered at x0,y0.
        """
        x = np.arange(-xdim//2 + 1, xdim//2 + 1) - x0
        y = np.arange(-ydim//2 + 1, ydim//2 + 1) - y0

        r = np.sqrt((x[None, :])**2 + (y[:, None])**2)
        psf = np.exp(-r**2./2./sigma**2.) / (2*np.pi*sigma**2.)

        return psf

    def makeTrailedSourceImage(self, params, xdim, ydim, npts=10, sigma=1.5):
        """
        Generate trailed source image.
        Integrate a psf over path.
        """
        F, xt, yt, xh, yh = params
        image = np.zeros((ydim, xdim))

        xs = np.linspace(xt, xh, npts)
        ys = np.linspace(yt, yh, npts)

        for i in range(xs.size):
            image += self.psf(ys[i], xs[i], xdim, ydim, sigma=sigma)

        image = F * (image / image.sum())
        return image

    def chiSquared(self, parameters, data, error, sigma, npts):
        xdim, ydim = data.shape
        model = self.makeTrailedSourceImage(parameters, xdim, ydim, sigma=sigma, npts=npts)
        return np.sum(((data - model)**2.0 / error))

    def optimize(self, measRecord, exposure):
        F0 = measRecord.getApInstFlux()
        x0, y0, x1, y1 = self.getEndPoints(measRecord)
        p0 = np.array([F0, x0, y0, x1, y1])

        data = exposure.image.array
        error = exposure.getVariance().array
        sigma = self.config.psfSigma
        npts = self.config.modelNumPoints
        method = self.config.optimizerMethod

        res = minimize(self.chiSquared, p0, args=(data, error, sigma, npts),
                       method=method, options={"maxiter": self.config.optimizerMaxIter})
        return res
