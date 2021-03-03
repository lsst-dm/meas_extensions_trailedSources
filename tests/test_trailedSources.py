import numpy as np
import unittest
import lsst.utils.tests
import lsst.meas.extensions.trailedSources
from lsst.meas.base.tests import AlgorithmTestCase
from lsst.utils.tests import methodParameters

# Trailed-source length, angle, and centroid.
np.random.seed(432)
nTrails = 50
Ls = np.random.uniform(2, 20, nTrails)
thetas = np.random.uniform(0, 2*np.pi, nTrails)
xcs = np.random.uniform(0, 100, nTrails)
ycs = np.random.uniform(0, 100, nTrails)


class TrailedSource:

    def __init__(self, instFlux, length, angle, xc, yc):
        self.instFlux = instFlux
        self.length = length
        self.angle = angle
        self.center = lsst.geom.Point2D(xc, yc)
        self.x0 = xc - length/2 * np.cos(angle)
        self.y0 = yc - length/2 * np.sin(angle)
        self.x1 = xc + length/2 * np.cos(angle)
        self.y1 = yc + length/2 * np.sin(angle)


# "Extend" meas.base.tests.TestDataset
class TrailedTestDataset(lsst.meas.base.tests.TestDataset):

    def __init__(self, bbox, threshold=10.0, exposure=None, **kwds):
        super().__init__(bbox, threshold, exposure, **kwds)

    def addToSchema(self, schema):
        # Add Naive keys to schema
        name = "ext_trailedSources_Naive"
        self.keys["flux"] = schema.addField(
            name + "_flux", type="D", doc="Trailed source flux.", units="count")
        self.keys["L"] = schema.addField(
            name + "_length", type="D", doc="Trail length.", units="pixel")
        self.keys["theta"] = schema.addField(
            name + "_angle", type="D", doc="Trail angle from +x-axis.")

        return schema

    def addToCatalog(self, catalog, record):
        # Get Naive measurements and add to catalog
        catalog[0].set(self.keys["flux"], record.get("ext_trailedSources_Naive_flux"))
        catalog[0].set(self.keys["L"], record.get("ext_trailedSources_Naive_length"))
        catalog[0].set(self.keys["theta"], record.get("ext_trailedSources_Naive_angle"))

    def addTrailedSource(self, trail):
        """Add a trailed source to the simulation.
        """

        record = self.catalog.addNew()
        record.set(self.keys["centroid"], trail.center)
        covariance = np.random.normal(0, 0.1, 4).reshape(2, 2)
        covariance[0, 1] = covariance[1, 0]
        record.set(self.keys["centroid_sigma"], covariance.astype(np.float32))
        record.set(self.keys["shape"], self.psfShape)
        record.set(self.keys["isStar"], False)

        numIter = int(5*trail.length)
        xp = np.linspace(trail.x0, trail.x1, num=numIter)
        yp = np.linspace(trail.y0, trail.y1, num=numIter)
        for (x, y) in zip(xp, yp):
            pt = lsst.geom.Point2D(x, y)
            im = self.drawGaussian(self.exposure.getBBox(), trail.instFlux,
                                   lsst.afw.geom.Ellipse(self.psfShape, pt))
            self.exposure.getMaskedImage().getImage().getArray()[:, :] += im.getArray()

        record.set(self.keys["instFlux"], self.exposure.getImage().array.sum())
        self._installFootprint(record, self.exposure.getImage())

        return record, self.exposure.getImage()


class TrailedSourcesTestCase(AlgorithmTestCase, lsst.utils.tests.TestCase):
    # Following from meas_base/test_NaiveCentroid.py
    # Taken from NaiveCentroidTestCase
    @methodParameters(L=Ls, theta=thetas, xc=xcs, yc=ycs)
    def setUp(self, L, theta, xc, yc):
        self.center = lsst.geom.Point2D(50.1, 49.8)
        self.bbox = lsst.geom.Box2I(lsst.geom.Point2I(-20, -30),
                                    lsst.geom.Extent2I(140, 160))
        self.dataset = TrailedTestDataset(self.bbox)

        self.trail = TrailedSource(10000.0, L, theta, xc, yc)
        self.dataset.addTrailedSource(self.trail)

    def tearDown(self):
        del self.center
        del self.bbox
        del self.trail
        del self.dataset

    def makeTrailedSourceMeasurementTask(self, plugin=None, dependencies=(),
                                         config=None, schema=None, algMetadata=None):
        config = self.makeSingleFrameMeasurementConfig(plugin=plugin,
                                                       dependencies=dependencies)
        config.slots.shape = "base_SdssShape"
        return self.makeSingleFrameMeasurementTask(plugin=plugin,
                                                   dependencies=dependencies,
                                                   config=config, schema=schema,
                                                   algMetadata=algMetadata)

    def testNaivePlugin(self):
        task = self.makeTrailedSourceMeasurementTask(
            plugin="ext_trailedSources_Naive",
            dependencies=("base_SdssCentroid", "base_SdssShape")
        )
        exposure, catalog = self.dataset.realize(10.0, task.schema, randomSeed=0)
        task.run(catalog, exposure)
        record = catalog[0]

        # Compare true with measured length and angle
        L = record.get("ext_trailedSources_Naive_length")
        theta = record.get("ext_trailedSources_Naive_angle")
        self.assertFloatsAlmostEqual(L, self.trail.length, atol=None, rtol=.1)
        self.assertFloatsAlmostEqual(theta, self.trail.angle, atol=None, rtol=.05)

        # Check test setup
        self.assertNotEqual(L, self.trail.length)
        self.assertNotEqual(theta, self.trail.angle)

    def testVeresPlugin(self):
        # First run the Naive method
        taskNaive = self.makeTrailedSourceMeasurementTask(
            plugin="ext_trailedSources_Naive", dependencies=("base_SdssShape",)
        )
        expNaive, catNaive = self.dataset.realize(10.0, taskNaive.schema, randomSeed=0)
        taskNaive.run(catNaive, expNaive)

        # Now save the Naive measurements to the new catalog
        taskVeres = self.makeTrailedSourceMeasurementTask(
            plugin="ext_trailedSources_Veres", dependencies=("base_SdssShape",)
        )
        schema = self.dataset.addToSchema(taskVeres.schema)
        expVeres, catVeres = self.dataset.realize(10.0, schema, randomSeed=0)
        self.dataset.addToCatalog(catVeres, catNaive[0])

        # And run the Veres measurement
        taskVeres.run(catVeres, expVeres)
        record = catVeres[0]

        # Compare measured trail length and angle to true values
        L = record.get("ext_trailedSources_Veres_length")
        theta = record.get("ext_trailedSources_Veres_angle")
        self.assertFloatsAlmostEqual(L, self.trail.length, atol=None, rtol=.05)
        self.assertFloatsAlmostEqual(theta, self.trail.angle, atol=None, rtol=.01)

        # Make sure test setup is working as expected
        self.assertNotEqual(L, self.trail.length)
        self.assertNotEqual(theta, self.trail.angle)

        # Test that reduced chi-squared is reasonable
        rChiSq = record.get("ext_trailedSources_Veres_rChiSq")
        self.assertGreater(rChiSq, 0.9)
        self.assertLess(rChiSq, 1.1)


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
