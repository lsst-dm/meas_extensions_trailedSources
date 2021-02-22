import numpy as np
import unittest
import lsst.utils.tests
import lsst.meas.extensions.trailedSources
from lsst.meas.base.tests import AlgorithmTestCase
from lsst.utils.tests import methodParameters

# Trailed source end points (chosen at random)
x0s = [-4.3561, 11.5712, 5.5694, -9.5507, 1.5940, 4.1805]
y0s = [-13.2116, 19.8783, 5.8372, -8.8392, 4.8078, -17.4946]
x1s = [2.1383, 2.9489, 13.6870, 18.0673, 13.2956, 13.2345]
y1s = [16.5114, 15.8209, 14.3787, 19.4292, 8.9292, 8.9737]


class TrailedSource:
    def __init__(self, instFlux, head, tail, center=None):
        self.instFlux = instFlux
        self.head = head
        self.tail = tail

        x0 = head.getX()
        y0 = head.getY()
        x1 = tail.getX()
        y1 = tail.getY()

        if center is None:
            center = lsst.geom.Point2D(0.5*(x1 + x0), 0.5*(y0 + y1))
        self.center = center

        self.length = np.sqrt((x1-x0)**2 + (y1-y0)**2)
        self.angle = np.arctan2(y1-y0, x1-x0)


# "Extend" meas.base.tests.TestDataset
class TrailedTestDataset(lsst.meas.base.tests.TestDataset):

    @classmethod
    def makeMinimalSchema(cls):
        schema = super().makeMinimalSchema()
        cls.keys["x0"] = schema.addField("x0", type="D", doc="true x0", units="pixel")
        cls.keys["y0"] = schema.addField("y0", type="D", doc="true y0", units="pixel")
        cls.keys["x1"] = schema.addField("x1", type="D", doc="true x1", units="pixel")
        cls.keys["y1"] = schema.addField("y1", type="D", doc="true y1", units="pixel")
        # Add head and tail errors later

        return schema

    def __init__(self, bbox, threshold=10.0, exposure=None, **kwds):
        super().__init__(bbox, threshold, exposure, **kwds)

    def addTrailedSource(self, trail):
        """Add a trailed source to the simulation.
        """
        x0 = trail.head.getX()
        y0 = trail.head.getY()
        x1 = trail.tail.getX()
        y1 = trail.tail.getY()

        record = self.catalog.addNew()
        record.set(self.keys["instFlux"], trail.instFlux)
        record.set(self.keys["centroid"], trail.center)
        record.set(self.keys["x0"], x0)
        record.set(self.keys["y0"], y0)
        record.set(self.keys["x1"], x1)
        record.set(self.keys["y1"], y1)
        covariance = np.random.normal(0, 0.1, 4).reshape(2, 2)
        covariance[0, 1] = covariance[1, 0]
        record.set(self.keys["centroid_sigma"], covariance.astype(np.float32))
        record.set(self.keys["shape"], self.psfShape)

        numIter = int(5*np.sqrt((x1 - x0)**2 + (y1 - y0)**2))
        xp = np.linspace(x0, x1, num=numIter)
        yp = np.linspace(y0, y1, num=numIter)
        for (x, y) in zip(xp, yp):
            pt = lsst.geom.Point2D(x, y)
            im = self.drawGaussian(self.exposure.getBBox(), trail.instFlux,
                                   lsst.afw.geom.Ellipse(self.psfShape, pt))
            self.exposure.getMaskedImage().getImage().getArray()[:, :] += im.getArray()

        self._installFootprint(record, self.exposure.getImage())

        return record, self.exposure.getImage()


class NaiveTrailTestCase(AlgorithmTestCase, lsst.utils.tests.TestCase):
    # Following from meas_base/test_NaiveCentroid.py
    # Taken from NaiveCentroidTestCase
    @methodParameters(x0=x0s, y0=y0s, x1=x1s, y1=y1s)
    def setUp(self, x0, y0, x1, y1):
        self.center = lsst.geom.Point2D(50.1, 49.8)
        self.bbox = lsst.geom.Box2I(lsst.geom.Point2I(-20, -30),
                                    lsst.geom.Extent2I(140, 160))
        self.dataset = TrailedTestDataset(self.bbox)

        head = lsst.geom.Point2D(x0, y0)
        tail = lsst.geom.Point2D(x1, y1)
        self.trail = TrailedSource(10000.0, head, tail)
        self.dataset.addTrailedSource(self.trail)

    def tearDown(self):
        del self.center
        del self.bbox
        del self.trail
        del self.dataset

    def testSingleFramePlugin(self):
        task = self.makeSingleFrameMeasurementTask("ext_trailedSources_Naive")
        exposure, catalog = self.dataset.realize(10.0, task.schema, randomSeed=0)
        task.run(catalog, exposure)
        record = catalog[0]

        x0 = record.get("ext_trailedSources_Naive_x0")
        y0 = record.get("ext_trailedSources_Naive_y0")
        x1 = record.get("ext_trailedSources_Naive_x1")
        y1 = record.get("ext_trailedSources_Naive_y1")
        self.assertFloatsAlmostEqual(x0, self.trail.head.getX(), atol=None, rtol=.02)
        self.assertFloatsAlmostEqual(y0, self.trail.head.getY(), atol=None, rtol=.02)
        self.assertFloatsAlmostEqual(x1, self.trail.tail.getX(), atol=None, rtol=.02)
        self.assertFloatsAlmostEqual(y1, self.trail.tail.getY(), atol=None, rtol=.02)


class VeresTrailTestCase(AlgorithmTestCase, lsst.utils.tests.TestCase):
    # Following from meas_base/test_NaiveCentroid.py
    # Taken from NaiveCentroidTestCase
    @methodParameters(x0=x0s, y0=y0s, x1=x1s, y1=y1s)
    def setUp(self, x0, y0, x1, y1):
        self.center = lsst.geom.Point2D(50.1, 49.8)
        self.bbox = lsst.geom.Box2I(lsst.geom.Point2I(-20, -30),
                                    lsst.geom.Extent2I(140, 160))
        self.dataset = TrailedTestDataset(self.bbox)

        head = lsst.geom.Point2D(x0, y0)
        tail = lsst.geom.Point2D(x1, y1)
        self.trail = TrailedSource(10000.0, head, tail)
        self.dataset.addTrailedSource(self.trail)

    def tearDown(self):
        del self.center
        del self.bbox
        del self.trail
        del self.dataset

    def addToSchema(self, schema):
        # Add Naive keys to schema
        name = "ext_trailedSources_Naive"
        self.keyX0 = schema.addField(name + "_x0", type="D", doc="Trail head X coordinate.", units="pixel")
        self.keyY0 = schema.addField(name + "_y0", type="D", doc="Trail head Y coordinate.", units="pixel")
        self.keyX1 = schema.addField(name + "_x1", type="D", doc="Trail tail X coordinate.", units="pixel")
        self.keyY1 = schema.addField(name + "_y1", type="D", doc="Trail tail Y coordinate.", units="pixel")
        self.keyFlux = schema.addField(name + "_flux", type="D", doc="Trailed source flux.", units="count")

        return schema

    def addToCatalog(self, catalog, record):
        # Get Naive measurements and add to catalog
        catalog[0].set(self.keyX0, record["ext_trailedSources_Naive_x0"])
        catalog[0].set(self.keyY0, record["ext_trailedSources_Naive_y0"])
        catalog[0].set(self.keyX1, record["ext_trailedSources_Naive_x1"])
        catalog[0].set(self.keyY1, record["ext_trailedSources_Naive_y1"])
        catalog[0].set(self.keyFlux, self.trail.instFlux)

    def testSingleFramePlugin(self):
        # First run the Naive method
        taskNaive = self.makeSingleFrameMeasurementTask("ext_trailedSources_Naive")
        expNaive, catNaive = self.dataset.realize(10.0, taskNaive.schema, randomSeed=0)
        taskNaive.run(catNaive, expNaive)

        # Now save the Naive measurements to the new catalog
        taskVeres = self.makeSingleFrameMeasurementTask("ext_trailedSources_Veres")
        schema = self.addToSchema(taskVeres.schema)
        expVeres, catVeres = self.dataset.realize(10.0, schema, randomSeed=0)
        self.addToCatalog(catVeres, catNaive[0])

        # And run the Veres measurement
        taskVeres.run(catVeres, expVeres)
        record = catVeres[0]

        # Test that end points are close to true value
        x0 = record.get("ext_trailedSources_Veres_x0")
        y0 = record.get("ext_trailedSources_Veres_y0")
        x1 = record.get("ext_trailedSources_Veres_x1")
        y1 = record.get("ext_trailedSources_Veres_y1")
        self.assertFloatsAlmostEqual(x0, self.trail.head.getX(), atol=None, rtol=.02)
        self.assertFloatsAlmostEqual(y0, self.trail.head.getY(), atol=None, rtol=.02)
        self.assertFloatsAlmostEqual(x1, self.trail.tail.getX(), atol=None, rtol=.02)
        self.assertFloatsAlmostEqual(y1, self.trail.tail.getY(), atol=None, rtol=.02)

        # Test that reduced chi-squared in [0.9,1.1]
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
