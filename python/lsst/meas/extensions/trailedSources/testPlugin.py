# Testing a new plugin

from lsst.meas.base import SingleFramePlugin, SingleFramePluginConfig
from lsst.meas.base.pluginRegistry import register

__all__ = ("SingleFrameTestConfig", "SingleFrameTestPlugin")


class SingleFrameTestConfig(SingleFramePluginConfig):
    """Config class for Test plugin
    """
    pass


@register("ext_trailedSources_Test")
class SingleFrameTestPlugin(SingleFramePlugin):
    """Test SingleFramePlugin.
    """

    ConfigClass = SingleFrameTestConfig

    @classmethod
    def getExecutionOrder(cls):
        return cls.CENTROID_ORDER  # change this later

    def __init__(self, config, name, schema, metadata):
        super().__init__(config, name, schema, metadata)
        self.keyTest = schema.addField(name + "_test", type="D", doc="TESTING PLUGIN", units="pixel")
        self.flag = schema.addField(name + "_flag", type="Flag", doc="Test Failed")

    def measure(self, measRecord, exposure):
        """Set self.keyTest = 1.0"""
        measRecord.set(self.keyTest, 1.)

    def fail(self, measRecord, error=None):
        """Record failure"""
        measRecord.set(self.flag, True)
