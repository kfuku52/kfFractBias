"""kfFractBias: offline fractionation-bias analysis."""

from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version("kffractbias")
except PackageNotFoundError:
    __version__ = "0+unknown"
