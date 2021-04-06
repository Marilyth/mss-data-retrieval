from juregrid3d import *
from . import misc
try:
    from .version import VERSION as __version__
    from .version import REVISION as __revision__
except ImportError:
    __version__ = "unbuilt-dev"
    __revision__ = "unbuilt-dev"
