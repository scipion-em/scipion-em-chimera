# vim: set expandtab shiftwidth=4 softtabstop=4:

from .constants import *
try:
    from .api import _MyAPI
    # Create the ``bundle_api`` object that ChimeraX expects.
    bundle_api = _MyAPI()
except:
    pass

