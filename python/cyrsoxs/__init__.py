from importlib import import_module as _import_module

# Load the compiled extension built by CMake
_ext = _import_module("CyRSoXS")

# Re-export everything from the extension
globals().update({name: getattr(_ext, name) for name in dir(_ext) if not name.startswith("_")})
