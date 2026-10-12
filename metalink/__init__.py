try:
    from importlib import metadata
except ImportError as e:
    import importlib_metadata as metadata

__version__ = metadata.version('Metalink')
