import pytest


def test_import_cyrsoxs_extension():
    CyRSoXS = pytest.importorskip("CyRSoXS")
    assert CyRSoXS.__doc__ is None or isinstance(CyRSoXS.__doc__, str)
