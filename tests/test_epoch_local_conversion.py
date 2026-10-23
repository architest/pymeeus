import datetime

import pytest

from pymeeus.Epoch import DAY2SEC, Epoch


@pytest.mark.parametrize("offset", [-14400, 20700, 50400, 0])
@pytest.mark.parametrize("utc", [
    datetime.datetime(2021, 7, 2, 0, 57, 11),
    datetime.datetime(2021, 12, 31, 23, 30),
    datetime.datetime(2022, 1, 1, 0, 30),
])
def test_local_epoch_conversion(monkeypatch, offset, utc):
    # Isolate conversion direction from how the current offset is obtained.
    monkeypatch.setattr(Epoch, "utc2local", staticmethod(lambda: offset))
    local = utc + datetime.timedelta(seconds=offset)
    local_args = (local.year, local.month, local.day, local.hour,
                  local.minute, local.second)
    utc_args = (utc.year, utc.month, utc.day, utc.hour, utc.minute,
                utc.second)
    reference = Epoch(*utc_args, utc=True)

    for options in [{"local": True}, {"local": True, "leap_seconds": 27}]:
        epoch = Epoch(*local_args, **options)
        assert abs(epoch - reference) * DAY2SEC < 0.0001

    for options in [{"local": True}, {"local": True, "leap_seconds": 27}]:
        actual = reference.get_full_date(**options)
        assert actual[:5] == local_args[:5]
        assert abs(actual[5] - local_args[5]) < 0.0001


def test_local_false_keeps_other_conversion_options(monkeypatch):
    monkeypatch.setattr(Epoch, "utc2local", staticmethod(lambda: 20700))
    epoch = Epoch(2021, 7, 2, 0, 57, 11, utc=True)
    assert epoch.get_full_date(local=False) == epoch.get_full_date()
    assert epoch.get_full_date(local=False, utc=True) == \
        epoch.get_full_date(utc=True)
