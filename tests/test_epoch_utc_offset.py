import calendar
import datetime
import importlib

import pytest


@pytest.mark.parametrize("utc, offset", [
    (datetime.datetime(2021, 7, 8, 0), -4 * 3600),
    (datetime.datetime(2021, 7, 7, 23, 30), 20700),
    (datetime.datetime(2021, 12, 31, 23, 45), 14 * 3600),
    (datetime.datetime(2021, 1, 1, 2), -12600),
    (datetime.datetime(2021, 7, 7, 12), 20700),
    (datetime.datetime(2021, 7, 7, 12), 0),
])
def test_current_utc_offset(monkeypatch, utc, offset):
    epoch_module = importlib.import_module("pymeeus.Epoch")
    local = utc + datetime.timedelta(seconds=offset)
    timestamp = calendar.timegm(utc.timetuple())

    class Clock:
        @staticmethod
        def now():
            return local

        @staticmethod
        def utcnow():
            return utc

    class CalendarClock:
        datetime = Clock

    class SystemClock:
        @staticmethod
        def time():
            return timestamp

        @staticmethod
        def localtime(value):
            return local.timetuple()

        @staticmethod
        def gmtime(value):
            return utc.timetuple()

    monkeypatch.setattr(epoch_module, "datetime", CalendarClock)
    monkeypatch.setattr(epoch_module, "time", SystemClock, raising=False)
    assert epoch_module.Epoch.utc2local() == offset
