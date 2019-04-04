Epoch examples
**************

Let's define a small helper function::

    def print_me(msg, val):

        print("{}: {}".format(msg, val))

Let's start creating and Epoch object, and printing it::

    e = Epoch(1987, 6, 19.5)

    print_me("JDE for 1987/6/19.5", e)

    # JDE for 1987/6/19.5: 2446966.00064

Redefine the Epoch object::

    e.set(333, 'Jan', 27, 12)

    print_me("JDE for 333/1/27.5", e)

    # JDE for 333/1/27.5: 1842713.0

We can create an Epoch from a ``date`` or ``datetime`` object::

    d = datetime.datetime(837, 4, 10, 7, 12, 0, 0)

    f = Epoch(d)

    print_me("JDE for 837/4/10.3", f)

    # JDE for 837/4/10.3: 2026871.8

Let's check if a given date belong to the Julian or the Gregorian calendar::

    print_me("Is 1590/4/21.4 a Julian date?", Epoch.is_julian(1590, 4, 21.4))

    # Is 1590/4/21.4 a Julian date?: False

We can also check if a given year is leap or not::

    print_me("Is -1000 a leap year?", Epoch.is_leap(-1000))

    # Is -1000 a leap year?: True

    print_me("Is 1800 a leap year?", Epoch.is_leap(1800))

    # Is 1800 a leap year?: False

    print_me("Is 2012 a leap year?", Epoch.is_leap(2012))

    # Is 2012 a leap year?: True

Get the Day Of Year (DOY) corresponding to a given date::

    print_me("Day Of Year (DOY) of 1978/11/14", Epoch.get_doy(1978, 11, 14))

    # Day Of Year (DOY) of 1978/11/14: 318.0

    print_me("Day Of Year (DOY) of -400/2/29.9", Epoch.get_doy(-400, 2, 29.9))

    # Day Of Year (DOY) of -400/2/29.9: 60.9

Now the opposite: Get a date from a DOY::

    t = Epoch.doy2date(2017, 365.7)

    s = str(t[0]) + "/" + str(t[1]) + "/" + str(round(t[2], 2))

    print_me("Date from DOY 2017:365.7", s)

    # Date from DOY 2017:365.7: 2017/12/31.7

    t = Epoch.doy2date(-4, 60)

    s = str(t[0]) + "/" + str(t[1]) + "/" + str(round(t[2], 2))

    print_me("Date from DOY -4:60", s)

    # Date from DOY -4:60: -4/2/29.0


There is an internal table which we can use to get the leap seconds::

    print_me("Number of leap seconds applied up to July 1983", Epoch.leap_seconds(1983, 7))

    # Number of leap seconds applied up to July 1983: 12

We can convert the internal JDE value back to a date::

    e = Epoch(2436116.31)

    y, m, d = e.get_date()

    s = str(y) + "/" + str(m) + "/" + str(round(d, 2))

    print_me("Date from JDE 2436116.31", s)

    # Date from JDE 2436116.31: 1957/10/4.81

It is possible to get the day of the week corresponding to a given date::

    e = Epoch(2018, 'Feb', 15)

    print_me("The day of week of 2018/2/15 is", e.dow(as_string=True))

    # The day of week of 2018/2/15 is: Thursday

In some cases it is useful to get the Modified Julian Day (MJD)::

    e = Epoch(1923, 'August', 23)

    print_me("Modified Julian Day for 1923/8/23", round(e.mjd(), 2))

    # Modified Julian Day for 1923/8/23: 23654.0

If your system is appropriately configured, you can get the difference in seconds between your local time and UTC::

    print_me("From local system time to UTC you must add/subtract" +

             " this amount of seconds", Epoch.utc2local())

    # From local system time to UTC you must add/subtract this amount of seconds: 7200.0

Compute DeltaT = TT - UT differences for various dates::

    print_me("DeltaT (TT - UT) for Feb/333", round(Epoch.tt2ut(333, 2), 1))

    # DeltaT (TT - UT) for Feb/333: 7358.5

    print_me("DeltaT (TT - UT) for Jan/1642", round(Epoch.tt2ut(1642, 1), 1))

    # DeltaT (TT - UT) for Jan/1642: 62.1

    print_me("DeltaT (TT - UT) for Feb/1928", round(Epoch.tt2ut(1928, 1), 1))

    # DeltaT (TT - UT) for Feb/1928: 24.2

    print_me("DeltaT (TT - UT) for Feb/1977", round(Epoch.tt2ut(1977, 2), 1))

    # DeltaT (TT - UT) for Feb/1977: 47.7

    print_me("DeltaT (TT - UT) for Jan/1998", round(Epoch.tt2ut(1998, 1), 1))

    # DeltaT (TT - UT) for Jan/1998: 63.0

The difference between civil day and sidereal day is almost 4 minutes::

    e = Epoch(1987, 4, 10)

    st1 = round(e.mean_sidereal_time(), 9)

    e = Epoch(1987, 4, 11)

    st2 = round(e.mean_sidereal_time(), 9)

    ds = (st2 - st1)*DAY2MIN

    msg = "{}m {}s".format(INT(ds), (ds % 1)*60.0)

    print_me("Difference between sidereal time 1987/4/11 and 1987/4/10", msg)

    # Difference between sidereal time 1987/4/11 and 1987/4/10: 3m 56.555424s

When correcting for nutation-related effects, we get the **apparent** sidereal time::

    e = Epoch(1987, 4, 10)

    print_me("e.apparent_sidereal_time(23.44357, (-3.788)/3600.0)",

             e.apparent_sidereal_time(23.44357, (-3.788)/3600.0))

    # e.apparent_sidereal_time(23.44357, (-3.788)/3600.0): 0.549145082637

Epoch class can also provide the date of Easter for a given year. Let's spice up the output a little bit, calling ``dow()`` and ``get_month()``::

    month, day = Epoch.easter(2019)

    e = Epoch(2019, month, day)

    s = e.dow(as_string=True) + ", " + str(day) + get_ordinal_suffix(day) + \

        " of " + Epoch.get_month(month, as_string=True)

    print_me("Easter day for 2019", s)

    # Easter day for 2019: Sunday, 21st of April

Compute the date of the Jewish Easter (Pesach) for a given year::

    month, day = Epoch.jewish_pesach(1990)

    s = str(day) + get_ordinal_suffix(day) + " of " + Epoch.get_month(month, as_string=True)

    print_me("Jewish Pesach day for 1990", s)

    # Jewish Pesach day for 1990: 10th of April

Now, let's convert a date in the Moslem calendar to the Gregorian calendar::

    y, m, d = Epoch.moslem2gregorian(1421, 1, 1)

    print_me("The date 1421/1/1 in the Moslem calendar is, in Gregorian " +

             "calendar", "{}/{}/{}".format(y, m, d))

    # The date 1421/1/1 in the Moslem calendar is, in Gregorian calendar: 2000/4/6

    y, m, d = Epoch.moslem2gregorian(1439, 9, 1)

    print_me("The start of Ramadan month (9/1) for Gregorian year 2018 is",

             "{}/{}/{}".format(y, m, d))

    # The start of Ramadan month (9/1) for Gregorian year 2018 is: 2018/5/16

We can go from the Gregorian calendar back to the Moslem calendar too::

    print_me("Date 1991/8/13 in Gregorian calendar is, in Moslem calendar",

             "{}/{}/{}".format(*Epoch.gregorian2moslem(1991, 8, 13)))

    # Date 1991/8/13 in Gregorian calendar is, in Moslem calendar: 1412/2/2

.. note:: The ``*`` before ``Epoch`` will **unpack** the tuple into components


It is possible to carry out some algebraic operations with Epochs.

- Add 10000 days to a given date::

    a = Epoch(1991, 7, 11)

    b = a + 10000

    y, m, d = b.get_date()

    s = str(y) + "/" + str(m) + "/" + str(round(d, 2))

    print_me("1991/7/11 plus 10000 days is", s)

    # 1991/7/11 plus 10000 days is: 2018/11/26.0

- Subtract two Epochs to find the number of days between them::

    a = Epoch(1986, 2, 9.0)

    b = Epoch(1910, 4, 20.0)

    print_me("The number of days between 1986/2/9 and 1910/4/20 is", round(a - b, 2))

    # The number of days between 1986/2/9 and 1910/4/20 is: 27689.0

- We can also subtract a given amount of days from an Epoch::

    a = Epoch(2003, 12, 31.0)

    b = a - 365.5

    y, m, d = b.get_date()

    s = str(y) + "/" + str(m) + "/" + str(round(d, 2))

    print_me("2003/12/31 minus 365.5 days is", s)

    # 2003/12/31 minus 365.5 days is: 2002/12/30.5

- Accumulative addition and subtraction of days is also allowed::

    a = Epoch(2003, 12, 31.0)

    a += 32.5

    y, m, d = a.get_date()

    s = str(y) + "/" + str(m) + "/" + str(round(d, 2))

    print_me("2003/12/31 plus 32.5 days is", s)

    # 2003/12/31 plus 32.5 days is: 2004/2/1.5



    a = Epoch(2001, 12, 31.0)

    a -= 2*365

    y, m, d = a.get_date()

    s = str(y) + "/" + str(m) + "/" + str(round(d, 2))

    print_me("2001/12/31 minus 2*365 days is", s)

    # 2001/12/31 minus 2*365 days is: 2000/1/1.0

- It is also possible to add days from the right::

    a = Epoch(2004, 2, 27.8)

    b = 2.2 + a

    y, m, d = b.get_date()

    s = str(y) + "/" + str(m) + "/" + str(round(d, 2))

    print_me("2.2 days plus 2004/2/27.8 is", s)

    # 2.2 days plus 2004/2/27.8 is: 2004/3/1.0

- Comparison operadors between epochs are also defined::

    a = Epoch(2007, 5, 20.0)

    b = Epoch(2007, 5, 20.000001)

    print_me("2007/5/20.0 == 2007/5/20.000001", a == b)

    # 2007/5/20.0 == 2007/5/20.000001: False

    print_me("2007/5/20.0 != 2007/5/20.000001", a != b)

    # 2007/5/20.0 != 2007/5/20.000001: True

    print_me("2007/5/20.0 > 2007/5/20.000001", a > b)

    # 2007/5/20.0 > 2007/5/20.000001: False

    print_me("2007/5/20.0 <= 2007/5/20.000001", a <= b)

    # 2007/5/20.0 <= 2007/5/20.000001: True

- Compute the time of rise and setting of the Sun in a given day::

    e = Epoch(2018, 5, 2)

    print("On May 2nd, 2018, Sun rising/setting times in Munich were (UTC):")

    latitude = Angle(48, 8, 0)

    longitude = Angle(11, 34, 0)

    altitude = 520.0

    rising, setting = e.rise_set(latitude, longitude, altitude)

    y, m, d, h, mi, s = rising.get_full_date()

    print("Rising time: {}:{}".format(h, mi))

    # Rising time: 3:50

    y, m, d, h, mi, s = setting.get_full_date()

    print("Setting time: {}:{}".format(h, mi))

    # Setting time: 18:33
