Base examples
*************

Let's define a small helper function::

    def print_me(msg, val):

        print("{}: {}".format(msg, val))

Let's print the tolerance::

    print_me("The default value for the tolerance is", TOL)

    The default value for the tolerance is: 1e-10

Find the accuracy of this computer::

    j, d = machine_accuracy()


    print_me("Number of significant BITS in the mantissa\t", j)

    Number of significant BITS in the mantissa	: 52.0


    print_me("Number of significant DIGITS in a decimal number", d)

    Number of significant DIGITS in a decimal number: 15

Print the suffixes for some ordinal numbers::

    print_me("The suffix for ordinal 2 is", get_ordinal_suffix(2))

    The suffix for ordinal 2 is: nd


    print_me("The suffix for ordinal 16 is", get_ordinal_suffix(16))

    The suffix for ordinal 16 is: th

