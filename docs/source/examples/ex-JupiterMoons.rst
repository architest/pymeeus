JupiterMoons examples
**********************

Let's define a small helper function::

    def print_me(msg, val):
        print("{}: {}".format(msg, val))


Lets compute the ascending node of Jupiter as well as the longitude of the node of
the equator of Jupiter on the ecliptic (psi)::

    utc_1992_12_16_00_00_00 = Epoch(1992, 12, 16, utc=True)
    psi_corrected, OMEGA_ascending_node_jupiter = JupiterMoons.jupiter_system_angles(utc_1992_12_16_00_00_00)
    print("Ascending node of Jupiter: ", OMEGA_ascending_node_jupiter)
    #100.39249942976576
    print("Longitude of the node of the eauator of Jupiter on the ecliptic (psi):", psi_corrected)
    #317.1058009213959

Lets compute the corrected rectangular geocentric position of Jupiter's satellites
for a given epoch, using the E5-theory::

    utc_1992_12_16_00_00_00 = Epoch(1992, 12, 16, utc=True)

    io, europa, ganymede, callisto = JupiterMoons.rectangular_positions_jovian_equatorial(utc_1992_12_16_00_00_00)

    print("Corrected rectangular geocentric position of Io [X, Y , Z]: ", io)
    #(-3.450168811390241, 0.21370246960509387, -4.818966623735296)

    print("Corrected rectangular geocentric position of Europa [X, Y , Z]: ", europa)
    #(7.441869121153001, 0.27524463479625677, -5.747104399729193)

    print("Corrected rectangular geocentric position of Ganymede [X, Y , Z]: ", ganymede)
    #(1.201111684800708, 0.5899903274317162, -14.940581367576527)

    print("Corrected rectangular geocentric position of Callisto [X, Y , Z]: ", callisto)
    #(7.071943240286434, 1.0289562923230684, -25.224137724734955)

Lets compute the uncorrected rectangular geocentric position of Jupiter's satellites for a given epoch,
using the E5-theory::

    #So the effects of different light-time and perspective described in Pymeeus page 313 - 314 are neglected

    utc_1992_12_16_00_00_00 = Epoch(1992, 12, 16, utc=True)
    io_uncorrected, europa_uncorrected, ganymede_uncorrected, callisto_uncorrected = \
        JupiterMoons.rectangular_positions_jovian_equatorial(utc_1992_12_16_00_00_00, do_correction=False)

    print("Uncorrected rectangular geocentric position of Io [X, Y , Z]: ", io_uncorrected)
    # (-3.4489935969836503, 0.21361563816963675, -4.818966623735296)

    print("Uncorrected rectangular geocentric position of Europa [X, Y , Z]: ", europa_uncorrected)
    # (7.438101803124541, 0.2751112576349763, -5.747104399729193)

    print("Uncorrected rectangular geocentric position of Ganymede [X, Y , Z]: ", ganymede_uncorrected)
    # (1.1990581804888616, 0.589247092847632, -14.940581367576527)

    print("Uncorrected rectangular geocentric position of Callisto [X, Y , Z]: ", callisto_uncorrected)
    # (7.056237832405445, 1.0267678919629089, -25.224137724734955)

Lets calculate the distance between Earth and Jupiter (DELTA) for a given epoch::

    utc_1992_12_16_00_00_00 = Epoch(1992, 12, 16, utc=True)
    delta, tau, l, b, r = JupiterMoons.calculate_delta(utc_1992_12_16_00_00_00)

    print("Distance between Earth and Jupiter in AU: ", delta)
    #5.6611211815432645

    print("Light-time from Earth to Jupiter in d (day): ", tau)
    #0.03269590898252075

Lets calculate the perspective distance in Jupiter radii of all satellites
for an eclipse of Io::

    io_ecc_start_2021_02_12_14_19_14 = Epoch(2021, 2, 12.5966898148148)

    result_matrix = JupiterMoons.check_phenomena(io_ecc_start_2021_02_12_14_19_14)

    #structure of result_matrix
    # Row 0: Io          Column 0: perspective distance as seen from the Earth
    # Row 1: Europa      Column 1: perspective distance as seen from the Sun
    # Row 2: Ganymede    Column 2: No use
    # Row 3: Callisto

    # print Row 0
    print("(perspective distance of Io (Earth View), perspective distance of Io (Sun view), No use): ")
    print(result_matrix[0])
    #[1.1926058680144362, 0.856027716233023, 0.0]

    # print Row 1
    print("(perspective distance of Europa (Earth View), perspective distance of Europa (Sun view), No use): ")
    print(result_matrix[1])
    #[-8.739720236890856, -8.893094092124032, 0.0]

    # print Row 2
    print("(perspective distance of Ganymede (Earth View), perspective distance of Ganymede (Sun view), No use): ")
    print(result_matrix[2])
    #[14.069121992481382, 13.8323491767871, 0.0]

    # print Row 3
    print("(perspective distance of Callisto (Earth View), perspective distance of Callisto (Sun view), No use): ")
    print(result_matrix[3])
    #[-2.934134686233644, -3.9904786452498144, 0.0]

Lets check if an eclipse or\and occultation for any of the four Galilean satellites is detected for a given epoch::

    io_ecc_start_2021_02_12_14_19_14 = Epoch(2021, 2, 12.5966898148148)

    #Structure of result matrix
    # Row 0: Io          Column 0: Occultation True\False
    # Row 1: Europa      Column 1: Eclipse True\False
    # Row 2: Ganymede    Column 2: No use
    # Row 3: Callisto

    result_matrix = JupiterMoons.is_phenomena(io_ecc_start_2021_02_12_14_19_14)

    #print Row 0
    print("(Occultation of Io, Eclipse of Io, No use): ")
    print(result_matrix[0])
    #[False, True, False]

    # print Row 1
    print(" (Occultation of Europa, Eclipse of Europa, No use): ")
    print(result_matrix[1])
    #[False, False, False]

    # print Row 2
    print(" (Occultation of Ganymede, Eclipse of Gaymede, No use): ")
    print(result_matrix[2])
    #[False,False,False]

    # print Row 3
    print("(Occultation of Callisto, Eclipse of Callisto, No use): ")
    print(result_matrix[3])
    #[False,False,False]

Calculation of the perspective distance ot the planet Io to the center of Jupiter
for December 16 at 0h UTC as seen from the Sun::

    utc_1992_12_16_00_00_00 = Epoch(1992, 12, 16, utc=True)

    result_matrix = JupiterMoons.rectangular_positions_jovian_equatorial(utc_1992_12_16_00_00_00, solar=True)

    #Structure of result matrix
    # Row 0: Io          Column 0: X coordinate of satellite in Jupiter radii
    # Row 1: Europa      Column 1: Y coordinate of satellite in Jupiter radii
    # Row 2: Ganymede    Column 2: Z coordinate of satellite in Jupiter radii
    # Row 3: Callisto

    io_radius_to_center_of_jupiter_sun = JupiterMoons.check_coordinates(result_matrix[0][0], result_matrix[0][1])

    print("Perspective distance of Io as seen from the Sun in Jupiter radii: ", io_radius_to_center_of_jupiter_sun)
    #3.457757270630766

Calculation of the perspective distance ot the planet Io to the center of Jupiter
for December 16 at 0h UTC as seen from the Earth::

    utc_1992_12_16_00_00_00 = Epoch(1992, 12, 16, utc=True)
    result_matrix = JupiterMoons.rectangular_positions_jovian_equatorial(utc_1992_12_16_00_00_00, solar=False)

    #Structure of result matrix
    # Row 0: Io          Column 0: X coordinate of satellite in Jupiter radii
    # Row 1: Europa      Column 1: Y coordinate of satellite in Jupiter radii
    # Row 2: Ganymede    Column 2: Z coordinate of satellite in Jupiter radii
    # Row 3: Callisto

    io_radius_to_center_of_jupiter_earth = JupiterMoons.check_coordinates(result_matrix[0][0], result_matrix[0][1])

    print("Perspective distance of Io as seen from the Earth in Jupiter radii: ", io_radius_to_center_of_jupiter_earth)
    # 2.553301264153796

