# from milo_1_0_3 import program_state as ps
# from milo_1_0_3 import containers
# from milo_1_0_3.tools import parse_frequencies as pf


# def test_parse_gaussian_frequency_data():
#     input_iterable = [
#         "               Frequencies ---  1682.1354 3524.4296 3668.7401\n",
#         "            Reduced masses ---     1.0895    1.0389    1.0827\n",
#         "           Force constants ---     1.8163    7.6032    8.5864\n",
#         "            IR Intensities ---    52.8486    4.2243    0.3831\n",
#         "         Coord Atom Element:\n",
#         "           1     1     8         -0.00000   0.00000  -0.00000\n",
#         "           2     1     8          0.00000  -0.00000  -0.07070\n",
#         "           3     1     8         -0.07382   0.04553  -0.00000\n",
#         "           1     2     1          0.00000   0.00000   0.00000\n",
#         "           2     2     1          0.39258   0.60700   0.56106\n",
#         "           3     2     1          0.58580  -0.36126  -0.42745\n",
#         "           1     3     1          0.00000  -0.00000   0.00000\n",
#         "           2     3     1         -0.39258  -0.60700   0.56106\n",
#         "           3     3     1          0.58580  -0.36126   0.42745\n",
#         "           Harmonic frequencies (cm**-1), IR intensities (KM/Mole), Raman scatt\n",
#         "           activities (A**4/AMU), depolarization ratios for plane and unpolariz\n",
#     ]
#     program_state = ps.ProgramState()
#     pf.parse_gaussian_frequency_data(
#         input_iterable, program_state)
#     expected_frequencies = containers.Frequencies(
#         [1682.1354, 3524.4296, 3668.7401])
#     assert program_state.frequencies == expected_frequencies
#     assert program_state.reduced_masses == containers.Masses(
#         [1.0895, 1.0389, 1.0827])
#     assert program_state.force_constants == [1.8163, 7.6032, 8.5864]
#     assert program_state.mode_displacements == [
#         containers.Positions(
#             [(0.0, 0.0, 0.0), (0.0, -0.0, -0.0707), (-0.07382, 0.04553, -0.0)]
#         ),
#         containers.Positions(
#             [(0.0, 0.0, 0.0), (0.39258, 0.607, 0.56106),
#              (0.5858, -0.36126, -0.42745)]
#         ),
#         containers.Positions(
#             [(0.0, -0.0, 0.0), (-0.39258, -0.607, 0.56106),
#              (0.5858, -0.36126, 0.42745)]
#         ),
#     ]
