# from milo_1_0_3.program_state import ProgramState, enums
# from milo_1_0_3.input_parser import parse_input


# def test_parse_input():
#     program_state = ProgramState()  # Assuming ProgramState class exists
#     input_iterable = [
#         "$job",
#         "  current_step 100",
#         "  energy_boost on 0.1 0.5",
#         "  fixed_mode_direction 3 1",
#         "  fixed_vibrational_quanta 2 3",
#         "$molecule",
#         "  0 1",
#         "  C 0.0 0.0 0.0",
#         "  H 1.0 0.0 0.0",
#         "$end",
#     ]
#     parse_input(input_iterable, program_state)

#     # Add assertions here to verify the expected behavior of parse_input
#     assert program_state.current_step == 100
#     assert program_state.energy_boost == enums.EnergyBoost.ON
#     assert program_state.energy_boost_min == 0.1
#     assert program_state.energy_boost_max == 0.5
#     assert program_state.fixed_mode_directions == {3: 1}
#     assert program_state.fixed_vibrational_quanta == {2: 3}
#     assert program_state.charge == 0
#     assert program_state.spin == 1
#     assert len(program_state.atoms) == 2
#     assert program_state.atoms[0].symbol == "C"
#     assert program_state.atoms[1].symbol == "H"
#     # Add more assertions to verify the rest of the program_state object
