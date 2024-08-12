import milo_1_0_3.program_state as ps
from milo_1_0_3.atom import Atom
from milo_1_0_3.main import print_accelerations


def test_print_accelerations(capsys):
    program_state = ps.ProgramState()
    atom1 = Atom("H", 0, 0, 0)
    atom2 = Atom("O", 1, 1, 1)
    program_state.atoms = [atom1, atom2]
    #program_state.accelerations = [[0.1, 0.2, 0.3], [0.4, 0.5, 0.6]]
    accelerations = container.

    print_accelerations(program_state)

    captured = capsys.readouterr()
    assert captured.out == (
        "  Accelerations:\n"
        "    H  1.000000e-01  2.000000e-01  3.000000e-01\n"
        "    O  4.000000e-01  5.000000e-01  6.000000e-01\n"
    )
