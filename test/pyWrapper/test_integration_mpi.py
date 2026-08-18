from test.utils.utils import *


@no_mpi_skip
@pytest.mark.mpi
@pytest.mark.wires
def test_airplane_case_with_mpi(tmp_path):
    fn = CASES_FOLDER + "airplane/airplane.fdtd.json"
    solver = FDTD(
        fn,
        path_to_exe=SEMBA_EXE,
        run_in_folder=tmp_path,
        flags=["-mapvtk"],
        mpi_command="mpirun -np 2",
    )
    solver.run()

    vtkmapfile = solver.getVTKMap()
    assert os.path.isfile(vtkmapfile)


@no_mpi_skip
@pytest.mark.mpi
@pytest.mark.wires
def test_simple_cabin_initialization_with_mpi(tmp_path):
    fn = CASES_FOLDER + "simple_cabin/simple_cabin.fdtd.json"
    solver = FDTD(
        fn,
        path_to_exe=SEMBA_EXE,
        run_in_folder=tmp_path,
        flags=["-mapvtk"],
        mpi_command="mpirun -np 2",
    )
    solver.run()

    vtkmapfile = solver.getVTKMap()
    assert os.path.isfile(vtkmapfile)
