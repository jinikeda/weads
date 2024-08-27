import os
import shutil
import pytest
import pandas as pd
import numpy as np
from unittest import mock
from io import StringIO
from click.testing import CliRunner
from numpy import dtype

from .general_functions import *
from .tidaldatums import *
from .tidaldatumsidw import *  # Import the interpolation function
# from src.CRMS_Discrete_Hydrographic2subsets import *
# from src.CRMS2Resample import *
# from src.CRMS2Plot import *
# from src.click_main import discrete_subcommand


# Sample fort.13 file
fort13_content = """\
fort13
6
2
primitive_weighting_in_continuity_equation
1
1
.030000
mannings_n_at_sea_floor
1
1
.025000
primitive_weighting_in_continuity_equation
2
4 0.020000
5 0.020000
mannings_n_at_sea_floor
4
1 0.030000
2 0.030000
3 0.030000
6 0.030000
"""


# Sample fort.14 file contents
fort14_content = """\
grid
4 6
1  0.0  0.0  6.27
2  1.0  0.0  7.13
3  2.0  0.0  8.11
4  0.0  1.0  9.09
5  1.0  1.0  9.05
6  2.0  1.0  0.00
1  3  1  2  5
2  3  2  3  6
3  3  1  5  4
4  3  2  6  5
"""


# Sample inundation.63 file contents
inundation63_content = """\
ESLR2021TCB ! 32 CHARACTER ALPH   AstronomicTides  ! 24 C  grid
2    6  0.1000000E+001        1     1 FileFmtVersion:    1050624
3.8880000000E+006        3888000
1     3.8880000000E+006
2     0.0000000000E+000
3     3.8880000000E+006
4     1.9440000000E+006
5     3.8880000000E+006
6     0.0000000000E+000
"""


# Sample inundation.63 file contents
fort53_content = """\
23
0.0000000000E+00  1.0000000   0.00000000  STEADY
0.2639203022E-05  1.0000000   0.00000000  MN
0.4925201824E-05  1.0000000   0.00000000  SM
0.6759774415E-04  1.0000000   0.00000000  O1
0.7292115836E-04  1.0000000   0.00000000  K1
0.1329544977E-03  1.0000000   0.00000000  MNS2
0.1355937007E-03  1.0000000   0.00000000  2MS2
0.1378796995E-03  1.0000000   0.00000000  N2
0.1405189025E-03  1.0000000   0.00000000  M2
0.1431581055E-03  1.0000000   0.00000000  2MN2
0.1454441043E-03  1.0000000   0.00000000  S2
0.1503693062E-03  1.0000000   0.00000000  2SM2
0.2783986020E-03  1.0000000   0.00000000  MN4
0.2810378050E-03  1.0000000   0.00000000  M4
0.2859630068E-03  1.0000000   0.00000000  MS4
0.4189175045E-03  1.0000000   0.00000000  2MN6
0.4215567075E-03  1.0000000   0.00000000  M6
0.4238427063E-03  1.0000000   0.00000000  MSN6
0.5620756100E-03  1.0000000   0.00000000  M8
0.7025945125E-03  1.0000000   0.00000000  M10
0.7252294600E-04  1.0000000   0.00000000  P1
0.1458423172E-03  1.0000000   0.00000000  K2
0.6495854113E-04  1.0000000   0.00000000  Q1
6
           1
   6.23506755E-001      0.0000
   3.90209773E-003     67.3278
   3.13572360E-003    249.2524
   8.37940447E-002    235.7702
   9.98960506E-002    247.7906
   5.39948435E-004    116.1835
   9.24451061E-004    331.5015
   1.10660759E-001    209.5442
   4.98798137E-001    222.1969
   2.90304162E-003    295.1856
   1.75485788E-001    249.6303
   1.22399709E-004     78.2958
   4.47184802E-005      3.9159
   4.71096089E-005    217.5453
   4.76668856E-005    230.1727
   2.65952175E-005      7.7592
   3.19476857E-005    220.5863
   3.20181990E-005     19.8499
   2.30957564E-005    215.5225
   1.76210579E-005    208.0938
   3.15420513E-002    243.6697
   4.27675661E-002    264.4085
   1.59801427E-002    213.2766
           2
   6.23343783E-001      0.0000
   4.50111197E-003     70.9706
   3.48568835E-003    248.3217
   8.41947911E-002    235.8382
   1.00458102E-001    247.8510
   5.22026318E-004    114.3054
   9.58132293E-004    330.8668
   1.11694362E-001    209.6204
   5.03635161E-001    222.3826
   2.91494465E-003    295.3532
   1.77335622E-001    249.7514
   1.34664929E-004     72.1287
   3.93887820E-004    265.6120
   1.26428573E-003    271.6857
   1.51736088E-003    303.4962
   5.42894827E-005     39.7217
   6.20443657E-005    111.6116
   2.64695332E-005     78.6074
   2.54215903E-005    232.2932
   3.17741424E-005    215.6186
   3.17422236E-002    243.7023
   4.33586826E-002    264.5548
   1.60530629E-002    213.4864
           3
   6.23462528E-001      0.0000
   4.15408450E-003     67.8958
   3.23836434E-003    248.1032
   8.44848306E-002    235.6671
   1.00928596E-001    247.6903
   5.21133240E-004    114.7142
   1.02198320E-003    331.9971
   1.12717153E-001    209.4443
   5.08714975E-001    222.3855
   2.94614257E-003    296.4891
   1.79405822E-001    249.5435
   1.58148883E-004     70.9147
   9.12367107E-004    255.1977
   2.97028168E-003    265.7956
   3.71271243E-003    299.7102
   1.35005130E-004     57.4882
   2.50259615E-004     97.4495
   9.21863407E-005    122.6333
   2.03245567E-005    218.0588
   5.58534060E-005    217.1492
   3.17789400E-002    243.4529
   4.36797953E-002    264.4543
   1.61372577E-002    213.3446
           4
   6.23353227E-001      0.0000
   4.50651787E-003     70.2248
   3.46904598E-003    247.6905
   8.48191717E-002    235.5985
   1.01433557E-001    247.6249
   5.08950120E-004    113.0635
   1.08989292E-003    331.6190
   1.13795258E-001    209.3772
   5.14135357E-001    222.4815
   2.97856562E-003    297.3081
   1.81467689E-001    249.4757
   1.81306254E-004     66.0492
   1.45603922E-003    253.2562
   4.69721712E-003    264.4313
   5.92357301E-003    298.0481
   2.11129493E-004     55.9216
   4.13248138E-004     90.2045
   1.54875213E-004    121.9304
   1.65945373E-005    230.5382
   7.87777309E-005    215.8962
   3.18756838E-002    243.3636
   4.41807388E-002    264.5494
   1.62070412E-002    213.3655
           5
   6.23442246E-001      0.0000
   4.23824990E-003     68.1869
   3.26390022E-003    247.7861
   8.50616877E-002    235.3844
   1.01873861E-001    247.3769
   5.11169499E-004    111.6853
   1.14660668E-003    334.1064
   1.15039482E-001    209.1047
   5.20647194E-001    222.4449
   2.98362256E-003    298.1522
   1.84033692E-001    249.1861
   2.21018887E-004     67.5600
   2.22724833E-003    250.5909
   7.22846907E-003    262.6860
   9.17801973E-003    297.0336
   3.35571561E-004     58.3435
   6.90737236E-004     88.7397
   2.58358349E-004    125.9564
   5.67858582E-006    106.1092
   1.11052501E-004    216.3046
   3.18891756E-002    243.0861
   4.46614268E-002    264.3456
   1.62884272E-002    213.1155
           6
   6.23350659E-001      0.0000
   4.48605874E-003     70.2029
   3.41682394E-003    248.0800
   8.54434972E-002    235.1952
   1.02517635E-001    247.2419
   4.80353607E-004    110.1139
   1.28374854E-003    333.2224
   1.16453781E-001    208.9890
   5.27754670E-001    222.5324
   3.03170929E-003    299.8265
   1.86793449E-001    249.0546
   2.65252409E-004     62.6164
   3.04172597E-003    250.0950
   9.86524180E-003    262.3900
   1.25991216E-002    296.6896
   4.58508506E-004     57.6805
   9.59349595E-004     86.6555
   3.65094238E-004    125.6835
   2.62826190E-005      4.9028
   1.44011838E-004    216.7807
   3.20356345E-002    242.9649
   4.52615479E-002    264.3687
   1.63854807E-002    213.0262
"""

@pytest.fixture
def mock_fort13_file(tmp_path):
    """Fixture to create a temporary fort.13 file."""
    # Create a temporary file path
    mock_file = tmp_path / "fort.13"
    # Write the sample content to the temporary file
    with open(mock_file, "w") as f:
        f.write(fort13_content)
    return mock_file


def test_read_fort13(mock_fort13_file):
    """Test the read_fort14 function with a mocked fort.14 file."""
    # Call the function with the mock file path
    mann, mann_indices, local_mann_indices, global_mann = read_fort13(mock_fort13_file)

    # Assert that the global manning's value is read correctly
    assert global_mann == 0.025

    # Assert that the mann_indices and local_mann_indices are read correctly
    expected_mann_indices = [7, 15]
    expected_local_mann_indices = [0, 1, 2, 5] # For Python index not ADCIRC

    assert mann_indices == expected_mann_indices
    assert local_mann_indices == expected_local_mann_indices

@pytest.fixture
def mock_fort14_file(tmp_path):
    """Fixture to create a temporary fort.14 file."""
    # Create a temporary file path
    mock_file = tmp_path / "fort.14"
    # Write the sample content to the temporary file
    with open(mock_file, "w") as f:
        f.write(fort14_content)
    return mock_file


def test_read_fort14(mock_fort14_file):
    """Test the read_fort14 function with a mocked fort.14 file."""
    # Call the function with the mock file path
    ADCIRC_nodes, nN, eN = read_fort14(mock_fort14_file)

    # Assert that the number of nodes and elements is correct
    assert nN == 6, f"Expected number of nodes to be 6, got {nN}"
    assert eN == 4, f"Expected number of elements to be 4, got {eN}"


@pytest.fixture
def mock_inundation63_file(tmp_path):
    """Fixture to create a temporary inundation.63 file."""
    # Create a temporary file path
    mock_file = tmp_path / "inundation.63"
    # Write the sample content to the temporary file
    with open(mock_file, "w") as f:
        f.write(inundation63_content)
    return mock_file


def test_read_inundationtime63(mock_inundation63_file):
    """Test the read_inundationtime63 function with a mocked inundation.63 file."""
    # Call the function with the mock file path
    inundationtime, nN, max_time = read_inundationtime63(mock_inundation63_file)

    # Assert that the number of nodes and max_time are read correctly
    assert nN == 6
    assert max_time == 3.8880000000E+006

    # Expected results for inundation time
    expected_inundationtime = np.array([
        (1, 1.0),
        (2, 0.0),
        (3, 1.0),
        (4, 0.5),
        (5, 1.0),
        (6, 0.0)
    ], dtype=[('nodeNum', int), ('time', float)])

    # Use numpy.isclose to compare the 'time' field within a tolerance
    assert np.isclose(inundationtime['time'], expected_inundationtime['time'], atol=1e-6).all()

    # Expected hydro class values
    expected_hydro_class = {
        1: int(2),
        2: int(0),
        3: int(2),
        4: int(1),
        5: int(2),
        6: int(0)
    }

    accuracy = 1.0 * 10 ** -6

    """ inunTBN = -99999.0: nodata,
                0: land, 1: intertidal, 2: subtidal """  # , 3: pond/lake"""

    # Create the masks
    inunT = inundationtime['time']
    mask_land = (0 - accuracy < inunT) & (inunT < 0 + accuracy)  # land mask
    mask_intertidal = (accuracy < inunT) & (inunT < 1 - accuracy)  # intertidal mask
    mask_water = (1 - accuracy < inunT) & (inunT < 1 + accuracy)  # water mask

    hydro_class = np.zeros([nN, 1], dtype=int)  # initialize the hydro_class
    hydro_class[mask_land] = int(0)  # fully dried (land) region
    hydro_class[mask_intertidal] = int(1)  # intertidal region
    hydro_class[mask_water] = int(2)  # water region

    # Flatten the hydro_class for easy comparison
    hydro_class = hydro_class.flatten()

    # Assert the values
    for node, expected_class in expected_hydro_class.items():  # node is ADCIRC index but hydro_class is using Python index
        assert hydro_class[
                   node - 1] == expected_class, f"Node {node} should be in class {expected_class}, but got {hydro_class[node - 1]}"


@pytest.fixture
def mock_fort53_file(tmp_path):
    """Fixture to create a temporary fort.53 file."""
    # Create a temporary file path
    mock_file = tmp_path / "fort.53"
    # Write the sample content to the temporary file
    with open(mock_file, "w") as f:
        f.write(fort53_content)
    return mock_file


def test_read_fort53(mock_fort53_file):
    """Test the read_fort53 function with a mocked fort.53 file."""
    # Call the function with the mock file path
    Harmonics_nodes, nN, numHarm, tidal_frequencies, tidal_constituents = read_fort53(mock_fort53_file)

    # Assert the number of harmonics and nodes are read correctly (the numHarm may be differed by User input of fort.15)
    assert numHarm == 23
    assert nN == 6

    # Assert that the tidal frequencies and constituents are read correctly  (harmonic_Freq.txt)
    expected_tidal_frequencies = [
        0.0, 0.2639203022E-05, 0.4925201824E-05, 0.6759774415E-04, 0.7292115836E-04,
        0.1329544977E-03, 0.1355937007E-03, 0.1378796995E-03, 0.1405189025E-03,
        0.1431581055E-03, 0.1454441043E-03, 0.1503693062E-03, 0.2783986020E-03,
        0.2810378050E-03, 0.2859630068E-03, 0.4189175045E-03, 0.4215567075E-03,
        0.4238427063E-03, 0.5620756100E-03, 0.7025945125E-03, 0.7252294600E-04,
        0.1458423172E-03, 0.6495854113E-04
    ]

    expected_tidal_constituents = [
        "STEADY", "MN", "SM", "O1", "K1", "MNS2", "2MS2", "N2", "M2", "2MN2",
        "S2", "2SM2", "MN4", "M4", "MS4", "2MN6", "M6", "MSN6", "M8", "M10",
        "P1", "K2", "Q1"
    ]

    assert tidal_frequencies == expected_tidal_frequencies
    assert tidal_constituents == expected_tidal_constituents

########################################################################################################################
### Tidal Calculation
########################################################################################################################
# Fixture to generate test data
@pytest.fixture
def mock_water_level_data():
    """Fixture to create mock water level data using a sine wave."""
    t = np.arange(0, 24 * 3 * 3600, 1800)  # 3 days simulation with 30 min intervals
    wl = np.sin(t / (24 * 3600) * 2 * np.pi)  # Full cycle (2*pi) for a sine wave
    return wl

# Test functions using the fixture
def test_mean_water_levels(mock_water_level_data):

    """Test the mean_high_water function with mocked water level data."""
    expected_MHW = 1  # Expected max value for a sine wave
    calculated_MHW = mean_high_water(mock_water_level_data)
    assert np.isclose(calculated_MHW, expected_MHW), f"Expected {expected_MHW}, but got {calculated_MHW}"

    """Test the mean_low_water function with mocked water level data."""
    expected_MLW = -1  # Expected min value for a sine wave
    calculated_MLW = mean_low_water(mock_water_level_data)
    assert np.isclose(calculated_MLW, expected_MLW), f"Expected {expected_MLW}, but got {calculated_MLW}"

    """Test the mean_sea_water level with mocked water level data."""
    expected_MSL = 0.0 # Expected average value for a sine wave
    calculated_MSL = np.average(mock_water_level_data)
    assert np.isclose(calculated_MSL,expected_MSL), f"Expected {expected_MSL}, but got {calculated_MSL}"


@pytest.fixture
def mock_data():
    """Fixture to create a small dataset for testing."""
    ndv = -99999.0  # No data value (ndv) using ADCIRC convention
    data = {
        'node_id': [1, 2, 3, 4, 5, 6],
        'x_prj': [0.0, 0.5, 2.0, 0.0, 1.0, 2.0],
        'y_prj': [0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
        'msl': [2.27, ndv, 3.11, ndv, 0.50, ndv],
        'mlw': [2.00, ndv, 3.00, ndv, 1.00, ndv],
        'mhw': [2.50, ndv, 3.50, ndv, 1.50, ndv],
        'HydroClass': [2, 0, 2, 1, 2, 0]

    }
    df = pd.DataFrame(data)
    return df


def test_idw_interpolation(mock_data):
    """Test IDW interpolation function with mock data."""

    # Prepare the data
    tidal_prj = mock_data
    HC = tidal_prj['HydroClass'].values  # Sample HC data

    # The indices of interpolated areas (not fully wetted areas)
    indices = np.where(HC < 2)
    # The indices of fully wetted areas (base points)
    inverse_indices = np.where(~(HC < 2))
    print(
        'Nodes of interpolation:\t', len(
            indices[0]), ', Nodes of references:\t', len(
            inverse_indices[0]))

    # Use a projected system
    scale_factor = 1  # scale factor for KDTree unit:
    xy_base = tidal_prj[['x_prj', 'y_prj']
                        ].iloc[inverse_indices].to_numpy() / scale_factor
    xy_interp = tidal_prj[['x_prj', 'y_prj']
                          ].iloc[indices].to_numpy() / scale_factor


    leafsize = 1  # leafsize for KDTree
    power = 2  # the power of the inverse distance
    numNeighbors = 2  # Number of neighbors to use in interpolation

    # Perform interpolation for 'mlw'
    z = tidal_prj['mlw'].values
    invdisttree = Invdisttree(xy_base, z[inverse_indices], leafsize=leafsize, stat=0)
    interpol, weight_factors = invdisttree(xy_interp, nnear=numNeighbors, eps=0.0, p=power)

    # Assert that interpolated values are close to expected
    expected_interpol = [1.833, 1.50, 2.00]
    assert np.allclose(interpol, expected_interpol, atol=0.01)

    print(f"Interpolated values for MLW: {interpol}")

########################################################################################################################
### Vegetation Calculation
########################################################################################################################


# Example usage
if __name__ == "__main__":
    pytest.main(["-v", __file__])