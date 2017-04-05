from shifty import mc_rms
from shifty import snr_to_unc


# This test is going fail
def test_mc_rms_1():
    assert mc_rms(1, 1, 1, 1, 1) == 1


def test_snr_to_unc():
    assert snr_to_unc(1) == 1
