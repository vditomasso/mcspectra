from shifty import mc_rms
import pytest

# This test is going fail
def test_mc_rms_1():
	assert mc_rms(1,1,1,1,1) == 1