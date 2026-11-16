import numpy as np
from TADA import run_tada

def test_tada_aligns_head_tail_equivalent_directors():
    positions = np.array([[0,0,0], [1,0,0], [2,0,0]], float)
    directors = np.array([[1,0,0], [-1,0,0], [1,0,0]], float)
    aligned, order = run_tada(positions, directors, 1.5, seed=0)
    assert np.allclose(aligned, [[1,0,0], [1,0,0], [1,0,0]])
    assert sorted(order) == [1, 2, 3]
