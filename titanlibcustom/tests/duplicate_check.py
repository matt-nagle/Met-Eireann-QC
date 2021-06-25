from __future__ import print_function
import unittest
import titanlibcustom
import numpy as np
import time


class Test(unittest.TestCase):
    def test_single(self):
        points = titanlibcustom.Points([0], [0])
        flags = titanlibcustom.duplicate_check(points, 100)
        np.testing.assert_array_almost_equal(flags, [0])

    def test_multiple_identical(self):
        points = titanlibcustom.Points([0, 0, 1, 1], [0, 1, 1, 1])
        flags = titanlibcustom.duplicate_check(points, 10000)
        np.testing.assert_array_almost_equal(flags, [0, 0, 0, 1])
        flags = titanlibcustom.duplicate_check(points, 0)
        np.testing.assert_array_almost_equal(flags, [0, 0, 0, 1])

    def test_multiple_radius(self):
        points = titanlibcustom.Points([0, 0, 1, 1.1], [0, 1, 1, 1])
        flags = titanlibcustom.duplicate_check(points, 12000)
        np.testing.assert_array_almost_equal(flags, [0, 0, 0, 1])
        flags = titanlibcustom.duplicate_check(points, 10000)
        np.testing.assert_array_almost_equal(flags, [0, 0, 0, 0])


if __name__ == '__main__':
    unittest.main()
