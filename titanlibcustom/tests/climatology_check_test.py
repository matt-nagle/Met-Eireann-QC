from __future__ import print_function
import unittest
import titanlibcustom


"""Convenient vectors used as inputs"""
len3 = [1, 2, 3]
len2 = [1, 2]
len1 = [1]
points3 = titanlibcustom.Points(len3, len3)
points2 = titanlibcustom.Points(len2, len2)
points1 = titanlibcustom.Points(len1, len1)
ut = 1546300800


class ClimatologyCheckTest(unittest.TestCase):
    def test_valid_input(self):
        """Check that the test doesn't fail"""
        # TODO: Reenable these
        #self.assertEqual(titanlibcustom.range_check_climatology(len1, len1, len1, len1, ut, len1, len1)[0],0)
        #self.assertEqual(titanlibcustom.range_check_climatology(len3, len3, len3, len3, ut, len1, len1)[0],0)
        #self.assertEqual(titanlibcustom.range_check_climatology(len3, len3, len3, len3, ut, len1, len3)[0],0)
        #self.assertEqual(titanlibcustom.range_check_climatology(len3, len3, len3, len3, ut, len1, len1)[0],0)
        pass

    def test_invalid_input(self):
        """Check that the test fails"""
        # TODO: Reenable these
        with self.assertRaises(RuntimeError):
            # Min/max inconsistent with inputs
            flags = titanlibcustom.range_check_climatology(points3, len3, ut, len2, len2)


if __name__ == '__main__':
    unittest.main()
