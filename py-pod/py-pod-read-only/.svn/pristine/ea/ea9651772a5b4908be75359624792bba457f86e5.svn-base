"""
Test cases for the linalg package.

@author: Christophe Alexandre <ch.alexandre at bluewin dot ch>

@license:
Copyright(C) 2010 Christophe Alexandre

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/lgpl.txt>.
"""

import unittest
import math
import logging

from pod import linalg

class VecSpaceTest(unittest.TestCase):
  
  def test_vector(self):
    v = linalg.Vector(4).set_values(1.0, -2.0, -1.0)
    self.assertEqual(v.get_component(1), -2.0)
    self.assertEqual(v.get_component(3), 0.0)
    try:
      x = v.get_component(4)
    except AssertionError:
      return
    self.fail()
    
  def test_matrix(self):
    m = linalg.Matrix(3, 4)
    m.set_value(0, 0, 1.0)
    m.set_value(1, 1, 1.0)
    m.set_value(2, 2, 1.0)
    m.set_value(0, 2, -1.0)
    self.assertEqual(m.get_dim_row(), 3)
    self.assertEqual(m.get_dim_col(), 4)
    m = m.transpose()
    m.set_value(3, 1, -0.6)
    self.assertEqual(m.get_dim_row(), 4)
    self.assertEqual(m.get_dim_col(), 3)
    self.assertEqual(m.get_value(3, 1), -0.6)
    
  def test_product_matrix(self):
    m1 = linalg.Matrix(2, 3)
    m1.set_value(0, 0, 0.0)
    m1.set_value(1, 0, 4.0)
    m1.set_value(0, 1, -1.0)
    m1.set_value(1, 1, 11.0)
    m1.set_value(1, 2, 2.0)
    m1.set_value(0, 2, 2.0)
    m2 = linalg.Matrix(3, 2)
    m2.set_value(0, 0, 3.0)
    m2.set_value(1, 0, 1.0)
    m2.set_value(2, 0, 6.0)
    m2.set_value(0, 1, -1.0)
    m2.set_value(1, 1, 2.0)
    m2.set_value(2, 1, 1.0)
    m = linalg.prod_matrix(m1, m2)
    self.assertEqual(m.get_value(0, 0), 11.0)
    self.assertEqual(m.get_value(1, 0), 35.0)
    self.assertEqual(m.get_value(0, 1), 0.0)
    self.assertEqual(m.get_value(1, 1), 20.0)
    
  def test_minor(self):
    m1 = linalg.Matrix(2, 3)
    m1.set_value(0, 0, 0.0)
    m1.set_value(1, 0, 4.0)
    m1.set_value(0, 1, -1.0)
    m1.set_value(1, 1, 11.0)
    m1.set_value(1, 2, 2.0)
    m1.set_value(0, 2, 2.0)
    self.assertEqual(m1.minor(1, 1).get_value(0, 1), 2.0)
    self.assertEqual(m1.minor(0, 2).get_value(0, 1), 11.0)
    
  def test_determinant1(self):
    m = linalg.Matrix(3)
    m.set_value(0, 0, 3.0)
    m.set_value(0, 1, 1.0)
    m.set_value(0, 2, 8.0)
    m.set_value(1, 0, 2.0)
    m.set_value(1, 1, -5.0)
    m.set_value(1, 2, 4.0)
    m.set_value(2, 0, -1.0)
    m.set_value(2, 1, 6.0)
    m.set_value(2, 2, -2.0)
    self.assertEqual(m.determinant(), 14)
    
  def test_determinant2(self):
    m = linalg.Matrix(2)
    m.set_value(0, 0, 3.0)
    m.set_value(0, 1, 1.0)
    m.set_value(1, 0, 2.0)
    m.set_value(1, 1, -5.0)
    self.assertEqual(m.determinant(), -17)
    
  def test_invert_matrix1(self):
    m = linalg.Matrix(4)
    t = [
      [1.0, 2.0, 4.0, 5.0],
      [-1.0, 0.0, 3.0, 4.0],
      [3.0, 4.0, 4.0, 4.0],
      [-2.0, -1.0, 4.0, 5.0]
    ]
    expected = [
      [2.6667e+00, -4.0000e+00, -1.0000e+00, 1.3333e+00],
      [-2.3333e+00, 4.0000e+00, 1.0000e+00, -1.6667e+00],
      [-1.3333e+00, -1.0000e+00, 1.0000e+00, 1.3333e+00],
      [1.6667e+00, 8.8832e-16, -1.0000e+00, -6.6667e-01]
    ]
    m.set_table(t)
    inv = m.inverse()
    for i in [0, 1, 2]:
      for j in [0, 1, 2]:
        self.assertAlmostEqual(inv.get_table()[i][j], expected[i][j], 4)
    
  def test_invert_matrix2(self):
    m = linalg.Matrix(2)
    m.set_value(0, 0, 2.0)
    m.set_value(0, 1, -3.0)
    m.set_value(1, 0, 4.0)
    m.set_value(1, 1, -7.0)
    expected = [
      [3.5, -1.5],
      [2.0, -1.0]
    ]
    inv = m.inverse()
    self.assertEqual(inv.get_table(), expected)
    
  def test_invert_matrix3(self):
    m = linalg.Matrix(3)
    t = [
      [1, 2, 3],
      [0, 4, 5],
      [1, 0, 6]
    ]
    expected = [
      [12.0/11.0, -6.0/11.0, -1.0/11.0],
      [5.0/22.0, 3.0/22.0, -5.0/22.0],
      [-2/11.0, 1.0/11.0, 2.0/11.0]
    ]
    m.set_table(t)
    inv = m.inverse()
    for i in [0, 1, 2]:
      for j in [0, 1, 2]:
        self.assertAlmostEqual(inv.get_table()[i][j], expected[i][j], 6)
    
if __name__ == '__main__':
    unittest.main()
