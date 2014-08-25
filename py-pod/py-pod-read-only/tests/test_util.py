"""
Test cases for the util package.

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

from pod import util

class UtilTest(unittest.TestCase):
  
  def test_prod_scalar(self):
    v1 = [1.0, 2.0, 2.0]
    v2 = [1.0, -2.0, 1.0]
    self.assertEqual(util.prod_scalar(v1, v2), -1.0)
    
  def test_scalar_fail(self):
    v1 = [1.0, 2.0, 2.0]
    v2 = [1.0, -2.0, 1.0, 0.5]
    try:
      x = util.prod_scalar(v1, v2)
    except AssertionError:
      return
    self.fail()
    
  def test_norm(self):
    v = [1.0, 2.0, 2.0]
    self.assertEqual(util.norm(v), 3.0)
    
if __name__ == '__main__':
    unittest.main()
