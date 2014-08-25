"""
Test cases for the vecspace package.

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

from pod import vecspace

class VecSpaceTest(unittest.TestCase):

  def test_projection_trivial1(self):
    space = vecspace.VectorSpace3D()
    p1 = space.define_point(0.0, 0.0, 1.0)
    p2 = space.define_point(0.0, 1.0, 0.0)
    p3 = space.define_point(0.0, 1.0, 1.0)
    plane = space.define_plane(p1, p2, p3)
    p = space.define_point(1.0, 1.0, 1.0)
    proj = plane.project(p)
    expected = [ 0.0, 1.0, 1.0 ]
    self.assertEqual(proj.projected.get_data(), expected)
    
  def test_projection_trivial2(self):
    space = vecspace.VectorSpace3D()
    p1 = space.define_point(1.0, 0.0, 1.0)
    p2 = space.define_point(1.0, 1.0, 0.0)
    p3 = space.define_point(1.0, 1.0, 1.0)
    plane = space.define_plane(p1, p2, p3)
    p = space.define_point(0.0, 0.0, 1.0)
    proj = plane.project(p)
    expected = [ 1.0, 0.0, 1.0 ]
    self.assertEqual(proj.projected.get_data(), expected)

  def test_projection1(self):
    space = vecspace.VectorSpace3D()
    p1 = space.define_point(0.0, 0.0, 1.0)
    p2 = space.define_point(0.0, 1.0, 0.0)
    p3 = space.define_point(1.0, 0.0, 0.0)
    plane = space.define_plane(p1, p2, p3)
    p = space.define_point(0.0, 0.0, 0.0)
    proj = plane.project(p)
    expected = [ 0.333, 0.333, 0.333 ]
    for i in [0, 1, 2]:
      self.assertAlmostEqual(proj.projected.get_data()[i], expected[i], 3)
            
  def test_projection_idem(self):
    space = vecspace.VectorSpace3D()
    p1 = space.define_point(1.0, 0.0, 1.0)
    p2 = space.define_point(1.0, 1.0, 0.0)
    p3 = space.define_point(1.0, 1.0, 1.0)
    plane = space.define_plane(p1, p2, p3)
    p = space.define_point(1.0, 0.0, 1.0)
    proj = plane.project(p)
    expected = [ 1.0, 0.0, 1.0 ]
    self.assertEqual(proj.projected.get_data(), expected)
    
  def test_projection2(self):
    space = vecspace.VectorSpace3D()
    p1 = space.define_point(0.0, 0.0, -1.6)
    p2 = space.define_point(-4.0, 0.0, 0.0)
    p3 = space.define_point(0.0, 4.0, 0.0)
    plane = space.define_plane(p1, p2, p3)
    p = space.define_point(4.0, -4.0, 3.0)
    proj = plane.project(p)
    self.assertAlmostEqual(proj.projector.norm(), 39.0 / math.sqrt(33.0), 8)
    
  def test_projection3(self):
    space = vecspace.VectorSpace3D()
    p1 = space.define_point(0.0, 0.0, 5.0)
    p2 = space.define_point(-2.5, 0.0, 0.0)
    p3 = space.define_point(0.0, 5.0/3.0, 0.0)
    plane = space.define_plane(p1, p2, p3)
    p = space.define_point(2.0, 10.0, -7.0)
    proj = plane.project(p)
    expected = [ 4.0, 7.0, -8.0 ]
    for i in [0, 1, 2]:
      self.assertAlmostEqual(proj.projected.get_data()[i], expected[i], 3)
      
  def test_projection_line1(self):
    """
    Projecting on a 3D line contained in a subspace
    """
    space = vecspace.VectorSpace3D()
    p1 = space.define_point(0.0, 0.0, 5.0)
    p2 = space.define_point(-2.5, 0.0, 0.0)
    line = space.define_line(p1, p2)
    p = space.define_point(2.0, 10.0, -7.0)
    proj = line.project(p)
    expected = [ -4.4, 0.0, -3.8 ]
    for i in [0, 1, 2]:
      self.assertAlmostEqual(proj.projected.get_data()[i], expected[i], 3)
    
  def test_projection_line2(self):
    """
    Projecting on a 3D line
    """
    space = vecspace.VectorSpace3D()
    p1 = space.define_point(0.0, 0.1, 5.0)
    p2 = space.define_point(-2.5, 0.0, 0.0)
    line = space.define_line(p1, p2)
    p = space.define_point(2.0, 10.0, -7.0)
    proj = line.project(p)
    expected = [ -4.319, -0.073, -3.639 ]
    for i in [0, 1, 2]:
      self.assertAlmostEqual(proj.projected.get_data()[i], expected[i], 3)
    
  def test_projection_001(self):
    """
    Projecting on a 2D line
    """
    space = vecspace.VectorSpace(2)
    p1 = space.define_point(0.01, 0.02)
    p2 = space.define_point(0.0, 0.0)
    p = space.define_point(0.02, 0.03)
    p_expected = space.define_point(0.08/5.0, 0.16/5.0)
    line = space.define_line(p1, p2)
    proj = line.project(p)
    self.assertEqual(proj.projected, p_expected)
    
if __name__ == '__main__':
    unittest.main()
