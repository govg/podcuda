"""
Basic linear algebra components and functions.

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

import math
import os
import logging

from util import NullHandler
from util import numbering
from util import prod_scalar

_h = NullHandler()
_logger = logging.getLogger('linalg')
_logger.addHandler(_h)

class Matrix(object):
  
  def __init__(self, dim_row, dim_col=None):
    if dim_col is None:
      dim_col = dim_row
    self.dim_row = dim_row
    self.dim_col = dim_col
    self._vectors = dict()
    
  def transpose(self):
    m = Matrix(self.get_dim_col(), self.get_dim_row())
    for i in range(self.get_dim_row()):
      for j in range(self.get_dim_col()):
        m.set_value(j, i, self.get_value(i, j))
    return m
    
  def get_value(self, i, j):
    assert i < self.get_dim_row(), 'row %d exceeding dimension %d' % (i, self.get_dim_row())
    assert j < self.get_dim_col(), 'column %d exceeding dimension %d' % (j, self.get_dim_col())
    if self._vectors.has_key(i):
      v = self._vectors[i]
    else:
      v = self._create_vector()
    return v.get_component(j)
  
  def set_table(self, table):
    assert len(table) == self.get_dim_row(), 'expected %d rows instead of %d'% (self.get_dim_row(), len(table))
    for i in range(self.get_dim_row()):
      assert len(table[i]) == self.get_dim_col(), 'expected %d columns instead of %d'% (self.get_dim_col(), len(table[i]))
      for j in range(self.get_dim_col()):
        self.set_value(i, j, table[i][j])
        
  def get_table(self):
    table = list()
    for i in range(self.get_dim_row()):
      row = list()
      for j in range(self.get_dim_col()):
        row.append(self.get_value(i, j))
      table.append(row)
    return table
  
  def get_dimension(self):
     return (self.dim_row, self.dim_col)
  
  def get_dim_row(self):
    return self.get_dimension()[0]
    
  def get_dim_col(self):
    return self.get_dimension()[1]
  
  def _create_vector(self):
    return Vector(self.get_dimension()[1])
    
  def set_value(self, i, j, val):
    msg = str(self.get_dimension())
    assert i < self.get_dim_row(), 'row %d exceeding dimension %s' % (i, msg)
    assert j < self.get_dim_col(), 'column %d exceeding dimension %s' % (j, msg)
    if self._vectors.has_key(i):
      v = self._vectors[i]
    else:
      v = self._create_vector()
      self._vectors[i] = v
    try:
      v.set_component(j, float(val))
    except TypeError, e:
      logging.error('not a number: %s' % str(val))
      raise e

  def minor(self, r, c):
    """
    Sub-matrix excluding the specified row and column.
    """
    m = Matrix(self.get_dim_row() - 1, self.get_dim_col() - 1)
    dr = self.get_dim_row()
    dc = self.get_dim_col()
    for k, i in numbering(range(0, r) + range(r+1, dr)):
      for l, j in numbering(range(0, c) + range(c+1, dc)):
        v = self.get_value(i, j)
        m.set_value(k, l, v)
    return m

  def determinant(self):
    size = self.get_dim_row()
    det = 0.0
    # stopping condition
    if size == 1:
      det = self.get_value(0, 0)
    else:
      for i in range(size):
        minor = self.minor(0, i)
        if i % 2 == 0:
          sign = 1.0
        else:
          sign = -1.0
        det += sign * self.get_value(0, i) * minor.determinant()
    return det

  def cofactors(self):
    m = Matrix(self.get_dim_row())
    for row in range(self.get_dim_row()):
      for col in range(self.get_dim_col()):
        if row % 2 == col % 2:
          sign = 1.0
        else:
          sign = -1.0
        v = sign * self.minor(row, col).determinant()
        m.set_value(row, col, v)
    return m
  
  def multiply(self, m):
    return prod_matrix(self, m)
    
  def sub(self, m):
    msg = 'matrices with different dimensions (%s - %s)' % (str(self.get_dimension()), str(m.get_dimension()))
    assert m.get_dimension() == self.get_dimension(), msg
    n = Matrix(self.get_dim_row(), self.get_dim_col())
    for row in range(self.get_dim_row()):
      for col in range(self.get_dim_col()):
        v = self.get_value(row, col) - m.get_value(row, col)
        n.set_value(row, col, v)
    return n
    
  def add(self, m):
    assert m.get_dimension() == self.get_dimension(), 'incompatible dimensions'
    n = Matrix(self.get_dim_row(), self.get_dim_col())
    for row in range(self.get_dim_row()):
      for col in range(self.get_dim_col()):
        v = self.get_value(row, col) + m.get_value(row, col)
        n.set_value(row, col, v)
    return n
      
  def scale(self, factor):
    m = Matrix(self.get_dim_row(), self.get_dim_col())
    for row in range(self.get_dim_row()):
      for col in range(self.get_dim_col()):
        v = self.get_value(row, col) * factor
        m.set_value(row, col, v)
    return m
    
  def inverse(self):
    if self.get_dimension() == (1,1):
      m = Matrix(1)
      m.set_value(0, 0, 1.0 / self.get_value(0, 0))
      return m
    det_inv = 1.0 / self.determinant()
    cofactors = self.cofactors()
    cofactors = cofactors.transpose()
    i = cofactors.scale(det_inv)
    return i
  
  def pseudo_inverse(self):
    """
    Full row rank pseudo inverse.
    
    A+ = transp(A).inv(A.transp(A))
    """
    a = self
    t_a = a.transpose()
    a_t_a = a.multiply(t_a)
    inv_a_transp_a = a_t_a.inverse()
    return t_a.multiply(inv_a_transp_a)
    
  def copy(self):
    c = Matrix(self.get_dim_row(), self.get_dim_col())
    for row in range(self.get_dim_row()):
      for col in range(self.get_dim_col()):
        c.set_value(row, col, self.get_value(row, col))
    return c

  def __repr__(self):
    out = '(M%dx%d)' % (self.get_dim_row(), self.get_dim_col()) + os.linesep
    for row in range(self.get_dim_row()):
      line = []
      for col in range(self.get_dim_col()):
        line.append(self.get_value(row, col))
      out += ', '.join(map(str, line)) + os.linesep
    return out

class VectorMatrix(Matrix):
  """
  Makes a vector compatible with Matrix operations.
  """
  
  def __init__(self, vector):
    """
    Column vector.
    
    Use transpose() to turn it into a row vector.
    """
    Matrix.__init__(self, vector.get_length(), 1)
    for row in range(vector.get_length()):
      self.set_value(row, 0, vector.get_component(row))
      
class Vector(object):
  """
  @todo: should somehow be merged with VectorMatrix.
  """
  
  def __init__(self, dim):
    self._dim = dim
    self._values = dict()
  
  def get_length(self):
    return self._dim
  
  def set_values(self, *values):
    """ Using specified values to initialize the vector. """
    for n, v in numbering(values):
      self.set_component(n, v)
    return self
  
  def set_data(self, data):
    """ Using raw data (python list) to initialize the vector. """
    self.set_values(*data)
  
  def get_component(self, i):
    assert i < self.get_length(), 'index %d exceeding dimension %d' % (i, self.get_length())
    assert i >= 0, 'non positive index %d' % i
    if not self._values.has_key(i):
      return 0.0
    else:
      return self._values[i]
      
  def set_component(self, i, v):
    assert i < self.get_length(), 'index %d exceeding dimension %d' % (i, self.get_length())
    assert i >= 0, 'non positive index %d' % i
    self._values[i] = v

  def get_data(self):
    """ Raw data as built-in python list."""
    l = list()
    for i in range(self.get_length()):
      l.append(self.get_component(i))
    return l

  def sub(self, vector):
    result = Vector(self.get_length())
    for i in range(self.get_length()):
      result.set_component(i, self.get_component(i) - vector.get_component(i))
    return result

  def add(self, vector):
    result = Vector(self.get_length())
    for i in range(self.get_length()):
      result.set_component(i, self.get_component(i) + vector.get_component(i))
    return result
    
  def scale(self, a):
    result = Vector(self.get_length())
    for i in range(self.get_length()):
      result.set_component(i, a * self.get_component(i))
    return result

  def product(self, vector):
    return prod_scalar(vector.get_data(), self.get_data())

  def norm(self):
    return math.sqrt(self.product(self))

  def symmetric(self):
    null_vector = Vector(self.get_length())
    return null_vector.sub(self)
    
  def units(self, unit_vector):
    """
    How many times the current vector fits in the specified units (in norm).
    
    The sign has a meaning only if both vectors are colinear.
    """
    ratio = self.norm() / unit_vector.norm()
    if self.sub(unit_vector).norm() > unit_vector.norm():
      # vectors are in opposite directions
      ratio = -ratio
    return ratio
      
  def __eq__(self, other):
    return self._values == other._values
  
  def __ne__(self, other):
    return self._values != other._values
  
  def __repr__(self):
    line = []
    for col in range(self.get_length()):
      line.append(self.get_component(col))
    out = '(V)[' + ', '.join(map(str, line)) + ']'
    return out
    
  def __hash__(self):
    return hash(tuple(self.get_data()))

class Point(Vector):
  """
  Equivalent to a Vector
  """
  def __init__(self, *coordinates):
    Vector.__init__(self, len(coordinates))
    self.set_values(*coordinates)
  
  def __repr__(self):
    line = []
    for col in range(self.get_length()):
      line.append(self.get_component(col))
    out = '(P)[' + ', '.join(map(str, line)) + ']'
    return out

def identity(dimension):
  matrix = Matrix(dimension)
  for i in range(dimension):
    matrix.set_value(i, i, 1.0)
  return matrix
  
def zero(dimension):
  zeros = [0.0] * dimension
  return Point(*zeros)

def prod_matrix(m1, m2):
  msg = 'incompatible dimensions: %s x %s' % (str(m1.get_dimension()), str(m2.get_dimension()))
  assert m1.get_dim_col() == m2.get_dim_row(), msg
  result = Matrix(m1.get_dim_row(), m2.get_dim_col())
  for i in range(m1.get_dim_row()):
    for j in range(m2.get_dim_col()):
      v = 0.0
      for k in range(m1.get_dim_col()):
        v += m1.get_value(i, k) * m2.get_value(k, j)
      result.set_value(i, j, v)
  return result

