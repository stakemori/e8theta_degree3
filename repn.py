# -*- coding: utf-8 -*-
from abc import ABCMeta, abstractmethod
from sage.modules.free_module_element import vector


class ReplSpaceElement(object):
    __metaclass__ = ABCMeta

    @property
    def vector(self):
        return self._v

    def __init__(self, v):
        self._v = vector((a for a in v))

    def __add__(self, other):
        if other == 0:
            return self
        elif isinstance(other, ReplSpaceElement):
            return ReplSpaceElement(self.vector + other.vector)
        else:
            raise NotImplementedError

    def __neg__(self):
        return ReplSpaceElement(-self.vector)

    def __radd__(self, other):
        return self.__add__(other)

    def __mul__(self, other):
        if other in self.vector.base_ring():
            return ReplSpaceElement(self * other)
        else:
            raise ValueError

    @abstractmethod
    def right_action(self, g):
        r'''
        :param g: an element of GL_n
        :returns g \dot self
        :rtype ReplSpaceElement
        '''
        pass
