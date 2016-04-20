# -*- coding: utf-8 -*-
from abc import ABCMeta, abstractmethod, abstractproperty


class ReplSpaceElement(object):
    __metaclass__ = ABCMeta

    @abstractproperty
    def vector(self):
        pass
