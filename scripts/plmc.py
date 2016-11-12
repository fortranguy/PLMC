import abc
import numpy as np

class PeriodicBox(metaclass=abc.ABCMeta):
    size = np.zeros(3)
    def __init__(self, size):
        self.size = size
    @abc.abstractmethod
    def folded(self, position):
        pass
    def vector(self, position_1, position_2):
        return self.folded(position_2 - position_1)


class XYZperiodicBox(PeriodicBox):
    def folded(self, position):
        folded = np.mod(position, self.size)
        return np.where(folded > 0.5*self.size, folded - self.size, folded)

class XYperiodicBox(PeriodicBox):
    def folded(self, position):
        folded = np.hstack([np.mod(position[0:2], self.size[0:2]), position[2]])
        folded[0:2] = np.where(folded[0:2] > 0.5*self.size[0:2], folded[0:2] - self.size[0:2], \
            folded[0:2])
        return folded

def newBox(i_box, generatingData):
    periodicity = generatingData["Environment"]["Boxes"]["periodicity"]
    boxSize = np.array(generatingData["Environment"]["Boxes"]["initial size"][i_box])
    if periodicity == "XYZ":
        return XYZperiodicBox(boxSize)
    elif periodicity == "XY":
        return XYperiodicBox(boxSize)
    else:
        raise NameError("Box periodicity unknown.")


class Component:
    num = 0
    positions = np.zeros([0, 3])
    orientations = np.zeros([0, 3])
    isDipolar = False
    def __init__(self, num, positions, orientations, isDipolar):
        self.num = num
        self.positions = positions
        self.orientations = orientations
        self.isDipolar = isDipolar

def newComponents(i_box, generatingData):
    numComponents = generatingData["Mixture"]["number of components"]
    if numComponents == 0:
        raise NameError("There is no components.")
    components = []
    for i_component in range(1, numComponents+1):
        initialNumber = generatingData["Mixture"]["Component %i"%i_component]["initial number"][i_box]
        isDipolar = generatingData["Mixture"]["Component %i"%i_component]["is dipolar"]
        components.append(Component(initialNumber, np.array([0, 3]), np.array([0, 3]), isDipolar))
    return components

def dipolarEnergyIsNegative(vector_ij, orientation_i, orientation_j):
    distance_ij = np.linalg.norm(vector_ij)
    return np.dot(orientation_i, orientation_j) / distance_ij**3 - \
        3. * np.dot(orientation_i, vector_ij)*np.dot(orientation_j, vector_ij) / distance_ij**5 < 0