import numpy as np
from ..forces import force

from numba import jitclass          # import the decorator
from numba import int32, float32, int64, float64    # import the types

spec = [
     ('nbodies', int64),
     ('child', int32[:]),
     ('bmin', float64[:]),
     ('bmax', float64[:]),
     ('center', float64[:]),
     ('box_size', float64[:]),
     ('ncell', int64),
     ('cell_center', float64[:, :]),
     ('cell_radius', float64[:, :]),
     ('center_of_mass', float64[:, :]),
     ('mass', float64[:])]

@jitclass(spec)
class quadArray:
    def __init__(self, bmin, bmax, size):
        self.nbodies = int(size)
        self.child = -np.ones(4*(2*size+1), dtype=np.int32)
        self.bmin = bmin
        self.bmax = bmax
        self.center = .5 * (self.bmin + self.bmax)
        self.box_size = (self.bmax - self.bmin)
        self.ncell = 0
        self.cell_center = np.zeros((2 * size + 1, 2))
        self.cell_radius = np.zeros((2 * size + 1, 2))
        self.cell_center[0] = self.center
        self.cell_radius[0] = self.box_size
        self.mass = np.zeros(1)
        self.center_of_mass = np.zeros((1, 2))

    def buildTree(self, particles):
        for ip in range(len(particles)):
            x = particles[ip, 0]
            y = particles[ip, 1]
            center = self.center.copy()
            box_size = self.box_size.copy()

            cell = 0

            childPath = 0
            if x > center[0]:
                childPath += 1
            if y > center[1]:
                childPath += 2

            childIndex = self.nbodies + childPath

            while (self.child[childIndex] > self.nbodies):
                cell = self.child[childIndex] - self.nbodies
                center[:] = self.cell_center[cell]
                childPath = 0
                if x > center[0]:
                    childPath += 1
                if y > center[1]:
                    childPath += 2
                childIndex = self.nbodies + 4 * cell + childPath
            # no particle on this cell, just add it
            if (self.child[childIndex] == -1):
                self.child[childIndex] = ip
                self.child[ip] = cell
            # this cell already has a particle
            # subdivide and set the two particles
            elif (self.child[childIndex] < self.nbodies):
                npart = self.child[childIndex]

                oldchildPath = newchildPath = childPath
                while (oldchildPath == newchildPath):
                    self.ncell += 1
                    self.child[childIndex] = self.nbodies + self.ncell
                    center[:] = self.cell_center[cell]
                    box_size = .5 * self.cell_radius[cell]
                    if (oldchildPath & 1):
                        center[0] += box_size[0]
                    else:
                        center[0] -= box_size[0]
                    if ((oldchildPath >> 1) & 1):
                        center[1] += box_size[1]
                    else:
                        center[1] -= box_size[1]

                    oldchildPath = 0
                    if particles[npart][0] > center[0]:
                        oldchildPath += 1
                    if particles[npart][1] > center[1]:
                        oldchildPath += 2

                    newchildPath = 0
                    if x > center[0]:
                        newchildPath += 1
                    if y > center[1]:
                        newchildPath += 2

                    cell = self.ncell

                    self.cell_center[self.ncell] = center
                    self.cell_radius[self.ncell] = box_size

                    childIndex = self.nbodies + 4 * self.ncell + oldchildPath

                self.child[childIndex] = npart
                self.child[npart] = self.ncell

                childIndex = self.nbodies + 4 * self.ncell + newchildPath
                self.child[childIndex] = ip
                self.child[ip] = self.ncell

    def computeMassDistribution(self, particles, mass):
        self.mass = np.zeros(self.nbodies + self.ncell + 1)
        self.mass[:self.nbodies] = mass
        self.center_of_mass = np.zeros((self.nbodies + self.ncell + 1, 2))
        self.center_of_mass[:self.nbodies] = particles[:, :2]
        for i in range(self.ncell, -1, -1):
            sum_com_1 = 0.0
            sum_com_0 = 0.0
            sum_mass = 0.0
            for j in range(self.nbodies + 4 * i, self.nbodies + 4 * i + 4):
                    element = self.child[j]
                    if element >= 0:
                        sum_mass += self.mass[element]
                        sum_com_0 += self.mass[element] * self.center_of_mass[element, 0]
                        sum_com_1 += self.mass[element] * self.center_of_mass[element, 1]
            self.mass[self.nbodies + i] = sum_mass
            self.center_of_mass[self.nbodies + i, 0] = sum_com_0 / sum_mass
            self.center_of_mass[self.nbodies + i, 1] = sum_com_1 / sum_mass

#             elements = self.child[self.nbodies + 4 * i:self.nbodies + 4 * i + 4]
#             # print('elements', i, elements, self.center_of_mass[elements[elements>=0]]*self.mass[elements[elements>=0]])
#             masked_elements = elements[elements >= 0]
#             self.mass[self.nbodies + i] = np.sum(self.mass[masked_elements])
#             self.center_of_mass[self.nbodies + i] = np.sum(self.center_of_mass[masked_elements] * np.atleast_2d(self.mass[masked_elements]), axis=0)
#             self.center_of_mass[self.nbodies + i] /= self.mass[self.nbodies + i]
        # print('mass', self.mass)
        # print('center_of_mass', self.center_of_mass)

    def computeForce(self, p):
        depth = 0
        localPos = np.zeros(2 * self.nbodies, dtype=np.int32)
        localNode = np.zeros(2 * self.nbodies, dtype=np.int32)
        localNode[0] = self.nbodies

#         pos = p[:2]
#         acc = np.zeros(2)
        x = p[0]
        y = p[1]
        acc0 = 0.0
        acc1 = 0.0

        while depth >= 0:
            while localPos[depth] < 4:
                child = self.child[localNode[depth] + localPos[depth]]
                # print('child 1', child, localNode[depth] + localPos[depth])
                localPos[depth] += 1
                if child >= 0:
                    if child < self.nbodies:
                        F = force((x, y), self.center_of_mass[child], self.mass[child])
                        acc0 += F[0]
                        acc1 += F[1]
                    else:
                        dx = self.center_of_mass[child, 0] - x
                        dy = self.center_of_mass[child, 1] - y
                        dist = np.sqrt(dx ** 2 + dy ** 2)
                        if dist != 0 and self.cell_radius[child - self.nbodies][0] / dist < .5:
                            F = force((x, y), self.center_of_mass[child], self.mass[child])
                            acc0 += F[0]
                            acc1 += F[1]
                        else:
                            depth += 1
                            localNode[depth] = self.nbodies + 4 * (child - self.nbodies)
                            localPos[depth] = 0
            depth -= 1
        return (acc0, acc1)

    def __str__(self):
        indent = ' ' * 2
        s = 'Tree :\n'
        for i in range(self.ncell + 1):
            s += indent + 'cell {i}\n'.format(i=i)
            cellElements = self.child[self.nbodies + 4 * i:self.nbodies + 4 * i + 4]
            s += 2 * indent + 'box: {min} {max} \n'.format(min=self.cell_center[i] - self.cell_radius[i], max=self.cell_center[i] + self.cell_radius[i])
            s += 2 * indent + 'particules: {p}\n'.format(p=cellElements[np.logical_and(0 <= cellElements, cellElements < self.nbodies)])
            s += 2 * indent + 'cells: {c}\n'.format(c=cellElements[cellElements >= self.nbodies] - self.nbodies)
        return s
