import numpy as np
cimport numpy as cnp
from .. cimport forces 


cdef class quadArray:
    cdef:
        public int nbodies, ncell
        public cnp.int32_t[::1] child
        public cnp.float64_t[::1] bmin, bmax, center, box_size, mass
        public cnp.float64_t[:, ::1] cell_center, cell_radius, center_of_mass

    def __cinit__(self, bmin, bmax, size):
        self.nbodies = size
        self.child = -np.ones(4*(2*size+1), dtype=np.int32)
        self.bmin = np.ascontiguousarray(bmin, dtype=np.float64)
        self.bmax = np.ascontiguousarray(bmax, dtype=np.float64)
        self.center = np.empty_like(self.bmin)
        self.box_size = np.empty_like(self.bmin)
        self.center[0] = .5*(bmin[0] + bmax[0])
        self.center[1] = .5*(bmin[1] + bmax[1])
        self.box_size[0] = (bmax[0] - bmin[0])
        self.box_size[1] = (bmax[1] - bmin[1])
        self.ncell = 0
        self.cell_center = np.zeros((2*size+1, 2))
        self.cell_radius = np.zeros((2*size+1, 2))
        self.cell_center[0, 0] = self.center[0]
        self.cell_center[0, 1] = self.center[1]
        self.cell_radius[0, 0] = self.box_size[0]
        self.cell_radius[0, 1] = self.box_size[1]

    def __dealloc__(self):
        self.bmin = None
        self.bmax = None
        self.cell_center = None
        self.cell_radius = None
        self.child = None
        self.mass = None
        self.center_of_mass = None

    def buildTree(self, particles):
        cdef:
            int childIndex, cell, ip, npart
            cnp.float64_t[::1] box_size, center

        for ip, p in enumerate(particles):
            center = self.center.copy()
            box_size = self.box_size.copy()
            x, y = p[:2]
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
                childIndex = self.nbodies + 4*cell + childPath
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
                    box_size = .5*np.asarray(self.cell_radius[cell])
                    if (oldchildPath&1):
                        center[0] += box_size[0]
                    else:
                        center[0] -= box_size[0]
                    if ((oldchildPath>>1)&1):
                        center[1] += box_size[1]
                    else:
                        center[1] -= box_size[1]

                    oldchildPath = 0
                    if particles[npart][0] > center[0]:
                        oldchildPath += 1
                    if particles[npart][1] > center[1]:
                        oldchildPath += 2

                    newchildPath = 0
                    if p[0] > center[0]:
                        newchildPath += 1
                    if p[1] > center[1]:
                        newchildPath += 2

                    cell = self.ncell

                    self.cell_center[self.ncell, 0] = center[0]
                    self.cell_center[self.ncell, 1] = center[1]
                    self.cell_radius[self.ncell, 0] = box_size[0]
                    self.cell_radius[self.ncell, 1] = box_size[1]

                    childIndex = self.nbodies + 4*self.ncell + oldchildPath

                self.child[childIndex] = npart
                self.child[npart] = self.ncell

                childIndex = self.nbodies + 4*self.ncell + newchildPath
                self.child[childIndex] = ip
                self.child[ip] = self.ncell

    def computeMassDistribution(self, cnp.float64_t[:, ::1] particles,
                                      cnp.float64_t[::1] mass):
        cdef:
            int i, j, element
            bint mask
            cnp.float64_t m, sum_mass, sum_com_0, sum_com_1
        i = self.nbodies + self.ncell + 1
        self.mass = np.zeros(i, dtype=np.float64)
        self.center_of_mass = np.zeros((i, 2), dtype=np.float64)
        for i in range(self.nbodies):
            self.mass[i] = mass[i]
            self.center_of_mass[i, 0] = particles[i, 0]
            self.center_of_mass[i, 1] = particles[i, 1]
        for i in range(self.ncell, -1, -1):
            sum_com_1 = sum_com_0 = sum_mass = 0.0
            for j in range(self.nbodies + 4*i, self.nbodies + 4*i + 4):
                    element = self.child[j]
                    if element >= 0:
                        sum_mass += self.mass[element]
                        sum_com_0 += self.mass[element] * self.center_of_mass[element, 0]
                        sum_com_1 += self.mass[element] * self.center_of_mass[element, 1]
            self.mass[self.nbodies + i] = sum_mass
            self.center_of_mass[self.nbodies + i, 0 ] = sum_com_0 / sum_mass
            self.center_of_mass[self.nbodies + i, 1 ] = sum_com_1 / sum_mass

    def computeForce(self, cnp.float64_t[::1] p):
        cdef:
            int depth, child
            cnp.int32_t[::1] localPos, localNode
            cnp.float64_t[::1] pos, acc, F
        depth = 0
        localPos = np.zeros(2*self.nbodies, dtype=np.int32)
        localNode = np.zeros(2*self.nbodies, dtype=np.int32)
        localNode[0] = self.nbodies

        pos = p[:2]
        acc = np.zeros(2)

        while depth >= 0:
            while localPos[depth] < 4:
                child = self.child[localNode[depth] + localPos[depth]]
                # print('child 1', child, localNode[depth] + localPos[depth])
                localPos[depth] += 1
                if child >= 0:
                    if child < self.nbodies:
                        F = forces.force(pos, self.center_of_mass[child], self.mass[child])
                        acc[0] = acc[0] + F[0]
                        acc[1] = acc[1] + F[1]
                    else:
                        dx = self.center_of_mass[child, 0] - pos[0]
                        dy = self.center_of_mass[child, 1] - pos[1]
                        dist = np.sqrt(dx**2 + dy**2)
                        if dist != 0 and self.cell_radius[child - self.nbodies][0]/dist <.5:
                            F = forces.force(pos, self.center_of_mass[child], self.mass[child])
                            acc[0] = acc[0] + F[0]
                            acc[1] = acc[1] + F[1]
                        else:
                            depth += 1
                            localNode[depth] = self.nbodies + 4*(child - self.nbodies)
                            localPos[depth] = 0
            depth -= 1
        return acc

    def __str__(self):
        cdef int i
        indent = ' '*2
        s = 'Tree :\n'
        for i in range(self.ncell+1):
            s += indent + 'cell {i}\n'.format(i=i)
            cellElements = self.child[self.nbodies + 4*i:self.nbodies + 4*i+4]
            s += 2*indent + 'box: {min} {max} \n'.format(min=np.asarray(self.cell_center[i])-np.asarray(self.cell_radius[i]), 
                                                         max=np.asarray(self.cell_center[i])+np.asarray(self.cell_radius[i]))
            s += 2*indent + 'particules: {p}\n'.format(p=cellElements[np.logical_and(0<=cellElements, cellElements<self.nbodies)])
            s += 2*indent + 'cells: {c}\n'.format(c=cellElements[cellElements>=self.nbodies]-self.nbodies)
        return s
