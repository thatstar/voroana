from libcpp.vector cimport vector
from libcpp cimport bool as cbool
from cython.operator cimport dereference


cdef extern from "voro++.hh" namespace "voro":
    cppclass container_base:
        pass

    cppclass container:
        container(double, double, double, double, double, double,
                  int, int, int, cbool, cbool, cbool, int) except +
        cbool compute_cell(voronoicell_neighbor &, c_loop_all &)
        cbool point_inside(double, double, double)
        void put(int, double, double, double)
        int total_particles()
        inline void add_wall(wall &)

    cppclass container_poly:
        container_poly(double, double, double, double, double, double,
                       int, int, int, cbool, cbool, cbool, int) except +
        cbool compute_cell(voronoicell_neighbor &, c_loop_all &)
        cbool point_inside(double, double, double)
        void put(int, double, double, double, double)
        int total_particles()
        inline void add_wall(wall &)

    cppclass container_periodic_base:
        pass

    cppclass container_periodic:
        container_periodic(double, double, double, double, double, double,
                           int, int, int, int) except +
        cbool compute_cell(voronoicell_neighbor &, c_loop_all_periodic &)
        void put(int, double, double, double)

    cppclass container_periodic_poly:
        container_periodic_poly(double, double, double, double, double, double,
                                int, int, int, int) except +
        cbool compute_cell(voronoicell_neighbor &, c_loop_all_periodic &)
        void put(int, double, double, double, double)

    cppclass voronoicell_neighbor:
        voronoicell()
        void init(double, double, double, double, double, double)
        cbool nplane(double, double, double, int)
        void centroid(double &, double &, double &)
        double volume()
        double max_radius_squared()
        double total_edge_distance()
        double surface_area()
        double number_of_faces()
        double number_of_edges()
        void vertex_orders(vector[int]&)
        void vertices(double, double, double, vector[double] &)
        void face_areas(vector[double] &)
        void face_orders(vector[int] &)
        void face_freq_table(vector[int] &)
        void face_vertices(vector[int] &)
        void face_perimeters(vector[double] &)
        void normals(vector[double] &)
        void neighbors(vector[int] &)

    cppclass c_loop_all:
        c_loop_all(container_base &)
        cbool start()
        cbool inc()
        int pid()
        void pos(double &, double &, double &)
        void pos(int &, double &, double &, double &, double &)

    cppclass c_loop_all_periodic:
        c_loop_all_periodic(container_periodic_base &)
        cbool start()
        cbool inc()
        int pid()
        void pos(double &, double &, double &)
        void pos(int &, double &, double &, double &, double &)

    cppclass wall:
        pass

    cppclass wall_plane:
        wall_plane(double, double, double, double)


cdef class Cell:
    cdef voronoicell_neighbor *thisptr
    cdef int _id
    cdef double x, y, z
    cdef double r

    def __cinit__(self):
        self.thisptr = new voronoicell_neighbor()

    def __dealloc__(self):
        del self.thisptr

    @property
    def pos(self):
        return self.x, self.y, self.z

    @property
    def radius(self):
        return self.r

    @property
    def id(self):
        return self._id

    def volume(self):
        return self.thisptr.volume()

    def max_radius_squared(self):
        return self.thisptr.max_radius_squared()

    def surface_area(self):
        return self.thisptr.surface_area()

    def number_of_edges(self):
        return self.thisptr.number_of_edges()

    def number_of_faces(self):
        return self.thisptr.number_of_faces()

    def centroid(self, cbool shift=False):
        cdef double cx = 0.0, cy = 0.0, cz = 0.0
        self.thisptr.centroid(cx, cy, cz)
        if shift:
            x, y, z = self.pos
            return cx + x, cy + y, cz + z
        else:
            return cx, cy, cz

    def vertex_orders(self):
        cdef vector[int] v
        self.thisptr.vertex_orders(v)
        return v

    def vertices(self):
        cdef vector[double] v
        self.thisptr.vertices(self.x, self.y, self.z, v)
        return list(zip(v[0::3], v[1::3], v[2::3]))

    def face_areas(self):
        cdef vector[double] v
        self.thisptr.face_areas(v)
        return v

    def face_orders(self):
        cdef vector[int] v
        self.thisptr.face_orders(v)
        return v

    def face_freq_table(self):
        cdef vector[int] v
        self.thisptr.face_freq_table(v)
        return v

    def face_vertices(self):
        cdef vector[int] v
        self.thisptr.face_vertices(v)

        mylist = []
        it = iter(v)
        while True:
            try:
                n = next(it)
            except StopIteration:
                break
            mylist.append([next(it) for _ in range(n)])
        return mylist

    def face_perimeters(self):
        cdef vector[double] v
        self.thisptr.face_perimeters(v)
        return v

    def normals(self):
        cdef vector[double] v
        self.thisptr.normals(v)
        return list(zip(v[0::3], v[1::3], v[2::3]))

    def neighbors(self):
        cdef vector[int] v
        self.thisptr.neighbors(v)
        return v

    def init(self, double xmin, double xmax, double ymin, double ymax, double zmin, double zmax):
        self.thisptr.init(xmin, xmax, ymin, ymax, zmin, zmax)

    def nplane(self, double x, double y, double z, int p_id):
        return bool(self.thisptr.nplane(x, y, z, p_id))

    def __str__(self):
        return '<Cell {0}>'.format(self._id)

    def __repr__(self):
        return '<Cell {0}>'.format(self._id)


cdef class Container:
    cdef container *thisptr

    def __cinit__(self, double ax_, double bx_, double ay_, double by_,
                  double az_, double bz_, int nx_, int ny_, int nz_,
                  cbool xperiodic_, cbool yperiodic_, cbool zperiodic_, int init_mem):
        self.thisptr = new container(ax_, bx_, ay_, by_, az_, bz_, nx_, ny_, nz_,
                                     xperiodic_, yperiodic_, zperiodic_, init_mem)

    def __dealloc__(self):
        del self.thisptr

    def put(self, int n, double x, double y, double z):
        assert self.thisptr.point_inside(x, y, z)
        self.thisptr.put(n, x, y, z)

    def get_cells(self):
        cdef container_base *baseptr = (<container_base *>(self.thisptr))
        cdef c_loop_all *vl = new c_loop_all(dereference(baseptr))

        cell = Cell()

        cdef int vcells_left = self.thisptr.total_particles()
        cdef int id

        mylist = [None for _ in range(vcells_left)]

        if not vl.start():
            del vl
            raise RuntimeError("Failed to start loop!")

        while True:
            if (self.thisptr.compute_cell(dereference(cell.thisptr), dereference(vl))):
                cell._id = vl.pid()
                assert cell._id < self.thisptr.total_particles(), \
                    (
                        "Cell id %s larger than total %s!" % (cell._id, self.thisptr.total_particles())
                    )
                vl.pos(cell.x, cell.y, cell.z)
                cell.r = 0
                mylist[cell._id] = cell
                vcells_left -= 1
                cell = Cell()

            if not vl.inc():
                break

        del vl

        if vcells_left != 0:
            raise RuntimeError("Voronoi computation failed!")
        return mylist

    def add_wall_plane(self, double x, double y, double z, double a):
        cdef wall_plane *wp = new wall_plane(x, y, z, a)
        cdef wall *w = (<wall *>(wp))
        self.thisptr.add_wall(dereference(w))


cdef class ContainerPoly:
    cdef container_poly *thisptr

    def __cinit__(self, double ax_, double bx_, double ay_, double by_,
                  double az_, double bz_, int nx_, int ny_, int nz_,
                  cbool xperiodic_, cbool yperiodic_, cbool zperiodic_, int init_mem):
        self.thisptr = new container_poly(ax_, bx_, ay_, by_, az_, bz_, nx_, ny_, nz_,
                                          xperiodic_, yperiodic_, zperiodic_, init_mem)

    def __dealloc__(self):
        del self.thisptr

    def put(self, int n, double x, double y, double z, double r):
        assert self.thisptr.point_inside(x, y, z)
        self.thisptr.put(n, x, y, z, r)

    def get_cells(self):
        cdef container_base *baseptr = (<container_base *>(self.thisptr))
        cdef c_loop_all *vl = new c_loop_all(dereference(baseptr))

        cell = Cell()

        cdef int vcells_left = self.thisptr.total_particles()
        cdef int id

        mylist = [None for _ in range(vcells_left)]

        if not vl.start():
            del vl
            raise RuntimeError("Failed to start loop!")

        while True:
            if self.thisptr.compute_cell(dereference(cell.thisptr), dereference(vl)):
                cell._id = vl.pid()
                assert cell._id < self.thisptr.total_particles(), \
                    (
                        "Cell id %s larger than total %s!" % (cell._id, self.thisptr.total_particles())
                    )
                mylist[cell._id] = cell
                vcells_left -= 1
                cell = Cell()

            if not vl.inc():
                break

        del vl

        if vcells_left != 0:
            raise RuntimeError("Voronoi computation failed!")
        return mylist

    def add_wall_plane(self, double x, double y, double z, double a):
        cdef wall_plane *wp = new wall_plane(x, y, z, a)
        cdef wall *w = (<wall *>(wp))
        self.thisptr.add_wall(dereference(w))


cdef class ContainerPeriodic:
    cdef container_periodic *thisptr
    cdef int nparticles

    def __cinit__(self,
                  double bx_, double bxy_, double by_, double bxz_, double byz_, double bz_,
                  int nx_, int ny_, int nz_, int init_mem):
        self.thisptr = new container_periodic(bx_, bxy_, by_, bxz_, byz_, bz_,
                                              nx_, ny_, nz_, init_mem)
        nparticles = 0

    def __dealloc__(self):
        del self.thisptr

    def put(self, int n, double x, double y, double z):
        self.thisptr.put(n, x, y, z)
        self.nparticles += 1

    def get_cells(self):
        cdef container_periodic_base *baseptr = (<container_periodic_base *>(self.thisptr))
        cdef c_loop_all_periodic *vl = new c_loop_all_periodic(dereference(baseptr))

        cell = Cell()

        cdef int vcells_left = self.nparticles
        cdef int id

        mylist = [None for _ in range(vcells_left)]

        if not vl.start():
            del vl
            raise RuntimeError("Failed to start loop!")

        while True:
            if (self.thisptr.compute_cell(dereference(cell.thisptr), dereference(vl))):
                cell._id = vl.pid()
                assert cell._id < self.nparticles, \
                    (
                        "Cell id %s larger than total {}!".format(cell._id, self.nparticles)
                    )
                vl.pos(cell.x, cell.y, cell.z)
                cell.r = 0
                mylist[cell._id] = cell
                vcells_left -= 1
                cell = Cell()

            if not vl.inc():
                break

        del vl

        if vcells_left != 0:
            raise RuntimeError("Voronoi computation failed!")
        return mylist


cdef class ContainerPeriodicPoly:
    cdef container_periodic_poly *thisptr
    cdef int nparticles

    def __cinit__(self,
                  double bx_, double bxy_, double by_, double bxz_, double byz_, double bz_,
                  int nx_, int ny_, int nz_, int init_mem):
        self.thisptr = new container_periodic_poly(bx_, bxy_, by_, bxz_, byz_, bz_,
                                                   nx_, ny_, nz_, init_mem)
        self.nparticles = 0

    def __dealloc__(self):
        del self.thisptr

    def put(self, int n, double x, double y, double z, double r):
        self.thisptr.put(n, x, y, z, r)
        self.nparticles += 1

    def get_cells(self):
        cdef container_periodic_base *baseptr = (<container_periodic_base *>(self.thisptr))
        cdef c_loop_all_periodic *vl = new c_loop_all_periodic(dereference(baseptr))

        cell = Cell()

        cdef int vcells_left = self.nparticles
        cdef int id

        mylist = [None for _ in range(vcells_left)]

        if not vl.start():
            del vl
            raise RuntimeError("Failed to start loop!")

        while True:
            if self.thisptr.compute_cell(dereference(cell.thisptr), dereference(vl)):
                cell._id = vl.pid()
                assert cell._id < self.nparticles, \
                    (
                        "Cell id %s larger than total {}!".format(cell._id, self.nparticles)
                    )
                mylist[cell._id] = cell
                vcells_left -= 1
                cell = Cell()

            if not vl.inc():
                break

        del vl

        if vcells_left != 0:
            raise RuntimeError("Voronoi computation failed!")
        return mylist
