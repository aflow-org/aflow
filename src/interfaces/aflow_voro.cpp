// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2021           *
// *                                                                         *
// ***************************************************************************

// aflow_voro.h and aflow_voro.cpp are based on the Voro++ created by Chris H. Rycroft
// licence details can be found in aflow_voro.h
#include "aflow_voro.h"

#include <cstdio>
#include <vector>

#include "AUROSTD/aurostd.h"
#include "AUROSTD/aurostd_xmatrix.h"
#include "AUROSTD/aurostd_xvector.h"

// Voro++, a 3D cell-based Voronoi library
//
// Author   : Chris H. Rycroft (LBL / UC Berkeley)
// Email    : chr@alum.mit.edu
// Date     : August 30th 2011

/** \file common.cc
 * \brief Implementations of the small helper functions. */

using aurostd::xmatrix;
using aurostd::xvector;
using std::vector;

namespace voro {

/** \brief Prints a vector of integers.
   *
   * Prints a vector of integers.
   * \param[in] v the vector to print.
   * \param[in] fp the file stream to print to. */
  void voro_print_vector(std::vector<int> &v, FILE *fp) {
    int k = 0;
    const int s = v.size();
    while (k + 4 < s) {
      fprintf(fp, "%d %d %d %d ", v[k], v[k + 1], v[k + 2], v[k + 3]);
      k += 4;
    }
    if (k + 3 <= s) {
      if (k + 4 == s) {
        fprintf(fp, "%d %d %d %d", v[k], v[k + 1], v[k + 2], v[k + 3]);
      } else {
        fprintf(fp, "%d %d %d", v[k], v[k + 1], v[k + 2]);
      }
    } else {
      if (k + 2 == s) {
        fprintf(fp, "%d %d", v[k], v[k + 1]);
      } else {
        fprintf(fp, "%d", v[k]);
      }
    }
  }

/** \brief Prints a vector of doubles.
   *
   * Prints a vector of doubles.
   * \param[in] v the vector to print.
   * \param[in] fp the file stream to print to. */
  void voro_print_vector(std::vector<double> &v, FILE *fp) {
    int k = 0;
    const int s = v.size();
    while (k + 4 < s) {
      fprintf(fp, "%g %g %g %g ", v[k], v[k + 1], v[k + 2], v[k + 3]);
      k += 4;
    }
    if (k + 3 <= s) {
      if (k + 4 == s) {
        fprintf(fp, "%g %g %g %g", v[k], v[k + 1], v[k + 2], v[k + 3]);
      } else {
        fprintf(fp, "%g %g %g", v[k], v[k + 1], v[k + 2]);
      }
    } else {
      if (k + 2 == s) {
        fprintf(fp, "%g %g", v[k], v[k + 1]);
      } else {
        fprintf(fp, "%g", v[k]);
      }
    }
  }

/** \brief Prints a vector a face vertex information.
   *
   * Prints a vector of face vertex information. A value is read, which
   * corresponds to the number of vertices in the next face. The routine reads
   * this number of values and prints them as a bracked list. This is repeated
   * until the end of the vector is reached.
   * \param[in] v the vector to interpret and print.
   * \param[in] fp the file stream to print to. */
  void voro_print_face_vertices(std::vector<int> &v, FILE *fp) {
    int j;
    int k = 0;
    int l;
    if (!v.empty()) {
      l = v[k++];
      if (l <= 1) {
        if (l == 1) {
          fprintf(fp, "(%d)", v[k++]);
        } else {
          fputs("()", fp);
        }
      } else {
        j = k + l;
        fprintf(fp, "(%d", v[k++]);
        while (k < j) {
          fprintf(fp, ",%d", v[k++]);
        }
        fputs(")", fp);
      }
      while ((unsigned int) k < v.size()) {
        l = v[k++];
        if (l <= 1) {
          if (l == 1) {
            fprintf(fp, " (%d)", v[k++]);
          } else {
            fputs(" ()", fp);
          }
        } else {
          j = k + l;
          fprintf(fp, " (%d", v[k++]);
          while (k < j) {
            fprintf(fp, ",%d", v[k++]);
          }
          fputs(")", fp);
        }
      }
    }
  }

} // namespace voro
// Voro++, a 3D cell-based Voronoi library
//
// Author   : Chris H. Rycroft (LBL / UC Berkeley)
// Email    : chr@alum.mit.edu
// Date     : August 30th 2011

/** \file cell.cc
 * \brief Function implementations for the voronoicell and related classes. */

#include <cmath>
#include <cstring>

namespace voro {

/** Constructs a Voronoi cell and sets up the initial memory. */
  voronoicell_base::voronoicell_base() :
      current_vertices(init_vertices),
      current_vertex_order(init_vertex_order),
      current_delete_size(init_delete_size),
      current_delete2_size(init_delete2_size),
      ed(new int *[current_vertices]),
      nu(new int[current_vertices]),
      pts(new double[3 * current_vertices]),
      mem(new int[current_vertex_order]),
      mec(new int[current_vertex_order]),
      mep(new int *[current_vertex_order]),
      ds(new int[current_delete_size]),
      stacke(ds + current_delete_size),
      ds2(new int[current_delete2_size]),
      stacke2(ds2 + current_delete_size),
      current_marginal(init_marginal),
      marg(new int[current_marginal]) {
    int i;
    for (i = 0; i < 3; i++) {
      mem[i] = init_n_vertices;
      mec[i] = 0;
      mep[i] = new int[init_n_vertices * ((i << 1) + 1)];
    }
    mem[3] = init_3_vertices;
    mec[3] = 0;
    mep[3] = new int[init_3_vertices * 7];
    for (i = 4; i < current_vertex_order; i++) {
      mem[i] = init_n_vertices;
      mec[i] = 0;
      mep[i] = new int[init_n_vertices * ((i << 1) + 1)];
    }
  }

/** The voronoicell destructor deallocates all the dynamic memory. */
  voronoicell_base::~voronoicell_base() {
    for (int i = current_vertex_order - 1; i >= 0; i--) {
      if (mem[i] > 0) {
        delete[] mep[i];
      }
    }
    delete[] marg;
    delete[] ds2;
    delete[] ds;
    delete[] mep;
    delete[] mec;
    delete[] mem;
    delete[] pts;
    delete[] nu;
    delete[] ed;
  }

/** Ensures that enough memory is allocated prior to carrying out a copy.
   * \param[in] vc a reference to the specialized version of the calling class.
   * \param[in] vb a pointered to the class to be copied. */
  template <class vc_class> void voronoicell_base::check_memory_for_copy(vc_class &vc, voronoicell_base *vb) {
    while (current_vertex_order < vb->current_vertex_order) {
      add_memory_vorder(vc);
    }
    for (int i = 0; i < current_vertex_order; i++) {
      while (mem[i] < vb->mec[i]) {
        add_memory(vc, i, ds2);
      }
    }
    while (current_vertices < vb->p) {
      add_memory_vertices(vc);
    }
  }

/** Copies the vertex and edge information from another class. The routine
   * assumes that enough memory is available for the copy.
   * \param[in] vb a pointer to the class to copy. */
  void voronoicell_base::copy(voronoicell_base *vb) {
    int i;
    int j;
    p = vb->p;
    up = 0;
    for (i = 0; i < current_vertex_order; i++) {
      mec[i] = vb->mec[i];
      for (j = 0; j < mec[i] * (2 * i + 1); j++) {
        mep[i][j] = vb->mep[i][j];
      }
      for (j = 0; j < mec[i] * (2 * i + 1); j += 2 * i + 1) {
        ed[mep[i][j + 2 * i]] = mep[i] + j;
      }
    }
    for (i = 0; i < p; i++) {
      nu[i] = vb->nu[i];
    }
    for (i = 0; i < 3 * p; i++) {
      pts[i] = vb->pts[i];
    }
  }

/** Copies the information from another voronoicell class into this
   * class, extending memory allocation if necessary.
   * \param[in] c the class to copy. */
  void voronoicell_neighbor::operator=(voronoicell &c) {
    voronoicell_base *vb = ((voronoicell_base *) &c);
    check_memory_for_copy(*this, vb);
    copy(vb);
    int i;
    int j;
    for (i = 0; i < c.current_vertex_order; i++) {
      for (j = 0; j < c.mec[i] * i; j++) {
        mne[i][j] = 0;
      }
      for (j = 0; j < c.mec[i]; j++) {
        ne[c.mep[i][(2 * i + 1) * j + 2 * i]] = mne[i] + (j * i);
      }
    }
  }

/** Copies the information from another voronoicell_neighbor class into this
   * class, extending memory allocation if necessary.
   * \param[in] c the class to copy. */
  void voronoicell_neighbor::operator=(voronoicell_neighbor &c) {
    voronoicell_base *vb = ((voronoicell_base *) &c);
    check_memory_for_copy(*this, vb);
    copy(vb);
    int i;
    int j;
    for (i = 0; i < c.current_vertex_order; i++) {
      for (j = 0; j < c.mec[i] * i; j++) {
        mne[i][j] = c.mne[i][j];
      }
      for (j = 0; j < c.mec[i]; j++) {
        ne[c.mep[i][(2 * i + 1) * j + 2 * i]] = mne[i] + (j * i);
      }
    }
  }

/** Translates the vertices of the Voronoi cell by a given vector.
   * \param[in] (x,y,z) the coordinates of the vector. */
  void voronoicell_base::translate(double x, double y, double z) const {
    x *= 2;
    y *= 2;
    z *= 2;
    double *ptsp = pts;
    while (ptsp < pts + 3 * p) {
      *(ptsp++) = x;
      *(ptsp++) = y;
      *(ptsp++) = z;
    }
  }

/** Increases the memory storage for a particular vertex order, by increasing
   * the size of the of the corresponding mep array. If the arrays already exist,
   * their size is doubled; if they don't exist, then new ones of size
   * init_n_vertices are allocated. The routine also ensures that the pointers in
   * the ed array are updated, by making use of the back pointers. For the cases
   * where the back pointer has been temporarily overwritten in the marginal
   * vertex code, the auxiliary delete stack is scanned to find out how to update
   * the ed value. If the template has been instantiated with the neighbor
   * tracking turned on, then the routine also reallocates the corresponding mne
   * array.
   * \param[in] i the order of the vertex memory to be increased. */
  template <class vc_class> void voronoicell_base::add_memory(vc_class &vc, int i, int *stackp2) {
    const int s = (i << 1) + 1;
    if (mem[i] == 0) {
      vc.n_allocate(i, init_n_vertices);
      mep[i] = new int[init_n_vertices * s];
      mem[i] = init_n_vertices;
#if VOROPP_VERBOSE >= 2
      fprintf(stderr, "Order %d vertex memory created\n", i);
#endif
    } else {
      int j = 0;
      int k;
      int *l;
      mem[i] <<= 1;
      if (mem[i] > max_n_vertices) {
        voro_fatal_error("Point memory allocation exceeded absolute maximum", VOROPP_MEMORY_ERROR);
      }
#if VOROPP_VERBOSE >= 2
      fprintf(stderr, "Order %d vertex memory scaled up to %d\n", i, mem[i]);
#endif
      l = new int[s * mem[i]];
      int m = 0;
      vc.n_allocate_aux1(i);
      while (j < s * mec[i]) {
        k = mep[i][j + (i << 1)];
        if (k >= 0) {
          ed[k] = l + j;
          vc.n_set_to_aux1_offset(k, m);
        } else {
          int *dsp;
          for (dsp = ds2; dsp < stackp2; dsp++) {
            if (ed[*dsp] == mep[i] + j) {
              ed[*dsp] = l + j;
              vc.n_set_to_aux1_offset(*dsp, m);
              break;
            }
          }
          if (dsp == stackp2) {
            voro_fatal_error("Couldn't relocate dangling pointer", VOROPP_INTERNAL_ERROR);
          }
#if VOROPP_VERBOSE >= 3
          fputs("Relocated dangling pointer", stderr);
#endif
        }
        for (k = 0; k < s; k++, j++) {
          l[j] = mep[i][j];
        }
        for (k = 0; k < i; k++, m++) {
          vc.n_copy_to_aux1(i, m);
        }
      }
      delete[] mep[i];
      mep[i] = l;
      vc.n_switch_to_aux1(i);
    }
  }

/** Doubles the maximum number of vertices allowed, by reallocating the ed, nu,
   * and pts arrays. If the allocation exceeds the absolute maximum set in
   * max_vertices, then the routine exits with a fatal error. If the template has
   * been instantiated with the neighbor tracking turned on, then the routine
   * also reallocates the ne array. */
  template <class vc_class> void voronoicell_base::add_memory_vertices(vc_class &vc) {
    const int i = (current_vertices << 1);
    int j;
    int **pp;
    int *pnu;
    if (i > max_vertices) {
      voro_fatal_error("Vertex memory allocation exceeded absolute maximum", VOROPP_MEMORY_ERROR);
    }
#if VOROPP_VERBOSE >= 2
    fprintf(stderr, "Vertex memory scaled up to %d\n", i);
#endif
    double *ppts;
    pp = new int *[i];
    for (j = 0; j < current_vertices; j++) {
      pp[j] = ed[j];
    }
    delete[] ed;
    ed = pp;
    vc.n_add_memory_vertices(i);
    pnu = new int[i];
    for (j = 0; j < current_vertices; j++) {
      pnu[j] = nu[j];
    }
    delete[] nu;
    nu = pnu;
    ppts = new double[3 * i];
    for (j = 0; j < 3 * current_vertices; j++) {
      ppts[j] = pts[j];
    }
    delete[] pts;
    pts = ppts;
    current_vertices = i;
  }

/** Doubles the maximum allowed vertex order, by reallocating mem, mep, and mec
   * arrays. If the allocation exceeds the absolute maximum set in
   * max_vertex_order, then the routine causes a fatal error. If the template has
   * been instantiated with the neighbor tracking turned on, then the routine
   * also reallocates the mne array. */
  template <class vc_class> void voronoicell_base::add_memory_vorder(vc_class &vc) {
    const int i = (current_vertex_order << 1);
    int j;
    int *p1;
    int **p2;
    if (i > max_vertex_order) {
      voro_fatal_error("Vertex order memory allocation exceeded absolute maximum", VOROPP_MEMORY_ERROR);
    }
#if VOROPP_VERBOSE >= 2
    fprintf(stderr, "Vertex order memory scaled up to %d\n", i);
#endif
    p1 = new int[i];
    for (j = 0; j < current_vertex_order; j++) {
      p1[j] = mem[j];
    }
    while (j < i) {
      p1[j++] = 0;
    }
    delete[] mem;
    mem = p1;
    p2 = new int *[i];
    for (j = 0; j < current_vertex_order; j++) {
      p2[j] = mep[j];
    }
    delete[] mep;
    mep = p2;
    p1 = new int[i];
    for (j = 0; j < current_vertex_order; j++) {
      p1[j] = mec[j];
    }
    while (j < i) {
      p1[j++] = 0;
    }
    delete[] mec;
    mec = p1;
    vc.n_add_memory_vorder(i);
    current_vertex_order = i;
  }

/** Doubles the size allocation of the main delete stack. If the allocation
   * exceeds the absolute maximum set in max_delete_size, then routine causes a
   * fatal error. */
  void voronoicell_base::add_memory_ds(int *&stackp) {
    current_delete_size <<= 1;
    if (current_delete_size > max_delete_size) {
      voro_fatal_error("Delete stack 1 memory allocation exceeded absolute maximum", VOROPP_MEMORY_ERROR);
    }
#if VOROPP_VERBOSE >= 2
    fprintf(stderr, "Delete stack 1 memory scaled up to %d\n", current_delete_size);
#endif
    int *dsn = new int[current_delete_size];
    int *dsnp = dsn;
    int *dsp = ds;
    while (dsp < stackp) {
      *(dsnp++) = *(dsp++);
    }
    delete[] ds;
    ds = dsn;
    stackp = dsnp;
    stacke = ds + current_delete_size;
  }

/** Doubles the size allocation of the auxiliary delete stack. If the
   * allocation exceeds the absolute maximum set in max_delete2_size, then the
   * routine causes a fatal error. */
  void voronoicell_base::add_memory_ds2(int *&stackp2) {
    current_delete2_size <<= 1;
    if (current_delete2_size > max_delete2_size) {
      voro_fatal_error("Delete stack 2 memory allocation exceeded absolute maximum", VOROPP_MEMORY_ERROR);
    }
#if VOROPP_VERBOSE >= 2
    fprintf(stderr, "Delete stack 2 memory scaled up to %d\n", current_delete2_size);
#endif
    int *dsn = new int[current_delete2_size];
    int *dsnp = dsn;
    int *dsp = ds2;
    while (dsp < stackp2) {
      *(dsnp++) = *(dsp++);
    }
    delete[] ds2;
    ds2 = dsn;
    stackp2 = dsnp;
    stacke2 = ds2 + current_delete2_size;
  }

/** Initializes a Voronoi cell as a rectangular box with the given dimensions.
   * \param[in] (xmin,xmax) the minimum and maximum x coordinates.
   * \param[in] (ymin,ymax) the minimum and maximum y coordinates.
   * \param[in] (zmin,zmax) the minimum and maximum z coordinates. */
  void voronoicell_base::init_base(double xmin, double xmax, double ymin, double ymax, double zmin, double zmax) {
    for (int i = 0; i < current_vertex_order; i++) {
      mec[i] = 0;
    }
    up = 0;
    mec[3] = p = 8;
    xmin *= 2;
    xmax *= 2;
    ymin *= 2;
    ymax *= 2;
    zmin *= 2;
    zmax *= 2;
    *pts = xmin;
    pts[1] = ymin;
    pts[2] = zmin;
    pts[3] = xmax;
    pts[4] = ymin;
    pts[5] = zmin;
    pts[6] = xmin;
    pts[7] = ymax;
    pts[8] = zmin;
    pts[9] = xmax;
    pts[10] = ymax;
    pts[11] = zmin;
    pts[12] = xmin;
    pts[13] = ymin;
    pts[14] = zmax;
    pts[15] = xmax;
    pts[16] = ymin;
    pts[17] = zmax;
    pts[18] = xmin;
    pts[19] = ymax;
    pts[20] = zmax;
    pts[21] = xmax;
    pts[22] = ymax;
    pts[23] = zmax;
    int *q = mep[3];
    *q = 1;
    q[1] = 4;
    q[2] = 2;
    q[3] = 2;
    q[4] = 1;
    q[5] = 0;
    q[6] = 0;
    q[7] = 3;
    q[8] = 5;
    q[9] = 0;
    q[10] = 2;
    q[11] = 1;
    q[12] = 0;
    q[13] = 1;
    q[14] = 0;
    q[15] = 6;
    q[16] = 3;
    q[17] = 2;
    q[18] = 1;
    q[19] = 0;
    q[20] = 2;
    q[21] = 2;
    q[22] = 7;
    q[23] = 1;
    q[24] = 2;
    q[25] = 1;
    q[26] = 0;
    q[27] = 3;
    q[28] = 6;
    q[29] = 0;
    q[30] = 5;
    q[31] = 2;
    q[32] = 1;
    q[33] = 0;
    q[34] = 4;
    q[35] = 4;
    q[36] = 1;
    q[37] = 7;
    q[38] = 2;
    q[39] = 1;
    q[40] = 0;
    q[41] = 5;
    q[42] = 7;
    q[43] = 2;
    q[44] = 4;
    q[45] = 2;
    q[46] = 1;
    q[47] = 0;
    q[48] = 6;
    q[49] = 5;
    q[50] = 3;
    q[51] = 6;
    q[52] = 2;
    q[53] = 1;
    q[54] = 0;
    q[55] = 7;
    *ed = q;
    ed[1] = q + 7;
    ed[2] = q + 14;
    ed[3] = q + 21;
    ed[4] = q + 28;
    ed[5] = q + 35;
    ed[6] = q + 42;
    ed[7] = q + 49;
    *nu = nu[1] = nu[2] = nu[3] = nu[4] = nu[5] = nu[6] = nu[7] = 3;
  }

/** Initializes a Voronoi cell as a regular octahedron.
   * \param[in] l The distance from the octahedron center to a vertex. Six
   *              vertices are initialized at (-l,0,0), (l,0,0), (0,-l,0),
   *              (0,l,0), (0,0,-l), and (0,0,l). */
  void voronoicell_base::init_octahedron_base(double l) {
    for (int i = 0; i < current_vertex_order; i++) {
      mec[i] = 0;
    }
    up = 0;
    mec[4] = p = 6;
    l *= 2;
    *pts = -l;
    pts[1] = 0;
    pts[2] = 0;
    pts[3] = l;
    pts[4] = 0;
    pts[5] = 0;
    pts[6] = 0;
    pts[7] = -l;
    pts[8] = 0;
    pts[9] = 0;
    pts[10] = l;
    pts[11] = 0;
    pts[12] = 0;
    pts[13] = 0;
    pts[14] = -l;
    pts[15] = 0;
    pts[16] = 0;
    pts[17] = l;
    int *q = mep[4];
    *q = 2;
    q[1] = 5;
    q[2] = 3;
    q[3] = 4;
    q[4] = 0;
    q[5] = 0;
    q[6] = 0;
    q[7] = 0;
    q[8] = 0;
    q[9] = 2;
    q[10] = 4;
    q[11] = 3;
    q[12] = 5;
    q[13] = 2;
    q[14] = 2;
    q[15] = 2;
    q[16] = 2;
    q[17] = 1;
    q[18] = 0;
    q[19] = 4;
    q[20] = 1;
    q[21] = 5;
    q[22] = 0;
    q[23] = 3;
    q[24] = 0;
    q[25] = 1;
    q[26] = 2;
    q[27] = 0;
    q[28] = 5;
    q[29] = 1;
    q[30] = 4;
    q[31] = 2;
    q[32] = 3;
    q[33] = 2;
    q[34] = 1;
    q[35] = 3;
    q[36] = 0;
    q[37] = 3;
    q[38] = 1;
    q[39] = 2;
    q[40] = 3;
    q[41] = 3;
    q[42] = 1;
    q[43] = 1;
    q[44] = 4;
    q[45] = 0;
    q[46] = 2;
    q[47] = 1;
    q[48] = 3;
    q[49] = 1;
    q[50] = 3;
    q[51] = 3;
    q[52] = 1;
    q[53] = 5;
    *ed = q;
    ed[1] = q + 9;
    ed[2] = q + 18;
    ed[3] = q + 27;
    ed[4] = q + 36;
    ed[5] = q + 45;
    *nu = nu[1] = nu[2] = nu[3] = nu[4] = nu[5] = 4;
  }

/** Initializes a Voronoi cell as a tetrahedron. It assumes that the normal to
   * the face for the first three vertices points inside.
   * \param (x0,y0,z0) a position vector for the first vertex.
   * \param (x1,y1,z1) a position vector for the second vertex.
   * \param (x2,y2,z2) a position vector for the third vertex.
   * \param (x3,y3,z3) a position vector for the fourth vertex. */
  void voronoicell_base::init_tetrahedron_base(double x0, double y0, double z0, double x1, double y1, double z1, double x2, double y2, double z2, double x3, double y3, double z3) {
    for (int i = 0; i < current_vertex_order; i++) {
      mec[i] = 0;
    }
    up = 0;
    mec[3] = p = 4;
    *pts = x0 * 2;
    pts[1] = y0 * 2;
    pts[2] = z0 * 2;
    pts[3] = x1 * 2;
    pts[4] = y1 * 2;
    pts[5] = z1 * 2;
    pts[6] = x2 * 2;
    pts[7] = y2 * 2;
    pts[8] = z2 * 2;
    pts[9] = x3 * 2;
    pts[10] = y3 * 2;
    pts[11] = z3 * 2;
    int *q = mep[3];
    *q = 1;
    q[1] = 3;
    q[2] = 2;
    q[3] = 0;
    q[4] = 0;
    q[5] = 0;
    q[6] = 0;
    q[7] = 0;
    q[8] = 2;
    q[9] = 3;
    q[10] = 0;
    q[11] = 2;
    q[12] = 1;
    q[13] = 1;
    q[14] = 0;
    q[15] = 3;
    q[16] = 1;
    q[17] = 2;
    q[18] = 2;
    q[19] = 1;
    q[20] = 2;
    q[21] = 0;
    q[22] = 1;
    q[23] = 2;
    q[24] = 1;
    q[25] = 2;
    q[26] = 1;
    q[27] = 3;
    *ed = q;
    ed[1] = q + 7;
    ed[2] = q + 14;
    ed[3] = q + 21;
    *nu = nu[1] = nu[2] = nu[3] = 3;
  }

/** Checks that the relational table of the Voronoi cell is accurate, and
   * prints out any errors. This algorithm is O(p), so running it every time the
   * plane routine is called will result in a significant slowdown. */
  void voronoicell_base::check_relations() const {
    int i;
    int j;
    for (i = 0; i < p; i++) {
      for (j = 0; j < nu[i]; j++) {
        if (ed[ed[i][j]][ed[i][nu[i] + j]] != i) {
          printf("Relational error at point %d, edge %d.\n", i, j);
        }
      }
    }
  }

/** This routine checks for any two vertices that are connected by more than
   * one edge. The plane algorithm is designed so that this should not happen, so
   * any occurrences are most likely errors. Note that the routine is O(p), so
   * running it every time the plane routine is called will result in a
   * significant slowdown. */
  void voronoicell_base::check_duplicates() const {
    int i;
    int j;
    int k;
    for (i = 0; i < p; i++) {
      for (j = 1; j < nu[i]; j++) {
        for (k = 0; k < j; k++) {
          if (ed[i][j] == ed[i][k]) {
            printf("Duplicate edges: (%d,%d) and (%d,%d) [%d]\n", i, j, i, k, ed[i][j]);
          }
        }
      }
    }
  }

/** Constructs the relational table if the edges have been specified. */
  void voronoicell_base::construct_relations() const {
    int i;
    int j;
    int k;
    int l;
    for (i = 0; i < p; i++) {
      for (j = 0; j < nu[i]; j++) {
        k = ed[i][j];
        l = 0;
        while (ed[k][l] != i) {
          l++;
          if (l == nu[k]) {
            voro_fatal_error("Relation table construction failed", VOROPP_INTERNAL_ERROR);
          }
        }
        ed[i][nu[i] + j] = l;
      }
    }
  }

/** Starting from a point within the current cutting plane, this routine attempts
   * to find an edge to a point outside the cutting plane. This prevents the plane
   * routine from .
   * \param[in] vc a reference to the specialized version of the calling class.
   * \param[in,out] up */
  template <class vc_class> inline bool voronoicell_base::search_for_outside_edge(vc_class &vc, int &up) {
    int i;
    int lp;
    int lw;
    int *j(ds2);
    int *stackp2(ds2);
    double l;
    *(stackp2++) = up;
    while (j < stackp2) {
      up = *(j++);
      for (i = 0; i < nu[up]; i++) {
        lp = ed[up][i];
        lw = m_test(lp, l);
        if (lw == -1) {
          return true;
        } else if (lw == 0) {
          add_to_stack(vc, lp, stackp2);
        }
      }
    }
    return false;
  }

/** Adds a point to the auxiliary delete stack if it is not already there.
   * \param[in] vc a reference to the specialized version of the calling class.
   * \param[in] lp the index of the point to add.
   * \param[in,out] stackp2 a pointer to the end of the stack entries. */
  template <class vc_class> inline void voronoicell_base::add_to_stack(vc_class &vc, int lp, int *&stackp2) {
    for (int *k(ds2); k < stackp2; k++) {
      if (*k == lp) {
        return;
      }
    }
    if (stackp2 == stacke2) {
      add_memory_ds2(stackp2);
    }
    *(stackp2++) = lp;
  }

/** Cuts the Voronoi cell by a particle whose center is at a separation of
   * (x,y,z) from the cell center. The value of rsq should be initially set to
   * \f$x^2+y^2+z^2\f$.
   * \param[in] vc a reference to the specialized version of the calling class.
   * \param[in] (x,y,z) the normal vector to the plane.
   * \param[in] rsq the distance along this vector of the plane.
   * \param[in] p_id the plane ID (for neighbor tracking only).
   * \return False if the plane cut deleted the cell entirely, true otherwise. */
  template <class vc_class> bool voronoicell_base::nplane(vc_class &vc, double x, double y, double z, double rsq, int p_id) {
    int count = 0;
    int i;
    int j;
    int k;
    int lp = up;
    int cp;
    int qp;
    int rp;
    int *stackp(ds);
    int *stackp2(ds2);
    int *dsp;
    int us = 0;
    int ls = 0;
    int qs;
    int iqs;
    int cs;
    int uw;
    int qw;
    int lw;
    int *edp;
    int *edd;
    double u;
    double l;
    double r;
    double q;
    bool complicated_setup = false;
    bool new_double_edge = false;
    bool double_edge = false;

    // Initialize the safe testing routine
    n_marg = 0;
    px = x;
    py = y;
    pz = z;
    prsq = rsq;

    // Test approximately sqrt(n)/4 points for their proximity to the plane
    // and keep the one which is closest
    uw = m_test(up, u);

    // Starting from an initial guess, we now move from vertex to vertex,
    // to try and find an edge which intersects the cutting plane,
    // or a vertex which is on the plane
    try {
      if (uw == 1) {
        // The test point is inside the cutting plane.
        us = 0;
        do {
          lp = ed[up][us];
          lw = m_test(lp, l);
          if (l < u) {
            break;
          }
          us++;
        } while (us < nu[up]);

        if (us == nu[up]) {
          return false;
        }

        ls = ed[up][nu[up] + us];
        while (lw == 1) {
          if (++count >= p) {
            throw true;
          }
          u = l;
          up = lp;
          for (us = 0; us < ls; us++) {
            lp = ed[up][us];
            lw = m_test(lp, l);
            if (l < u) {
              break;
            }
          }
          if (us == ls) {
            us++;
            while (us < nu[up]) {
              lp = ed[up][us];
              lw = m_test(lp, l);
              if (l < u) {
                break;
              }
              us++;
            }
            if (us == nu[up]) {
              return false;
            }
          }
          ls = ed[up][nu[up] + us];
        }

        // If the last point in the iteration is within the
        // plane, we need to do the complicated setup
        // routine. Otherwise, we use the regular iteration.
        if (lw == 0) {
          up = lp;
          complicated_setup = true;
        } else {
          complicated_setup = false;
        }
      } else if (uw == -1) {
        us = 0;
        do {
          qp = ed[up][us];
          qw = m_test(qp, q);
          if (u < q) {
            break;
          }
          us++;
        } while (us < nu[up]);
        if (us == nu[up]) {
          return true;
        }

        while (qw == -1) {
          qs = ed[up][nu[up] + us];
          if (++count >= p) {
            throw true;
          }
          u = q;
          up = qp;
          for (us = 0; us < qs; us++) {
            qp = ed[up][us];
            qw = m_test(qp, q);
            if (u < q) {
              break;
            }
          }
          if (us == qs) {
            us++;
            while (us < nu[up]) {
              qp = ed[up][us];
              qw = m_test(qp, q);
              if (u < q) {
                break;
              }
              us++;
            }
            if (us == nu[up]) {
              return true;
            }
          }
        }
        if (qw == 1) {
          lp = up;
          ls = us;
          l = u;
          up = qp;
          us = ed[lp][nu[lp] + ls];
          u = q;
          complicated_setup = false;
        } else {
          up = qp;
          complicated_setup = true;
        }
      } else {
        // Our original test point was on the plane, so we
        // automatically head for the complicated setup
        // routine
        complicated_setup = true;
      }
    } catch (bool except) {
      // This routine is a fall-back, in case floating point errors
      // cause the usual search routine to fail. In the fall-back
      // routine, we just test every edge to find one straddling
      // the plane.
#if VOROPP_VERBOSE >= 1
      fputs("Bailed out of convex calculation\n", stderr);
#endif
      qw = 1;
      lw = 0;
      for (qp = 0; qp < p; qp++) {
        qw = m_test(qp, q);
        if (qw == 1) {
          // The point is inside the cutting space. Now
          // see if we can find a neighbor which isn't.
          for (us = 0; us < nu[qp]; us++) {
            lp = ed[qp][us];
            if (lp < qp) {
              lw = m_test(lp, l);
              if (lw != 1) {
                break;
              }
            }
          }
          if (us < nu[qp]) {
            up = qp;
            if (lw == 0) {
              complicated_setup = true;
            } else {
              complicated_setup = false;
              u = q;
              ls = ed[up][nu[up] + us];
            }
            break;
          }
        } else if (qw == -1) {
          // The point is outside the cutting space. See
          // if we can find a neighbor which isn't.
          for (ls = 0; ls < nu[qp]; ls++) {
            up = ed[qp][ls];
            if (up < qp) {
              uw = m_test(up, u);
              if (uw != -1) {
                break;
              }
            }
          }
          if (ls < nu[qp]) {
            if (uw == 0) {
              up = qp;
              complicated_setup = true;
            } else {
              complicated_setup = false;
              lp = qp;
              l = q;
              us = ed[lp][nu[lp] + ls];
            }
            break;
          }
        } else {
          // The point is in the plane, so we just
          // proceed with the complicated setup routine
          up = qp;
          complicated_setup = true;
          break;
        }
      }
      if (qp == p) {
        return qw == -1 ? true : false;
      }
    }

    // We're about to add the first point of the new facet. In either
    // routine, we have to add a point, so first check there's space for
    // it.
    if (p == current_vertices) {
      add_memory_vertices(vc);
    }

    if (complicated_setup) {
      // We want to be strict about reaching the conclusion that the
      // cell is entirely within the cutting plane. It's not enough
      // to find a vertex that has edges which are all inside or on
      // the plane. If the vertex has neighbors that are also on the
      // plane, we should check those too.
      if (!search_for_outside_edge(vc, up)) {
        return false;
      }

      // The search algorithm found a point which is on the cutting
      // plane. We leave that point in place, and create a new one at
      // the same location.
      pts[3 * p] = pts[3 * up];
      pts[3 * p + 1] = pts[3 * up + 1];
      pts[3 * p + 2] = pts[3 * up + 2];

      // Search for a collection of edges of the test vertex which
      // are outside of the cutting space. Begin by testing the
      // zeroth edge.
      i = 0;
      lp = *ed[up];
      lw = m_test(lp, l);
      if (lw != -1) {
        // The first edge is either inside the cutting space,
        // or lies within the cutting plane. Test the edges
        // sequentially until we find one that is outside.
        rp = lw;
        do {
          i++;

          // If we reached the last edge with no luck
          // then all of the vertices are inside
          // or on the plane, so the cell is completely
          // deleted
          if (i == nu[up]) {
            return false;
          }
          lp = ed[up][i];
          lw = m_test(lp, l);
        } while (lw != -1);
        j = i + 1;

        // We found an edge outside the cutting space. Keep
        // moving through these edges until we find one that's
        // inside or on the plane.
        while (j < nu[up]) {
          lp = ed[up][j];
          lw = m_test(lp, l);
          if (lw != -1) {
            break;
          }
          j++;
        }

        // Compute the number of edges for the new vertex. In
        // general it will be the number of outside edges
        // found, plus two. But we need to recognize the
        // special case when all but one edge is outside, and
        // the remaining one is on the plane. For that case we
        // have to reduce the edge count by one to prevent
        // doubling up.
        if (j == nu[up] && i == 1 && rp == 0) {
          nu[p] = nu[up];
          double_edge = true;
        } else {
          nu[p] = j - i + 2;
        }
        k = 1;

        // Add memory for the new vertex if needed, and
        // initialize
        while (nu[p] >= current_vertex_order) {
          add_memory_vorder(vc);
        }
        if (mec[nu[p]] == mem[nu[p]]) {
          add_memory(vc, nu[p], stackp2);
        }
        vc.n_set_pointer(p, nu[p]);
        ed[p] = mep[nu[p]] + ((nu[p] << 1) + 1) * mec[nu[p]]++;
        ed[p][nu[p] << 1] = p;

        // Copy the edges of the original vertex into the new
        // one. Delete the edges of the original vertex, and
        // update the relational table.
        us = cycle_down(i, up);
        while (i < j) {
          qp = ed[up][i];
          qs = ed[up][nu[up] + i];
          vc.n_copy(p, k, up, i);
          ed[p][k] = qp;
          ed[p][nu[p] + k] = qs;
          ed[qp][qs] = p;
          ed[qp][nu[qp] + qs] = k;
          ed[up][i] = -1;
          i++;
          k++;
        }
        qs = i == nu[up] ? 0 : i;
      } else {
        // In this case, the zeroth edge is outside the cutting
        // plane. Begin by searching backwards from the last
        // edge until we find an edge which isn't outside.
        i = nu[up] - 1;
        lp = ed[up][i];
        lw = m_test(lp, l);
        while (lw == -1) {
          i--;

          // If i reaches zero, then we have a point in
          // the plane all of whose edges are outside
          // the cutting space, so we just exit
          if (i == 0) {
            return true;
          }
          lp = ed[up][i];
          lw = m_test(lp, l);
        }

        // Now search forwards from zero
        j = 1;
        qp = ed[up][j];
        qw = m_test(qp, q);
        while (qw == -1) {
          j++;
          qp = ed[up][j];
          qw = m_test(qp, l);
        }

        // Compute the number of edges for the new vertex. In
        // general it will be the number of outside edges
        // found, plus two. But we need to recognize the
        // special case when all but one edge is outside, and
        // the remaining one is on the plane. For that case we
        // have to reduce the edge count by one to prevent
        // doubling up.
        if (i == j && qw == 0) {
          double_edge = true;
          nu[p] = nu[up];
        } else {
          nu[p] = nu[up] - i + j + 1;
        }

        // Add memory to store the vertex if it doesn't exist
        // already
        k = 1;
        while (nu[p] >= current_vertex_order) {
          add_memory_vorder(vc);
        }
        if (mec[nu[p]] == mem[nu[p]]) {
          add_memory(vc, nu[p], stackp2);
        }

        // Copy the edges of the original vertex into the new
        // one. Delete the edges of the original vertex, and
        // update the relational table.
        vc.n_set_pointer(p, nu[p]);
        ed[p] = mep[nu[p]] + ((nu[p] << 1) + 1) * mec[nu[p]]++;
        ed[p][nu[p] << 1] = p;
        us = i++;
        while (i < nu[up]) {
          qp = ed[up][i];
          qs = ed[up][nu[up] + i];
          vc.n_copy(p, k, up, i);
          ed[p][k] = qp;
          ed[p][nu[p] + k] = qs;
          ed[qp][qs] = p;
          ed[qp][nu[qp] + qs] = k;
          ed[up][i] = -1;
          i++;
          k++;
        }
        i = 0;
        while (i < j) {
          qp = ed[up][i];
          qs = ed[up][nu[up] + i];
          vc.n_copy(p, k, up, i);
          ed[p][k] = qp;
          ed[p][nu[p] + k] = qs;
          ed[qp][qs] = p;
          ed[qp][nu[qp] + qs] = k;
          ed[up][i] = -1;
          i++;
          k++;
        }
        qs = j;
      }
      if (!double_edge) {
        vc.n_copy(p, k, up, qs);
        vc.n_set(p, 0, p_id);
      } else {
        vc.n_copy(p, 0, up, qs);
      }

      // Add this point to the auxiliary delete stack
      if (stackp2 == stacke2) {
        add_memory_ds2(stackp2);
      }
      *(stackp2++) = up;

      // Look at the edges on either side of the group that was
      // detected. We're going to commence facet computation by
      // moving along one of them. We are going to end up coming back
      // along the other one.
      cs = k;
      qp = up;
      q = u;
      i = ed[up][us];
      us = ed[up][nu[up] + us];
      up = i;
      ed[qp][nu[qp] << 1] = -p;

    } else {
      // The search algorithm found an intersected edge between the
      // points lp and up. Create a new vertex between them which
      // lies on the cutting plane. Since u and l differ by at least
      // the tolerance, this division should never screw up.
      if (stackp == stacke) {
        add_memory_ds(stackp);
      }
      *(stackp++) = up;
      r = u / (u - l);
      l = 1 - r;
      pts[3 * p] = pts[3 * lp] * r + pts[3 * up] * l;
      pts[3 * p + 1] = pts[3 * lp + 1] * r + pts[3 * up + 1] * l;
      pts[3 * p + 2] = pts[3 * lp + 2] * r + pts[3 * up + 2] * l;

      // This point will always have three edges. Connect one of them
      // to lp.
      nu[p] = 3;
      if (mec[3] == mem[3]) {
        add_memory(vc, 3, stackp2);
      }
      vc.n_set_pointer(p, 3);
      vc.n_set(p, 0, p_id);
      vc.n_copy(p, 1, up, us);
      vc.n_copy(p, 2, lp, ls);
      ed[p] = mep[3] + 7 * mec[3]++;
      ed[p][6] = p;
      ed[up][us] = -1;
      ed[lp][ls] = p;
      ed[lp][nu[lp] + ls] = 1;
      ed[p][1] = lp;
      ed[p][nu[p] + 1] = ls;
      cs = 2;

      // Set the direction to move in
      qs = cycle_up(us, up);
      qp = up;
      q = u;
    }

    // When the code reaches here, we have initialized the first point, and
    // we have a direction for moving it to construct the rest of the facet
    cp = p;
    rp = p;
    p++;
    while (qp != up || qs != us) {
      // We're currently tracing round an intersected facet. Keep
      // moving around it until we find a point or edge which
      // intersects the plane.
      lp = ed[qp][qs];
      lw = m_test(lp, l);

      if (lw == 1) {
        // The point is still in the cutting space. Just add it
        // to the delete stack and keep moving.
        qs = cycle_up(ed[qp][nu[qp] + qs], lp);
        qp = lp;
        q = l;
        if (stackp == stacke) {
          add_memory_ds(stackp);
        }
        *(stackp++) = qp;

      } else if (lw == -1) {
        // The point is outside of the cutting space, so we've
        // found an intersected edge. Introduce a regular point
        // at the point of intersection. Connect it to the
        // point we just tested. Also connect it to the previous
        // new point in the facet we're constructing.
        if (p == current_vertices) {
          add_memory_vertices(vc);
        }
        r = q / (q - l);
        l = 1 - r;
        pts[3 * p] = pts[3 * lp] * r + pts[3 * qp] * l;
        pts[3 * p + 1] = pts[3 * lp + 1] * r + pts[3 * qp + 1] * l;
        pts[3 * p + 2] = pts[3 * lp + 2] * r + pts[3 * qp + 2] * l;
        nu[p] = 3;
        if (mec[3] == mem[3]) {
          add_memory(vc, 3, stackp2);
        }
        ls = ed[qp][qs + nu[qp]];
        vc.n_set_pointer(p, 3);
        vc.n_set(p, 0, p_id);
        vc.n_copy(p, 1, qp, qs);
        vc.n_copy(p, 2, lp, ls);
        ed[p] = mep[3] + 7 * mec[3]++;
        *ed[p] = cp;
        ed[p][1] = lp;
        ed[p][3] = cs;
        ed[p][4] = ls;
        ed[p][6] = p;
        ed[lp][ls] = p;
        ed[lp][nu[lp] + ls] = 1;
        ed[cp][cs] = p;
        ed[cp][nu[cp] + cs] = 0;
        ed[qp][qs] = -1;
        qs = cycle_up(qs, qp);
        cp = p++;
        cs = 2;
      } else {
        // We've found a point which is on the cutting plane.
        // We're going to introduce a new point right here, but
        // first we need to figure out the number of edges it
        // has.
        if (p == current_vertices) {
          add_memory_vertices(vc);
        }

        // If the previous vertex detected a double edge, our
        // new vertex will have one less edge.
        k = double_edge ? 0 : 1;
        qs = ed[qp][nu[qp] + qs];
        qp = lp;
        iqs = qs;

        // Start testing the edges of the current point until
        // we find one which isn't outside the cutting space
        do {
          k++;
          qs = cycle_up(qs, qp);
          lp = ed[qp][qs];
          lw = m_test(lp, l);
        } while (lw == -1);

        // Now we need to find out whether this marginal vertex
        // we are on has been visited before, because if that's
        // the case, we need to add vertices to the existing
        // new vertex, rather than creating a fresh one. We also
        // need to figure out whether we're in a case where we
        // might be creating a duplicate edge.
        j = -ed[qp][nu[qp] << 1];
        if (qp == up && qs == us) {
          // If we're heading into the final part of the
          // new facet, then we never worry about the
          // duplicate edge calculation.
          new_double_edge = false;
          if (j > 0) {
            k += nu[j];
          }
        } else {
          if (j > 0) {
            // This vertex was visited before, so
            // count those vertices to the ones we
            // already have.
            k += nu[j];

            // The only time when we might make a
            // duplicate edge is if the point we're
            // going to move to next is also a
            // marginal point, so test for that
            // first.
            if (lw == 0) {
              // Now see whether this marginal point
              // has been visited before.
              i = -ed[lp][nu[lp] << 1];
              if (i > 0) {
                // Now see if the last edge of that other
                // marginal point actually ends up here.
                if (ed[i][nu[i] - 1] == j) {
                  new_double_edge = true;
                  k -= 1;
                } else {
                  new_double_edge = false;
                }
              } else {
                // That marginal point hasn't been visited
                // before, so we probably don't have to worry
                // about duplicate edges, except in the
                // case when that's the way into the end
                // of the facet, because that way always creates
                // an edge.
                if (j == rp && lp == up && ed[qp][nu[qp] + qs] == us) {
                  new_double_edge = true;
                  k -= 1;
                } else {
                  new_double_edge = false;
                }
              }
            } else {
              new_double_edge = false;
            }
          } else {
            // The vertex hasn't been visited
            // before, but let's see if it's
            // marginal
            if (lw == 0) {
              // If it is, we need to check
              // for the case that it's a
              // small branch, and that we're
              // heading right back to where
              // we came from
              i = -ed[lp][nu[lp] << 1];
              if (i == cp) {
                new_double_edge = true;
                k -= 1;
              } else {
                new_double_edge = false;
              }
            } else {
              new_double_edge = false;
            }
          }
        }

        // k now holds the number of edges of the new vertex
        // we are forming. Add memory for it if it doesn't exist
        // already.
        while (k >= current_vertex_order) {
          add_memory_vorder(vc);
        }
        if (mec[k] == mem[k]) {
          add_memory(vc, k, stackp2);
        }

        // Now create a new vertex with order k, or augment
        // the existing one
        if (j > 0) {
          // If we're augmenting a vertex but we don't
          // actually need any more edges, just skip this
          // routine to avoid memory confusion
          if (nu[j] != k) {
            // Allocate memory and copy the edges
            // of the previous instance into it
            vc.n_set_aux1(k);
            edp = mep[k] + ((k << 1) + 1) * mec[k]++;
            i = 0;
            while (i < nu[j]) {
              vc.n_copy_aux1(j, i);
              edp[i] = ed[j][i];
              edp[k + i] = ed[j][nu[j] + i];
              i++;
            }
            edp[k << 1] = j;

            // Remove the previous instance with
            // fewer vertices from the memory
            // structure
            edd = mep[nu[j]] + ((nu[j] << 1) + 1) * --mec[nu[j]];
            if (edd != ed[j]) {
              for (lw = 0; lw <= (nu[j] << 1); lw++) {
                ed[j][lw] = edd[lw];
              }
              vc.n_set_aux2_copy(j, nu[j]);
              vc.n_copy_pointer(edd[nu[j] << 1], j);
              ed[edd[nu[j] << 1]] = ed[j];
            }
            vc.n_set_to_aux1(j);
            ed[j] = edp;
          } else {
            i = nu[j];
          }
        } else {
          // Allocate a new vertex of order k
          vc.n_set_pointer(p, k);
          ed[p] = mep[k] + ((k << 1) + 1) * mec[k]++;
          ed[p][k << 1] = p;
          if (stackp2 == stacke2) {
            add_memory_ds2(stackp2);
          }
          *(stackp2++) = qp;
          pts[3 * p] = pts[3 * qp];
          pts[3 * p + 1] = pts[3 * qp + 1];
          pts[3 * p + 2] = pts[3 * qp + 2];
          ed[qp][nu[qp] << 1] = -p;
          j = p++;
          i = 0;
        }
        nu[j] = k;

        // Unless the previous case was a double edge, connect
        // the first available edge of the new vertex to the
        // last one in the facet
        if (!double_edge) {
          ed[j][i] = cp;
          ed[j][nu[j] + i] = cs;
          vc.n_set(j, i, p_id);
          ed[cp][cs] = j;
          ed[cp][nu[cp] + cs] = i;
          i++;
        }

        // Copy in the edges of the underlying vertex,
        // and do one less if this was a double edge
        qs = iqs;
        while (i < (new_double_edge ? k : k - 1)) {
          qs = cycle_up(qs, qp);
          lp = ed[qp][qs];
          ls = ed[qp][nu[qp] + qs];
          vc.n_copy(j, i, qp, qs);
          ed[j][i] = lp;
          ed[j][nu[j] + i] = ls;
          ed[lp][ls] = j;
          ed[lp][nu[lp] + ls] = i;
          ed[qp][qs] = -1;
          i++;
        }
        qs = cycle_up(qs, qp);
        cs = i;
        cp = j;
        vc.n_copy(j, new_double_edge ? 0 : cs, qp, qs);

        // Update the double_edge flag, to pass it
        // to the next instance of this routine
        double_edge = new_double_edge;
      }
    }

    // Connect the final created vertex to the initial one
    ed[cp][cs] = rp;
    *ed[rp] = cp;
    ed[cp][nu[cp] + cs] = 0;
    ed[rp][nu[rp]] = cs;

    // Delete points: first, remove any duplicates
    dsp = ds;
    while (dsp < stackp) {
      j = *dsp;
      if (ed[j][nu[j]] != -1) {
        ed[j][nu[j]] = -1;
        dsp++;
      } else {
        *dsp = *(--stackp);
      }
    }

    // Add the points in the auxiliary delete stack,
    // and reset their back pointers
    for (dsp = ds2; dsp < stackp2; dsp++) {
      j = *dsp;
      ed[j][nu[j] << 1] = j;
      if (ed[j][nu[j]] != -1) {
        ed[j][nu[j]] = -1;
        if (stackp == stacke) {
          add_memory_ds(stackp);
        }
        *(stackp++) = j;
      }
    }

    // Scan connections and add in extras
    for (dsp = ds; dsp < stackp; dsp++) {
      cp = *dsp;
      for (edp = ed[cp]; edp < ed[cp] + nu[cp]; edp++) {
        qp = *edp;
        if (qp != -1 && ed[qp][nu[qp]] != -1) {
          if (stackp == stacke) {
            const int dis = stackp - dsp;
            add_memory_ds(stackp);
            dsp = ds + dis;
          }
          *(stackp++) = qp;
          ed[qp][nu[qp]] = -1;
        }
      }
    }
    up = 0;

    // Delete them from the array structure
    while (stackp > ds) {
      --p;
      while (ed[p][nu[p]] == -1) {
        j = nu[p];
        edp = ed[p];
        edd = (mep[j] + ((j << 1) + 1) * --mec[j]);
        while (edp < ed[p] + (j << 1) + 1) {
          *(edp++) = *(edd++);
        }
        vc.n_set_aux2_copy(p, j);
        vc.n_copy_pointer(ed[p][(j << 1)], p);
        ed[ed[p][(j << 1)]] = ed[p];
        --p;
      }
      up = *(--stackp);
      if (up < p) {
        // Vertex management
        pts[3 * up] = pts[3 * p];
        pts[3 * up + 1] = pts[3 * p + 1];
        pts[3 * up + 2] = pts[3 * p + 2];

        // Memory management
        j = nu[up];
        edp = ed[up];
        edd = (mep[j] + ((j << 1) + 1) * --mec[j]);
        while (edp < ed[up] + (j << 1) + 1) {
          *(edp++) = *(edd++);
        }
        vc.n_set_aux2_copy(up, j);
        vc.n_copy_pointer(ed[up][j << 1], up);
        vc.n_copy_pointer(up, p);
        ed[ed[up][j << 1]] = ed[up];

        // Edge management
        ed[up] = ed[p];
        nu[up] = nu[p];
        for (i = 0; i < nu[up]; i++) {
          ed[ed[up][i]][ed[up][nu[up] + i]] = up;
        }
        ed[up][nu[up] << 1] = up;
      } else {
        up = p++;
      }
    }

    // Check for any vertices of zero order
    if (*mec > 0) {
      voro_fatal_error("Zero order vertex formed", VOROPP_INTERNAL_ERROR);
    }

    // Collapse any order 2 vertices and exit
    return collapse_order2(vc);
  }

/** During the creation of a new facet in the plane routine, it is possible
   * that some order two vertices may arise. This routine removes them.
   * Suppose an order two vertex joins c and d. If there's a edge between
   * c and d already, then the order two vertex is just removed; otherwise,
   * the order two vertex is removed and c and d are joined together directly.
   * It is possible this process will create order two or order one vertices,
   * and the routine is continually run until all of them are removed.
   * \return False if the vertex removal was unsuccessful, indicative of the cell
   *         reducing to zero volume and disappearing; true if the vertex removal
   *         was successful. */
  template <class vc_class> inline bool voronoicell_base::collapse_order2(vc_class &vc) {
    if (!collapse_order1(vc)) {
      return false;
    }
    int a;
    int b;
    int i;
    int j;
    int k;
    int l;
    while (mec[2] > 0) {
      // Pick a order 2 vertex and read in its edges
      i = --mec[2];
      j = mep[2][5 * i];
      k = mep[2][5 * i + 1];
      if (j == k) {
#if VOROPP_VERBOSE >= 1
        fputs("Order two vertex joins itself", stderr);
#endif
        return false;
      }

      // Scan the edges of j to see if joins k
      for (l = 0; l < nu[j]; l++) {
        if (ed[j][l] == k) {
          break;
        }
      }

      // If j doesn't already join k, join them together.
      // Otherwise delete the connection to the current
      // vertex from j and k.
      a = mep[2][5 * i + 2];
      b = mep[2][5 * i + 3];
      i = mep[2][5 * i + 4];
      if (l == nu[j]) {
        ed[j][a] = k;
        ed[k][b] = j;
        ed[j][nu[j] + a] = b;
        ed[k][nu[k] + b] = a;
      } else {
        if (!delete_connection(vc, j, a, false)) {
          return false;
        }
        if (!delete_connection(vc, k, b, true)) {
          return false;
        }
      }

      // Compact the memory
      --p;
      if (up == i) {
        up = 0;
      }
      if (p != i) {
        if (up == p) {
          up = i;
        }
        pts[3 * i] = pts[3 * p];
        pts[3 * i + 1] = pts[3 * p + 1];
        pts[3 * i + 2] = pts[3 * p + 2];
        for (k = 0; k < nu[p]; k++) {
          ed[ed[p][k]][ed[p][nu[p] + k]] = i;
        }
        vc.n_copy_pointer(i, p);
        ed[i] = ed[p];
        nu[i] = nu[p];
        ed[i][nu[i] << 1] = i;
      }

      // Collapse any order 1 vertices if they were created
      if (!collapse_order1(vc)) {
        return false;
      }
    }
    return true;
  }

/** Order one vertices can potentially be created during the order two collapse
   * routine. This routine keeps removing them until there are none left.
   * \return False if the vertex removal was unsuccessful, indicative of the cell
   *         having zero volume and disappearing; true if the vertex removal was
   *         successful. */
  template <class vc_class> inline bool voronoicell_base::collapse_order1(vc_class &vc) {
    int i;
    int j;
    int k;
    while (mec[1] > 0) {
      up = 0;
#if VOROPP_VERBOSE >= 1
      fputs("Order one collapse\n", stderr);
#endif
      i = --mec[1];
      j = mep[1][3 * i];
      k = mep[1][3 * i + 1];
      i = mep[1][3 * i + 2];
      if (!delete_connection(vc, j, k, false)) {
        return false;
      }
      --p;
      if (up == i) {
        up = 0;
      }
      if (p != i) {
        if (up == p) {
          up = i;
        }
        pts[3 * i] = pts[3 * p];
        pts[3 * i + 1] = pts[3 * p + 1];
        pts[3 * i + 2] = pts[3 * p + 2];
        for (k = 0; k < nu[p]; k++) {
          ed[ed[p][k]][ed[p][nu[p] + k]] = i;
        }
        vc.n_copy_pointer(i, p);
        ed[i] = ed[p];
        nu[i] = nu[p];
        ed[i][nu[i] << 1] = i;
      }
    }
    return true;
  }

/** This routine deletes the kth edge of vertex j and reorganizes the memory.
   * If the neighbor computation is enabled, we also have to supply an handedness
   * flag to decide whether to preserve the plane on the left or right of the
   * connection.
   * \return False if a zero order vertex was formed, indicative of the cell
   *         disappearing; true if the vertex removal was successful. */
  template <class vc_class> inline bool voronoicell_base::delete_connection(vc_class &vc, int j, int k, bool hand) {
    const int q = hand ? k : cycle_up(k, j);
    int i = nu[j] - 1, l, *edp, *edd, m;
#if VOROPP_VERBOSE >= 1
    if (i < 1) {
      fputs("Zero order vertex formed\n", stderr);
      return false;
    }
#endif
    if (mec[i] == mem[i]) {
      add_memory(vc, i, ds2);
    }
    vc.n_set_aux1(i);
    for (l = 0; l < q; l++) {
      vc.n_copy_aux1(j, l);
    }
    while (l < i) {
      vc.n_copy_aux1_shift(j, l);
      l++;
    }
    edp = mep[i] + ((i << 1) + 1) * mec[i]++;
    edp[i << 1] = j;
    for (l = 0; l < k; l++) {
      edp[l] = ed[j][l];
      edp[l + i] = ed[j][l + nu[j]];
    }
    while (l < i) {
      m = ed[j][l + 1];
      edp[l] = m;
      k = ed[j][l + nu[j] + 1];
      edp[l + i] = k;
      ed[m][nu[m] + k]--;
      l++;
    }

    edd = mep[nu[j]] + ((nu[j] << 1) + 1) * --mec[nu[j]];
    for (l = 0; l <= (nu[j] << 1); l++) {
      ed[j][l] = edd[l];
    }
    vc.n_set_aux2_copy(j, nu[j]);
    vc.n_set_to_aux2(edd[nu[j] << 1]);
    vc.n_set_to_aux1(j);
    ed[edd[nu[j] << 1]] = edd;
    ed[j] = edp;
    nu[j] = i;
    return true;
  }

/** Calculates the volume of the Voronoi cell, by decomposing the cell into
   * tetrahedra extending outward from the zeroth vertex, whose volumes are
   * evaluated using a scalar triple product.
   * \return A floating point number holding the calculated volume. */
  double voronoicell_base::volume() {
    const double fe = 1 / 48.0;
    double vol = 0;
    int i;
    int j;
    int k;
    int l;
    int m;
    int n;
    double ux;
    double uy;
    double uz;
    double vx;
    double vy;
    double vz;
    double wx;
    double wy;
    double wz;
    for (i = 1; i < p; i++) {
      ux = *pts - pts[3 * i];
      uy = pts[1] - pts[3 * i + 1];
      uz = pts[2] - pts[3 * i + 2];
      for (j = 0; j < nu[i]; j++) {
        k = ed[i][j];
        if (k >= 0) {
          ed[i][j] = -1 - k;
          l = cycle_up(ed[i][nu[i] + j], k);
          vx = pts[3 * k] - *pts;
          vy = pts[3 * k + 1] - pts[1];
          vz = pts[3 * k + 2] - pts[2];
          m = ed[k][l];
          ed[k][l] = -1 - m;
          while (m != i) {
            n = cycle_up(ed[k][nu[k] + l], m);
            wx = pts[3 * m] - *pts;
            wy = pts[3 * m + 1] - pts[1];
            wz = pts[3 * m + 2] - pts[2];
            vol += ux * vy * wz + uy * vz * wx + uz * vx * wy - uz * vy * wx - uy * vx * wz - ux * vz * wy;
            k = m;
            l = n;
            vx = wx;
            vy = wy;
            vz = wz;
            m = ed[k][l];
            ed[k][l] = -1 - m;
          }
        }
      }
    }
    reset_edges();
    return vol * fe;
  }

/** Calculates the areas of each face of the Voronoi cell and prints the
   * results to an output stream.
   * \param[out] v the vector to store the results in. */
  void voronoicell_base::face_areas(std::vector<double> &v) {
    double area;
    v.clear();
    int i;
    int j;
    int k;
    int l;
    int m;
    int n;
    double ux;
    double uy;
    double uz;
    double vx;
    double vy;
    double vz;
    double wx;
    double wy;
    double wz;
    for (i = 1; i < p; i++) {
      for (j = 0; j < nu[i]; j++) {
        k = ed[i][j];
        if (k >= 0) {
          area = 0;
          ed[i][j] = -1 - k;
          l = cycle_up(ed[i][nu[i] + j], k);
          m = ed[k][l];
          ed[k][l] = -1 - m;
          while (m != i) {
            n = cycle_up(ed[k][nu[k] + l], m);
            ux = pts[3 * k] - pts[3 * i];
            uy = pts[3 * k + 1] - pts[3 * i + 1];
            uz = pts[3 * k + 2] - pts[3 * i + 2];
            vx = pts[3 * m] - pts[3 * i];
            vy = pts[3 * m + 1] - pts[3 * i + 1];
            vz = pts[3 * m + 2] - pts[3 * i + 2];
            wx = uy * vz - uz * vy;
            wy = uz * vx - ux * vz;
            wz = ux * vy - uy * vx;
            area += sqrt(wx * wx + wy * wy + wz * wz);
            k = m;
            l = n;
            m = ed[k][l];
            ed[k][l] = -1 - m;
          }
          v.push_back(0.125 * area);
        }
      }
    }
    reset_edges();
  }

/** Calculates the total surface area of the Voronoi cell.
   * \return The computed area. */
  double voronoicell_base::surface_area() {
    double area = 0;
    int i;
    int j;
    int k;
    int l;
    int m;
    int n;
    double ux;
    double uy;
    double uz;
    double vx;
    double vy;
    double vz;
    double wx;
    double wy;
    double wz;
    for (i = 1; i < p; i++) {
      for (j = 0; j < nu[i]; j++) {
        k = ed[i][j];
        if (k >= 0) {
          ed[i][j] = -1 - k;
          l = cycle_up(ed[i][nu[i] + j], k);
          m = ed[k][l];
          ed[k][l] = -1 - m;
          while (m != i) {
            n = cycle_up(ed[k][nu[k] + l], m);
            ux = pts[3 * k] - pts[3 * i];
            uy = pts[3 * k + 1] - pts[3 * i + 1];
            uz = pts[3 * k + 2] - pts[3 * i + 2];
            vx = pts[3 * m] - pts[3 * i];
            vy = pts[3 * m + 1] - pts[3 * i + 1];
            vz = pts[3 * m + 2] - pts[3 * i + 2];
            wx = uy * vz - uz * vy;
            wy = uz * vx - ux * vz;
            wz = ux * vy - uy * vx;
            area += sqrt(wx * wx + wy * wy + wz * wz);
            k = m;
            l = n;
            m = ed[k][l];
            ed[k][l] = -1 - m;
          }
        }
      }
    }
    reset_edges();
    return 0.125 * area;
  }

/** Calculates the centroid of the Voronoi cell, by decomposing the cell into
   * tetrahedra extending outward from the zeroth vertex.
   * \param[out] (cx,cy,cz) references to floating point numbers in which to
   *                        pass back the centroid vector. */
  void voronoicell_base::centroid(double &cx, double &cy, double &cz) {
    double tvol;
    double vol = 0;
    cx = cy = cz = 0;
    int i;
    int j;
    int k;
    int l;
    int m;
    int n;
    double ux;
    double uy;
    double uz;
    double vx;
    double vy;
    double vz;
    double wx;
    double wy;
    double wz;
    for (i = 1; i < p; i++) {
      ux = *pts - pts[3 * i];
      uy = pts[1] - pts[3 * i + 1];
      uz = pts[2] - pts[3 * i + 2];
      for (j = 0; j < nu[i]; j++) {
        k = ed[i][j];
        if (k >= 0) {
          ed[i][j] = -1 - k;
          l = cycle_up(ed[i][nu[i] + j], k);
          vx = pts[3 * k] - *pts;
          vy = pts[3 * k + 1] - pts[1];
          vz = pts[3 * k + 2] - pts[2];
          m = ed[k][l];
          ed[k][l] = -1 - m;
          while (m != i) {
            n = cycle_up(ed[k][nu[k] + l], m);
            wx = pts[3 * m] - *pts;
            wy = pts[3 * m + 1] - pts[1];
            wz = pts[3 * m + 2] - pts[2];
            tvol = ux * vy * wz + uy * vz * wx + uz * vx * wy - uz * vy * wx - uy * vx * wz - ux * vz * wy;
            vol += tvol;
            cx += (wx + vx - ux) * tvol;
            cy += (wy + vy - uy) * tvol;
            cz += (wz + vz - uz) * tvol;
            k = m;
            l = n;
            vx = wx;
            vy = wy;
            vz = wz;
            m = ed[k][l];
            ed[k][l] = -1 - m;
          }
        }
      }
    }
    reset_edges();
    if (vol > tolerance_sq) {
      vol = 0.125 / vol;
      cx = cx * vol + 0.5 * (*pts);
      cy = cy * vol + 0.5 * pts[1];
      cz = cz * vol + 0.5 * pts[2];
    } else {
      cx = cy = cz = 0;
    }
  }

/** Computes the maximum radius squared of a vertex from the center of the
   * cell. It can be used to determine when enough particles have been testing an
   * all planes that could cut the cell have been considered.
   * \return The maximum radius squared of a vertex.*/
  double voronoicell_base::max_radius_squared() const {
    double r;
    double s;
    double *ptsp = pts + 3;
    double *ptse = pts + 3 * p;
    r = *pts * (*pts) + pts[1] * pts[1] + pts[2] * pts[2];
    while (ptsp < ptse) {
      s = *ptsp * (*ptsp);
      ptsp++;
      s += *ptsp * (*ptsp);
      ptsp++;
      s += *ptsp * (*ptsp);
      ptsp++;
      if (s > r) {
        r = s;
      }
    }
    return r;
  }

/** Calculates the total edge distance of the Voronoi cell.
   * \return A floating point number holding the calculated distance. */
  double voronoicell_base::total_edge_distance() const {
    int i;
    int j;
    int k;
    double dis = 0;
    double dx;
    double dy;
    double dz;
    for (i = 0; i < p - 1; i++) {
      for (j = 0; j < nu[i]; j++) {
        k = ed[i][j];
        if (k > i) {
          dx = pts[3 * k] - pts[3 * i];
          dy = pts[3 * k + 1] - pts[3 * i + 1];
          dz = pts[3 * k + 2] - pts[3 * i + 2];
          dis += sqrt(dx * dx + dy * dy + dz * dz);
        }
      }
    }
    return 0.5 * dis;
  }

/** Outputs the edges of the Voronoi cell in POV-Ray format to an open file
   * stream, displacing the cell by given vector.
   * \param[in] (x,y,z) a displacement vector to be added to the cell's position.
   * \param[in] fp a file handle to write to. */
  void voronoicell_base::draw_pov(double x, double y, double z, FILE *fp) const {
    int i;
    int j;
    int k;
    double *ptsp = pts;
    double *pt2;
    char posbuf1[128];
    char posbuf2[128];
    for (i = 0; i < p; i++, ptsp += 3) {
      snprintf(posbuf1, sizeof(posbuf1), "%g,%g,%g", x + *ptsp * 0.5, y + ptsp[1] * 0.5, z + ptsp[2] * 0.5);
      fprintf(fp, "sphere{<%s>,r}\n", posbuf1);
      for (j = 0; j < nu[i]; j++) {
        k = ed[i][j];
        if (k < i) {
          pt2 = pts + 3 * k;
          snprintf(posbuf2, sizeof(posbuf2), "%g,%g,%g", x + *pt2 * 0.5, y + 0.5 * pt2[1], z + 0.5 * pt2[2]);
          if (strcmp(posbuf1, posbuf2) != 0) {
            fprintf(fp, "cylinder{<%s>,<%s>,r}\n", posbuf1, posbuf2);
          }
        }
      }
    }
  }

/** Outputs the edges of the Voronoi cell in gnuplot format to an output stream.
   * \param[in] (x,y,z) a displacement vector to be added to the cell's position.
   * \param[in] fp a file handle to write to. */
  void voronoicell_base::draw_gnuplot(double x, double y, double z, FILE *fp) {
    int i;
    int j;
    int k;
    int l;
    int m;
    for (i = 1; i < p; i++) {
      for (j = 0; j < nu[i]; j++) {
        k = ed[i][j];
        if (k >= 0) {
          fprintf(fp, "%g %g %g\n", x + 0.5 * pts[3 * i], y + 0.5 * pts[3 * i + 1], z + 0.5 * pts[3 * i + 2]);
          l = i;
          m = j;
          do {
            ed[k][ed[l][nu[l] + m]] = -1 - l;
            ed[l][m] = -1 - k;
            l = k;
            fprintf(fp, "%g %g %g\n", x + 0.5 * pts[3 * k], y + 0.5 * pts[3 * k + 1], z + 0.5 * pts[3 * k + 2]);
          } while (search_edge(l, m, k));
          fputs("\n\n", fp);
        }
      }
    }
    reset_edges();
  }

  inline bool voronoicell_base::search_edge(int l, int &m, int &k) const {
    for (m = 0; m < nu[l]; m++) {
      k = ed[l][m];
      if (k >= 0) {
        return true;
      }
    }
    return false;
  }

/// AFLOW //HE20220908
  void voronoicell_base::get_vertex_facets(double x, double y, double z, vector<xvector<double>> &points, vector<vector<uint>> &facets) {
    int i;
    int j;
    int k;
    int l;
    int m;
    int n;
    double *ptsp = pts;
    points.clear();
    facets.clear();
    for (i = 0; i < p; i++, ptsp += 3) {
      points.push_back({x + *ptsp * 0.5, y + ptsp[1] * 0.5, z + ptsp[2] * 0.5});
    }
//    fprintf(fp,"}\nface_indices {\n%d\n",(p-2)<<1);
    for (i = 1; i < p; i++) {
      for (j = 0; j < nu[i]; j++) {
        k = ed[i][j];
        if (k >= 0) {
          ed[i][j] = -1 - k;
          l = cycle_up(ed[i][nu[i] + j], k);
          m = ed[k][l];
          ed[k][l] = -1 - m;
          while (m != i) {
            n = cycle_up(ed[k][nu[k] + l], m);
            facets.push_back({(uint) i, (uint) k, (uint) m});
//            fprintf(fp,",<%d,%d,%d>\n",i,k,m);
            k = m;
            l = n;
            m = ed[k][l];
            ed[k][l] = -1 - m;
          }
        }
      }
    }
    reset_edges();
  }

  /// Aflow //HE20220912
  double voronoicell_base::pbc_distance_min(const xmatrix<double> &lattice, const xvector<double> &p1, const xvector<double> &p2) {
    xvector<double> p_select;
    xvector<double> distance(27, 1);
    uint ni = 1;
    for (const int x : {-1, 0, 1}) {
      for (const int y : {-1, 0, 1}) {
        for (const int z : {-1, 0, 1}) {
          p_select = p2 + x * lattice(1) + y * lattice(2) + z * lattice(3);
          distance[ni] = aurostd::modulus(p1 - p_select);
          ni++;
        }
      }
    }
    return aurostd::min(distance);
  }

  double voronoicell_base::pbc_distance(const xvector<double> &normal, const xmatrix<double> &lattice, const xvector<double> &p1, const xvector<double> &p2) {
    xvector<double> p_select;
    xvector<double> fit(27, 1);
    xvector<double> distance(27, 1);
    vector<double> possible_distances;
    xvector<double> compare_normal;
    const double tolerance = 1E-6;
    uint ni = 1;
    for (const int x : {-1, 0, 1}) {
      for (const int y : {-1, 0, 1}) {
        for (const int z : {-1, 0, 1}) {
          p_select = p2 + x * lattice(1) + y * lattice(2) + z * lattice(3);
          distance[ni] = aurostd::modulus(p1 - p_select);
          compare_normal = (p_select - p1) / distance[ni];
          fit[ni] = aurostd::modulus(normal - compare_normal);
          ni++;
        }
      }
    }
    aurostd::quicksort2(27, fit, distance);
    const uint final_ni = 0;
    const double scale = min(distance);
    const double base_fit = fit[1] / scale;
    for (ni = 1; ni <= 27; ni++) {
      fit[ni] = fit[ni] / scale;
      if (std::abs(fit[ni] - base_fit) < tolerance) {
        possible_distances.push_back(distance[ni]);
      } else {
        break;
      }
    }
//    cout << "fit " << fit[1] << ", " << fit[2] << endl;
//
//    if (distance[1] > 20){
//      aurostd::x3DWriter w;
//      w.addLatticeBox(lattice, 0.1);
//      w.addSphere(p1, 0.3, "bm_green");
//      w.addSphere(p2, 0.3, "bm_green");
//      w.addOpenCylinder(p1, p1+normal*30, 0.05);
//      uint ni = 1;
//      for (int x : {-1, 0, 1}) {
//        for (int y: {-1, 0, 1}) {
//          for (int z: {-1, 0, 1}) {
//            p_select = p2 + x * lattice(1) + y * lattice(2) + z * lattice(3);
//            distance[ni] =aurostd::modulus(p1 - p_select);
//            compare_normal = (p_select-p1)/distance[ni];
//            w.addSphere(p_select, 0.3, "bm_blue");
//            fit[ni] = round(aurostd::modulus(normal-compare_normal),8);
//            ni++;
//          }
//        }
//      }
//      aurostd::quicksort2(27, fit, distance);
//      cout << "fit *dis " << fit[1] << ", " << fit[27] << endl;
//      cout << "dis " << distance[1] << endl;
//
//
//      aurostd::string2file(w.toHTML(),
//                           "../testing/soliquidy/debug.html");
//      exit(0);
//    }
    return aurostd::min(possible_distances);
  }
/// AFLOW //HE20220910
  void voronoicell_base::fill_laplacian(const int p_id, const vector<xvector<double>> &vertices, const xmatrix<double> &lattice, xmatrix<double> &laplacian) {
    vector<double> areas;
    vector<double> normal_list;
    xvector<double> normal;
    vector<int> neighbors_ids;
    double dist = 0.0;
    const double dist_min = 0.0;
    const int k = 0;
    double temp = 0.0;
    int n_id = 0;
//    cout << p_id << endl;
    face_areas(areas);
    neighbors(neighbors_ids);
    normals(normal_list);
    for (size_t i = 0; i < neighbors_ids.size(); i++) {
      n_id = neighbors_ids[i];
      if (n_id != p_id and n_id >= 0) {
//        cout << p_id  << "|" << n_id << endl;
        normal = {normal_list[i * 3], normal_list[i * 3 + 1], normal_list[i * 3 + 2]};
//        cout << normal << endl;
//        cout << (vertices[p_id] - vertices[n_id])/aurostd::modulus((vertices[p_id] - vertices[n_id])) << endl;

        dist = pbc_distance(normal, lattice, vertices[p_id], vertices[n_id]);
//          dist_min = pbc_distance_min(lattice, vertices[p_id], vertices[n_id]);
//          cout << dist << " | " << areas[i] << " | " << areas[i] / (2*dist) << endl;
//          if (abs(dist-dist_min) > 0.1) cout << "============================" << endl;
//        dist = aurostd::distance(vertices[p_id], vertices[n_id]);
//        dist = aurostd::modulus(vertices[p_id] - vertices[n_id]); //<< "| pbc: " << dist << endl;
//        cout << "-----" << endl;
        temp = areas[i] / (2 * dist);
        laplacian[p_id][n_id] += -temp;
        laplacian[p_id][p_id] += temp;
      }
    }
  }

/** Outputs the Voronoi cell in the POV mesh2 format, described in section
   * 1.3.2.2 of the POV-Ray documentation. The mesh2 output consists of a list of
   * vertex vectors, followed by a list of triangular faces. The routine also
   * makes use of the optional inside_vector specification, which makes the mesh
   * object solid, so the the POV-Ray Constructive Solid Geometry (CSG) can be
   * applied.
   * \param[in] (x,y,z) a displacement vector to be added to the cell's position.
   * \param[in] fp a file handle to write to. */
  void voronoicell_base::draw_pov_mesh(double x, double y, double z, FILE *fp) {
    int i;
    int j;
    int k;
    int l;
    int m;
    int n;
    double *ptsp = pts;
    fprintf(fp, "mesh2 {\nvertex_vectors {\n%d\n", p);
    for (i = 0; i < p; i++, ptsp += 3) {
      fprintf(fp, ",<%g,%g,%g>\n", x + *ptsp * 0.5, y + ptsp[1] * 0.5, z + ptsp[2] * 0.5);
    }
    fprintf(fp, "}\nface_indices {\n%d\n", (p - 2) << 1);
    for (i = 1; i < p; i++) {
      for (j = 0; j < nu[i]; j++) {
        k = ed[i][j];
        if (k >= 0) {
          ed[i][j] = -1 - k;
          l = cycle_up(ed[i][nu[i] + j], k);
          m = ed[k][l];
          ed[k][l] = -1 - m;
          while (m != i) {
            n = cycle_up(ed[k][nu[k] + l], m);
            fprintf(fp, ",<%d,%d,%d>\n", i, k, m);
            k = m;
            l = n;
            m = ed[k][l];
            ed[k][l] = -1 - m;
          }
        }
      }
    }
    fputs("}\ninside_vector <0,0,1>\n}\n", fp);
    reset_edges();
  }

/** Several routines in the class that gather cell-based statistics internally
   * track their progress by flipping edges to negative so that they know what
   * parts of the cell have already been tested. This function resets them back
   * to positive. When it is called, it assumes that every edge in the routine
   * should have already been flipped to negative, and it bails out with an
   * internal error if it encounters a positive edge. */
  inline void voronoicell_base::reset_edges() const {
    int i;
    int j;
    for (i = 0; i < p; i++) {
      for (j = 0; j < nu[i]; j++) {
        if (ed[i][j] >= 0) {
          voro_fatal_error("Edge reset routine found a previously untested edge", VOROPP_INTERNAL_ERROR);
        }
        ed[i][j] = -1 - ed[i][j];
      }
    }
  }

/** Checks to see if a given vertex is inside, outside or within the test
   * plane. If the point is far away from the test plane, the routine immediately
   * returns whether it is inside or outside. If the routine is close the the
   * plane and within the specified tolerance, then the special check_marginal()
   * routine is called.
   * \param[in] n the vertex to test.
   * \param[out] ans the result of the scalar product used in evaluating the
   *                 location of the point.
   * \return -1 if the point is inside the plane, 1 if the point is outside the
   *         plane, or 0 if the point is within the plane. */
  inline int voronoicell_base::m_test(int n, double &ans) {
    double *pp = pts + n + (n << 1);
    ans = *(pp++) * px;
    ans += *(pp++) * py;
    ans += *pp * pz - prsq;
    if (ans < -tolerance2) {
      return -1;
    } else if (ans > tolerance2) {
      return 1;
    }
    return check_marginal(n, ans);
  }

/** Checks to see if a given vertex is inside, outside or within the test
   * plane, for the case when the point has been detected to be very close to the
   * plane. The routine ensures that the returned results are always consistent
   * with previous tests, by keeping a table of any marginal results. The routine
   * first sees if the vertex is in the table, and if it finds a previously
   * computed result it uses that. Otherwise, it computes a result for this
   * vertex and adds it the table.
   * \param[in] n the vertex to test.
   * \param[in] ans the result of the scalar product used in evaluating
   *                the location of the point.
   * \return -1 if the point is inside the plane, 1 if the point is outside the
   *         plane, or 0 if the point is within the plane. */
  int voronoicell_base::check_marginal(int n, double &ans) {
    int i;
    for (i = 0; i < n_marg; i += 2) {
      if (marg[i] == n) {
        return marg[i + 1];
      }
    }
    if (n_marg == current_marginal) {
      current_marginal <<= 1;
      if (current_marginal > max_marginal) {
        voro_fatal_error("Marginal case buffer allocation exceeded absolute maximum", VOROPP_MEMORY_ERROR);
      }
#if VOROPP_VERBOSE >= 2
      fprintf(stderr, "Marginal cases buffer scaled up to %d\n", i);
#endif
      int *pmarg = new int[current_marginal];
      for (int j = 0; j < n_marg; j++) {
        pmarg[j] = marg[j];
      }
      delete[] marg;
      marg = pmarg;
    }
    marg[n_marg++] = n;
    marg[n_marg++] = ans > tolerance ? 1 : (ans < -tolerance ? -1 : 0);
    return marg[n_marg - 1];
  }

/** For each face of the Voronoi cell, this routine prints the out the normal
   * vector of the face, and scales it to the distance from the cell center to
   * that plane.
   * \param[out] v the vector to store the results in. */
  void voronoicell_base::normals(std::vector<double> &v) {
    int i;
    int j;
    int k;
    v.clear();
    for (i = 1; i < p; i++) {
      for (j = 0; j < nu[i]; j++) {
        k = ed[i][j];
        if (k >= 0) {
          normals_search(v, i, j, k);
        }
      }
    }
    reset_edges();
  }

/** This inline routine is called by normals(). It attempts to construct a
   * single normal vector that is associated with a particular face. It first
   * traces around the face, trying to find two vectors along the face edges
   * whose vector product is above the numerical tolerance. It then constructs
   * the normal vector using this product. If the face is too small, and none of
   * the vector products are large enough, the routine may return (0,0,0) as the
   * normal vector.
   * \param[in] v the vector to store the results in.
   * \param[in] i the initial vertex of the face to test.
   * \param[in] j the index of an edge of the vertex.
   * \param[in] k the neighboring vertex of i, set to ed[i][j]. */
  inline void voronoicell_base::normals_search(std::vector<double> &v, int i, int j, int k) {
    ed[i][j] = -1 - k;
    int l = cycle_up(ed[i][nu[i] + j], k);
    int m;
    double ux;
    double uy;
    double uz;
    double vx;
    double vy;
    double vz;
    double wx;
    double wy;
    double wz;
    double wmag;
    do {
      m = ed[k][l];
      ed[k][l] = -1 - m;
      ux = pts[3 * m] - pts[3 * k];
      uy = pts[3 * m + 1] - pts[3 * k + 1];
      uz = pts[3 * m + 2] - pts[3 * k + 2];

      // Test to see if the length of this edge is above the tolerance
      if (ux * ux + uy * uy + uz * uz > tolerance_sq) {
        while (m != i) {
          l = cycle_up(ed[k][nu[k] + l], m);
          k = m;
          m = ed[k][l];
          ed[k][l] = -1 - m;
          vx = pts[3 * m] - pts[3 * k];
          vy = pts[3 * m + 1] - pts[3 * k + 1];
          vz = pts[3 * m + 2] - pts[3 * k + 2];

          // Construct the vector product of this edge with
          // the previous one
          wx = uz * vy - uy * vz;
          wy = ux * vz - uz * vx;
          wz = uy * vx - ux * vy;
          wmag = wx * wx + wy * wy + wz * wz;

          // Test to see if this vector product of the
          // two edges is above the tolerance
          if (wmag > tolerance_sq) {
            // Construct the normal vector and print it
            wmag = 1 / sqrt(wmag);
            v.push_back(wx * wmag);
            v.push_back(wy * wmag);
            v.push_back(wz * wmag);

            // Mark all of the remaining edges of this
            // face and exit
            while (m != i) {
              l = cycle_up(ed[k][nu[k] + l], m);
              k = m;
              m = ed[k][l];
              ed[k][l] = -1 - m;
            }
            return;
          }
        }
        v.push_back(0);
        v.push_back(0);
        v.push_back(0);
        return;
      }
      l = cycle_up(ed[k][nu[k] + l], m);
      k = m;
    } while (k != i);
    v.push_back(0);
    v.push_back(0);
    v.push_back(0);
  }

/** Returns the number of faces of a computed Voronoi cell.
   * \return The number of faces. */
  int voronoicell_base::number_of_faces() {
    int i;
    int j;
    int k;
    int l;
    int m;
    int s = 0;
    for (i = 1; i < p; i++) {
      for (j = 0; j < nu[i]; j++) {
        k = ed[i][j];
        if (k >= 0) {
          s++;
          ed[i][j] = -1 - k;
          l = cycle_up(ed[i][nu[i] + j], k);
          do {
            m = ed[k][l];
            ed[k][l] = -1 - m;
            l = cycle_up(ed[k][nu[k] + l], m);
            k = m;
          } while (k != i);
        }
      }
    }
    reset_edges();
    return s;
  }

/** Returns a vector of the vertex orders.
   * \param[out] v the vector to store the results in. */
  void voronoicell_base::vertex_orders(std::vector<int> &v) const {
    v.resize(p);
    for (int i = 0; i < p; i++) {
      v[i] = nu[i];
    }
  }

/** Outputs the vertex orders.
   * \param[out] fp the file handle to write to. */
  void voronoicell_base::output_vertex_orders(FILE *fp) const {
    if (p > 0) {
      fprintf(fp, "%d", *nu);
      for (int *nup = nu + 1; nup < nu + p; nup++) {
        fprintf(fp, " %d", *nup);
      }
    }
  }

/** Returns a vector of the vertex vectors using the local coordinate system.
   * \param[out] v the vector to store the results in. */
  void voronoicell_base::vertices(std::vector<double> &v) const {
    v.resize(3 * p);
    double *ptsp = pts;
    for (int i = 0; i < 3 * p; i += 3) {
      v[i] = *(ptsp++) * 0.5;
      v[i + 1] = *(ptsp++) * 0.5;
      v[i + 2] = *(ptsp++) * 0.5;
    }
  }

/** Outputs the vertex vectors using the local coordinate system.
   * \param[out] fp the file handle to write to. */
  void voronoicell_base::output_vertices(FILE *fp) const {
    if (p > 0) {
      fprintf(fp, "(%g,%g,%g)", *pts * 0.5, pts[1] * 0.5, pts[2] * 0.5);
      for (double *ptsp = pts + 3; ptsp < pts + 3 * p; ptsp += 3) {
        fprintf(fp, " (%g,%g,%g)", *ptsp * 0.5, ptsp[1] * 0.5, ptsp[2] * 0.5);
      }
    }
  }

/** Returns a vector of the vertex vectors in the global coordinate system.
   * \param[out] v the vector to store the results in.
   * \param[in] (x,y,z) the position vector of the particle in the global
   *                    coordinate system. */
  void voronoicell_base::vertices(double x, double y, double z, std::vector<double> &v) const {
    v.resize(3 * p);
    double *ptsp = pts;
    for (int i = 0; i < 3 * p; i += 3) {
      v[i] = x + *(ptsp++) * 0.5;
      v[i + 1] = y + *(ptsp++) * 0.5;
      v[i + 2] = z + *(ptsp++) * 0.5;
    }
  }

/** Outputs the vertex vectors using the global coordinate system.
   * \param[out] fp the file handle to write to.
   * \param[in] (x,y,z) the position vector of the particle in the global
   *                    coordinate system. */
  void voronoicell_base::output_vertices(double x, double y, double z, FILE *fp) const {
    if (p > 0) {
      fprintf(fp, "(%g,%g,%g)", x + *pts * 0.5, y + pts[1] * 0.5, z + pts[2] * 0.5);
      for (double *ptsp = pts + 3; ptsp < pts + 3 * p; ptsp += 3) {
        fprintf(fp, " (%g,%g,%g)", x + *ptsp * 0.5, y + ptsp[1] * 0.5, z + ptsp[2] * 0.5);
      }
    }
  }

/** This routine returns the perimeters of each face.
   * \param[out] v the vector to store the results in. */
  void voronoicell_base::face_perimeters(std::vector<double> &v) {
    v.clear();
    int i;
    int j;
    int k;
    int l;
    int m;
    double dx;
    double dy;
    double dz;
    double perim;
    for (i = 1; i < p; i++) {
      for (j = 0; j < nu[i]; j++) {
        k = ed[i][j];
        if (k >= 0) {
          dx = pts[3 * k] - pts[3 * i];
          dy = pts[3 * k + 1] - pts[3 * i + 1];
          dz = pts[3 * k + 2] - pts[3 * i + 2];
          perim = sqrt(dx * dx + dy * dy + dz * dz);
          ed[i][j] = -1 - k;
          l = cycle_up(ed[i][nu[i] + j], k);
          do {
            m = ed[k][l];
            dx = pts[3 * m] - pts[3 * k];
            dy = pts[3 * m + 1] - pts[3 * k + 1];
            dz = pts[3 * m + 2] - pts[3 * k + 2];
            perim += sqrt(dx * dx + dy * dy + dz * dz);
            ed[k][l] = -1 - m;
            l = cycle_up(ed[k][nu[k] + l], m);
            k = m;
          } while (k != i);
          v.push_back(0.5 * perim);
        }
      }
    }
    reset_edges();
  }

/** For each face, this routine outputs a bracketed sequence of numbers
   * containing a list of all the vertices that make up that face.
   * \param[out] v the vector to store the results in. */
  void voronoicell_base::face_vertices(std::vector<int> &v) {
    int i;
    int j;
    int k;
    int l;
    int m;
    int vp(0);
    int vn;
    v.clear();
    for (i = 1; i < p; i++) {
      for (j = 0; j < nu[i]; j++) {
        k = ed[i][j];
        if (k >= 0) {
          v.push_back(0);
          v.push_back(i);
          ed[i][j] = -1 - k;
          l = cycle_up(ed[i][nu[i] + j], k);
          do {
            v.push_back(k);
            m = ed[k][l];
            ed[k][l] = -1 - m;
            l = cycle_up(ed[k][nu[k] + l], m);
            k = m;
          } while (k != i);
          vn = v.size();
          v[vp] = vn - vp - 1;
          vp = vn;
        }
      }
    }
    reset_edges();
  }

/** Outputs a list of the number of edges in each face.
   * \param[out] v the vector to store the results in. */
  void voronoicell_base::face_orders(std::vector<int> &v) {
    int i;
    int j;
    int k;
    int l;
    int m;
    int q;
    v.clear();
    for (i = 1; i < p; i++) {
      for (j = 0; j < nu[i]; j++) {
        k = ed[i][j];
        if (k >= 0) {
          q = 1;
          ed[i][j] = -1 - k;
          l = cycle_up(ed[i][nu[i] + j], k);
          do {
            q++;
            m = ed[k][l];
            ed[k][l] = -1 - m;
            l = cycle_up(ed[k][nu[k] + l], m);
            k = m;
          } while (k != i);
          v.push_back(q);
          ;
        }
      }
    }
    reset_edges();
  }

/** Computes the number of edges that each face has and outputs a frequency
   * table of the results.
   * \param[out] v the vector to store the results in. */
  void voronoicell_base::face_freq_table(std::vector<int> &v) {
    int i;
    int j;
    int k;
    int l;
    int m;
    int q;
    v.clear();
    for (i = 1; i < p; i++) {
      for (j = 0; j < nu[i]; j++) {
        k = ed[i][j];
        if (k >= 0) {
          q = 1;
          ed[i][j] = -1 - k;
          l = cycle_up(ed[i][nu[i] + j], k);
          do {
            q++;
            m = ed[k][l];
            ed[k][l] = -1 - m;
            l = cycle_up(ed[k][nu[k] + l], m);
            k = m;
          } while (k != i);
          if ((unsigned int) q >= v.size()) {
            v.resize(q + 1, 0);
          }
          v[q]++;
        }
      }
    }
    reset_edges();
  }

/** This routine tests to see whether the cell intersects a plane by starting
   * from the guess point up. If up intersects, then it immediately returns true.
   * Otherwise, it calls the plane_intersects_track() routine.
   * \param[in] (x,y,z) the normal vector to the plane.
   * \param[in] rsq the distance along this vector of the plane.
   * \return False if the plane does not intersect the plane, true if it does. */
  bool voronoicell_base::plane_intersects(double x, double y, double z, double rsq) {
    const double g = x * pts[3 * up] + y * pts[3 * up + 1] + z * pts[3 * up + 2];
    if (g < rsq) {
      return plane_intersects_track(x, y, z, rsq, g);
    }
    return true;
  }

/** This routine tests to see if a cell intersects a plane. It first tests a
   * random sample of approximately sqrt(p)/4 points. If any of those are
   * intersect, then it immediately returns true. Otherwise, it takes the closest
   * point and passes that to plane_intersect_track() routine.
   * \param[in] (x,y,z) the normal vector to the plane.
   * \param[in] rsq the distance along this vector of the plane.
   * \return False if the plane does not intersect the plane, true if it does. */
  bool voronoicell_base::plane_intersects_guess(double x, double y, double z, double rsq) {
    up = 0;
    double g = x * pts[3 * up] + y * pts[3 * up + 1] + z * pts[3 * up + 2];
    if (g < rsq) {
      int ca = 1;
      const int cc = p >> 3;
      int mp = 1;
      double m;
      while (ca < cc) {
        m = x * pts[3 * mp] + y * pts[3 * mp + 1] + z * pts[3 * mp + 2];
        if (m > g) {
          if (m > rsq) {
            return true;
          }
          g = m;
          up = mp;
        }
        ca += mp++;
      }
      return plane_intersects_track(x, y, z, rsq, g);
    }
    return true;
  }

/* This routine tests to see if a cell intersects a plane, by tracing over the cell from
   * vertex to vertex, starting at up. It is meant to be called either by plane_intersects()
   * or plane_intersects_track(), when those routines cannot immediately resolve the case.
   * \param[in] (x,y,z) the normal vector to the plane.
   * \param[in] rsq the distance along this vector of the plane.
   * \param[in] g the distance of up from the plane.
   * \return False if the plane does not intersect the plane, true if it does. */
  inline bool voronoicell_base::plane_intersects_track(double x, double y, double z, double rsq, double g) {
    int count = 0;
    int ls;
    int us;
    int tp;
    double t;

    // The test point is outside of the cutting space
    for (us = 0; us < nu[up]; us++) {
      tp = ed[up][us];
      t = x * pts[3 * tp] + y * pts[3 * tp + 1] + z * pts[3 * tp + 2];
      if (t > g) {
        ls = ed[up][nu[up] + us];
        up = tp;
        while (t < rsq) {
          if (++count >= p) {
#if VOROPP_VERBOSE >= 1
            fputs("Bailed out of convex calculation", stderr);
#endif
            for (tp = 0; tp < p; tp++) {
              if (x * pts[3 * tp] + y * pts[3 * tp + 1] + z * pts[3 * tp + 2] > rsq) {
                return true;
              }
            }
            return false;
          }

          // Test all the neighbors of the current point
          // and find the one which is closest to the
          // plane
          for (us = 0; us < ls; us++) {
            tp = ed[up][us];
            g = x * pts[3 * tp] + y * pts[3 * tp + 1] + z * pts[3 * tp + 2];
            if (g > t) {
              break;
            }
          }
          if (us == ls) {
            us++;
            while (us < nu[up]) {
              tp = ed[up][us];
              g = x * pts[3 * tp] + y * pts[3 * tp + 1] + z * pts[3 * tp + 2];
              if (g > t) {
                break;
              }
              us++;
            }
            if (us == nu[up]) {
              return false;
            }
          }
          ls = ed[up][nu[up] + us];
          up = tp;
          t = g;
        }
        return true;
      }
    }
    return false;
  }

/** Counts the number of edges of the Voronoi cell.
   * \return the number of edges. */
  int voronoicell_base::number_of_edges() const {
    int edges = 0;
    int *nup = nu;
    while (nup < nu + p) {
      edges += *(nup++);
    }
    return edges >> 1;
  }

/** Outputs a custom string of information about the Voronoi cell. The string
   * of information follows a similar style as the C printf command, and detailed
   * information about its format is available at
   * http://math.lbl.gov/voro++/doc/custom.html.
   * \param[in] format the custom string to print.
   * \param[in] i the ID of the particle associated with this Voronoi cell.
   * \param[in] (x,y,z) the position of the particle associated with this Voronoi
   *                    cell.
   * \param[in] r a radius associated with the particle.
   * \param[in] fp the file handle to write to. */
  void voronoicell_base::output_custom(const char *format, int i, double x, double y, double z, double r, FILE *fp) {
    char *fmp = (const_cast<char *>(format));
    std::vector<int> vi;
    std::vector<double> vd;
    while (*fmp != 0) {
      if (*fmp == '%') {
        fmp++;
        switch (*fmp) {
          // Particle-related output
          case 'i': fprintf(fp, "%d", i); break;
          case 'x': fprintf(fp, "%g", x); break;
          case 'y': fprintf(fp, "%g", y); break;
          case 'z': fprintf(fp, "%g", z); break;
          case 'q': fprintf(fp, "%g %g %g", x, y, z); break;
          case 'r':
            fprintf(fp, "%g", r);
            break;

            // Vertex-related output
          case 'w': fprintf(fp, "%d", p); break;
          case 'p': output_vertices(fp); break;
          case 'P': output_vertices(x, y, z, fp); break;
          case 'o': output_vertex_orders(fp); break;
          case 'm':
            fprintf(fp, "%g", 0.25 * max_radius_squared());
            break;

            // Edge-related output
          case 'g': fprintf(fp, "%d", number_of_edges()); break;
          case 'E': fprintf(fp, "%g", total_edge_distance()); break;
          case 'e':
            face_perimeters(vd);
            voro_print_vector(vd, fp);
            break;

            // Face-related output
          case 's': fprintf(fp, "%d", number_of_faces()); break;
          case 'F': fprintf(fp, "%g", surface_area()); break;
          case 'A': {
            face_freq_table(vi);
            voro_print_vector(vi, fp);
          } break;
          case 'a':
            face_orders(vi);
            voro_print_vector(vi, fp);
            break;
          case 'f':
            face_areas(vd);
            voro_print_vector(vd, fp);
            break;
          case 't': {
            face_vertices(vi);
            voro_print_face_vertices(vi, fp);
          } break;
          case 'l':
            normals(vd);
            voro_print_positions(vd, fp);
            break;
          case 'n':
            neighbors(vi);
            voro_print_vector(vi, fp);
            break;

            // Volume-related output
          case 'v': fprintf(fp, "%g", volume()); break;
          case 'c': {
            double cx;
            double cy;
            double cz;
            centroid(cx, cy, cz);
            fprintf(fp, "%g %g %g", cx, cy, cz);
          } break;
          case 'C': {
            double cx;
            double cy;
            double cz;
            centroid(cx, cy, cz);
            fprintf(fp, "%g %g %g", x + cx, y + cy, z + cz);
          } break;

            // End-of-string reached
          case 0:
            fmp--;
            break;

            // The percent sign is not part of a
            // control sequence
          default: putc('%', fp); putc(*fmp, fp);
        }
      } else {
        putc(*fmp, fp);
      }
      fmp++;
    }
    fputs("\n", fp);
  }

/** This initializes the class to be a rectangular box. It calls the base class
   * initialization routine to set up the edge and vertex information, and then
   * sets up the neighbor information, with initial faces being assigned ID
   * numbers from -1 to -6.
   * \param[in] (xmin,xmax) the minimum and maximum x coordinates.
   * \param[in] (ymin,ymax) the minimum and maximum y coordinates.
   * \param[in] (zmin,zmax) the minimum and maximum z coordinates. */
  void voronoicell_neighbor::init(double xmin, double xmax, double ymin, double ymax, double zmin, double zmax) {
    init_base(xmin, xmax, ymin, ymax, zmin, zmax);
    int *q = mne[3];
    *q = -5;
    q[1] = -3;
    q[2] = -1;
    q[3] = -5;
    q[4] = -2;
    q[5] = -3;
    q[6] = -5;
    q[7] = -1;
    q[8] = -4;
    q[9] = -5;
    q[10] = -4;
    q[11] = -2;
    q[12] = -6;
    q[13] = -1;
    q[14] = -3;
    q[15] = -6;
    q[16] = -3;
    q[17] = -2;
    q[18] = -6;
    q[19] = -4;
    q[20] = -1;
    q[21] = -6;
    q[22] = -2;
    q[23] = -4;
    *ne = q;
    ne[1] = q + 3;
    ne[2] = q + 6;
    ne[3] = q + 9;
    ne[4] = q + 12;
    ne[5] = q + 15;
    ne[6] = q + 18;
    ne[7] = q + 21;
  }

/** This initializes the class to be an octahedron. It calls the base class
   * initialization routine to set up the edge and vertex information, and then
   * sets up the neighbor information, with the initial faces being assigned ID
   * numbers from -1 to -8.
   * \param[in] l The distance from the octahedron center to a vertex. Six
   *              vertices are initialized at (-l,0,0), (l,0,0), (0,-l,0),
   *              (0,l,0), (0,0,-l), and (0,0,l). */
  void voronoicell_neighbor::init_octahedron(double l) {
    init_octahedron_base(l);
    int *q = mne[4];
    *q = -5;
    q[1] = -6;
    q[2] = -7;
    q[3] = -8;
    q[4] = -1;
    q[5] = -2;
    q[6] = -3;
    q[7] = -4;
    q[8] = -6;
    q[9] = -5;
    q[10] = -2;
    q[11] = -1;
    q[12] = -8;
    q[13] = -7;
    q[14] = -4;
    q[15] = -3;
    q[16] = -5;
    q[17] = -8;
    q[18] = -3;
    q[19] = -2;
    q[20] = -7;
    q[21] = -6;
    q[22] = -1;
    q[23] = -4;
    *ne = q;
    ne[1] = q + 4;
    ne[2] = q + 8;
    ne[3] = q + 12;
    ne[4] = q + 16;
    ne[5] = q + 20;
  }

/** This initializes the class to be a tetrahedron. It calls the base class
   * initialization routine to set up the edge and vertex information, and then
   * sets up the neighbor information, with the initial faces being assigned ID
   * numbers from -1 to -4.
   * \param (x0,y0,z0) a position vector for the first vertex.
   * \param (x1,y1,z1) a position vector for the second vertex.
   * \param (x2,y2,z2) a position vector for the third vertex.
   * \param (x3,y3,z3) a position vector for the fourth vertex. */
  void voronoicell_neighbor::init_tetrahedron(double x0, double y0, double z0, double x1, double y1, double z1, double x2, double y2, double z2, double x3, double y3, double z3) {
    init_tetrahedron_base(x0, y0, z0, x1, y1, z1, x2, y2, z2, x3, y3, z3);
    int *q = mne[3];
    *q = -4;
    q[1] = -3;
    q[2] = -2;
    q[3] = -3;
    q[4] = -4;
    q[5] = -1;
    q[6] = -4;
    q[7] = -2;
    q[8] = -1;
    q[9] = -2;
    q[10] = -3;
    q[11] = -1;
    *ne = q;
    ne[1] = q + 3;
    ne[2] = q + 6;
    ne[3] = q + 9;
  }

/** This routine checks to make sure the neighbor information of each face is
   * consistent. */
  void voronoicell_neighbor::check_facets() {
    int i;
    int j;
    int k;
    int l;
    int m;
    int q;
    for (i = 1; i < p; i++) {
      for (j = 0; j < nu[i]; j++) {
        k = ed[i][j];
        if (k >= 0) {
          ed[i][j] = -1 - k;
          q = ne[i][j];
          l = cycle_up(ed[i][nu[i] + j], k);
          do {
            m = ed[k][l];
            ed[k][l] = -1 - m;
            if (ne[k][l] != q) {
              fprintf(stderr, "Facet error at (%d,%d)=%d, started from (%d,%d)=%d\n", k, l, ne[k][l], i, j, q);
            }
            l = cycle_up(ed[k][nu[k] + l], m);
            k = m;
          } while (k != i);
        }
      }
    }
    reset_edges();
  }

/** The class constructor allocates memory for storing neighbor information. */
  voronoicell_neighbor::voronoicell_neighbor() {
    int i;
    mne = new int *[current_vertex_order];
    ne = new int *[current_vertices];
    for (i = 0; i < 3; i++) {
      mne[i] = new int[init_n_vertices * i];
    }
    mne[3] = new int[init_3_vertices * 3];
    for (i = 4; i < current_vertex_order; i++) {
      mne[i] = new int[init_n_vertices * i];
    }
  }

/** The class destructor frees the dynamically allocated memory for storing
   * neighbor information. */
  voronoicell_neighbor::~voronoicell_neighbor() {
    for (int i = current_vertex_order - 1; i >= 0; i--) {
      if (mem[i] > 0) {
        delete[] mne[i];
      }
    }
    delete[] mne;
    delete[] ne;
  }

/** Computes a vector list of neighbors. */
  void voronoicell_neighbor::neighbors(std::vector<int> &v) {
    v.clear();
    int i;
    int j;
    int k;
    int l;
    int m;
    for (i = 1; i < p; i++) {
      for (j = 0; j < nu[i]; j++) {
        k = ed[i][j];
        if (k >= 0) {
          v.push_back(ne[i][j]);
          ed[i][j] = -1 - k;
          l = cycle_up(ed[i][nu[i] + j], k);
          do {
            m = ed[k][l];
            ed[k][l] = -1 - m;
            l = cycle_up(ed[k][nu[k] + l], m);
            k = m;
          } while (k != i);
        }
      }
    }
    reset_edges();
  }

/** Prints the vertices, their edges, the relation table, and also notifies if
   * any memory errors are visible. */
  void voronoicell_base::print_edges() {
    int j;
    double *ptsp = pts;
    for (int i = 0; i < p; i++, ptsp += 3) {
      printf("%d %d  ", i, nu[i]);
      for (j = 0; j < nu[i]; j++) {
        printf(" %d", ed[i][j]);
      }
      printf("  ");
      while (j < (nu[i] << 1)) {
        printf(" %d", ed[i][j]);
      }
      printf("   %d", ed[i][j]);
      print_edges_neighbors(i);
      printf("  %g %g %g %p", *ptsp, ptsp[1], ptsp[2], (void *) ed[i]);
      if (ed[i] >= mep[nu[i]] + mec[nu[i]] * ((nu[i] << 1) + 1)) {
        puts(" Memory error");
      } else {
        puts("");
      }
    }
  }

/** This prints out the neighbor information for vertex i. */
  void voronoicell_neighbor::print_edges_neighbors(int i) {
    if (nu[i] > 0) {
      int j = 0;
      printf("     (");
      while (j < nu[i] - 1) {
        printf("%d,", ne[i][j++]);
      }
      printf("%d)", ne[i][j]);
    } else {
      printf("     ()");
    }
  }

// Explicit instantiation
  template bool voronoicell_base::nplane(voronoicell &, double, double, double, double, int);
  template bool voronoicell_base::nplane(voronoicell_neighbor &, double, double, double, double, int);
  template void voronoicell_base::check_memory_for_copy(voronoicell &, voronoicell_base *);
  template void voronoicell_base::check_memory_for_copy(voronoicell_neighbor &, voronoicell_base *);

} // namespace voro
// Voro++, a 3D cell-based Voronoi library
//
// Author   : Chris H. Rycroft (LBL / UC Berkeley)
// Email    : chr@alum.mit.edu
// Date     : August 30th 2011

/** \file v_base.cc
 * \brief Function implementations for the base Voronoi container class. */

namespace voro {

/** This function is called during container construction. The routine scans
   * all of the worklists in the wl[] array. For a given worklist of blocks
   * labeled \f$w_1\f$ to \f$w_n\f$, it computes a sequence \f$r_0\f$ to
   * \f$r_n\f$ so that $r_i$ is the minimum distance to all the blocks
   * \f$w_{j}\f$ where \f$j>i\f$ and all blocks outside the worklist. The values
   * of \f$r_n\f$ is calculated first, as the minimum distance to any block in
   * the shell surrounding the worklist. The \f$r_i\f$ are then computed in
   * reverse order by considering the distance to \f$w_{i+1}\f$. */
  voro_base::voro_base(int nx_, int ny_, int nz_, double boxx_, double boxy_, double boxz_) :
      nx(nx_), ny(ny_), nz(nz_), nxy(nx_ * ny_), nxyz(nxy * nz_), boxx(boxx_), boxy(boxy_), boxz(boxz_), xsp(1 / boxx_), ysp(1 / boxy_), zsp(1 / boxz_), mrad(new double[wl_hgridcu * wl_seq_length]) {
    const unsigned int b1 = 1 << 21;
    const unsigned int b2 = 1 << 22;
    const unsigned int b3 = 1 << 24;
    const unsigned int b4 = 1 << 25;
    const unsigned int b5 = 1 << 27;
    const unsigned int b6 = 1 << 28;
    const double xstep = boxx / wl_fgrid;
    const double ystep = boxy / wl_fgrid;
    const double zstep = boxz / wl_fgrid;
    int i;
    int j;
    int k;
    int lx;
    int ly;
    int lz;
    int q;
    unsigned int f;
    unsigned int *e = const_cast<unsigned int *>(wl);
    double xlo;
    double ylo;
    double zlo;
    double xhi;
    double yhi;
    double zhi;
    double minr;
    double *radp = mrad;
    for (zlo = 0, zhi = zstep, lz = 0; lz < wl_hgrid; zlo = zhi, zhi += zstep, lz++) {
      for (ylo = 0, yhi = ystep, ly = 0; ly < wl_hgrid; ylo = yhi, yhi += ystep, ly++) {
        for (xlo = 0, xhi = xstep, lx = 0; lx < wl_hgrid; xlo = xhi, xhi += xstep, lx++) {
          minr = large_number;
          for (q = e[0] + 1; q < wl_seq_length; q++) {
            f = e[q];
            i = (f & 127) - 64;
            j = (f >> 7 & 127) - 64;
            k = (f >> 14 & 127) - 64;
            if ((f & b2) == b2) {
              compute_minimum(minr, xlo, xhi, ylo, yhi, zlo, zhi, i - 1, j, k);
              if ((f & b1) == 0) {
                compute_minimum(minr, xlo, xhi, ylo, yhi, zlo, zhi, i + 1, j, k);
              }
            } else if ((f & b1) == b1) {
              compute_minimum(minr, xlo, xhi, ylo, yhi, zlo, zhi, i + 1, j, k);
            }
            if ((f & b4) == b4) {
              compute_minimum(minr, xlo, xhi, ylo, yhi, zlo, zhi, i, j - 1, k);
              if ((f & b3) == 0) {
                compute_minimum(minr, xlo, xhi, ylo, yhi, zlo, zhi, i, j + 1, k);
              }
            } else if ((f & b3) == b3) {
              compute_minimum(minr, xlo, xhi, ylo, yhi, zlo, zhi, i, j + 1, k);
            }
            if ((f & b6) == b6) {
              compute_minimum(minr, xlo, xhi, ylo, yhi, zlo, zhi, i, j, k - 1);
              if ((f & b5) == 0) {
                compute_minimum(minr, xlo, xhi, ylo, yhi, zlo, zhi, i, j, k + 1);
              }
            } else if ((f & b5) == b5) {
              compute_minimum(minr, xlo, xhi, ylo, yhi, zlo, zhi, i, j, k + 1);
            }
          }
          q--;
          while (q > 0) {
            radp[q] = minr;
            f = e[q];
            i = (f & 127) - 64;
            j = (f >> 7 & 127) - 64;
            k = (f >> 14 & 127) - 64;
            compute_minimum(minr, xlo, xhi, ylo, yhi, zlo, zhi, i, j, k);
            q--;
          }
          *radp = minr;
          e += wl_seq_length;
          radp += wl_seq_length;
        }
      }
    }
  }

/** Computes the minimum distance from a subregion to a given block. If this distance
   * is smaller than the value of minr, then it passes
   * \param[in,out] minr a pointer to the current minimum distance. If the distance
   *                     computed in this function is smaller, then this distance is
   *                     set to the new one.
   * \param[out] (xlo,ylo,zlo) the lower coordinates of the subregion being
   *                           considered.
   * \param[out] (xhi,yhi,zhi) the upper coordinates of the subregion being
   *                           considered.
   * \param[in] (ti,tj,tk) the coordinates of the block. */
  void voro_base::compute_minimum(double &minr, double &xlo, double &xhi, double &ylo, double &yhi, double &zlo, double &zhi, int ti, int tj, int tk) const {
    double radsq;
    double temp;
    if (ti > 0) {
      temp = boxx * ti - xhi;
      radsq = temp * temp;
    } else if (ti < 0) {
      temp = xlo - boxx * (1 + ti);
      radsq = temp * temp;
    } else {
      radsq = 0;
    }

    if (tj > 0) {
      temp = boxy * tj - yhi;
      radsq += temp * temp;
    } else if (tj < 0) {
      temp = ylo - boxy * (1 + tj);
      radsq += temp * temp;
    }

    if (tk > 0) {
      temp = boxz * tk - zhi;
      radsq += temp * temp;
    } else if (tk < 0) {
      temp = zlo - boxz * (1 + tk);
      radsq += temp * temp;
    }

    if (radsq < minr) {
      minr = radsq;
    }
  }

/** Checks to see whether "%n" appears in a format sequence to determine
   * whether neighbor information is required or not.
   * \param[in] format the format string to check.
   * \return True if a "%n" is found, false otherwise. */
  bool voro_base::contains_neighbor(const char *format) {
    char *fmp = (const_cast<char *>(format));

    // Check to see if "%n" appears in the format sequence
    while (*fmp != 0) {
      if (*fmp == '%') {
        fmp++;
        if (*fmp == 'n') {
          return true;
        } else if (*fmp == 0) {
          return false;
        }
      }
      fmp++;
    }

    return false;
  }

// Voro++, a 3D cell-based Voronoi library
//
// Author   : Chris H. Rycroft (LBL / UC Berkeley)
// Email    : chr@alum.mit.edu
// Date     : August 30th 2011

/** \file v_base_wl.cc
   * \brief The table of block worklists that are used during the cell
   * computation, which is part of the voro_base class.
   *
   * This file is automatically generated by worklist_gen.pl and it is not
   * intended to be edited by hand. */

  const unsigned int voro_base::wl[wl_seq_length * wl_hgridcu] = {
      7,          0x10203f,   0x101fc0,   0xfe040,    0xfe03f,    0x101fbf,   0xfdfc0,    0xfdfbf,    0x10fe0bf,  0x11020bf,  0x11020c0,  0x10fe0c0,  0x2fe041,   0x302041,   0x301fc1,   0x2fdfc1,   0x8105fc0,
      0x8106040,  0x810603f,  0x8105fbf,  0x701fbe,   0x70203e,   0x6fe03e,   0x6fdfbe,   0x30fdf3f,  0x3101f3f,  0x3101f40,  0x30fdf40,  0x180f9fc0, 0x180fa040, 0x180fa03f, 0x180f9fbf, 0x12fe0c1,  0x13020c1,
      0x91060c0,  0x91060bf,  0x8306041,  0x8305fc1,  0x3301f41,  0x32fdf41,  0x182f9fc1, 0x182fa041, 0x190fa0c0, 0x190fa0bf, 0x16fe0be,  0x17020be,  0x870603e,  0x8705fbe,  0xb105f3f,  0xb105f40,  0x3701f3e,
      0x36fdf3e,  0x186f9fbe, 0x186fa03e, 0x1b0f9f3f, 0x1b0f9f40, 0x93060c1,  0x192fa0c1, 0x97060be,  0xb305f41,  0x1b2f9f41, 0x196fa0be, 0xb705f3e,  0x1b6f9f3e, 11,         0x101fc0,   0xfe040,    0xfdfc0,
      0x10203f,   0x101fbf,   0xfe03f,    0xfdfbf,    0xfdfc1,    0x101fc1,   0x102041,   0xfe041,    0x10fe0c0,  0x11020c0,  0x8106040,  0x8105fc0,  0x8105fbf,  0x810603f,  0x11020bf,  0x10fe0bf,  0x180fa040,
      0x180f9fc0, 0x30fdf40,  0x3101f40,  0x3101f3f,  0x30fdf3f,  0x180f9fbf, 0x180fa03f, 0x6fe03e,   0x70203e,   0x701fbe,   0x6fdfbe,   0x8105fc1,  0x8106041,  0x11020c1,  0x10fe0c1,  0x180fa041, 0x180f9fc1,
      0x30fdf41,  0x3101f41,  0x91060c0,  0x91060bf,  0x190fa0c0, 0x190fa0bf, 0xb105f40,  0xb105f3f,  0x8705fbe,  0x870603e,  0x97020be,  0x16fe0be,  0x1b0f9f40, 0x1b0f9f3f, 0x36fdf3e,  0xb701f3e,  0x1b6f9fbe,
      0x196fa03e, 0x93060c1,  0xb305f41,  0x192fa0c1, 0x1b2f9f41, 0x1b2fdfc2, 0xb301fc2,  0x9302042,  0x192fe042, 11,         0x101fc0,   0xfe040,    0xfdfc0,    0xfdfbf,    0x101fbf,   0x10203f,   0xfe03f,
      0xfe041,    0x102041,   0x101fc1,   0xfdfc1,    0x8105fc0,  0x8106040,  0x11020c0,  0x10fe0c0,  0x10fe0bf,  0x11020bf,  0x810603f,  0x8105fbf,  0x3101f40,  0x30fdf40,  0x180f9fc0, 0x180fa040, 0x180fa03f,
      0x180f9fbf, 0x30fdf3f,  0x3101f3f,  0x8105fc1,  0x8106041,  0x11020c1,  0x10fe0c1,  0x180fa041, 0x180f9fc1, 0x30fdf41,  0x3101f41,  0x701fbe,   0x70203e,   0x6fe03e,   0x6fdfbe,   0x91060c0,  0x91060bf,
      0xb105f40,  0xb105f3f,  0x190fa0c0, 0x190fa0bf, 0x93060c1,  0x1b0f9f40, 0x1b0f9f3f, 0xb305f41,  0x192fa0c1, 0x16fe0be,  0x17020be,  0x970603e,  0x8705fbe,  0xb701f3e,  0x36fdf3e,  0x1b2f9f41, 0x1b2fdfc2,
      0x192fe042, 0x9302042,  0xb301fc2,  0x1b6f9fbe, 0x196fa03e, 11,         0x101fc0,   0xfe040,    0xfdfc0,    0xfdfbf,    0x101fbf,   0x10203f,   0xfe03f,    0xfe041,    0x102041,   0x101fc1,   0xfdfc1,
      0x8105fc0,  0x8106040,  0x11020c0,  0x10fe0c0,  0x10fe0bf,  0x11020bf,  0x810603f,  0x8105fbf,  0x3101f40,  0x30fdf40,  0x180f9fc0, 0x180fa040, 0x10fe0c1,  0x11020c1,  0x8106041,  0x8105fc1,  0x3101f3f,
      0x30fdf3f,  0x180f9fbf, 0x180fa03f, 0x180fa041, 0x180f9fc1, 0x30fdf41,  0x3101f41,  0x91060c0,  0x91060bf,  0x70203e,   0x701fbe,   0x6fdfbe,   0x6fe03e,   0x190fa0c0, 0xb105f40,  0xb105f3f,  0x93060c1,
      0x190fa0bf, 0x192fa0c1, 0x1b0f9f40, 0x1b0f9f3f, 0xb305f41,  0xb301fc2,  0x9302042,  0x192fe042, 0x1b2fdfc2, 0x1b2f9f41, 0x16fe0be,  0x17020be,  0x970603e,  0x8705fbe,  0xb701f3e,  0x36fdf3e,  0x1b6f9fbe,
      0x196fa03e, 11,         0x10203f,   0xfe040,    0xfe03f,    0x101fc0,   0x101fbf,   0xfdfc0,    0xfdfbf,    0xfe0bf,    0x1020bf,   0x1020c0,   0xfe0c0,    0x2fe041,   0x302041,   0x8106040,  0x810603f,
      0x8105fbf,  0x8105fc0,  0x301fc1,   0x2fdfc1,   0x180fa040, 0x180fa03f, 0x6fe03e,   0x70203e,   0x701fbe,   0x6fdfbe,   0x180f9fbf, 0x180f9fc0, 0x30fdf40,  0x3101f40,  0x3101f3f,  0x30fdf3f,  0x81060bf,
      0x81060c0,  0x3020c1,   0x2fe0c1,   0x180fa0c0, 0x180fa0bf, 0x6fe0be,   0x7020be,   0x8306041,  0x8305fc1,  0x182fa041, 0x182f9fc1, 0x870603e,  0x8705fbe,  0xb105f3f,  0xb105f40,  0xb301f41,  0x32fdf41,
      0x186fa03e, 0x186f9fbe, 0x36fdf3e,  0xb701f3e,  0x1b6f9f3f, 0x1b2f9f40, 0x93060c1,  0x97060be,  0x192fa0c1, 0x196fa0be, 0x196fe13f, 0x970213f,  0x9302140,  0x192fe140, 9,          0xfe040,    0x10203f,
      0x101fc0,   0xfdfc0,    0xfe03f,    0xfdfbf,    0x101fbf,   0x102041,   0xfe041,    0x10fe0c0,  0x11020c0,  0x11020bf,  0x10fe0bf,  0xfdfc1,    0x101fc1,   0x8106040,  0x810603f,  0x8105fc0,  0x8105fbf,
      0x180fa040, 0x180f9fc0, 0x180fa03f, 0x180f9fbf, 0x10fe0c1,  0x11020c1,  0x70203e,   0x6fe03e,   0x6fdfbe,   0x701fbe,   0x3101f3f,  0x3101f40,  0x30fdf40,  0x30fdf3f,  0x8106041,  0x91060c0,  0x91060bf,
      0x8305fc1,  0x180fa041, 0x190fa0c0, 0x190fa0bf, 0x180f9fc1, 0x30fdf41,  0x3301f41,  0x17020be,  0x16fe0be,  0x93060c1,  0x870603e,  0x8705fbe,  0xb105f3f,  0xb105f40,  0x192fa0c1, 0x1b0f9f40, 0x1b0f9f3f,
      0x186f9fbe, 0x196fa03e, 0x1b6fdf3e, 0xb701f3e,  0x97060be,  0xb305f41,  0x1b2f9f41, 0x1b2fdfc2, 0x192fe042, 0xa302042,  11,         0xfe040,    0x101fc0,   0xfdfc0,    0xfe03f,    0x10203f,   0x101fbf,
      0xfdfbf,    0xfe041,    0x102041,   0x101fc1,   0xfdfc1,    0x10fe0c0,  0x11020c0,  0x11020bf,  0x10fe0bf,  0x8106040,  0x8105fc0,  0x810603f,  0x8105fbf,  0x11020c1,  0x10fe0c1,  0x180fa040, 0x180f9fc0,
      0x180fa03f, 0x180f9fbf, 0x30fdf40,  0x3101f40,  0x8106041,  0x8105fc1,  0x3101f3f,  0x30fdf3f,  0x91060c0,  0x91060bf,  0x180fa041, 0x180f9fc1, 0x6fe03e,   0x70203e,   0x701fbe,   0x6fdfbe,   0x190fa0c0,
      0x190fa0bf, 0x30fdf41,  0x3101f41,  0x93060c1,  0x192fa0c1, 0xb105f40,  0xb105f3f,  0x17020be,  0x16fe0be,  0x970603e,  0x8705fbe,  0x1b0f9f40, 0x1b0f9f3f, 0x186f9fbe, 0x196fa03e, 0xb305f41,  0xb301fc2,
      0x9302042,  0x192fe042, 0x1b2fdfc2, 0x1b2f9f41, 0x1b6fdf3e, 0xb701f3e,  11,         0xfe040,    0x101fc0,   0xfdfc0,    0xfe03f,    0x10203f,   0x101fbf,   0xfdfbf,    0xfe041,    0x102041,   0x101fc1,
      0xfdfc1,    0x10fe0c0,  0x11020c0,  0x11020bf,  0x10fe0bf,  0x8106040,  0x8105fc0,  0x11020c1,  0x10fe0c1,  0x810603f,  0x8105fbf,  0x180fa040, 0x180f9fc0, 0x8106041,  0x8105fc1,  0x3101f40,  0x180fa03f,
      0x180f9fbf, 0x30fdf40,  0x180fa041, 0x180f9fc1, 0x91060c0,  0x3101f3f,  0x30fdf3f,  0x30fdf41,  0x3101f41,  0x91060bf,  0x190fa0c0, 0x91060c1,  0x190fa0bf, 0x186fe03e, 0x70203e,   0x3701fbe,  0x1b6fdfbe,
      0x190fa0c1, 0x182fe042, 0xb105f40,  0xb105f3f,  0x302042,   0x3301fc2,  0x1b2fdfc2, 0x1b0f9f40, 0x1b6f9f3f, 0xb105f41,  0x17020be,  0x196fe0be, 0x970603e,  0xb705fbe,  0x1b2f9f41, 0x192fe0c2, 0x13020c2,
      0x9306042,  0xb305fc2,  11,         0x10203f,   0xfe040,    0xfe03f,    0xfdfbf,    0x101fbf,   0x101fc0,   0xfdfc0,    0xfe0c0,    0x1020c0,   0x1020bf,   0xfe0bf,    0x810603f,  0x8106040,  0x302041,
      0x2fe041,   0x2fdfc1,   0x301fc1,   0x8105fc0,  0x8105fbf,  0x70203e,   0x6fe03e,   0x180fa03f, 0x180fa040, 0x180f9fc0, 0x180f9fbf, 0x6fdfbe,   0x701fbe,   0x81060bf,  0x81060c0,  0x3020c1,   0x2fe0c1,
      0x180fa0c0, 0x180fa0bf, 0x6fe0be,   0x7020be,   0x3101f3f,  0x3101f40,  0x30fdf40,  0x30fdf3f,  0x8306041,  0x8305fc1,  0x870603e,  0x8705fbe,  0x182fa041, 0x182f9fc1, 0x93060c1,  0x186fa03e, 0x186f9fbe,
      0x97060be,  0x192fa0c1, 0x32fdf41,  0x3301f41,  0xb305f40,  0xb105f3f,  0xb701f3e,  0x36fdf3e,  0x196fa0be, 0x196fe13f, 0x192fe140, 0x9302140,  0x970213f,  0x1b6f9f3f, 0x1b2f9f40, 11,         0xfe040,
      0x10203f,   0xfe03f,    0xfdfc0,    0x101fc0,   0x101fbf,   0xfdfbf,    0xfe0c0,    0x1020c0,   0x1020bf,   0xfe0bf,    0x2fe041,   0x302041,   0x301fc1,   0x2fdfc1,   0x8106040,  0x810603f,  0x8105fc0,
      0x8105fbf,  0x3020c1,   0x2fe0c1,   0x180fa040, 0x180fa03f, 0x180f9fc0, 0x180f9fbf, 0x6fe03e,   0x70203e,   0x81060c0,  0x81060bf,  0x701fbe,   0x6fdfbe,   0x8306041,  0x8305fc1,  0x180fa0c0, 0x180fa0bf,
      0x30fdf40,  0x3101f40,  0x3101f3f,  0x30fdf3f,  0x182fa041, 0x182f9fc1, 0x6fe0be,   0x7020be,   0x93060c1,  0x192fa0c1, 0x870603e,  0x8705fbe,  0x3301f41,  0x32fdf41,  0xb305f40,  0xb105f3f,  0x186fa03e,
      0x186f9fbe, 0x1b0f9f3f, 0x1b2f9f40, 0x97060be,  0x970213f,  0x9302140,  0x192fe140, 0x196fe13f, 0x196fa0be, 0x1b6fdf3e, 0xb701f3e,  15,         0xfe040,    0xfe03f,    0x10203f,   0x101fc0,   0xfdfc0,
      0xfdfbf,    0x101fbf,   0x102041,   0xfe041,    0xfe0c0,    0x1020c0,   0x1020bf,   0xfe0bf,    0xfdfc1,    0x101fc1,   0x8106040,  0x1020c1,   0xfe0c1,    0x810603f,  0x8105fc0,  0x8105fbf,  0x180fa040,
      0x180fa03f, 0x180f9fc0, 0x180f9fbf, 0x8106041,  0x81060c0,  0x81060bf,  0x8105fc1,  0x180fa041, 0x180fa0c0, 0x180fa0bf, 0x6fe03e,   0x70203e,   0x3101f40,  0x30fdf40,  0x180f9fc1, 0x30fdf3f,  0x3101f3f,
      0x3701fbe,  0x36fdfbe,  0x93060c1,  0x192fa0c1, 0x6fe0be,   0x7020be,   0x3101f41,  0x30fdf41,  0xb305f40,  0xb105f3f,  0x970603e,  0xb705fbe,  0x196fa03e, 0x186f9fbe, 0x1b2f9f40, 0x1b6f9f3f, 0x192fe042,
      0x9302042,  0xb301fc2,  0x1b2fdfc2, 0x192fe140, 0x9302140,  0x970213f,  0x196fe13f, 15,         0xfe040,    0xfdfc0,    0x101fc0,   0x10203f,   0xfe03f,    0xfdfbf,    0x101fbf,   0x102041,   0xfe041,
      0xfdfc1,    0x101fc1,   0x1020c0,   0xfe0c0,    0xfe0bf,    0x1020bf,   0x1020c1,   0xfe0c1,    0x8106040,  0x8105fc0,  0x810603f,  0x8105fbf,  0x180fa040, 0x180f9fc0, 0x8106041,  0x8105fc1,  0x81060c0,
      0x180fa03f, 0x180f9fbf, 0x180fa041, 0x180f9fc1, 0x180fa0c0, 0x81060bf,  0x91060c1,  0x3101f40,  0x30fdf40,  0x30fdf3f,  0x180fa0bf, 0x190fa0c1, 0x3101f3f,  0x3101f41,  0x30fdf41,  0x186fe03e, 0x70203e,
      0x3701fbe,  0x1b6fdfbe, 0x186fe0be, 0x7020be,   0x8302042,  0x182fe042, 0x1b2fdfc2, 0xb301fc2,  0xb105f40,  0xb705f3f,  0xb305f41,  0x1b2f9f40, 0x1b6f9f3f, 0x192fe140, 0x9302140,  0x93020c2,  0x192fe0c2,
      0x196fe13f, 0x970213f,  0xa70603e,  11,         0x10203f,   0xfe040,    0xfe03f,    0xfdfbf,    0x101fbf,   0x101fc0,   0xfdfc0,    0xfe0c0,    0x1020c0,   0x1020bf,   0xfe0bf,    0x810603f,  0x8106040,
      0x302041,   0x2fe041,   0x2fdfc1,   0x301fc1,   0x8105fc0,  0x8105fbf,  0x70203e,   0x6fe03e,   0x180fa03f, 0x180fa040, 0x2fe0c1,   0x3020c1,   0x81060c0,  0x81060bf,  0x701fbe,   0x6fdfbe,   0x180f9fbf,
      0x180f9fc0, 0x180fa0c0, 0x180fa0bf, 0x6fe0be,   0x7020be,   0x8306041,  0x8305fc1,  0x3101f40,  0x3101f3f,  0x30fdf3f,  0x30fdf40,  0x182fa041, 0x870603e,  0x8705fbe,  0x93060c1,  0x182f9fc1, 0x192fa0c1,
      0x186fa03e, 0x186f9fbe, 0x97060be,  0x970213f,  0x9302140,  0x192fe140, 0x196fe13f, 0x196fa0be, 0x32fdf41,  0x3301f41,  0xb305f40,  0xb105f3f,  0xb701f3e,  0x36fdf3e,  0x1b6f9f3f, 0x1b2f9f40, 11,
      0xfe040,    0x10203f,   0xfe03f,    0xfdfc0,    0x101fc0,   0x101fbf,   0xfdfbf,    0xfe0c0,    0x1020c0,   0x1020bf,   0xfe0bf,    0x2fe041,   0x302041,   0x301fc1,   0x2fdfc1,   0x8106040,  0x810603f,
      0x3020c1,   0x2fe0c1,   0x8105fc0,  0x8105fbf,  0x180fa040, 0x180fa03f, 0x81060c0,  0x81060bf,  0x70203e,   0x180f9fc0, 0x180f9fbf, 0x6fe03e,   0x180fa0c0, 0x180fa0bf, 0x8306041,  0x701fbe,   0x6fdfbe,
      0x6fe0be,   0x7020be,   0x8305fc1,  0x182fa041, 0x83060c1,  0x182f9fc1, 0x1b0fdf40, 0x3101f40,  0x3701f3f,  0x1b6fdf3f, 0x182fa0c1, 0x190fe140, 0x870603e,  0x8705fbe,  0x1102140,  0x170213f,  0x196fe13f,
      0x186fa03e, 0x1b6f9fbe, 0x87060be,  0x3301f41,  0x1b2fdf41, 0xb305f40,  0xb705f3f,  0x196fa0be, 0x192fe141, 0x1302141,  0x9306140,  0x970613f,  15,         0xfe040,    0xfe03f,    0x10203f,   0x101fc0,
      0xfdfc0,    0xfdfbf,    0x101fbf,   0x1020c0,   0xfe0c0,    0xfe0bf,    0x1020bf,   0x102041,   0xfe041,    0xfdfc1,    0x101fc1,   0x1020c1,   0xfe0c1,    0x8106040,  0x810603f,  0x8105fc0,  0x8105fbf,
      0x180fa040, 0x180fa03f, 0x81060c0,  0x81060bf,  0x8106041,  0x180f9fc0, 0x180f9fbf, 0x180fa0c0, 0x180fa0bf, 0x180fa041, 0x8105fc1,  0x83060c1,  0x70203e,   0x6fe03e,   0x6fdfbe,   0x180f9fc1, 0x182fa0c1,
      0x701fbe,   0x7020be,   0x6fe0be,   0x1b0fdf40, 0x3101f40,  0x3701f3f,  0x1b6fdf3f, 0x1b0fdf41, 0x3101f41,  0x9102140,  0x190fe140, 0x196fe13f, 0x970213f,  0x870603e,  0xb705fbe,  0x97060be,  0x196fa03e,
      0x1b6f9fbe, 0x192fe042, 0x9302042,  0x9302141,  0x192fe141, 0x1b2fdfc2, 0xb301fc2,  0xb505f40,  17,         0xfe040,    0xfe03f,    0x10203f,   0x101fc0,   0xfdfc0,    0xfe041,    0x102041,   0x1020c0,
      0xfe0c0,    0xfdfbf,    0x101fbf,   0x1020bf,   0xfe0bf,    0xfdfc1,    0x101fc1,   0x1020c1,   0xfe0c1,    0x8106040,  0x810603f,  0x8105fc0,  0x8106041,  0x81060c0,  0x180fa040, 0x180fa03f, 0x180f9fc0,
      0x8105fbf,  0x81060bf,  0x8105fc1,  0x180fa041, 0x180fa0c0, 0x180f9fbf, 0x81060c1,  0x180fa0bf, 0x180f9fc1, 0x180fa0c1, 0x186fe03e, 0x70203e,   0x3101f40,  0x1b0fdf40, 0x1b0fdf3f, 0x3101f3f,  0x3701fbe,
      0x1b6fdfbe, 0x186fe0be, 0x7020be,   0x9102140,  0x190fe140, 0x182fe042, 0x8302042,  0x3101f41,  0x1b0fdf41, 0x1b2fdfc2, 0xb301fc2,  0x83020c2,  0x182fe0c2, 0x196fe13f, 0x970213f,  0x9302141,  0x192fe141,
      0x970603e,  0xb305f40,  0xb105f3f,  0xb705fbe,  11,         0x10203f,   0x101fc0,   0x101fbf,   0xfe040,    0xfe03f,    0xfdfc0,    0xfdfbf,    0x105fbf,   0x10603f,   0x106040,   0x105fc0,   0x301fc1,
      0x302041,   0x11020c0,  0x11020bf,  0x10fe0bf,  0x10fe0c0,  0x2fe041,   0x2fdfc1,   0x3101f40,  0x3101f3f,  0x701fbe,   0x70203e,   0x6fe03e,   0x6fdfbe,   0x30fdf3f,  0x30fdf40,  0x180f9fc0, 0x180fa040,
      0x180fa03f, 0x180f9fbf, 0x11060bf,  0x11060c0,  0x306041,   0x305fc1,   0x3105f40,  0x3105f3f,  0x705fbe,   0x70603e,   0x13020c1,  0x12fe0c1,  0x3301f41,  0x32fdf41,  0x17020be,  0x16fe0be,  0x190fa0bf,
      0x190fa0c0, 0x192fa041, 0x182f9fc1, 0x3701f3e,  0x36fdf3e,  0x186f9fbe, 0x196fa03e, 0x1b6f9f3f, 0x1b2f9f40, 0x93060c1,  0x97060be,  0xb305f41,  0xb705f3e,  0xb709fbf,  0x970a03f,  0x930a040,  0xb309fc0,
      9,          0x101fc0,   0x10203f,   0xfe040,    0xfdfc0,    0x101fbf,   0xfdfbf,    0xfe03f,    0x102041,   0x101fc1,   0x8105fc0,  0x8106040,  0x810603f,  0x8105fbf,  0xfdfc1,    0xfe041,    0x11020c0,
      0x11020bf,  0x10fe0c0,  0x10fe0bf,  0x3101f40,  0x30fdf40,  0x3101f3f,  0x30fdf3f,  0x8105fc1,  0x8106041,  0x70203e,   0x701fbe,   0x6fdfbe,   0x6fe03e,   0x180fa03f, 0x180fa040, 0x180f9fc0, 0x180f9fbf,
      0x11020c1,  0x91060c0,  0x91060bf,  0x12fe0c1,  0x3101f41,  0xb105f40,  0xb105f3f,  0x30fdf41,  0x180f9fc1, 0x182fa041, 0x870603e,  0x8705fbe,  0x93060c1,  0x17020be,  0x16fe0be,  0x190fa0bf, 0x190fa0c0,
      0xb305f41,  0x1b0f9f40, 0x1b0f9f3f, 0x36fdf3e,  0xb701f3e,  0x1b6f9fbe, 0x196fa03e, 0x97060be,  0x192fa0c1, 0x1b2f9f41, 0x1b2fdfc2, 0xb301fc2,  0x11302042, 11,         0x101fc0,   0xfe040,    0xfdfc0,
      0x101fbf,   0x10203f,   0xfe03f,    0xfdfbf,    0x101fc1,   0x102041,   0xfe041,    0xfdfc1,    0x8105fc0,  0x8106040,  0x810603f,  0x8105fbf,  0x11020c0,  0x10fe0c0,  0x11020bf,  0x10fe0bf,  0x8106041,
      0x8105fc1,  0x3101f40,  0x30fdf40,  0x3101f3f,  0x30fdf3f,  0x180f9fc0, 0x180fa040, 0x11020c1,  0x10fe0c1,  0x180fa03f, 0x180f9fbf, 0x91060c0,  0x91060bf,  0x3101f41,  0x30fdf41,  0x701fbe,   0x70203e,
      0x6fe03e,   0x6fdfbe,   0xb105f40,  0xb105f3f,  0x180f9fc1, 0x180fa041, 0x93060c1,  0xb305f41,  0x190fa0c0, 0x190fa0bf, 0x870603e,  0x8705fbe,  0x97020be,  0x16fe0be,  0x1b0f9f40, 0x1b0f9f3f, 0x36fdf3e,
      0xb701f3e,  0x192fa0c1, 0x192fe042, 0x9302042,  0xb301fc2,  0x1b2fdfc2, 0x1b2f9f41, 0x1b6f9fbe, 0x196fa03e, 11,         0x101fc0,   0xfe040,    0xfdfc0,    0x101fbf,   0x10203f,   0xfe03f,    0xfdfbf,
      0x101fc1,   0x102041,   0xfe041,    0xfdfc1,    0x8105fc0,  0x8106040,  0x810603f,  0x8105fbf,  0x11020c0,  0x10fe0c0,  0x8106041,  0x8105fc1,  0x11020bf,  0x10fe0bf,  0x3101f40,  0x30fdf40,  0x11020c1,
      0x10fe0c1,  0x180fa040, 0x3101f3f,  0x30fdf3f,  0x180f9fc0, 0x3101f41,  0x30fdf41,  0x91060c0,  0x180fa03f, 0x180f9fbf, 0x180f9fc1, 0x180fa041, 0x91060bf,  0xb105f40,  0x91060c1,  0xb105f3f,  0x3701fbe,
      0x70203e,   0x186fe03e, 0x1b6fdfbe, 0xb105f41,  0x3301fc2,  0x190fa0c0, 0x190fa0bf, 0x302042,   0x182fe042, 0x1b2fdfc2, 0x1b0f9f40, 0x1b6f9f3f, 0x190fa0c1, 0x870603e,  0xb705fbe,  0x97020be,  0x196fe0be,
      0x1b2f9f41, 0xb305fc2,  0x8306042,  0x93020c2,  0x192fe0c2, 9,          0x10203f,   0x101fc0,   0xfe040,    0xfe03f,    0x101fbf,   0xfdfbf,    0xfdfc0,    0x1020c0,   0x1020bf,   0x810603f,  0x8106040,
      0x8105fc0,  0x8105fbf,  0xfe0bf,    0xfe0c0,    0x302041,   0x301fc1,   0x2fe041,   0x2fdfc1,   0x70203e,   0x6fe03e,   0x701fbe,   0x6fdfbe,   0x81060bf,  0x81060c0,  0x3101f40,  0x3101f3f,  0x30fdf3f,
      0x30fdf40,  0x180f9fc0, 0x180fa040, 0x180fa03f, 0x180f9fbf, 0x3020c1,   0x8306041,  0x8305fc1,  0x12fe0c1,  0x7020be,   0x870603e,  0x8705fbe,  0x6fe0be,   0x180fa0bf, 0x190fa0c0, 0xb105f40,  0xb105f3f,
      0x93060c1,  0x3301f41,  0x32fdf41,  0x182f9fc1, 0x182fa041, 0x97060be,  0x186fa03e, 0x186f9fbe, 0x36fdf3e,  0xb701f3e,  0x1b6f9f3f, 0x1b2f9f40, 0xb305f41,  0x192fa0c1, 0x196fa0be, 0x196fe13f, 0x970213f,
      0x11302140, 7,          0x10203f,   0x101fc0,   0xfe040,    0xfe03f,    0x101fbf,   0xfdfc0,    0xfdfbf,    0x302041,   0x1020c0,   0x8106040,  0x810603f,  0x1020bf,   0xfe0c0,    0x2fe041,   0x301fc1,
      0x8105fc0,  0x8105fbf,  0xfe0bf,    0x2fdfc1,   0x3020c1,   0x8306041,  0x81060c0,  0x81060bf,  0x2fe0c1,   0x8305fc1,  0x3101f40,  0x3101f3f,  0x701fbe,   0x70203e,   0x6fe03e,   0x180fa03f, 0x180fa040,
      0x180f9fc0, 0x30fdf40,  0x30fdf3f,  0x6fdfbe,   0x180f9fbf, 0x180fa0c0, 0x182fa041, 0x93060c1,  0x870603e,  0x7020be,   0x6fe0be,   0x180fa0bf, 0x182f9fc1, 0x32fdf41,  0x3301f41,  0xb105f40,  0xb105f3f,
      0x8705fbe,  0x97060be,  0x192fa0c1, 0xb305f41,  0xb701f3e,  0x36fdf3e,  0x1b0f9f3f, 0x1b2f9f40, 0x1b6f9fbe, 0x196fa03e, 0x196fe13f, 0x9302140,  0x970213f,  0x192fe140, 11,         0x101fc0,   0xfe040,
      0xfdfc0,    0x10203f,   0x101fbf,   0xfe03f,    0xfdfbf,    0x102041,   0x101fc1,   0xfe041,    0xfdfc1,    0x11020c0,  0x8106040,  0x8105fc0,  0x10fe0c0,  0x11020bf,  0x810603f,  0x8105fbf,  0x10fe0bf,
      0x11020c1,  0x8106041,  0x8105fc1,  0x10fe0c1,  0x91060c0,  0x91060bf,  0x3101f40,  0x30fdf40,  0x180f9fc0, 0x180fa040, 0x180fa03f, 0x180f9fbf, 0x30fdf3f,  0x3101f3f,  0x701fbe,   0x70203e,   0x6fe03e,
      0x6fdfbe,   0x93060c1,  0x3101f41,  0x30fdf41,  0x180f9fc1, 0x180fa041, 0x190fa0c0, 0x190fa0bf, 0xb105f40,  0xb105f3f,  0x8705fbe,  0x870603e,  0x17020be,  0x16fe0be,  0x192fa0c1, 0xb305f41,  0xb301fc2,
      0x9302042,  0x192fe042, 0x1b2fdfc2, 0x1b2f9f40, 0x1b0f9f3f, 0x186f9fbe, 0x196fa03e, 0x97060be,  0xb701f3e,  0x1b6fdf3e, 11,         0x101fc0,   0xfe040,    0xfdfc0,    0x10203f,   0x101fbf,   0xfe03f,
      0xfdfbf,    0x102041,   0x101fc1,   0xfe041,    0xfdfc1,    0x11020c0,  0x8106040,  0x8105fc0,  0x10fe0c0,  0x11020bf,  0x810603f,  0x8105fbf,  0x10fe0bf,  0x11020c1,  0x8106041,  0x8105fc1,  0x10fe0c1,
      0x91060c0,  0x91060bf,  0x3101f40,  0x30fdf40,  0x180f9fc0, 0x180fa040, 0x180fa03f, 0x180f9fbf, 0x30fdf3f,  0x3101f3f,  0x3101f41,  0x91060c1,  0x180fa041, 0x180f9fc1, 0x30fdf41,  0xb105f40,  0x70203e,
      0x3701fbe,  0x186fe03e, 0x1b6fdfbe, 0x190fa0c0, 0x190fa0bf, 0x190fa0c1, 0xb105f3f,  0xb105f41,  0x3301fc2,  0x302042,   0x182fe042, 0x1b2fdfc2, 0x1b0f9f40, 0x17020be,  0x970603e,  0xb705fbe,  0x196fe0be,
      0x1b6f9f3f, 0x1b2f9f41, 0x13020c2,  0x9306042,  0xb305fc2,  0x192fe0c2, 11,         0x10203f,   0xfe040,    0xfe03f,    0x101fbf,   0x101fc0,   0xfdfc0,    0xfdfbf,    0x1020bf,   0x1020c0,   0xfe0c0,
      0xfe0bf,    0x810603f,  0x8106040,  0x8105fc0,  0x8105fbf,  0x302041,   0x2fe041,   0x301fc1,   0x2fdfc1,   0x81060c0,  0x81060bf,  0x70203e,   0x6fe03e,   0x701fbe,   0x6fdfbe,   0x180fa03f, 0x180fa040,
      0x3020c1,   0x2fe0c1,   0x180f9fc0, 0x180f9fbf, 0x8306041,  0x8305fc1,  0x7020be,   0x6fe0be,   0x3101f3f,  0x3101f40,  0x30fdf40,  0x30fdf3f,  0x870603e,  0x8705fbe,  0x180fa0bf, 0x180fa0c0, 0x93060c1,
      0x97060be,  0x182fa041, 0x182f9fc1, 0xb105f40,  0xb105f3f,  0xb301f41,  0x32fdf41,  0x186fa03e, 0x186f9fbe, 0x36fdf3e,  0xb701f3e,  0x192fa0c1, 0x192fe140, 0x9302140,  0x970213f,  0x196fe13f, 0x196fa0be,
      0x1b6f9f3f, 0x1b2f9f40, 11,         0x10203f,   0xfe040,    0xfe03f,    0x101fc0,   0x101fbf,   0xfdfc0,    0xfdfbf,    0x1020c0,   0x1020bf,   0xfe0c0,    0xfe0bf,    0x302041,   0x8106040,  0x810603f,
      0x2fe041,   0x301fc1,   0x8105fc0,  0x8105fbf,  0x2fdfc1,   0x3020c1,   0x81060c0,  0x81060bf,  0x2fe0c1,   0x8306041,  0x8305fc1,  0x70203e,   0x6fe03e,   0x180fa03f, 0x180fa040, 0x180f9fc0, 0x180f9fbf,
      0x6fdfbe,   0x701fbe,   0x3101f3f,  0x3101f40,  0x30fdf40,  0x30fdf3f,  0x93060c1,  0x7020be,   0x6fe0be,   0x180fa0bf, 0x180fa0c0, 0x182fa041, 0x182f9fc1, 0x870603e,  0x8705fbe,  0xb105f3f,  0xb105f40,
      0x3301f41,  0x32fdf41,  0x192fa0c1, 0x97060be,  0x970213f,  0x9302140,  0x192fe140, 0x196fe13f, 0x196fa03e, 0x186f9fbe, 0x1b0f9f3f, 0x1b2f9f40, 0xb305f41,  0xb701f3e,  0x1b6fdf3e, 15,         0xfe040,
      0x10203f,   0x101fc0,   0xfdfc0,    0xfe03f,    0x101fbf,   0xfdfbf,    0x102041,   0x1020c0,   0xfe0c0,    0xfe041,    0x101fc1,   0xfdfc1,    0x1020bf,   0xfe0bf,    0x8106040,  0x810603f,  0x8105fc0,
      0x8105fbf,  0x1020c1,   0xfe0c1,    0x8106041,  0x81060c0,  0x81060bf,  0x8105fc1,  0x180fa040, 0x180f9fc0, 0x180fa03f, 0x180f9fbf, 0x93060c1,  0x70203e,   0x6fe03e,   0x3101f40,  0x1b0fdf40, 0x3101f3f,
      0x3701fbe,  0x6fdfbe,   0x1b6fdf3f, 0x180fa041, 0x180fa0c0, 0x180fa0bf, 0x180f9fc1, 0x1b0fdf41, 0x3101f41,  0x7020be,   0x6fe0be,   0x870603e,  0x8705fbe,  0xb105f40,  0xb705f3f,  0x192fa0c1, 0x192fe140,
      0x9302140,  0x970213f,  0x97060be,  0x9302042,  0x192fe042, 0xb301fc2,  0xb305f41,  0x1b2fdfc2, 0x196fe13f, 0x196fa03e, 0x1b6f9fbe, 14,         0xfe040,    0x101fc0,   0xfdfc0,    0x10203f,   0xfe03f,
      0x101fbf,   0xfdfbf,    0x102041,   0xfe041,    0x101fc1,   0xfdfc1,    0x1020c0,   0xfe0c0,    0x1020bf,   0x8106040,  0xfe0bf,    0x8105fc0,  0x1020c1,   0xfe0c1,    0x810603f,  0x8105fbf,  0x8106041,
      0x8105fc1,  0x81060c0,  0x81060bf,  0x91060c1,  0x180fa040, 0x180f9fc0, 0x180fa03f, 0x180f9fbf, 0x180fa041, 0x180f9fc1, 0x1b0fdf40, 0x3101f40,  0x3101f3f,  0x1b0fdf3f, 0x180fa0c0, 0x190fa0bf, 0x186fe03e,
      0x70203e,   0x3701fbe,  0x3101f41,  0x1b0fdf41, 0x1b6fdfbe, 0x190fa0c1, 0x182fe042, 0x302042,   0xb105f40,  0xb105f3f,  0x3301fc2,  0x1b2fdfc2, 0x7020be,   0x196fe0be, 0x970603e,  0xb705fbe,  0xb105f41,
      0x9302140,  0x192fe140, 0x13020c2,  0x192fe0c2, 0x9306042,  0xb305fc2,  0x1170213f, 11,         0x10203f,   0xfe040,    0xfe03f,    0x101fbf,   0x101fc0,   0xfdfc0,    0xfdfbf,    0x1020bf,   0x1020c0,
      0xfe0c0,    0xfe0bf,    0x810603f,  0x8106040,  0x8105fc0,  0x8105fbf,  0x302041,   0x2fe041,   0x81060c0,  0x81060bf,  0x301fc1,   0x2fdfc1,   0x70203e,   0x6fe03e,   0x3020c1,   0x2fe0c1,   0x180fa040,
      0x701fbe,   0x6fdfbe,   0x180fa03f, 0x7020be,   0x6fe0be,   0x8306041,  0x180f9fc0, 0x180f9fbf, 0x180fa0bf, 0x180fa0c0, 0x8305fc1,  0x870603e,  0x83060c1,  0x8705fbe,  0x3701f3f,  0x3101f40,  0x1b0fdf40,
      0x1b6fdf3f, 0x87060be,  0x170213f,  0x182fa041, 0x182f9fc1, 0x1102140,  0x190fe140, 0x196fe13f, 0x186fa03e, 0x1b6f9fbe, 0x182fa0c1, 0xb105f40,  0xb705f3f,  0xb301f41,  0x1b2fdf41, 0x196fa0be, 0x970613f,
      0x9106140,  0x9302141,  0x192fe141, 11,         0x10203f,   0xfe040,    0xfe03f,    0x101fc0,   0x101fbf,   0xfdfc0,    0xfdfbf,    0x1020c0,   0x1020bf,   0xfe0c0,    0xfe0bf,    0x302041,   0x8106040,
      0x810603f,  0x2fe041,   0x301fc1,   0x8105fc0,  0x8105fbf,  0x2fdfc1,   0x3020c1,   0x81060c0,  0x81060bf,  0x2fe0c1,   0x8306041,  0x8305fc1,  0x70203e,   0x6fe03e,   0x180fa03f, 0x180fa040, 0x180f9fc0,
      0x180f9fbf, 0x6fdfbe,   0x701fbe,   0x7020be,   0x83060c1,  0x180fa0c0, 0x180fa0bf, 0x6fe0be,   0x870603e,  0x3101f40,  0x3701f3f,  0x1b0fdf40, 0x1b6fdf3f, 0x182fa041, 0x182f9fc1, 0x182fa0c1, 0x8705fbe,
      0x87060be,  0x170213f,  0x1102140,  0x190fe140, 0x196fe13f, 0x186fa03e, 0x3301f41,  0xb305f40,  0xb705f3f,  0x1b2fdf41, 0x1b6f9fbe, 0x196fa0be, 0x1302141,  0x9306140,  0x970613f,  0x192fe141, 14,
      0xfe040,    0x10203f,   0xfe03f,    0x101fc0,   0xfdfc0,    0x101fbf,   0xfdfbf,    0x1020c0,   0xfe0c0,    0x1020bf,   0xfe0bf,    0x102041,   0xfe041,    0x101fc1,   0x8106040,  0xfdfc1,    0x810603f,
      0x1020c1,   0xfe0c1,    0x8105fc0,  0x8105fbf,  0x81060c0,  0x81060bf,  0x8106041,  0x8105fc1,  0x83060c1,  0x180fa040, 0x180fa03f, 0x180f9fc0, 0x180f9fbf, 0x180fa0c0, 0x180fa0bf, 0x186fe03e, 0x70203e,
      0x701fbe,   0x186fdfbe, 0x180fa041, 0x182f9fc1, 0x1b0fdf40, 0x3101f40,  0x3701f3f,  0x7020be,   0x186fe0be, 0x1b6fdf3f, 0x182fa0c1, 0x190fe140, 0x1102140,  0x870603e,  0x8705fbe,  0x170213f,  0x196fe13f,
      0x3101f41,  0x1b2fdf41, 0xb305f40,  0xb705f3f,  0x87060be,  0x9302042,  0x192fe042, 0x1302141,  0x192fe141, 0x9306140,  0x970613f,  0x13301fc2, 17,         0xfe040,    0x10203f,   0x101fc0,   0xfdfc0,
      0xfe03f,    0x1020c0,   0x102041,   0xfe041,    0xfe0c0,    0x101fbf,   0xfdfbf,    0x1020bf,   0xfe0bf,    0x101fc1,   0xfdfc1,    0x1020c1,   0xfe0c1,    0x8106040,  0x810603f,  0x8105fc0,  0x8106041,
      0x81060c0,  0x8105fbf,  0x81060bf,  0x8105fc1,  0x81060c1,  0x180fa040, 0x180f9fc0, 0x180fa03f, 0x180fa0c0, 0x180fa041, 0x180f9fbf, 0x180fa0bf, 0x180f9fc1, 0x180fa0c1, 0x70203e,   0x186fe03e, 0x3101f40,
      0x1b0fdf40, 0x3101f3f,  0x3701fbe,  0x186fdfbe, 0x1b6fdf3f, 0x3101f41,  0x1b0fdf41, 0x8302042,  0x182fe042, 0x9102140,  0x7020be,   0x186fe0be, 0x190fe140, 0x970213f,  0x196fe13f, 0x9102141,  0xb301fc2,
      0x1b2fdfc2, 0x93020c2,  0x182fe0c2, 0x192fe141, 0x970603e,  0xb305f40,  0xb105f3f,  0xb705fbe,  11,         0x10203f,   0x101fc0,   0x101fbf,   0xfdfbf,    0xfe03f,    0xfe040,    0xfdfc0,    0x105fc0,
      0x106040,   0x10603f,   0x105fbf,   0x11020bf,  0x11020c0,  0x302041,   0x301fc1,   0x2fdfc1,   0x2fe041,   0x10fe0c0,  0x10fe0bf,  0x70203e,   0x701fbe,   0x3101f3f,  0x3101f40,  0x30fdf40,  0x30fdf3f,
      0x6fdfbe,   0x6fe03e,   0x11060bf,  0x11060c0,  0x306041,   0x305fc1,   0x3105f40,  0x3105f3f,  0x705fbe,   0x70603e,   0x180fa03f, 0x180fa040, 0x180f9fc0, 0x180f9fbf, 0x13020c1,  0x12fe0c1,  0x17020be,
      0x16fe0be,  0x3301f41,  0x32fdf41,  0x93060c1,  0x3701f3e,  0x36fdf3e,  0x97060be,  0xb305f41,  0x182f9fc1, 0x182fa041, 0x192fa0c0, 0x190fa0bf, 0x196fa03e, 0x186f9fbe, 0xb705f3e,  0xb709fbf,  0xb309fc0,
      0x930a040,  0x970a03f,  0x1b6f9f3f, 0x1b2f9f40, 11,         0x101fc0,   0x10203f,   0x101fbf,   0xfdfc0,    0xfe040,    0xfe03f,    0xfdfbf,    0x105fc0,   0x106040,   0x10603f,   0x105fbf,   0x301fc1,
      0x302041,   0x2fe041,   0x2fdfc1,   0x11020c0,  0x11020bf,  0x10fe0c0,  0x10fe0bf,  0x306041,   0x305fc1,   0x3101f40,  0x3101f3f,  0x30fdf40,  0x30fdf3f,  0x701fbe,   0x70203e,   0x11060c0,  0x11060bf,
      0x6fe03e,   0x6fdfbe,   0x13020c1,  0x12fe0c1,  0x3105f40,  0x3105f3f,  0x180f9fc0, 0x180fa040, 0x180fa03f, 0x180f9fbf, 0x3301f41,  0x32fdf41,  0x705fbe,   0x70603e,   0x93060c1,  0xb305f41,  0x17020be,
      0x16fe0be,  0x182fa041, 0x182f9fc1, 0x192fa0c0, 0x190fa0bf, 0x3701f3e,  0x36fdf3e,  0x1b0f9f3f, 0x1b2f9f40, 0x97060be,  0x970a03f,  0x930a040,  0xb309fc0,  0xb709fbf,  0xb705f3e,  0x1b6f9fbe, 0x196fa03e,
      15,         0x101fc0,   0x101fbf,   0x10203f,   0xfe040,    0xfdfc0,    0xfdfbf,    0xfe03f,    0x102041,   0x101fc1,   0x105fc0,   0x106040,   0x10603f,   0x105fbf,   0xfdfc1,    0xfe041,    0x11020c0,
      0x106041,   0x105fc1,   0x11020bf,  0x10fe0c0,  0x10fe0bf,  0x3101f40,  0x3101f3f,  0x30fdf40,  0x30fdf3f,  0x11020c1,  0x11060c0,  0x11060bf,  0x10fe0c1,  0x3101f41,  0x3105f40,  0x3105f3f,  0x701fbe,
      0x70203e,   0x180fa040, 0x180f9fc0, 0x30fdf41,  0x180f9fbf, 0x180fa03f, 0x186fe03e, 0x186fdfbe, 0x93060c1,  0xb305f41,  0x705fbe,   0x70603e,   0x180fa041, 0x180f9fc1, 0x192fa0c0, 0x190fa0bf, 0x97020be,
      0x196fe0be, 0xb701f3e,  0x36fdf3e,  0x1b2f9f40, 0x1b6f9f3f, 0xb301fc2,  0x9302042,  0x192fe042, 0x1b2fdfc2, 0xb309fc0,  0x930a040,  0x970a03f,  0xb709fbf,  15,         0x101fc0,   0xfdfc0,    0xfe040,
      0x10203f,   0x101fbf,   0xfdfbf,    0xfe03f,    0x102041,   0x101fc1,   0xfdfc1,    0xfe041,    0x106040,   0x105fc0,   0x105fbf,   0x10603f,   0x106041,   0x105fc1,   0x11020c0,  0x10fe0c0,  0x11020bf,
      0x10fe0bf,  0x3101f40,  0x30fdf40,  0x11020c1,  0x10fe0c1,  0x11060c0,  0x3101f3f,  0x30fdf3f,  0x3101f41,  0x30fdf41,  0x3105f40,  0x11060bf,  0x91060c1,  0x180fa040, 0x180f9fc0, 0x180f9fbf, 0x3105f3f,
      0xb105f41,  0x180fa03f, 0x180fa041, 0x180f9fc1, 0x3701fbe,  0x70203e,   0x186fe03e, 0x1b6fdfbe, 0x3705fbe,  0x70603e,   0x1302042,  0x3301fc2,  0x1b2fdfc2, 0x192fe042, 0x190fa0c0, 0x196fa0bf, 0x192fa0c1,
      0x1b2f9f40, 0x1b6f9f3f, 0xb309fc0,  0x930a040,  0x9306042,  0xb305fc2,  0xb709fbf,  0x970a03f,  0x117020be, 11,         0x10203f,   0x101fc0,   0x101fbf,   0xfe03f,    0xfe040,    0xfdfc0,    0xfdfbf,
      0x10603f,   0x106040,   0x105fc0,   0x105fbf,   0x11020bf,  0x11020c0,  0x10fe0c0,  0x10fe0bf,  0x302041,   0x301fc1,   0x2fe041,   0x2fdfc1,   0x11060c0,  0x11060bf,  0x70203e,   0x701fbe,   0x6fe03e,
      0x6fdfbe,   0x3101f3f,  0x3101f40,  0x306041,   0x305fc1,   0x30fdf40,  0x30fdf3f,  0x13020c1,  0x12fe0c1,  0x70603e,   0x705fbe,   0x180fa03f, 0x180fa040, 0x180f9fc0, 0x180f9fbf, 0x17020be,  0x16fe0be,
      0x3105f3f,  0x3105f40,  0x93060c1,  0x97060be,  0x3301f41,  0x32fdf41,  0x190fa0c0, 0x190fa0bf, 0x192fa041, 0x182f9fc1, 0x3701f3e,  0x36fdf3e,  0x186f9fbe, 0x196fa03e, 0xb305f41,  0xb309fc0,  0x930a040,
      0x970a03f,  0xb709fbf,  0xb705f3e,  0x1b6f9f3f, 0x1b2f9f40, 11,         0x10203f,   0x101fc0,   0x101fbf,   0xfe040,    0xfe03f,    0xfdfc0,    0xfdfbf,    0x106040,   0x10603f,   0x105fc0,   0x105fbf,
      0x302041,   0x11020c0,  0x11020bf,  0x301fc1,   0x2fe041,   0x10fe0c0,  0x10fe0bf,  0x2fdfc1,   0x306041,   0x11060c0,  0x11060bf,  0x305fc1,   0x13020c1,  0x12fe0c1,  0x70203e,   0x701fbe,   0x3101f3f,
      0x3101f40,  0x30fdf40,  0x30fdf3f,  0x6fdfbe,   0x6fe03e,   0x180fa03f, 0x180fa040, 0x180f9fc0, 0x180f9fbf, 0x93060c1,  0x70603e,   0x705fbe,   0x3105f3f,  0x3105f40,  0x3301f41,  0x32fdf41,  0x17020be,
      0x16fe0be,  0x190fa0bf, 0x190fa0c0, 0x182fa041, 0x182f9fc1, 0xb305f41,  0x97060be,  0x970a03f,  0x930a040,  0xb309fc0,  0xb709fbf,  0xb701f3e,  0x36fdf3e,  0x1b0f9f3f, 0x1b2f9f40, 0x192fa0c1, 0x196fa03e,
      0x1b6f9fbe, 15,         0x101fc0,   0x10203f,   0xfe040,    0xfdfc0,    0x101fbf,   0xfe03f,    0xfdfbf,    0x102041,   0x106040,   0x105fc0,   0x101fc1,   0xfe041,    0xfdfc1,    0x10603f,   0x105fbf,
      0x11020c0,  0x11020bf,  0x10fe0c0,  0x10fe0bf,  0x106041,   0x105fc1,   0x11020c1,  0x11060c0,  0x11060bf,  0x10fe0c1,  0x3101f40,  0x30fdf40,  0x3101f3f,  0x30fdf3f,  0x93060c1,  0x70203e,   0x701fbe,
      0x180fa040, 0x1b0f9fc0, 0x180fa03f, 0x186fe03e, 0x6fdfbe,   0x1b6f9fbf, 0x3101f41,  0x3105f40,  0x3105f3f,  0x30fdf41,  0x1b0f9fc1, 0x180fa041, 0x70603e,   0x705fbe,   0x17020be,  0x16fe0be,  0x190fa0c0,
      0x196fa0bf, 0xb305f41,  0xb309fc0,  0x930a040,  0x970a03f,  0x97060be,  0x9302042,  0xb301fc2,  0x192fe042, 0x192fa0c1, 0x1b2fdfc2, 0xb709fbf,  0xb701f3e,  0x1b6fdf3e, 14,         0x101fc0,   0xfe040,
      0xfdfc0,    0x10203f,   0x101fbf,   0xfe03f,    0xfdfbf,    0x102041,   0x101fc1,   0xfe041,    0xfdfc1,    0x106040,   0x105fc0,   0x10603f,   0x11020c0,  0x105fbf,   0x10fe0c0,  0x106041,   0x105fc1,
      0x11020bf,  0x10fe0bf,  0x11020c1,  0x10fe0c1,  0x11060c0,  0x11060bf,  0x91060c1,  0x3101f40,  0x30fdf40,  0x3101f3f,  0x30fdf3f,  0x3101f41,  0x30fdf41,  0x1b0f9fc0, 0x180fa040, 0x180fa03f, 0x1b0f9fbf,
      0x3105f40,  0xb105f3f,  0x3701fbe,  0x70203e,   0x186fe03e, 0x180fa041, 0x1b0f9fc1, 0x1b6fdfbe, 0xb105f41,  0x3301fc2,  0x302042,   0x190fa0c0, 0x190fa0bf, 0x182fe042, 0x1b2fdfc2, 0x70603e,   0xb705fbe,
      0x97020be,  0x196fe0be, 0x190fa0c1, 0x930a040,  0xb309fc0,  0x8306042,  0xb305fc2,  0x93020c2,  0x192fe0c2, 0xa70a03f,  15,         0x10203f,   0x101fbf,   0x101fc0,   0xfe040,    0xfe03f,    0xfdfbf,
      0xfdfc0,    0x1020c0,   0x1020bf,   0x10603f,   0x106040,   0x105fc0,   0x105fbf,   0xfe0bf,    0xfe0c0,    0x302041,   0x1060c0,   0x1060bf,   0x301fc1,   0x2fe041,   0x2fdfc1,   0x70203e,   0x701fbe,
      0x6fe03e,   0x6fdfbe,   0x3020c1,   0x306041,   0x305fc1,   0x2fe0c1,   0x7020be,   0x70603e,   0x705fbe,   0x3101f3f,  0x3101f40,  0x180fa040, 0x180fa03f, 0x6fe0be,   0x180f9fbf, 0x180f9fc0, 0x1b0fdf40,
      0x1b0fdf3f, 0x93060c1,  0x97060be,  0x3105f3f,  0x3105f40,  0x180fa0c0, 0x180fa0bf, 0x192fa041, 0x182f9fc1, 0xb301f41,  0x1b2fdf41, 0xb701f3e,  0x36fdf3e,  0x196fa03e, 0x1b6f9fbe, 0x970213f,  0x9302140,
      0x192fe140, 0x196fe13f, 0x970a03f,  0x930a040,  0xb309fc0,  0xb709fbf,  15,         0x10203f,   0x101fc0,   0xfe040,    0xfe03f,    0x101fbf,   0xfdfc0,    0xfdfbf,    0x1020c0,   0x106040,   0x10603f,
      0x1020bf,   0xfe0c0,    0xfe0bf,    0x105fc0,   0x105fbf,   0x302041,   0x301fc1,   0x2fe041,   0x2fdfc1,   0x1060c0,   0x1060bf,   0x3020c1,   0x306041,   0x305fc1,   0x2fe0c1,   0x70203e,   0x6fe03e,
      0x701fbe,   0x6fdfbe,   0x93060c1,  0x3101f40,  0x3101f3f,  0x180fa040, 0x186fa03f, 0x180f9fc0, 0x1b0fdf40, 0x30fdf3f,  0x1b6f9fbf, 0x7020be,   0x70603e,   0x705fbe,   0x6fe0be,   0x186fa0bf, 0x180fa0c0,
      0x3105f40,  0x3105f3f,  0x3301f41,  0x32fdf41,  0x182fa041, 0x1b2f9fc1, 0x97060be,  0x970a03f,  0x930a040,  0xb309fc0,  0xb305f41,  0x9302140,  0x970213f,  0x192fe140, 0x192fa0c1, 0x196fe13f, 0xb709fbf,
      0xb701f3e,  0x1b6fdf3e, 16,         0x10203f,   0x101fc0,   0xfe040,    0xfe03f,    0x101fbf,   0xfdfc0,    0xfdfbf,    0x102041,   0x1020c0,   0x106040,   0x10603f,   0x1020bf,   0xfe0c0,    0xfe041,
      0x101fc1,   0x105fc0,   0x105fbf,   0xfe0bf,    0xfdfc1,    0x1020c1,   0x106041,   0x1060c0,   0x1060bf,   0xfe0c1,    0x105fc1,   0x93060c1,  0x70203e,   0x3101f40,  0x180fa040, 0x180fa03f, 0x186fe03e,
      0x701fbe,   0x3701f3f,  0x30fdf40,  0x1b0f9fc0, 0x180f9fbf, 0x186fdfbe, 0x1b6fdf3f, 0x3101f41,  0x3105f40,  0xb105f3f,  0x70603e,   0x7020be,   0x6fe0be,   0x180fa0c0, 0x180fa041, 0x182f9fc1, 0x1b2fdf41,
      0xb705fbe,  0x186fa0bf, 0x192fa0c1, 0x97060be,  0xb305f41,  0x9302042,  0x192fe042, 0x13301fc2, 0x930a040,  0x970a03f,  0xb509fc0,  0x9302140,  0x970213f,  0x192fe140, 0x196fe13f, 15,         0x101fc0,
      0xfe040,    0xfdfc0,    0x10203f,   0x101fbf,   0xfe03f,    0x102041,   0x101fc1,   0xfdfbf,    0xfe041,    0x1020c0,   0x106040,   0xfdfc1,    0x105fc0,   0xfe0c0,    0x1020bf,   0x10603f,   0x105fbf,
      0xfe0bf,    0x1020c1,   0x106041,   0x105fc1,   0xfe0c1,    0x1060c0,   0x11060bf,  0x91060c1,  0x3101f40,  0x180fa040, 0x180f9fc0, 0x1b0fdf40, 0x3101f3f,  0x30fdf3f,  0x180fa03f, 0x1b0f9fbf, 0x180fa041,
      0x180f9fc1, 0x3101f41,  0x1b0fdf41, 0x70203e,   0x3701fbe,  0x186fe03e, 0x1b6fdfbe, 0x3105f40,  0xb105f3f,  0x180fa0c0, 0x190fa0bf, 0x192fa0c1, 0x302042,   0x3301fc2,  0xb305f41,  0x182fe042, 0x1b2fdfc2,
      0x17020be,  0x170603e,  0xb705fbe,  0x196fe0be, 0x9502140,  0x194fe140, 0x193020c2, 0xa306042,  0x930a040,  0xb309fc0,  0xa70a03f,  15,         0x10203f,   0xfe03f,    0xfe040,    0x101fc0,   0x101fbf,
      0xfdfbf,    0xfdfc0,    0x1020c0,   0x1020bf,   0xfe0bf,    0xfe0c0,    0x106040,   0x10603f,   0x105fbf,   0x105fc0,   0x1060c0,   0x1060bf,   0x302041,   0x2fe041,   0x301fc1,   0x2fdfc1,   0x70203e,
      0x6fe03e,   0x3020c1,   0x2fe0c1,   0x306041,   0x701fbe,   0x6fdfbe,   0x7020be,   0x6fe0be,   0x70603e,   0x305fc1,   0x83060c1,  0x180fa040, 0x180fa03f, 0x180f9fbf, 0x705fbe,   0x87060be,  0x180f9fc0,
      0x180fa0c0, 0x180fa0bf, 0x3701f3f,  0x3101f40,  0x1b0fdf40, 0x1b6fdf3f, 0x3705f3f,  0x3105f40,  0x1302140,  0x170213f,  0x196fe13f, 0x192fe140, 0x182fa041, 0x1b2f9fc1, 0x192fa0c1, 0x196fa03e, 0x1b6f9fbe,
      0x970a03f,  0x930a040,  0x9306140,  0x970613f,  0xb709fbf,  0xb309fc0,  0x13301f41, 14,         0x10203f,   0xfe040,    0xfe03f,    0x101fc0,   0x101fbf,   0xfdfc0,    0xfdfbf,    0x1020c0,   0x1020bf,
      0xfe0c0,    0xfe0bf,    0x106040,   0x10603f,   0x105fc0,   0x302041,   0x105fbf,   0x2fe041,   0x1060c0,   0x1060bf,   0x301fc1,   0x2fdfc1,   0x3020c1,   0x2fe0c1,   0x306041,   0x305fc1,   0x83060c1,
      0x70203e,   0x6fe03e,   0x701fbe,   0x6fdfbe,   0x7020be,   0x6fe0be,   0x186fa03f, 0x180fa040, 0x180f9fc0, 0x186f9fbf, 0x70603e,   0x8705fbe,  0x3701f3f,  0x3101f40,  0x1b0fdf40, 0x180fa0c0, 0x186fa0bf,
      0x1b6fdf3f, 0x87060be,  0x170213f,  0x1102140,  0x182fa041, 0x182f9fc1, 0x190fe140, 0x196fe13f, 0x3105f40,  0xb705f3f,  0xb301f41,  0x1b2fdf41, 0x182fa0c1, 0x930a040,  0x970a03f,  0x9106140,  0x970613f,
      0x9302141,  0x192fe141, 0xb509fc0,  15,         0x10203f,   0xfe040,    0xfe03f,    0x101fc0,   0x101fbf,   0xfdfc0,    0x1020c0,   0x1020bf,   0xfdfbf,    0xfe0c0,    0x102041,   0x106040,   0xfe0bf,
      0x10603f,   0xfe041,    0x101fc1,   0x105fc0,   0x105fbf,   0xfdfc1,    0x1020c1,   0x1060c0,   0x1060bf,   0xfe0c1,    0x106041,   0x305fc1,   0x83060c1,  0x70203e,   0x180fa040, 0x180fa03f, 0x186fe03e,
      0x701fbe,   0x6fdfbe,   0x180f9fc0, 0x186f9fbf, 0x180fa0c0, 0x180fa0bf, 0x7020be,   0x186fe0be, 0x3101f40,  0x3701f3f,  0x1b0fdf40, 0x1b6fdf3f, 0x70603e,   0x8705fbe,  0x180fa041, 0x182f9fc1, 0x192fa0c1,
      0x1102140,  0x170213f,  0x97060be,  0x190fe140, 0x196fe13f, 0x3301f41,  0x3305f40,  0xb705f3f,  0x1b2fdf41, 0xa302042,  0x1a2fe042, 0x19302141, 0x9506140,  0x930a040,  0x970a03f,  0xb509fc0,  17,
      0xfe040,    0x10203f,   0x101fc0,   0xfdfc0,    0xfe03f,    0x1020c0,   0x102041,   0x101fbf,   0xfe041,    0xfe0c0,    0x106040,   0xfdfbf,    0x1020bf,   0x101fc1,   0xfdfc1,    0xfe0bf,    0x1020c1,
      0x10603f,   0x105fc0,   0xfe0c1,    0x1060c0,   0x106041,   0x8105fbf,  0x81060bf,  0x8105fc1,  0x81060c1,  0x180fa040, 0x180f9fc0, 0x180fa03f, 0x180fa0c0, 0x180fa041, 0x180f9fbf, 0x70203e,   0x186fe03e,
      0x3101f40,  0x1b0fdf40, 0x180f9fc1, 0x180fa0bf, 0x701fbe,   0x3701f3f,  0x1b0fdf3f, 0x1b6fdfbe, 0x7020be,   0x186fe0be, 0x192fa0c1, 0x9102140,  0x190fe140, 0x8302042,  0x3101f41,  0x1b0fdf41, 0x182fe042,
      0xb301fc2,  0xb305f40,  0x870603e,  0x970213f,  0x196fe13f, 0x11102141, 0x113020c2, 0x1b2fdfc2, 0xb105f3f,  0xb705fbe,  0x97060be,  0xa50a040,  11,         0x10203f,   0x101fc0,   0x101fbf,   0xfdfbf,
      0xfe03f,    0xfe040,    0xfdfc0,    0x105fc0,   0x106040,   0x10603f,   0x105fbf,   0x11020bf,  0x11020c0,  0x302041,   0x301fc1,   0x2fdfc1,   0x2fe041,   0x10fe0c0,  0x10fe0bf,  0x70203e,   0x701fbe,
      0x3101f3f,  0x3101f40,  0x305fc1,   0x306041,   0x11060c0,  0x11060bf,  0x6fe03e,   0x6fdfbe,   0x30fdf3f,  0x30fdf40,  0x3105f40,  0x3105f3f,  0x705fbe,   0x70603e,   0x13020c1,  0x12fe0c1,  0x180fa040,
      0x180fa03f, 0x180f9fbf, 0x180f9fc0, 0x3301f41,  0x17020be,  0x16fe0be,  0x93060c1,  0x32fdf41,  0xb305f41,  0x3701f3e,  0x36fdf3e,  0x97060be,  0x970a03f,  0x930a040,  0xb309fc0,  0xb709fbf,  0xb705f3e,
      0x182f9fc1, 0x182fa041, 0x192fa0c0, 0x190fa0bf, 0x196fa03e, 0x186f9fbe, 0x1b6f9f3f, 0x1b2f9f40, 11,         0x101fc0,   0x10203f,   0x101fbf,   0xfdfc0,    0xfe040,    0xfe03f,    0xfdfbf,    0x105fc0,
      0x106040,   0x10603f,   0x105fbf,   0x301fc1,   0x302041,   0x2fe041,   0x2fdfc1,   0x11020c0,  0x11020bf,  0x306041,   0x305fc1,   0x10fe0c0,  0x10fe0bf,  0x3101f40,  0x3101f3f,  0x11060c0,  0x11060bf,
      0x70203e,   0x30fdf40,  0x30fdf3f,  0x701fbe,   0x3105f40,  0x3105f3f,  0x13020c1,  0x6fe03e,   0x6fdfbe,   0x705fbe,   0x70603e,   0x12fe0c1,  0x3301f41,  0x13060c1,  0x32fdf41,  0x1b0f9fc0, 0x180fa040,
      0x186fa03f, 0x1b6f9fbf, 0x3305f41,  0xb109fc0,  0x17020be,  0x16fe0be,  0x810a040,  0x870a03f,  0xb709fbf,  0x3701f3e,  0x1b6fdf3e, 0x17060be,  0x182fa041, 0x1b2f9fc1, 0x192fa0c0, 0x196fa0bf, 0xb705f3e,
      0xb309fc1,  0x830a041,  0x930a0c0,  0x970a0bf,  15,         0x101fc0,   0x101fbf,   0x10203f,   0xfe040,    0xfdfc0,    0xfdfbf,    0xfe03f,    0x106040,   0x105fc0,   0x105fbf,   0x10603f,   0x102041,
      0x101fc1,   0xfdfc1,    0xfe041,    0x106041,   0x105fc1,   0x11020c0,  0x11020bf,  0x10fe0c0,  0x10fe0bf,  0x3101f40,  0x3101f3f,  0x11060c0,  0x11060bf,  0x11020c1,  0x30fdf40,  0x30fdf3f,  0x3105f40,
      0x3105f3f,  0x3101f41,  0x10fe0c1,  0x13060c1,  0x70203e,   0x701fbe,   0x6fdfbe,   0x30fdf41,  0x3305f41,  0x6fe03e,   0x70603e,   0x705fbe,   0x1b0f9fc0, 0x180fa040, 0x186fa03f, 0x1b6f9fbf, 0x1b0f9fc1,
      0x180fa041, 0x910a040,  0xb109fc0,  0xb709fbf,  0x970a03f,  0x17020be,  0x196fe0be, 0x97060be,  0xb701f3e,  0x1b6fdf3e, 0xb301fc2,  0x9302042,  0x930a041,  0xb309fc1,  0x1b2fdfc2, 0x192fe042, 0x194fa0c0,
      17,         0x101fc0,   0x101fbf,   0x10203f,   0xfe040,    0xfdfc0,    0x101fc1,   0x102041,   0x106040,   0x105fc0,   0xfdfbf,    0xfe03f,    0x10603f,   0x105fbf,   0xfdfc1,    0xfe041,    0x106041,
      0x105fc1,   0x11020c0,  0x11020bf,  0x10fe0c0,  0x11020c1,  0x11060c0,  0x3101f40,  0x3101f3f,  0x30fdf40,  0x10fe0bf,  0x11060bf,  0x10fe0c1,  0x3101f41,  0x3105f40,  0x30fdf3f,  0x11060c1,  0x3105f3f,
      0x30fdf41,  0x3105f41,  0x3701fbe,  0x70203e,   0x180fa040, 0x1b0f9fc0, 0x1b0f9fbf, 0x180fa03f, 0x186fe03e, 0x1b6fdfbe, 0x3705fbe,  0x70603e,   0x910a040,  0xb109fc0,  0x3301fc2,  0x1302042,  0x180fa041,
      0x1b0f9fc1, 0x1b2fdfc2, 0x192fe042, 0x1306042,  0x3305fc2,  0xb709fbf,  0x970a03f,  0x930a041,  0xb309fc1,  0x97020be,  0x192fa0c0, 0x190fa0bf, 0x196fe0be, 11,         0x10203f,   0x101fc0,   0x101fbf,
      0xfe03f,    0xfe040,    0xfdfc0,    0xfdfbf,    0x10603f,   0x106040,   0x105fc0,   0x105fbf,   0x11020bf,  0x11020c0,  0x10fe0c0,  0x10fe0bf,  0x302041,   0x301fc1,   0x11060c0,  0x11060bf,  0x2fe041,
      0x2fdfc1,   0x70203e,   0x701fbe,   0x306041,   0x305fc1,   0x3101f40,  0x6fe03e,   0x6fdfbe,   0x3101f3f,  0x70603e,   0x705fbe,   0x13020c1,  0x30fdf40,  0x30fdf3f,  0x3105f3f,  0x3105f40,  0x12fe0c1,
      0x17020be,  0x13060c1,  0x16fe0be,  0x186fa03f, 0x180fa040, 0x1b0f9fc0, 0x1b6f9fbf, 0x17060be,  0x870a03f,  0x3301f41,  0x32fdf41,  0x810a040,  0xb109fc0,  0xb709fbf,  0x3701f3e,  0x1b6fdf3e, 0x3305f41,
      0x190fa0c0, 0x196fa0bf, 0x192fa041, 0x1b2f9fc1, 0xb705f3e,  0x970a0bf,  0x910a0c0,  0x930a041,  0xb309fc1,  11,         0x10203f,   0x101fc0,   0x101fbf,   0xfe040,    0xfe03f,    0xfdfc0,    0xfdfbf,
      0x106040,   0x10603f,   0x105fc0,   0x105fbf,   0x302041,   0x11020c0,  0x11020bf,  0x301fc1,   0x2fe041,   0x10fe0c0,  0x10fe0bf,  0x2fdfc1,   0x306041,   0x11060c0,  0x11060bf,  0x305fc1,   0x13020c1,
      0x12fe0c1,  0x70203e,   0x701fbe,   0x3101f3f,  0x3101f40,  0x30fdf40,  0x30fdf3f,  0x6fdfbe,   0x6fe03e,   0x70603e,   0x13060c1,  0x3105f40,  0x3105f3f,  0x705fbe,   0x17020be,  0x180fa040, 0x186fa03f,
      0x1b0f9fc0, 0x1b6f9fbf, 0x3301f41,  0x32fdf41,  0x3305f41,  0x16fe0be,  0x17060be,  0x870a03f,  0x810a040,  0xb109fc0,  0xb709fbf,  0x3701f3e,  0x182fa041, 0x192fa0c0, 0x196fa0bf, 0x1b2f9fc1, 0x1b6fdf3e,
      0xb705f3e,  0x830a041,  0x930a0c0,  0x970a0bf,  0xb309fc1,  14,         0x101fc0,   0x10203f,   0x101fbf,   0xfe040,    0xfdfc0,    0xfe03f,    0xfdfbf,    0x106040,   0x105fc0,   0x10603f,   0x105fbf,
      0x102041,   0x101fc1,   0xfe041,    0x11020c0,  0xfdfc1,    0x11020bf,  0x106041,   0x105fc1,   0x10fe0c0,  0x10fe0bf,  0x11060c0,  0x11060bf,  0x11020c1,  0x10fe0c1,  0x13060c1,  0x3101f40,  0x3101f3f,
      0x30fdf40,  0x30fdf3f,  0x3105f40,  0x3105f3f,  0x3701fbe,  0x70203e,   0x6fe03e,   0x36fdfbe,  0x3101f41,  0x32fdf41,  0x1b0f9fc0, 0x180fa040, 0x186fa03f, 0x70603e,   0x3705fbe,  0x1b6f9fbf, 0x3305f41,
      0xb109fc0,  0x810a040,  0x17020be,  0x16fe0be,  0x870a03f,  0xb709fbf,  0x180fa041, 0x1b2f9fc1, 0x192fa0c0, 0x196fa0bf, 0x17060be,  0x9302042,  0xb301fc2,  0x830a041,  0xb309fc1,  0x930a0c0,  0x970a0bf,
      0x1a2fe042, 17,         0x101fc0,   0x10203f,   0xfe040,    0xfdfc0,    0x101fbf,   0x106040,   0x102041,   0x101fc1,   0x105fc0,   0xfe03f,    0xfdfbf,    0x10603f,   0x105fbf,   0xfe041,    0xfdfc1,
      0x106041,   0x105fc1,   0x11020c0,  0x11020bf,  0x10fe0c0,  0x11020c1,  0x11060c0,  0x10fe0bf,  0x11060bf,  0x10fe0c1,  0x11060c1,  0x3101f40,  0x30fdf40,  0x3101f3f,  0x3105f40,  0x3101f41,  0x30fdf3f,
      0x3105f3f,  0x30fdf41,  0x3105f41,  0x70203e,   0x3701fbe,  0x180fa040, 0x1b0f9fc0, 0x180fa03f, 0x186fe03e, 0x36fdfbe,  0x1b6f9fbf, 0x180fa041, 0x1b0f9fc1, 0x1302042,  0x3301fc2,  0x910a040,  0x70603e,
      0x3705fbe,  0xb109fc0,  0x970a03f,  0xb709fbf,  0x910a041,  0x192fe042, 0x1b2fdfc2, 0x9306042,  0x3305fc2,  0xb309fc1,  0x97020be,  0x192fa0c0, 0x190fa0bf, 0x196fe0be, 15,         0x10203f,   0x101fbf,
      0x101fc0,   0xfe040,    0xfe03f,    0xfdfbf,    0xfdfc0,    0x106040,   0x10603f,   0x105fbf,   0x105fc0,   0x1020c0,   0x1020bf,   0xfe0bf,    0xfe0c0,    0x1060c0,   0x1060bf,   0x302041,   0x301fc1,
      0x2fe041,   0x2fdfc1,   0x70203e,   0x701fbe,   0x306041,   0x305fc1,   0x3020c1,   0x6fe03e,   0x6fdfbe,   0x70603e,   0x705fbe,   0x7020be,   0x2fe0c1,   0x13060c1,  0x3101f40,  0x3101f3f,  0x30fdf3f,
      0x6fe0be,   0x17060be,  0x30fdf40,  0x3105f40,  0x3105f3f,  0x186fa03f, 0x180fa040, 0x1b0f9fc0, 0x1b6f9fbf, 0x186fa0bf, 0x180fa0c0, 0x830a040,  0x870a03f,  0xb709fbf,  0xb309fc0,  0x3301f41,  0x1b2fdf41,
      0xb305f41,  0xb701f3e,  0x1b6fdf3e, 0x970213f,  0x9302140,  0x930a0c0,  0x970a0bf,  0x196fe13f, 0x192fe140, 0x1a2fa041, 14,         0x10203f,   0x101fc0,   0x101fbf,   0xfe040,    0xfe03f,    0xfdfc0,
      0xfdfbf,    0x106040,   0x10603f,   0x105fc0,   0x105fbf,   0x1020c0,   0x1020bf,   0xfe0c0,    0x302041,   0xfe0bf,    0x301fc1,   0x1060c0,   0x1060bf,   0x2fe041,   0x2fdfc1,   0x306041,   0x305fc1,
      0x3020c1,   0x2fe0c1,   0x13060c1,  0x70203e,   0x701fbe,   0x6fe03e,   0x6fdfbe,   0x70603e,   0x705fbe,   0x3701f3f,  0x3101f40,  0x30fdf40,  0x36fdf3f,  0x7020be,   0x16fe0be,  0x186fa03f, 0x180fa040,
      0x1b0f9fc0, 0x3105f40,  0x3705f3f,  0x1b6f9fbf, 0x17060be,  0x870a03f,  0x810a040,  0x3301f41,  0x32fdf41,  0xb109fc0,  0xb709fbf,  0x180fa0c0, 0x196fa0bf, 0x192fa041, 0x1b2f9fc1, 0x3305f41,  0x9302140,
      0x970213f,  0x910a0c0,  0x970a0bf,  0x930a041,  0xb309fc1,  0x194fe140, 15,         0x10203f,   0x101fc0,   0x101fbf,   0xfe040,    0xfe03f,    0xfdfc0,    0x106040,   0x10603f,   0xfdfbf,    0x105fc0,
      0x102041,   0x1020c0,   0x105fbf,   0x1020bf,   0x101fc1,   0xfe041,    0xfe0c0,    0xfe0bf,    0xfdfc1,    0x106041,   0x1060c0,   0x1060bf,   0x105fc1,   0x1020c1,   0x2fe0c1,   0x13060c1,  0x70203e,
      0x3101f40,  0x3101f3f,  0x3701fbe,  0x6fe03e,   0x6fdfbe,   0x30fdf40,  0x36fdf3f,  0x3105f40,  0x3105f3f,  0x70603e,   0x3705fbe,  0x180fa040, 0x186fa03f, 0x1b0f9fc0, 0x1b6f9fbf, 0x7020be,   0x16fe0be,
      0x3101f41,  0x32fdf41,  0xb305f41,  0x810a040,  0x870a03f,  0x97060be,  0xb109fc0,  0xb709fbf,  0x182fa041, 0x182fa0c0, 0x196fa0bf, 0x1b2f9fc1, 0x11302042, 0x13301fc2, 0xb30a041,  0x950a0c0,  0x9302140,
      0x970213f,  0x194fe140, 17,         0x101fc0,   0x10203f,   0xfe040,    0xfdfc0,    0x101fbf,   0x106040,   0x102041,   0xfe03f,    0x101fc1,   0x105fc0,   0x1020c0,   0xfdfbf,    0x10603f,   0xfe041,
      0xfdfc1,    0x105fbf,   0x106041,   0x1020bf,   0xfe0c0,    0x105fc1,   0x1060c0,   0x1020c1,   0x10fe0bf,  0x11060bf,  0x10fe0c1,  0x11060c1,  0x3101f40,  0x30fdf40,  0x3101f3f,  0x3105f40,  0x3101f41,
      0x30fdf3f,  0x70203e,   0x3701fbe,  0x180fa040, 0x1b0f9fc0, 0x30fdf41,  0x3105f3f,  0x6fe03e,   0x186fa03f, 0x1b0f9fbf, 0x1b6fdfbe, 0x70603e,   0x3705fbe,  0xb305f41,  0x910a040,  0xb109fc0,  0x1302042,
      0x180fa041, 0x1b0f9fc1, 0x3301fc2,  0x192fe042, 0x192fa0c0, 0x17020be,  0x970a03f,  0xb709fbf,  0xa10a041,  0xa306042,  0x1b2fdfc2, 0x190fa0bf, 0x196fe0be, 0x97060be,  0x11502140, 17,         0x10203f,
      0x101fbf,   0x101fc0,   0xfe040,    0xfe03f,    0x1020bf,   0x1020c0,   0x106040,   0x10603f,   0xfdfbf,    0xfdfc0,    0x105fc0,   0x105fbf,   0xfe0bf,    0xfe0c0,    0x1060c0,   0x1060bf,   0x302041,
      0x301fc1,   0x2fe041,   0x3020c1,   0x306041,   0x70203e,   0x701fbe,   0x6fe03e,   0x2fdfc1,   0x305fc1,   0x2fe0c1,   0x7020be,   0x70603e,   0x6fdfbe,   0x3060c1,   0x705fbe,   0x6fe0be,   0x7060be,
      0x3701f3f,  0x3101f40,  0x180fa040, 0x186fa03f, 0x186f9fbf, 0x180f9fc0, 0x1b0fdf40, 0x1b6fdf3f, 0x3705f3f,  0x3105f40,  0x830a040,  0x870a03f,  0x170213f,  0x1302140,  0x180fa0c0, 0x186fa0bf, 0x196fe13f,
      0x192fe140, 0x1306140,  0x170613f,  0xb709fbf,  0xb309fc0,  0x930a0c0,  0x970a0bf,  0xb301f41,  0x192fa041, 0x182f9fc1, 0x1b2fdf41, 17,         0x10203f,   0x101fc0,   0xfe040,    0xfe03f,    0x101fbf,
      0x106040,   0x1020c0,   0x1020bf,   0x10603f,   0xfdfc0,    0xfdfbf,    0x105fc0,   0x105fbf,   0xfe0c0,    0xfe0bf,    0x1060c0,   0x1060bf,   0x302041,   0x301fc1,   0x2fe041,   0x3020c1,   0x306041,
      0x2fdfc1,   0x305fc1,   0x2fe0c1,   0x3060c1,   0x70203e,   0x6fe03e,   0x701fbe,   0x70603e,   0x7020be,   0x6fdfbe,   0x705fbe,   0x6fe0be,   0x7060be,   0x3101f40,  0x3701f3f,  0x180fa040, 0x186fa03f,
      0x180f9fc0, 0x1b0fdf40, 0x36fdf3f,  0x1b6f9fbf, 0x180fa0c0, 0x186fa0bf, 0x1302140,  0x170213f,  0x830a040,  0x3105f40,  0x3705f3f,  0x870a03f,  0xb309fc0,  0xb709fbf,  0x830a0c0,  0x192fe140, 0x196fe13f,
      0x9306140,  0x170613f,  0x970a0bf,  0xb301f41,  0x192fa041, 0x182f9fc1, 0x1b2fdf41, 17,         0x10203f,   0x101fc0,   0xfe040,    0xfe03f,    0x101fbf,   0x106040,   0x1020c0,   0xfdfc0,    0x1020bf,
      0x10603f,   0x102041,   0xfdfbf,    0x105fc0,   0xfe0c0,    0xfe0bf,    0x105fbf,   0x1060c0,   0x101fc1,   0xfe041,    0x1060bf,   0x106041,   0x1020c1,   0x2fdfc1,   0x305fc1,   0x2fe0c1,   0x3060c1,
      0x70203e,   0x6fe03e,   0x701fbe,   0x70603e,   0x7020be,   0x6fdfbe,   0x3101f40,  0x3701f3f,  0x180fa040, 0x186fa03f, 0x6fe0be,   0x705fbe,   0x30fdf40,  0x1b0f9fc0, 0x186f9fbf, 0x1b6fdf3f, 0x3105f40,
      0x3705f3f,  0x97060be,  0x830a040,  0x870a03f,  0x1302140,  0x180fa0c0, 0x186fa0bf, 0x170213f,  0x192fe140, 0x192fa041, 0x3301f41,  0xb309fc0,  0xb709fbf,  0x850a0c0,  0x9506140,  0x196fe13f, 0x182f9fc1,
      0x1b2fdf41, 0xb305f41,  0x12302042, 16,         0x10203f,   0x101fc0,   0xfe040,    0x102041,   0x1020c0,   0x106040,   0x101fbf,   0xfe03f,    0xfdfc0,    0x101fc1,   0x105fc0,   0x10603f,   0x1020bf,
      0xfe0c0,    0xfe041,    0xfdfbf,    0x1020c1,   0x106041,   0x1060c0,   0x105fbf,   0xfe0bf,    0xfdfc1,    0x105fc1,   0x1060bf,   0xfe0c1,    0x83060c1,  0x70203e,   0x3101f40,  0x180fa040, 0x180fa03f,
      0x186fe03e, 0x701fbe,   0x3701f3f,  0x30fdf40,  0x1b0f9fc0, 0x180fa041, 0x180fa0c0, 0x7020be,   0x70603e,   0x3105f40,  0xb101f41,  0x9302042,  0x1102140,  0x930a040,  0x6fdfbe,   0x36fdf3f,  0x1b6f9fbf,
      0x190fa0bf, 0x196fe0be, 0x170213f,  0x196fe140, 0x182f9fc1, 0x1b2fdf41, 0xb301fc2,  0x1a2fe042, 0x192fa0c1, 0x8705fbe,  0xb705f3f,  0xb309fc0,  0xa70a03f,  0x97060be,  0x9706140,  0x11302141};

} // namespace voro
// Voro++, a 3D cell-based Voronoi library
//
// Author   : Chris H. Rycroft (LBL / UC Berkeley)
// Email    : chr@alum.mit.edu
// Date     : August 30th 2011

/** \file c_loops.cc
 * \brief Function implementations for the loop classes. */

namespace voro {

/** Initializes a c_loop_subset object to scan over all particles within a
   * given sphere.
   * \param[in] (vx,vy,vz) the position vector of the center of the sphere.
   * \param[in] r the radius of the sphere.
   * \param[in] bounds_test whether to do detailed bounds checking. If this is
   *                        false then the class will loop over all particles in
   *                        blocks that overlap the given sphere. If it is true,
   *                        the particle will only loop over the particles which
   *                        actually lie within the sphere.
   * \return True if there is any valid point to loop over, false otherwise. */
  void c_loop_subset::setup_sphere(double vx, double vy, double vz, double r, bool bounds_test) {
    if (bounds_test) {
      mode = sphere;
      v0 = vx;
      v1 = vy;
      v2 = vz;
      v3 = r * r;
    } else {
      mode = no_check;
    }
    ai = step_int((vx - ax - r) * xsp);
    bi = step_int((vx - ax + r) * xsp);
    aj = step_int((vy - ay - r) * ysp);
    bj = step_int((vy - ay + r) * ysp);
    ak = step_int((vz - az - r) * zsp);
    bk = step_int((vz - az + r) * zsp);
    setup_common();
  }

/** Initializes the class to loop over all particles in a rectangular subgrid
   * of blocks.
   * \param[in] (ai_,bi_) the subgrid range in the x-direction, inclusive of both
   *                      ends.
   * \param[in] (aj_,bj_) the subgrid range in the y-direction, inclusive of both
   *                      ends.
   * \param[in] (ak_,bk_) the subgrid range in the z-direction, inclusive of both
   *                      ends.
   * \return True if there is any valid point to loop over, false otherwise. */
  void c_loop_subset::setup_intbox(int ai_, int bi_, int aj_, int bj_, int ak_, int bk_) {
    ai = ai_;
    bi = bi_;
    aj = aj_;
    bj = bj_;
    ak = ak_;
    bk = bk_;
    mode = no_check;
    setup_common();
  }

/** Sets up all of the common constants used for the loop.
   * \return True if there is any valid point to loop over, false otherwise. */
  void c_loop_subset::setup_common() {
    if (!xperiodic) {
      if (ai < 0) {
        ai = 0;
        if (bi < 0) {
          bi = 0;
        }
      }
      if (bi >= nx) {
        bi = nx - 1;
        if (ai >= nx) {
          ai = nx - 1;
        }
      }
    }
    if (!yperiodic) {
      if (aj < 0) {
        aj = 0;
        if (bj < 0) {
          bj = 0;
        }
      }
      if (bj >= ny) {
        bj = ny - 1;
        if (aj >= ny) {
          aj = ny - 1;
        }
      }
    }
    if (!zperiodic) {
      if (ak < 0) {
        ak = 0;
        if (bk < 0) {
          bk = 0;
        }
      }
      if (bk >= nz) {
        bk = nz - 1;
        if (ak >= nz) {
          ak = nz - 1;
        }
      }
    }
    ci = ai;
    cj = aj;
    ck = ak;
    di = i = step_mod(ci, nx);
    apx = px = step_div(ci, nx) * sx;
    dj = j = step_mod(cj, ny);
    apy = py = step_div(cj, ny) * sy;
    dk = k = step_mod(ck, nz);
    apz = pz = step_div(ck, nz) * sz;
    inc1 = di - step_mod(bi, nx);
    inc2 = nx * (ny + dj - step_mod(bj, ny)) + inc1;
    inc1 += nx;
    ijk = di + nx * (dj + ny * dk);
    q = 0;
  }

/** Starts the loop by finding the first particle within the container to
   * consider.
   * \return True if there is any particle to consider, false otherwise. */
  bool c_loop_subset::start() {
    while (co[ijk] == 0) {
      if (!next_block()) {
        return false;
      }
    }
    while (mode != no_check && out_of_bounds()) {
      q++;
      while (q >= co[ijk]) {
        q = 0;
        if (!next_block()) {
          return false;
        }
      }
    }
    return true;
  }

/** Initializes the class to loop over all particles in a rectangular box.
   * \param[in] (xmin,xmax) the minimum and maximum x coordinates of the box.
   * \param[in] (ymin,ymax) the minimum and maximum y coordinates of the box.
   * \param[in] (zmin,zmax) the minimum and maximum z coordinates of the box.
   * \param[in] bounds_test whether to do detailed bounds checking. If this is
   *                        false then the class will loop over all particles in
   *                        blocks that overlap the given box. If it is true, the
   *                        particle will only loop over the particles which
   *                        actually lie within the box.
   * \return True if there is any valid point to loop over, false otherwise. */
  void c_loop_subset::setup_box(double xmin, double xmax, double ymin, double ymax, double zmin, double zmax, bool bounds_test) {
    if (bounds_test) {
      mode = box;
      v0 = xmin;
      v1 = xmax;
      v2 = ymin;
      v3 = ymax;
      v4 = zmin;
      v5 = zmax;
    } else {
      mode = no_check;
    }
    ai = step_int((xmin - ax) * xsp);
    bi = step_int((xmax - ax) * xsp);
    aj = step_int((ymin - ay) * ysp);
    bj = step_int((ymax - ay) * ysp);
    ak = step_int((zmin - az) * zsp);
    bk = step_int((zmax - az) * zsp);
    setup_common();
  }

/** Computes whether the current point is out of bounds, relative to the
   * current loop setup.
   * \return True if the point is out of bounds, false otherwise. */
  bool c_loop_subset::out_of_bounds() {
    double *pp = p[ijk] + ps * q;
    if (mode == sphere) {
      const double fx(*pp + px - v0);
      const double fy(pp[1] + py - v1);
      const double fz(pp[2] + pz - v2);
      return fx * fx + fy * fy + fz * fz > v3;
    } else {
      double f(*pp + px);
      if (f < v0 || f > v1) {
        return true;
      }
      f = pp[1] + py;
      if (f < v2 || f > v3) {
        return true;
      }
      f = pp[2] + pz;
      return f < v4 || f > v5;
    }
  }

/** Returns the next block to be tested in a loop, and updates the periodicity
   * vector if necessary. */
  bool c_loop_subset::next_block() {
    if (i < bi) {
      i++;
      if (ci < nx - 1) {
        ci++;
        ijk++;
      } else {
        ci = 0;
        ijk += 1 - nx;
        px += sx;
      }
      return true;
    } else if (j < bj) {
      i = ai;
      ci = di;
      px = apx;
      j++;
      if (cj < ny - 1) {
        cj++;
        ijk += inc1;
      } else {
        cj = 0;
        ijk += inc1 - nxy;
        py += sy;
      }
      return true;
    } else if (k < bk) {
      i = ai;
      ci = di;
      j = aj;
      cj = dj;
      px = apx;
      py = apy;
      k++;
      if (ck < nz - 1) {
        ck++;
        ijk += inc2;
      } else {
        ck = 0;
        ijk += inc2 - nxyz;
        pz += sz;
      }
      return true;
    } else {
      return false;
    }
  }

/** Extends the memory available for storing the ordering. */
  void particle_order::add_ordering_memory() {
    int *no = new int[size << 2];
    int *nop = no;
    int *opp = o;
    while (opp < op) {
      *(nop++) = *(opp++);
    }
    delete[] o;
    size <<= 1;
    o = no;
    op = nop;
  }

} // namespace voro
// Voro++, a 3D cell-based Voronoi library
//
// Author   : Chris H. Rycroft (LBL / UC Berkeley)
// Email    : chr@alum.mit.edu
// Date     : August 30th 2011

/** \file v_compute.cc
 * \brief Function implementantions for the voro_compute template. */

namespace voro {

/** The class constructor initializes constants from the container class, and
   * sets up the mask and queue used for Voronoi computations.
   * \param[in] con_ a reference to the container class to use.
   * \param[in] (hx_,hy_,hz_) the size of the mask to use. */
  template <class c_class>
  voro_compute<c_class>::voro_compute(c_class &con_, int hx_, int hy_, int hz_) :
      con(con_),
      boxx(con_.boxx),
      boxy(con_.boxy),
      boxz(con_.boxz),
      xsp(con_.xsp),
      ysp(con_.ysp),
      zsp(con_.zsp),
      hx(hx_),
      hy(hy_),
      hz(hz_),
      hxy(hx_ * hy_),
      hxyz(hxy * hz_),
      ps(con_.ps),
      id(con_.id),
      p(con_.p),
      co(con_.co),
      bxsq(boxx * boxx + boxy * boxy + boxz * boxz),
      mv(0),
      qu_size(3 * (3 + hxy + hz * (hx + hy))),
      wl(con_.wl),
      mrad(con_.mrad),
      mask(new unsigned int[hxyz]),
      qu(new int[qu_size]),
      qu_l(qu + qu_size) {
    reset_mask();
  }

/** Scans all of the particles within a block to see if any of them have a
   * smaller distance to the given test vector. If one is found, the routine
   * updates the minimum distance and store information about this particle.
   * \param[in] ijk the index of the block.
   * \param[in] (x,y,z) the test vector to consider (which may have already had a
   *                    periodic displacement applied to it).
   * \param[in] (di,dj,dk) the coordinates of the current block, to store if the
   *			 particle record is updated.
   * \param[in,out] w a reference to a particle record in which to store
   *		    information about the particle whose Voronoi cell the
   *		    vector is within.
   * \param[in,out] mrs the current minimum distance, that may be updated if a
   * 		      closer particle is found. */
  template <class c_class> inline void voro_compute<c_class>::scan_all(int ijk, double x, double y, double z, int di, int dj, int dk, particle_record &w, double &mrs) {
    double x1;
    double y1;
    double z1;
    double rs;
    bool in_block = false;
    for (int l = 0; l < co[ijk]; l++) {
      x1 = p[ijk][ps * l] - x;
      y1 = p[ijk][ps * l + 1] - y;
      z1 = p[ijk][ps * l + 2] - z;
      rs = con.r_current_sub(x1 * x1 + y1 * y1 + z1 * z1, ijk, l);
      if (rs < mrs) {
        mrs = rs;
        w.l = l;
        in_block = true;
      }
    }
    if (in_block) {
      w.ijk = ijk;
      w.di = di;
      w.dj = dj, w.dk = dk;
    }
  }

/** Finds the Voronoi cell that given vector is within. For containers that are
   * not radially dependent, this corresponds to findig the particle that is
   * closest to the vector; for the radical tessellation containers, this
   * corresponds to a finding the minimum weighted distance.
   * \param[in] (x,y,z) the vector to consider.
   * \param[in] (ci,cj,ck) the coordinates of the block that the test particle is
   *                       in relative to the container data structure.
   * \param[in] ijk the index of the block that the test particle is in.
   * \param[out] w a reference to a particle record in which to store information
   * 		 about the particle whose Voronoi cell the vector is within.
   * \param[out] mrs the minimum computed distance. */
  template <class c_class> void voro_compute<c_class>::find_voronoi_cell(double x, double y, double z, int ci, int cj, int ck, int ijk, particle_record &w, double &mrs) {
    double qx = 0;
    double qy = 0;
    double qz = 0;
    double rs;
    int i;
    int j;
    int k;
    int di;
    int dj;
    int dk;
    int ei;
    int ej;
    int ek;
    int f;
    int g;
    int disp;
    double fx;
    double fy;
    double fz;
    double mxs;
    double mys;
    double mzs;
    double *radp;
    unsigned int q;
    unsigned int *e;
    unsigned int *mijk;

    // Init setup for parameters to return
    w.ijk = -1;
    mrs = large_number;

    con.initialize_search(ci, cj, ck, ijk, i, j, k, disp);

    // Test all particles in the particle's local region first
    scan_all(ijk, x, y, z, 0, 0, 0, w, mrs);

    // Now compute the fractional position of the particle within its
    // region and store it in (fx,fy,fz). We use this to compute an index
    // (di,dj,dk) of which subregion the particle is within.
    unsigned int m1;
    unsigned int m2;
    con.frac_pos(x, y, z, ci, cj, ck, fx, fy, fz);
    di = int(fx * xsp * wl_fgrid);
    dj = int(fy * ysp * wl_fgrid);
    dk = int(fz * zsp * wl_fgrid);

    // The indices (di,dj,dk) tell us which worklist to use, to test the
    // blocks in the optimal order. But we only store worklists for the
    // eighth of the region where di, dj, and dk are all less than half the
    // full grid. The rest of the cases are handled by symmetry. In this
    // section, we detect for these cases, by reflecting high values of di,
    // dj, and dk. For these cases, a mask is constructed in m1 and m2
    // which is used to flip the worklist information when it is loaded.
    if (di >= wl_hgrid) {
      mxs = boxx - fx;
      m1 = 127 + (3 << 21);
      m2 = 1 + (1 << 21);
      di = wl_fgrid - 1 - di;
      if (di < 0) {
        di = 0;
      }
    } else {
      m1 = m2 = 0;
      mxs = fx;
    }
    if (dj >= wl_hgrid) {
      mys = boxy - fy;
      m1 |= (127 << 7) + (3 << 24);
      m2 |= (1 << 7) + (1 << 24);
      dj = wl_fgrid - 1 - dj;
      if (dj < 0) {
        dj = 0;
      }
    } else {
      mys = fy;
    }
    if (dk >= wl_hgrid) {
      mzs = boxz - fz;
      m1 |= (127 << 14) + (3 << 27);
      m2 |= (1 << 14) + (1 << 27);
      dk = wl_fgrid - 1 - dk;
      if (dk < 0) {
        dk = 0;
      }
    } else {
      mzs = fz;
    }

    // Do a quick test to account for the case when the minimum radius is
    // small enought that no other blocks need to be considered
    rs = con.r_max_add(mrs);
    if (mxs * mxs > rs && mys * mys > rs && mzs * mzs > rs) {
      return;
    }

    // Now compute which worklist we are going to use, and set radp and e to
    // point at the right offsets
    ijk = di + wl_hgrid * (dj + wl_hgrid * dk);
    radp = mrad + ijk * wl_seq_length;
    e = (const_cast<unsigned int *>(wl)) + ijk * wl_seq_length;

    // Read in how many items in the worklist can be tested without having to
    // worry about writing to the mask
    f = e[0];
    g = 0;
    do {
      // If mrs is less than the minimum distance to any untested
      // block, then we are done
      if (con.r_max_add(mrs) < radp[g]) {
        return;
      }
      g++;

      // Load in a block off the worklist, permute it with the
      // symmetry mask, and decode its position. These are all
      // integer bit operations so they should run very fast.
      q = e[g];
      q ^= m1;
      q += m2;
      di = q & 127;
      di -= 64;
      dj = (q >> 7) & 127;
      dj -= 64;
      dk = (q >> 14) & 127;
      dk -= 64;

      // Check that the worklist position is in range
      ei = di + i;
      if (ei < 0 || ei >= hx) {
        continue;
      }
      ej = dj + j;
      if (ej < 0 || ej >= hy) {
        continue;
      }
      ek = dk + k;
      if (ek < 0 || ek >= hz) {
        continue;
      }

      // Call the compute_min_max_radius() function. This returns
      // true if the minimum distance to the block is bigger than the
      // current mrs, in which case we skip this block and move on.
      // Otherwise, it computes the maximum distance to the block and
      // returns it in crs.
      if (compute_min_radius(di, dj, dk, fx, fy, fz, mrs)) {
        continue;
      }

      // Now compute which region we are going to loop over, adding a
      // displacement for the periodic cases
      ijk = con.region_index(ci, cj, ck, ei, ej, ek, qx, qy, qz, disp);

      // If mrs is bigger than the maximum distance to the block,
      // then we have to test all particles in the block for
      // intersections. Otherwise, we do additional checks and skip
      // those particles which can't possibly intersect the block.
      scan_all(ijk, x - qx, y - qy, z - qz, di, dj, dk, w, mrs);
    } while (g < f);

    // Update mask value and initialize queue
    mv++;
    if (mv == 0) {
      reset_mask();
      mv = 1;
    }
    int *qu_s = qu;
    int *qu_e = qu;

    while (g < wl_seq_length - 1) {
      // If mrs is less than the minimum distance to any untested
      // block, then we are done
      if (con.r_max_add(mrs) < radp[g]) {
        return;
      }
      g++;

      // Load in a block off the worklist, permute it with the
      // symmetry mask, and decode its position. These are all
      // integer bit operations so they should run very fast.
      q = e[g];
      q ^= m1;
      q += m2;
      di = q & 127;
      di -= 64;
      dj = (q >> 7) & 127;
      dj -= 64;
      dk = (q >> 14) & 127;
      dk -= 64;

      // Compute the position in the mask of the current block. If
      // this lies outside the mask, then skip it. Otherwise, mark
      // it.
      ei = di + i;
      if (ei < 0 || ei >= hx) {
        continue;
      }
      ej = dj + j;
      if (ej < 0 || ej >= hy) {
        continue;
      }
      ek = dk + k;
      if (ek < 0 || ek >= hz) {
        continue;
      }
      mijk = mask + ei + hx * (ej + hy * ek);
      *mijk = mv;

      // Skip this block if it is further away than the current
      // minimum radius
      if (compute_min_radius(di, dj, dk, fx, fy, fz, mrs)) {
        continue;
      }

      // Now compute which region we are going to loop over, adding a
      // displacement for the periodic cases
      ijk = con.region_index(ci, cj, ck, ei, ej, ek, qx, qy, qz, disp);
      scan_all(ijk, x - qx, y - qy, z - qz, di, dj, dk, w, mrs);

      if (qu_e > qu_l - 18) {
        add_list_memory(qu_s, qu_e);
      }
      scan_bits_mask_add(q, mijk, ei, ej, ek, qu_e);
    }

    // Do a check to see if we've reached the radius cutoff
    if (con.r_max_add(mrs) < radp[g]) {
      return;
    }

    // We were unable to completely compute the cell based on the blocks in
    // the worklist, so now we have to go block by block, reading in items
    // off the list
    while (qu_s != qu_e) {
      // Read the next entry of the queue
      if (qu_s == qu_l) {
        qu_s = qu;
      }
      ei = *(qu_s++);
      ej = *(qu_s++);
      ek = *(qu_s++);
      di = ei - i;
      dj = ej - j;
      dk = ek - k;
      if (compute_min_radius(di, dj, dk, fx, fy, fz, mrs)) {
        continue;
      }

      ijk = con.region_index(ci, cj, ck, ei, ej, ek, qx, qy, qz, disp);
      scan_all(ijk, x - qx, y - qy, z - qz, di, dj, dk, w, mrs);

      // Test the neighbors of the current block, and add them to the
      // block list if they haven't already been tested
      if ((qu_s <= qu_e ? (qu_l - qu_e) + (qu_s - qu) : qu_s - qu_e) < 18) {
        add_list_memory(qu_s, qu_e);
      }
      add_to_mask(ei, ej, ek, qu_e);
    }
  }

/** Scans the six orthogonal neighbors of a given block and adds them to the
   * queue if they haven't been considered already. It assumes that the queue
   * will definitely have enough memory to add six entries at the end.
   * \param[in] (ei,ej,ek) the block to consider.
   * \param[in,out] qu_e a pointer to the end of the queue. */
  template <class c_class> inline void voro_compute<c_class>::add_to_mask(int ei, int ej, int ek, int *&qu_e) {
    unsigned int *mijk = mask + ei + hx * (ej + hy * ek);
    if (ek > 0) {
      if (*(mijk - hxy) != mv) {
        if (qu_e == qu_l) {
          qu_e = qu;
        }
        *(mijk - hxy) = mv;
        *(qu_e++) = ei;
        *(qu_e++) = ej;
        *(qu_e++) = ek - 1;
      }
    }
    if (ej > 0) {
      if (*(mijk - hx) != mv) {
        if (qu_e == qu_l) {
          qu_e = qu;
        }
        *(mijk - hx) = mv;
        *(qu_e++) = ei;
        *(qu_e++) = ej - 1;
        *(qu_e++) = ek;
      }
    }
    if (ei > 0) {
      if (*(mijk - 1) != mv) {
        if (qu_e == qu_l) {
          qu_e = qu;
        }
        *(mijk - 1) = mv;
        *(qu_e++) = ei - 1;
        *(qu_e++) = ej;
        *(qu_e++) = ek;
      }
    }
    if (ei < hx - 1) {
      if (*(mijk + 1) != mv) {
        if (qu_e == qu_l) {
          qu_e = qu;
        }
        *(mijk + 1) = mv;
        *(qu_e++) = ei + 1;
        *(qu_e++) = ej;
        *(qu_e++) = ek;
      }
    }
    if (ej < hy - 1) {
      if (*(mijk + hx) != mv) {
        if (qu_e == qu_l) {
          qu_e = qu;
        }
        *(mijk + hx) = mv;
        *(qu_e++) = ei;
        *(qu_e++) = ej + 1;
        *(qu_e++) = ek;
      }
    }
    if (ek < hz - 1) {
      if (*(mijk + hxy) != mv) {
        if (qu_e == qu_l) {
          qu_e = qu;
        }
        *(mijk + hxy) = mv;
        *(qu_e++) = ei;
        *(qu_e++) = ej;
        *(qu_e++) = ek + 1;
      }
    }
  }

/** Scans a worklist entry and adds any blocks to the queue
   * \param[in] (ei,ej,ek) the block to consider.
   * \param[in,out] qu_e a pointer to the end of the queue. */
  template <class c_class> inline void voro_compute<c_class>::scan_bits_mask_add(unsigned int q, unsigned int *mijk, int ei, int ej, int ek, int *&qu_e) {
    const unsigned int b1 = 1 << 21;
    const unsigned int b2 = 1 << 22;
    const unsigned int b3 = 1 << 24;
    const unsigned int b4 = 1 << 25;
    const unsigned int b5 = 1 << 27;
    const unsigned int b6 = 1 << 28;
    if ((q & b2) == b2) {
      if (ei > 0) {
        *(mijk - 1) = mv;
        *(qu_e++) = ei - 1;
        *(qu_e++) = ej;
        *(qu_e++) = ek;
      }
      if ((q & b1) == 0 && ei < hx - 1) {
        *(mijk + 1) = mv;
        *(qu_e++) = ei + 1;
        *(qu_e++) = ej;
        *(qu_e++) = ek;
      }
    } else if ((q & b1) == b1 && ei < hx - 1) {
      *(mijk + 1) = mv;
      *(qu_e++) = ei + 1;
      *(qu_e++) = ej;
      *(qu_e++) = ek;
    }
    if ((q & b4) == b4) {
      if (ej > 0) {
        *(mijk - hx) = mv;
        *(qu_e++) = ei;
        *(qu_e++) = ej - 1;
        *(qu_e++) = ek;
      }
      if ((q & b3) == 0 && ej < hy - 1) {
        *(mijk + hx) = mv;
        *(qu_e++) = ei;
        *(qu_e++) = ej + 1;
        *(qu_e++) = ek;
      }
    } else if ((q & b3) == b3 && ej < hy - 1) {
      *(mijk + hx) = mv;
      *(qu_e++) = ei;
      *(qu_e++) = ej + 1;
      *(qu_e++) = ek;
    }
    if ((q & b6) == b6) {
      if (ek > 0) {
        *(mijk - hxy) = mv;
        *(qu_e++) = ei;
        *(qu_e++) = ej;
        *(qu_e++) = ek - 1;
      }
      if ((q & b5) == 0 && ek < hz - 1) {
        *(mijk + hxy) = mv;
        *(qu_e++) = ei;
        *(qu_e++) = ej;
        *(qu_e++) = ek + 1;
      }
    } else if ((q & b5) == b5 && ek < hz - 1) {
      *(mijk + hxy) = mv;
      *(qu_e++) = ei;
      *(qu_e++) = ej;
      *(qu_e++) = ek + 1;
    }
  }

/** This routine computes a Voronoi cell for a single particle in the
   * container. It can be called by the user, but is also forms the core part of
   * several of the main functions, such as store_cell_volumes(), print_all(),
   * and the drawing routines. The algorithm constructs the cell by testing over
   * the neighbors of the particle, working outwards until it reaches those
   * particles which could not possibly intersect the cell. For maximum
   * efficiency, this algorithm is divided into three parts. In the first
   * section, the algorithm tests over the blocks which are in the immediate
   * vicinity of the particle, by making use of one of the precomputed worklists.
   * The code then continues to test blocks on the worklist, but also begins to
   * construct a list of neighboring blocks outside the worklist which may need
   * to be test. In the third section, the routine starts testing these
   * neighboring blocks, evaluating whether or not a particle in them could
   * possibly intersect the cell. For blocks that intersect the cell, it tests
   * the particles in that block, and then adds the block neighbors to the list
   * of potential places to consider.
   * \param[in,out] c a reference to a voronoicell object.
   * \param[in] ijk the index of the block that the test particle is in.
   * \param[in] s the index of the particle within the test block.
   * \param[in] (ci,cj,ck) the coordinates of the block that the test particle is
   *                       in relative to the container data structure.
   * \return False if the Voronoi cell was completely removed during the
   *         computation and has zero volume, true otherwise. */
  template <class c_class> template <class v_cell> bool voro_compute<c_class>::compute_cell(v_cell &c, int ijk, int s, int ci, int cj, int ck) {
    static const int count_list[8] = {7, 11, 15, 19, 26, 35, 45, 59};
    static const int *count_e = count_list + 8;
    double x;
    double y;
    double z;
    double x1;
    double y1;
    double z1;
    double qx = 0;
    double qy = 0;
    double qz = 0;
    double xlo;
    double ylo;
    double zlo;
    double xhi;
    double yhi;
    double zhi;
    double x2;
    double y2;
    double z2;
    double rs;
    int i;
    int j;
    int k;
    int di;
    int dj;
    int dk;
    int ei;
    int ej;
    int ek;
    int f;
    int g;
    int l;
    int disp;
    double fx;
    double fy;
    double fz;
    double gxs;
    double gys;
    double gzs;
    double *radp;
    unsigned int q;
    unsigned int *e;
    unsigned int *mijk;

    if (!con.initialize_voronoicell(c, ijk, s, ci, cj, ck, i, j, k, x, y, z, disp)) {
      return false;
    }
    con.r_init(ijk, s);

    // Initialize the Voronoi cell to fill the entire container
    double crs;
    double mrs;

    int next_count = 3;
    int *count_p = (const_cast<int *>(count_list));

    // Test all particles in the particle's local region first
    for (l = 0; l < s; l++) {
      x1 = p[ijk][ps * l] - x;
      y1 = p[ijk][ps * l + 1] - y;
      z1 = p[ijk][ps * l + 2] - z;
      rs = con.r_scale(x1 * x1 + y1 * y1 + z1 * z1, ijk, l);
      if (!c.nplane(x1, y1, z1, rs, id[ijk][l])) {
        return false;
      }
    }
    l++;
    while (l < co[ijk]) {
      x1 = p[ijk][ps * l] - x;
      y1 = p[ijk][ps * l + 1] - y;
      z1 = p[ijk][ps * l + 2] - z;
      rs = con.r_scale(x1 * x1 + y1 * y1 + z1 * z1, ijk, l);
      if (!c.nplane(x1, y1, z1, rs, id[ijk][l])) {
        return false;
      }
      l++;
    }

    // Now compute the maximum distance squared from the cell center to a
    // vertex. This is used to cut off the calculation since we only need
    // to test out to twice this range.
    mrs = c.max_radius_squared();

    // Now compute the fractional position of the particle within its
    // region and store it in (fx,fy,fz). We use this to compute an index
    // (di,dj,dk) of which subregion the particle is within.
    unsigned int m1;
    unsigned int m2;
    con.frac_pos(x, y, z, ci, cj, ck, fx, fy, fz);
    di = int(fx * xsp * wl_fgrid);
    dj = int(fy * ysp * wl_fgrid);
    dk = int(fz * zsp * wl_fgrid);

    // The indices (di,dj,dk) tell us which worklist to use, to test the
    // blocks in the optimal order. But we only store worklists for the
    // eighth of the region where di, dj, and dk are all less than half the
    // full grid. The rest of the cases are handled by symmetry. In this
    // section, we detect for these cases, by reflecting high values of di,
    // dj, and dk. For these cases, a mask is constructed in m1 and m2
    // which is used to flip the worklist information when it is loaded.
    if (di >= wl_hgrid) {
      gxs = fx;
      m1 = 127 + (3 << 21);
      m2 = 1 + (1 << 21);
      di = wl_fgrid - 1 - di;
      if (di < 0) {
        di = 0;
      }
    } else {
      m1 = m2 = 0;
      gxs = boxx - fx;
    }
    if (dj >= wl_hgrid) {
      gys = fy;
      m1 |= (127 << 7) + (3 << 24);
      m2 |= (1 << 7) + (1 << 24);
      dj = wl_fgrid - 1 - dj;
      if (dj < 0) {
        dj = 0;
      }
    } else {
      gys = boxy - fy;
    }
    if (dk >= wl_hgrid) {
      gzs = fz;
      m1 |= (127 << 14) + (3 << 27);
      m2 |= (1 << 14) + (1 << 27);
      dk = wl_fgrid - 1 - dk;
      if (dk < 0) {
        dk = 0;
      }
    } else {
      gzs = boxz - fz;
    }
    gxs *= gxs;
    gys *= gys;
    gzs *= gzs;

    // Now compute which worklist we are going to use, and set radp and e to
    // point at the right offsets
    ijk = di + wl_hgrid * (dj + wl_hgrid * dk);
    radp = mrad + ijk * wl_seq_length;
    e = (const_cast<unsigned int *>(wl)) + ijk * wl_seq_length;

    // Read in how many items in the worklist can be tested without having to
    // worry about writing to the mask
    f = e[0];
    g = 0;
    do {
      // At the intervals specified by count_list, we recompute the
      // maximum radius squared
      if (g == next_count) {
        mrs = c.max_radius_squared();
        if (count_p != count_e) {
          next_count = *(count_p++);
        }
      }

      // If mrs is less than the minimum distance to any untested
      // block, then we are done
      if (con.r_ctest(radp[g], mrs)) {
        return true;
      }
      g++;

      // Load in a block off the worklist, permute it with the
      // symmetry mask, and decode its position. These are all
      // integer bit operations so they should run very fast.
      q = e[g];
      q ^= m1;
      q += m2;
      di = q & 127;
      di -= 64;
      dj = (q >> 7) & 127;
      dj -= 64;
      dk = (q >> 14) & 127;
      dk -= 64;

      // Check that the worklist position is in range
      ei = di + i;
      if (ei < 0 || ei >= hx) {
        continue;
      }
      ej = dj + j;
      if (ej < 0 || ej >= hy) {
        continue;
      }
      ek = dk + k;
      if (ek < 0 || ek >= hz) {
        continue;
      }

      // Call the compute_min_max_radius() function. This returns
      // true if the minimum distance to the block is bigger than the
      // current mrs, in which case we skip this block and move on.
      // Otherwise, it computes the maximum distance to the block and
      // returns it in crs.
      if (compute_min_max_radius(di, dj, dk, fx, fy, fz, gxs, gys, gzs, crs, mrs)) {
        continue;
      }

      // Now compute which region we are going to loop over, adding a
      // displacement for the periodic cases
      ijk = con.region_index(ci, cj, ck, ei, ej, ek, qx, qy, qz, disp);

      // If mrs is bigger than the maximum distance to the block,
      // then we have to test all particles in the block for
      // intersections. Otherwise, we do additional checks and skip
      // those particles which can't possibly intersect the block.
      if (co[ijk] > 0) {
        l = 0;
        x2 = x - qx;
        y2 = y - qy;
        z2 = z - qz;
        if (!con.r_ctest(crs, mrs)) {
          do {
            x1 = p[ijk][ps * l] - x2;
            y1 = p[ijk][ps * l + 1] - y2;
            z1 = p[ijk][ps * l + 2] - z2;
            rs = con.r_scale(x1 * x1 + y1 * y1 + z1 * z1, ijk, l);
            if (!c.nplane(x1, y1, z1, rs, id[ijk][l])) {
              return false;
            }
            l++;
          } while (l < co[ijk]);
        } else {
          do {
            x1 = p[ijk][ps * l] - x2;
            y1 = p[ijk][ps * l + 1] - y2;
            z1 = p[ijk][ps * l + 2] - z2;
            rs = x1 * x1 + y1 * y1 + z1 * z1;
            if (con.r_scale_check(rs, mrs, ijk, l) && !c.nplane(x1, y1, z1, rs, id[ijk][l])) {
              return false;
            }
            l++;
          } while (l < co[ijk]);
        }
      }
    } while (g < f);

    // If we reach here, we were unable to compute the entire cell using
    // the first part of the worklist. This section of the algorithm
    // continues the worklist, but it now starts preparing the mask that we
    // need if we end up going block by block. We do the same as before,
    // but we put a mark down on the mask for every block that's tested.
    // The worklist also contains information about which neighbors of each
    // block are not also on the worklist, and we start storing those
    // points in a list in case we have to go block by block. Update the
    // mask counter, and if it wraps around then reset the whole mask; that
    // will only happen once every 2^32 tries.
    mv++;
    if (mv == 0) {
      reset_mask();
      mv = 1;
    }

    // Set the queue pointers
    int *qu_s = qu;
    int *qu_e = qu;

    while (g < wl_seq_length - 1) {
      // At the intervals specified by count_list, we recompute the
      // maximum radius squared
      if (g == next_count) {
        mrs = c.max_radius_squared();
        if (count_p != count_e) {
          next_count = *(count_p++);
        }
      }

      // If mrs is less than the minimum distance to any untested
      // block, then we are done
      if (con.r_ctest(radp[g], mrs)) {
        return true;
      }
      g++;

      // Load in a block off the worklist, permute it with the
      // symmetry mask, and decode its position. These are all
      // integer bit operations so they should run very fast.
      q = e[g];
      q ^= m1;
      q += m2;
      di = q & 127;
      di -= 64;
      dj = (q >> 7) & 127;
      dj -= 64;
      dk = (q >> 14) & 127;
      dk -= 64;

      // Compute the position in the mask of the current block. If
      // this lies outside the mask, then skip it. Otherwise, mark
      // it.
      ei = di + i;
      if (ei < 0 || ei >= hx) {
        continue;
      }
      ej = dj + j;
      if (ej < 0 || ej >= hy) {
        continue;
      }
      ek = dk + k;
      if (ek < 0 || ek >= hz) {
        continue;
      }
      mijk = mask + ei + hx * (ej + hy * ek);
      *mijk = mv;

      // Call the compute_min_max_radius() function. This returns
      // true if the minimum distance to the block is bigger than the
      // current mrs, in which case we skip this block and move on.
      // Otherwise, it computes the maximum distance to the block and
      // returns it in crs.
      if (compute_min_max_radius(di, dj, dk, fx, fy, fz, gxs, gys, gzs, crs, mrs)) {
        continue;
      }

      // Now compute which region we are going to loop over, adding a
      // displacement for the periodic cases
      ijk = con.region_index(ci, cj, ck, ei, ej, ek, qx, qy, qz, disp);

      // If mrs is bigger than the maximum distance to the block,
      // then we have to test all particles in the block for
      // intersections. Otherwise, we do additional checks and skip
      // those particles which can't possibly intersect the block.
      if (co[ijk] > 0) {
        l = 0;
        x2 = x - qx;
        y2 = y - qy;
        z2 = z - qz;
        if (!con.r_ctest(crs, mrs)) {
          do {
            x1 = p[ijk][ps * l] - x2;
            y1 = p[ijk][ps * l + 1] - y2;
            z1 = p[ijk][ps * l + 2] - z2;
            rs = con.r_scale(x1 * x1 + y1 * y1 + z1 * z1, ijk, l);
            if (!c.nplane(x1, y1, z1, rs, id[ijk][l])) {
              return false;
            }
            l++;
          } while (l < co[ijk]);
        } else {
          do {
            x1 = p[ijk][ps * l] - x2;
            y1 = p[ijk][ps * l + 1] - y2;
            z1 = p[ijk][ps * l + 2] - z2;
            rs = x1 * x1 + y1 * y1 + z1 * z1;
            if (con.r_scale_check(rs, mrs, ijk, l) && !c.nplane(x1, y1, z1, rs, id[ijk][l])) {
              return false;
            }
            l++;
          } while (l < co[ijk]);
        }
      }

      // If there might not be enough memory on the list for these
      // additions, then add more
      if (qu_e > qu_l - 18) {
        add_list_memory(qu_s, qu_e);
      }

      // Test the parts of the worklist element which tell us what
      // neighbors of this block are not on the worklist. Store them
      // on the block list, and mark the mask.
      scan_bits_mask_add(q, mijk, ei, ej, ek, qu_e);
    }

    // Do a check to see if we've reached the radius cutoff
    if (con.r_ctest(radp[g], mrs)) {
      return true;
    }

    // We were unable to completely compute the cell based on the blocks in
    // the worklist, so now we have to go block by block, reading in items
    // off the list
    while (qu_s != qu_e) {
      // If we reached the end of the list memory loop back to the
      // start
      if (qu_s == qu_l) {
        qu_s = qu;
      }

      // Read in a block off the list, and compute the upper and lower
      // coordinates in each of the three dimensions
      ei = *(qu_s++);
      ej = *(qu_s++);
      ek = *(qu_s++);
      xlo = (ei - i) * boxx - fx;
      xhi = xlo + boxx;
      ylo = (ej - j) * boxy - fy;
      yhi = ylo + boxy;
      zlo = (ek - k) * boxz - fz;
      zhi = zlo + boxz;

      // Carry out plane tests to see if any particle in this block
      // could possibly intersect the cell
      if (ei > i) {
        if (ej > j) {
          if (ek > k) {
            if (corner_test(c, xlo, ylo, zlo, xhi, yhi, zhi)) {
              continue;
            }
          } else if (ek < k) {
            if (corner_test(c, xlo, ylo, zhi, xhi, yhi, zlo)) {
              continue;
            }
          } else {
            if (edge_z_test(c, xlo, ylo, zlo, xhi, yhi, zhi)) {
              continue;
            }
          }
        } else if (ej < j) {
          if (ek > k) {
            if (corner_test(c, xlo, yhi, zlo, xhi, ylo, zhi)) {
              continue;
            }
          } else if (ek < k) {
            if (corner_test(c, xlo, yhi, zhi, xhi, ylo, zlo)) {
              continue;
            }
          } else {
            if (edge_z_test(c, xlo, yhi, zlo, xhi, ylo, zhi)) {
              continue;
            }
          }
        } else {
          if (ek > k) {
            if (edge_y_test(c, xlo, ylo, zlo, xhi, yhi, zhi)) {
              continue;
            }
          } else if (ek < k) {
            if (edge_y_test(c, xlo, ylo, zhi, xhi, yhi, zlo)) {
              continue;
            }
          } else {
            if (face_x_test(c, xlo, ylo, zlo, yhi, zhi)) {
              continue;
            }
          }
        }
      } else if (ei < i) {
        if (ej > j) {
          if (ek > k) {
            if (corner_test(c, xhi, ylo, zlo, xlo, yhi, zhi)) {
              continue;
            }
          } else if (ek < k) {
            if (corner_test(c, xhi, ylo, zhi, xlo, yhi, zlo)) {
              continue;
            }
          } else {
            if (edge_z_test(c, xhi, ylo, zlo, xlo, yhi, zhi)) {
              continue;
            }
          }
        } else if (ej < j) {
          if (ek > k) {
            if (corner_test(c, xhi, yhi, zlo, xlo, ylo, zhi)) {
              continue;
            }
          } else if (ek < k) {
            if (corner_test(c, xhi, yhi, zhi, xlo, ylo, zlo)) {
              continue;
            }
          } else {
            if (edge_z_test(c, xhi, yhi, zlo, xlo, ylo, zhi)) {
              continue;
            }
          }
        } else {
          if (ek > k) {
            if (edge_y_test(c, xhi, ylo, zlo, xlo, yhi, zhi)) {
              continue;
            }
          } else if (ek < k) {
            if (edge_y_test(c, xhi, ylo, zhi, xlo, yhi, zlo)) {
              continue;
            }
          } else {
            if (face_x_test(c, xhi, ylo, zlo, yhi, zhi)) {
              continue;
            }
          }
        }
      } else {
        if (ej > j) {
          if (ek > k) {
            if (edge_x_test(c, xlo, ylo, zlo, xhi, yhi, zhi)) {
              continue;
            }
          } else if (ek < k) {
            if (edge_x_test(c, xlo, ylo, zhi, xhi, yhi, zlo)) {
              continue;
            }
          } else {
            if (face_y_test(c, xlo, ylo, zlo, xhi, zhi)) {
              continue;
            }
          }
        } else if (ej < j) {
          if (ek > k) {
            if (edge_x_test(c, xlo, yhi, zlo, xhi, ylo, zhi)) {
              continue;
            }
          } else if (ek < k) {
            if (edge_x_test(c, xlo, yhi, zhi, xhi, ylo, zlo)) {
              continue;
            }
          } else {
            if (face_y_test(c, xlo, yhi, zlo, xhi, zhi)) {
              continue;
            }
          }
        } else {
          if (ek > k) {
            if (face_z_test(c, xlo, ylo, zlo, xhi, yhi)) {
              continue;
            }
          } else if (ek < k) {
            if (face_z_test(c, xlo, ylo, zhi, xhi, yhi)) {
              continue;
            }
          } else {
            voro_fatal_error("Compute cell routine revisiting central block, which should never\nhappen.", VOROPP_INTERNAL_ERROR);
          }
        }
      }

      // Now compute the region that we are going to test over, and
      // set a displacement vector for the periodic cases
      ijk = con.region_index(ci, cj, ck, ei, ej, ek, qx, qy, qz, disp);

      // Loop over all the elements in the block to test for cuts. It
      // would be possible to exclude some of these cases by testing
      // against mrs, but this will probably not save time.
      if (co[ijk] > 0) {
        l = 0;
        x2 = x - qx;
        y2 = y - qy;
        z2 = z - qz;
        do {
          x1 = p[ijk][ps * l] - x2;
          y1 = p[ijk][ps * l + 1] - y2;
          z1 = p[ijk][ps * l + 2] - z2;
          rs = con.r_scale(x1 * x1 + y1 * y1 + z1 * z1, ijk, l);
          if (!c.nplane(x1, y1, z1, rs, id[ijk][l])) {
            return false;
          }
          l++;
        } while (l < co[ijk]);
      }

      // If there's not much memory on the block list then add more
      if ((qu_s <= qu_e ? (qu_l - qu_e) + (qu_s - qu) : qu_s - qu_e) < 18) {
        add_list_memory(qu_s, qu_e);
      }

      // Test the neighbors of the current block, and add them to the
      // block list if they haven't already been tested
      add_to_mask(ei, ej, ek, qu_e);
    }

    return true;
  }

/** This function checks to see whether a particular block can possibly have
   * any intersection with a Voronoi cell, for the case when the closest point
   * from the cell center to the block is at a corner.
   * \param[in,out] c a reference to a Voronoi cell.
   * \param[in] (xl,yl,zl) the relative coordinates of the corner of the block
   *                       closest to the cell center.
   * \param[in] (xh,yh,zh) the relative coordinates of the corner of the block
   *                       furthest away from the cell center.
   * \return False if the block may intersect, true if does not. */
  template <class c_class> template <class v_cell> bool voro_compute<c_class>::corner_test(v_cell &c, double xl, double yl, double zl, double xh, double yh, double zh) {
    con.r_prime(xl * xl + yl * yl + zl * zl);
    if (c.plane_intersects_guess(xh, yl, zl, con.r_cutoff(xl * xh + yl * yl + zl * zl))) {
      return false;
    }
    if (c.plane_intersects(xh, yh, zl, con.r_cutoff(xl * xh + yl * yh + zl * zl))) {
      return false;
    }
    if (c.plane_intersects(xl, yh, zl, con.r_cutoff(xl * xl + yl * yh + zl * zl))) {
      return false;
    }
    if (c.plane_intersects(xl, yh, zh, con.r_cutoff(xl * xl + yl * yh + zl * zh))) {
      return false;
    }
    if (c.plane_intersects(xl, yl, zh, con.r_cutoff(xl * xl + yl * yl + zl * zh))) {
      return false;
    }
    if (c.plane_intersects(xh, yl, zh, con.r_cutoff(xl * xh + yl * yl + zl * zh))) {
      return false;
    }
    return true;
  }

/** This function checks to see whether a particular block can possibly have
   * any intersection with a Voronoi cell, for the case when the closest point
   * from the cell center to the block is on an edge which points along the x
   * direction.
   * \param[in,out] c a reference to a Voronoi cell.
   * \param[in] (x0,x1) the minimum and maximum relative x coordinates of the
   *                    block.
   * \param[in] (yl,zl) the relative y and z coordinates of the corner of the
   *                    block closest to the cell center.
   * \param[in] (yh,zh) the relative y and z coordinates of the corner of the
   *                    block furthest away from the cell center.
   * \return False if the block may intersect, true if does not. */
  template <class c_class> template <class v_cell> inline bool voro_compute<c_class>::edge_x_test(v_cell &c, double x0, double yl, double zl, double x1, double yh, double zh) {
    con.r_prime(yl * yl + zl * zl);
    if (c.plane_intersects_guess(x0, yl, zh, con.r_cutoff(yl * yl + zl * zh))) {
      return false;
    }
    if (c.plane_intersects(x1, yl, zh, con.r_cutoff(yl * yl + zl * zh))) {
      return false;
    }
    if (c.plane_intersects(x1, yl, zl, con.r_cutoff(yl * yl + zl * zl))) {
      return false;
    }
    if (c.plane_intersects(x0, yl, zl, con.r_cutoff(yl * yl + zl * zl))) {
      return false;
    }
    if (c.plane_intersects(x0, yh, zl, con.r_cutoff(yl * yh + zl * zl))) {
      return false;
    }
    if (c.plane_intersects(x1, yh, zl, con.r_cutoff(yl * yh + zl * zl))) {
      return false;
    }
    return true;
  }

/** This function checks to see whether a particular block can possibly have
   * any intersection with a Voronoi cell, for the case when the closest point
   * from the cell center to the block is on an edge which points along the y
   * direction.
   * \param[in,out] c a reference to a Voronoi cell.
   * \param[in] (y0,y1) the minimum and maximum relative y coordinates of the
   *                    block.
   * \param[in] (xl,zl) the relative x and z coordinates of the corner of the
   *                    block closest to the cell center.
   * \param[in] (xh,zh) the relative x and z coordinates of the corner of the
   *                    block furthest away from the cell center.
   * \return False if the block may intersect, true if does not. */
  template <class c_class> template <class v_cell> inline bool voro_compute<c_class>::edge_y_test(v_cell &c, double xl, double y0, double zl, double xh, double y1, double zh) {
    con.r_prime(xl * xl + zl * zl);
    if (c.plane_intersects_guess(xl, y0, zh, con.r_cutoff(xl * xl + zl * zh))) {
      return false;
    }
    if (c.plane_intersects(xl, y1, zh, con.r_cutoff(xl * xl + zl * zh))) {
      return false;
    }
    if (c.plane_intersects(xl, y1, zl, con.r_cutoff(xl * xl + zl * zl))) {
      return false;
    }
    if (c.plane_intersects(xl, y0, zl, con.r_cutoff(xl * xl + zl * zl))) {
      return false;
    }
    if (c.plane_intersects(xh, y0, zl, con.r_cutoff(xl * xh + zl * zl))) {
      return false;
    }
    if (c.plane_intersects(xh, y1, zl, con.r_cutoff(xl * xh + zl * zl))) {
      return false;
    }
    return true;
  }

/** This function checks to see whether a particular block can possibly have
   * any intersection with a Voronoi cell, for the case when the closest point
   * from the cell center to the block is on an edge which points along the z
   * direction.
   * \param[in,out] c a reference to a Voronoi cell.
   * \param[in] (z0,z1) the minimum and maximum relative z coordinates of the block.
   * \param[in] (xl,yl) the relative x and y coordinates of the corner of the
   *                    block closest to the cell center.
   * \param[in] (xh,yh) the relative x and y coordinates of the corner of the
   *                    block furthest away from the cell center.
   * \return False if the block may intersect, true if does not. */
  template <class c_class> template <class v_cell> inline bool voro_compute<c_class>::edge_z_test(v_cell &c, double xl, double yl, double z0, double xh, double yh, double z1) {
    con.r_prime(xl * xl + yl * yl);
    if (c.plane_intersects_guess(xl, yh, z0, con.r_cutoff(xl * xl + yl * yh))) {
      return false;
    }
    if (c.plane_intersects(xl, yh, z1, con.r_cutoff(xl * xl + yl * yh))) {
      return false;
    }
    if (c.plane_intersects(xl, yl, z1, con.r_cutoff(xl * xl + yl * yl))) {
      return false;
    }
    if (c.plane_intersects(xl, yl, z0, con.r_cutoff(xl * xl + yl * yl))) {
      return false;
    }
    if (c.plane_intersects(xh, yl, z0, con.r_cutoff(xl * xh + yl * yl))) {
      return false;
    }
    if (c.plane_intersects(xh, yl, z1, con.r_cutoff(xl * xh + yl * yl))) {
      return false;
    }
    return true;
  }

/** This function checks to see whether a particular block can possibly have
   * any intersection with a Voronoi cell, for the case when the closest point
   * from the cell center to the block is on a face aligned with the x direction.
   * \param[in,out] c a reference to a Voronoi cell.
   * \param[in] xl the minimum distance from the cell center to the face.
   * \param[in] (y0,y1) the minimum and maximum relative y coordinates of the
   *                    block.
   * \param[in] (z0,z1) the minimum and maximum relative z coordinates of the
   *                    block.
   * \return False if the block may intersect, true if does not. */
  template <class c_class> template <class v_cell> inline bool voro_compute<c_class>::face_x_test(v_cell &c, double xl, double y0, double z0, double y1, double z1) {
    con.r_prime(xl * xl);
    if (c.plane_intersects_guess(xl, y0, z0, con.r_cutoff(xl * xl))) {
      return false;
    }
    if (c.plane_intersects(xl, y0, z1, con.r_cutoff(xl * xl))) {
      return false;
    }
    if (c.plane_intersects(xl, y1, z1, con.r_cutoff(xl * xl))) {
      return false;
    }
    if (c.plane_intersects(xl, y1, z0, con.r_cutoff(xl * xl))) {
      return false;
    }
    return true;
  }

/** This function checks to see whether a particular block can possibly have
   * any intersection with a Voronoi cell, for the case when the closest point
   * from the cell center to the block is on a face aligned with the y direction.
   * \param[in,out] c a reference to a Voronoi cell.
   * \param[in] yl the minimum distance from the cell center to the face.
   * \param[in] (x0,x1) the minimum and maximum relative x coordinates of the
   *                    block.
   * \param[in] (z0,z1) the minimum and maximum relative z coordinates of the
   *                    block.
   * \return False if the block may intersect, true if does not. */
  template <class c_class> template <class v_cell> inline bool voro_compute<c_class>::face_y_test(v_cell &c, double x0, double yl, double z0, double x1, double z1) {
    con.r_prime(yl * yl);
    if (c.plane_intersects_guess(x0, yl, z0, con.r_cutoff(yl * yl))) {
      return false;
    }
    if (c.plane_intersects(x0, yl, z1, con.r_cutoff(yl * yl))) {
      return false;
    }
    if (c.plane_intersects(x1, yl, z1, con.r_cutoff(yl * yl))) {
      return false;
    }
    if (c.plane_intersects(x1, yl, z0, con.r_cutoff(yl * yl))) {
      return false;
    }
    return true;
  }

/** This function checks to see whether a particular block can possibly have
   * any intersection with a Voronoi cell, for the case when the closest point
   * from the cell center to the block is on a face aligned with the z direction.
   * \param[in,out] c a reference to a Voronoi cell.
   * \param[in] zl the minimum distance from the cell center to the face.
   * \param[in] (x0,x1) the minimum and maximum relative x coordinates of the
   *                    block.
   * \param[in] (y0,y1) the minimum and maximum relative y coordinates of the
   *                    block.
   * \return False if the block may intersect, true if does not. */
  template <class c_class> template <class v_cell> inline bool voro_compute<c_class>::face_z_test(v_cell &c, double x0, double y0, double zl, double x1, double y1) {
    con.r_prime(zl * zl);
    if (c.plane_intersects_guess(x0, y0, zl, con.r_cutoff(zl * zl))) {
      return false;
    }
    if (c.plane_intersects(x0, y1, zl, con.r_cutoff(zl * zl))) {
      return false;
    }
    if (c.plane_intersects(x1, y1, zl, con.r_cutoff(zl * zl))) {
      return false;
    }
    if (c.plane_intersects(x1, y0, zl, con.r_cutoff(zl * zl))) {
      return false;
    }
    return true;
  }

/** This routine checks to see whether a point is within a particular distance
   * of a nearby region. If the point is within the distance of the region, then
   * the routine returns true, and computes the maximum distance from the point
   * to the region. Otherwise, the routine returns false.
   * \param[in] (di,dj,dk) the position of the nearby region to be tested,
   *                       relative to the region that the point is in.
   * \param[in] (fx,fy,fz) the displacement of the point within its region.
   * \param[in] (gxs,gys,gzs) the maximum squared distances from the point to the
   *                          sides of its region.
   * \param[out] crs a reference in which to return the maximum distance to the
   *                 region (only computed if the routine returns false).
   * \param[in] mrs the distance to be tested.
   * \return True if the region is further away than mrs, false if the region in
   *         within mrs. */
  template <class c_class> bool voro_compute<c_class>::compute_min_max_radius(int di, int dj, int dk, double fx, double fy, double fz, double gxs, double gys, double gzs, double &crs, double mrs) {
    double xlo;
    double ylo;
    double zlo;
    if (di > 0) {
      xlo = di * boxx - fx;
      crs = xlo * xlo;
      if (dj > 0) {
        ylo = dj * boxy - fy;
        crs += ylo * ylo;
        if (dk > 0) {
          zlo = dk * boxz - fz;
          crs += zlo * zlo;
          if (con.r_ctest(crs, mrs)) {
            return true;
          }
          crs += bxsq + 2 * (boxx * xlo + boxy * ylo + boxz * zlo);
        } else if (dk < 0) {
          zlo = (dk + 1) * boxz - fz;
          crs += zlo * zlo;
          if (con.r_ctest(crs, mrs)) {
            return true;
          }
          crs += bxsq + 2 * (boxx * xlo + boxy * ylo - boxz * zlo);
        } else {
          if (con.r_ctest(crs, mrs)) {
            return true;
          }
          crs += boxx * (2 * xlo + boxx) + boxy * (2 * ylo + boxy) + gzs;
        }
      } else if (dj < 0) {
        ylo = (dj + 1) * boxy - fy;
        crs += ylo * ylo;
        if (dk > 0) {
          zlo = dk * boxz - fz;
          crs += zlo * zlo;
          if (con.r_ctest(crs, mrs)) {
            return true;
          }
          crs += bxsq + 2 * (boxx * xlo - boxy * ylo + boxz * zlo);
        } else if (dk < 0) {
          zlo = (dk + 1) * boxz - fz;
          crs += zlo * zlo;
          if (con.r_ctest(crs, mrs)) {
            return true;
          }
          crs += bxsq + 2 * (boxx * xlo - boxy * ylo - boxz * zlo);
        } else {
          if (con.r_ctest(crs, mrs)) {
            return true;
          }
          crs += boxx * (2 * xlo + boxx) + boxy * (-2 * ylo + boxy) + gzs;
        }
      } else {
        if (dk > 0) {
          zlo = dk * boxz - fz;
          crs += zlo * zlo;
          if (con.r_ctest(crs, mrs)) {
            return true;
          }
          crs += boxz * (2 * zlo + boxz);
        } else if (dk < 0) {
          zlo = (dk + 1) * boxz - fz;
          crs += zlo * zlo;
          if (con.r_ctest(crs, mrs)) {
            return true;
          }
          crs += boxz * (-2 * zlo + boxz);
        } else {
          if (con.r_ctest(crs, mrs)) {
            return true;
          }
          crs += gzs;
        }
        crs += gys + boxx * (2 * xlo + boxx);
      }
    } else if (di < 0) {
      xlo = (di + 1) * boxx - fx;
      crs = xlo * xlo;
      if (dj > 0) {
        ylo = dj * boxy - fy;
        crs += ylo * ylo;
        if (dk > 0) {
          zlo = dk * boxz - fz;
          crs += zlo * zlo;
          if (con.r_ctest(crs, mrs)) {
            return true;
          }
          crs += bxsq + 2 * (-boxx * xlo + boxy * ylo + boxz * zlo);
        } else if (dk < 0) {
          zlo = (dk + 1) * boxz - fz;
          crs += zlo * zlo;
          if (con.r_ctest(crs, mrs)) {
            return true;
          }
          crs += bxsq + 2 * (-boxx * xlo + boxy * ylo - boxz * zlo);
        } else {
          if (con.r_ctest(crs, mrs)) {
            return true;
          }
          crs += boxx * (-2 * xlo + boxx) + boxy * (2 * ylo + boxy) + gzs;
        }
      } else if (dj < 0) {
        ylo = (dj + 1) * boxy - fy;
        crs += ylo * ylo;
        if (dk > 0) {
          zlo = dk * boxz - fz;
          crs += zlo * zlo;
          if (con.r_ctest(crs, mrs)) {
            return true;
          }
          crs += bxsq + 2 * (-boxx * xlo - boxy * ylo + boxz * zlo);
        } else if (dk < 0) {
          zlo = (dk + 1) * boxz - fz;
          crs += zlo * zlo;
          if (con.r_ctest(crs, mrs)) {
            return true;
          }
          crs += bxsq + 2 * (-boxx * xlo - boxy * ylo - boxz * zlo);
        } else {
          if (con.r_ctest(crs, mrs)) {
            return true;
          }
          crs += boxx * (-2 * xlo + boxx) + boxy * (-2 * ylo + boxy) + gzs;
        }
      } else {
        if (dk > 0) {
          zlo = dk * boxz - fz;
          crs += zlo * zlo;
          if (con.r_ctest(crs, mrs)) {
            return true;
          }
          crs += boxz * (2 * zlo + boxz);
        } else if (dk < 0) {
          zlo = (dk + 1) * boxz - fz;
          crs += zlo * zlo;
          if (con.r_ctest(crs, mrs)) {
            return true;
          }
          crs += boxz * (-2 * zlo + boxz);
        } else {
          if (con.r_ctest(crs, mrs)) {
            return true;
          }
          crs += gzs;
        }
        crs += gys + boxx * (-2 * xlo + boxx);
      }
    } else {
      if (dj > 0) {
        ylo = dj * boxy - fy;
        crs = ylo * ylo;
        if (dk > 0) {
          zlo = dk * boxz - fz;
          crs += zlo * zlo;
          if (con.r_ctest(crs, mrs)) {
            return true;
          }
          crs += boxz * (2 * zlo + boxz);
        } else if (dk < 0) {
          zlo = (dk + 1) * boxz - fz;
          crs += zlo * zlo;
          if (con.r_ctest(crs, mrs)) {
            return true;
          }
          crs += boxz * (-2 * zlo + boxz);
        } else {
          if (con.r_ctest(crs, mrs)) {
            return true;
          }
          crs += gzs;
        }
        crs += boxy * (2 * ylo + boxy);
      } else if (dj < 0) {
        ylo = (dj + 1) * boxy - fy;
        crs = ylo * ylo;
        if (dk > 0) {
          zlo = dk * boxz - fz;
          crs += zlo * zlo;
          if (con.r_ctest(crs, mrs)) {
            return true;
          }
          crs += boxz * (2 * zlo + boxz);
        } else if (dk < 0) {
          zlo = (dk + 1) * boxz - fz;
          crs += zlo * zlo;
          if (con.r_ctest(crs, mrs)) {
            return true;
          }
          crs += boxz * (-2 * zlo + boxz);
        } else {
          if (con.r_ctest(crs, mrs)) {
            return true;
          }
          crs += gzs;
        }
        crs += boxy * (-2 * ylo + boxy);
      } else {
        if (dk > 0) {
          zlo = dk * boxz - fz;
          crs = zlo * zlo;
          if (con.r_ctest(crs, mrs)) {
            return true;
          }
          crs += boxz * (2 * zlo + boxz);
        } else if (dk < 0) {
          zlo = (dk + 1) * boxz - fz;
          crs = zlo * zlo;
          if (con.r_ctest(crs, mrs)) {
            return true;
          }
          crs += boxz * (-2 * zlo + boxz);
        } else {
          crs = 0;
          voro_fatal_error("Min/max radius function called for central block, which should never\nhappen.", VOROPP_INTERNAL_ERROR);
        }
        crs += gys;
      }
      crs += gxs;
    }
    return false;
  }

  template <class c_class> bool voro_compute<c_class>::compute_min_radius(int di, int dj, int dk, double fx, double fy, double fz, double mrs) {
    double t;
    double crs;

    if (di > 0) {
      t = di * boxx - fx;
      crs = t * t;
    } else if (di < 0) {
      t = (di + 1) * boxx - fx;
      crs = t * t;
    } else {
      crs = 0;
    }

    if (dj > 0) {
      t = dj * boxy - fy;
      crs += t * t;
    } else if (dj < 0) {
      t = (dj + 1) * boxy - fy;
      crs += t * t;
    }

    if (dk > 0) {
      t = dk * boxz - fz;
      crs += t * t;
    } else if (dk < 0) {
      t = (dk + 1) * boxz - fz;
      crs += t * t;
    }

    return crs > con.r_max_add(mrs);
  }

/** Adds memory to the queue.
   * \param[in,out] qu_s a reference to the queue start pointer.
   * \param[in,out] qu_e a reference to the queue end pointer. */
  template <class c_class> inline void voro_compute<c_class>::add_list_memory(int *&qu_s, int *&qu_e) {
    qu_size <<= 1;
    int *qu_n = new int[qu_size], *qu_c = qu_n;
#if VOROPP_VERBOSE >= 2
    fprintf(stderr, "List memory scaled up to %d\n", qu_size);
#endif
    if (qu_s <= qu_e) {
      while (qu_s < qu_e) {
        *(qu_c++) = *(qu_s++);
      }
    } else {
      while (qu_s < qu_l) {
        *(qu_c++) = *(qu_s++);
      }
      qu_s = qu;
      while (qu_s < qu_e) {
        *(qu_c++) = *(qu_s++);
      }
    }
    delete[] qu;
    qu_s = qu = qu_n;
    qu_l = qu + qu_size;
    qu_e = qu_c;
  }

// Explicit template instantiation
  template voro_compute<container>::voro_compute(container &, int, int, int);
  template voro_compute<container_poly>::voro_compute(container_poly &, int, int, int);
  template bool voro_compute<container>::compute_cell(voronoicell &, int, int, int, int, int);
  template bool voro_compute<container>::compute_cell(voronoicell_neighbor &, int, int, int, int, int);
  template void voro_compute<container>::find_voronoi_cell(double, double, double, int, int, int, int, particle_record &, double &);
  template bool voro_compute<container_poly>::compute_cell(voronoicell &, int, int, int, int, int);
  template bool voro_compute<container_poly>::compute_cell(voronoicell_neighbor &, int, int, int, int, int);
  template void voro_compute<container_poly>::find_voronoi_cell(double, double, double, int, int, int, int, particle_record &, double &);

// Explicit template instantiation
  template voro_compute<container_periodic>::voro_compute(container_periodic &, int, int, int);
  template voro_compute<container_periodic_poly>::voro_compute(container_periodic_poly &, int, int, int);
  template bool voro_compute<container_periodic>::compute_cell(voronoicell &, int, int, int, int, int);
  template bool voro_compute<container_periodic>::compute_cell(voronoicell_neighbor &, int, int, int, int, int);
  template void voro_compute<container_periodic>::find_voronoi_cell(double, double, double, int, int, int, int, particle_record &, double &);
  template bool voro_compute<container_periodic_poly>::compute_cell(voronoicell &, int, int, int, int, int);
  template bool voro_compute<container_periodic_poly>::compute_cell(voronoicell_neighbor &, int, int, int, int, int);
  template void voro_compute<container_periodic_poly>::find_voronoi_cell(double, double, double, int, int, int, int, particle_record &, double &);

} // namespace voro
// Voro++, a 3D cell-based Voronoi library
//
// Author   : Chris H. Rycroft (LBL / UC Berkeley)
// Email    : chr@alum.mit.edu
// Date     : August 30th 2011

/** \file container.cc
 * \brief Function implementations for the container and related classes. */

namespace voro {

/** The class constructor sets up the geometry of container, initializing the
   * minimum and maximum coordinates in each direction, and setting whether each
   * direction is periodic or not. It divides the container into a rectangular
   * grid of blocks, and allocates memory for each of these for storing particle
   * positions and IDs.
   * \param[in] (ax_,bx_) the minimum and maximum x coordinates.
   * \param[in] (ay_,by_) the minimum and maximum y coordinates.
   * \param[in] (az_,bz_) the minimum and maximum z coordinates.
   * \param[in] (nx_,ny_,nz_) the number of grid blocks in each of the three
   *			    coordinate directions.
   * \param[in] (xperiodic_,yperiodic_,zperiodic_) flags setting whether the
   *                                               container is periodic in each
   *                                               coordinate direction.
   * \param[in] init_mem the initial memory allocation for each block.
   * \param[in] ps_ the number of floating point entries to store for each
   *                particle. */
  container_base::container_base(double ax_, double bx_, double ay_, double by_, double az_, double bz_, int nx_, int ny_, int nz_, bool xperiodic_, bool yperiodic_, bool zperiodic_, int init_mem, int ps_) :
      voro_base(nx_, ny_, nz_, (bx_ - ax_) / nx_, (by_ - ay_) / ny_, (bz_ - az_) / nz_),
      ax(ax_),
      bx(bx_),
      ay(ay_),
      by(by_),
      az(az_),
      bz(bz_),
      xperiodic(xperiodic_),
      yperiodic(yperiodic_),
      zperiodic(zperiodic_),
      id(new int *[nxyz]),
      p(new double *[nxyz]),
      co(new int[nxyz]),
      mem(new int[nxyz]),
      ps(ps_) {
    int l;
    for (l = 0; l < nxyz; l++) {
      co[l] = 0;
    }
    for (l = 0; l < nxyz; l++) {
      mem[l] = init_mem;
    }
    for (l = 0; l < nxyz; l++) {
      id[l] = new int[init_mem];
    }
    for (l = 0; l < nxyz; l++) {
      p[l] = new double[ps * init_mem];
    }
  }

/** The container destructor frees the dynamically allocated memory. */
  container_base::~container_base() {
    int l;
    for (l = 0; l < nxyz; l++) {
      delete[] p[l];
    }
    for (l = 0; l < nxyz; l++) {
      delete[] id[l];
    }
    delete[] id;
    delete[] p;
    delete[] co;
    delete[] mem;
  }

/** The class constructor sets up the geometry of container.
   * \param[in] (ax_,bx_) the minimum and maximum x coordinates.
   * \param[in] (ay_,by_) the minimum and maximum y coordinates.
   * \param[in] (az_,bz_) the minimum and maximum z coordinates.
   * \param[in] (nx_,ny_,nz_) the number of grid blocks in each of the three
   *                       coordinate directions.
   * \param[in] (xperiodic_,yperiodic_,zperiodic_) flags setting whether the
   *                                               container is periodic in each
   *                                               coordinate direction.
   * \param[in] init_mem the initial memory allocation for each block. */
  container::container(double ax_, double bx_, double ay_, double by_, double az_, double bz_, int nx_, int ny_, int nz_, bool xperiodic_, bool yperiodic_, bool zperiodic_, int init_mem) :
      container_base(ax_, bx_, ay_, by_, az_, bz_, nx_, ny_, nz_, xperiodic_, yperiodic_, zperiodic_, init_mem, 3),
      vc(*this, xperiodic_ ? 2 * nx_ + 1 : nx_, yperiodic_ ? 2 * ny_ + 1 : ny_, zperiodic_ ? 2 * nz_ + 1 : nz_) {}

/** The class constructor sets up the geometry of container.
   * \param[in] (ax_,bx_) the minimum and maximum x coordinates.
   * \param[in] (ay_,by_) the minimum and maximum y coordinates.
   * \param[in] (az_,bz_) the minimum and maximum z coordinates.
   * \param[in] (nx_,ny_,nz_) the number of grid blocks in each of the three
   *                       coordinate directions.
   * \param[in] (xperiodic_,yperiodic_,zperiodic_) flags setting whether the
   *                                               container is periodic in each
   *                                               coordinate direction.
   * \param[in] init_mem the initial memory allocation for each block. */
  container_poly::container_poly(double ax_, double bx_, double ay_, double by_, double az_, double bz_, int nx_, int ny_, int nz_, bool xperiodic_, bool yperiodic_, bool zperiodic_, int init_mem) :
      container_base(ax_, bx_, ay_, by_, az_, bz_, nx_, ny_, nz_, xperiodic_, yperiodic_, zperiodic_, init_mem, 4),
      vc(*this, xperiodic_ ? 2 * nx_ + 1 : nx_, yperiodic_ ? 2 * ny_ + 1 : ny_, zperiodic_ ? 2 * nz_ + 1 : nz_) {
    ppr = p;
  }

/** Put a particle into the correct region of the container.
   * \param[in] n the numerical ID of the inserted particle.
   * \param[in] (x,y,z) the position vector of the inserted particle. */
  void container::put(int n, double x, double y, double z) {
    int ijk;
    if (put_locate_block(ijk, x, y, z)) {
      id[ijk][co[ijk]] = n;
      double *pp = p[ijk] + 3 * co[ijk]++;
      *(pp++) = x;
      *(pp++) = y;
      *pp = z;
    }
  }

/** Put a particle into the correct region of the container.
   * \param[in] n the numerical ID of the inserted particle.
   * \param[in] (x,y,z) the position vector of the inserted particle.
   * \param[in] r the radius of the particle. */
  void container_poly::put(int n, double x, double y, double z, double r) {
    int ijk;
    if (put_locate_block(ijk, x, y, z)) {
      id[ijk][co[ijk]] = n;
      double *pp = p[ijk] + 4 * co[ijk]++;
      *(pp++) = x;
      *(pp++) = y;
      *(pp++) = z;
      *pp = r;
      if (max_radius < r) {
        max_radius = r;
      }
    }
  }

/** Put a particle into the correct region of the container, also recording
   * into which region it was stored.
   * \param[in] vo the ordering class in which to record the region.
   * \param[in] n the numerical ID of the inserted particle.
   * \param[in] (x,y,z) the position vector of the inserted particle. */
  void container::put(particle_order &vo, int n, double x, double y, double z) {
    int ijk;
    if (put_locate_block(ijk, x, y, z)) {
      id[ijk][co[ijk]] = n;
      vo.add(ijk, co[ijk]);
      double *pp = p[ijk] + 3 * co[ijk]++;
      *(pp++) = x;
      *(pp++) = y;
      *pp = z;
    }
  }

/** Put a particle into the correct region of the container, also recording
   * into which region it was stored.
   * \param[in] vo the ordering class in which to record the region.
   * \param[in] n the numerical ID of the inserted particle.
   * \param[in] (x,y,z) the position vector of the inserted particle.
   * \param[in] r the radius of the particle. */
  void container_poly::put(particle_order &vo, int n, double x, double y, double z, double r) {
    int ijk;
    if (put_locate_block(ijk, x, y, z)) {
      id[ijk][co[ijk]] = n;
      vo.add(ijk, co[ijk]);
      double *pp = p[ijk] + 4 * co[ijk]++;
      *(pp++) = x;
      *(pp++) = y;
      *(pp++) = z;
      *pp = r;
      if (max_radius < r) {
        max_radius = r;
      }
    }
  }

/** This routine takes a particle position vector, tries to remap it into the
   * primary domain. If successful, it computes the region into which it can be
   * stored and checks that there is enough memory within this region to store
   * it.
   * \param[out] ijk the region index.
   * \param[in,out] (x,y,z) the particle position, remapped into the primary
   *                        domain if necessary.
   * \return True if the particle can be successfully placed into the container,
   * false otherwise. */
  bool container_base::put_locate_block(int &ijk, double &x, double &y, double &z) {
    if (put_remap(ijk, x, y, z)) {
      if (co[ijk] == mem[ijk]) {
        add_particle_memory(ijk);
      }
      return true;
    }
#if VOROPP_REPORT_OUT_OF_BOUNDS == 1
    fprintf(stderr, "Out of bounds: (x,y,z)=(%g,%g,%g)\n", x, y, z);
#endif
    return false;
  }

/** Takes a particle position vector and computes the region index into which
   * it should be stored. If the container is periodic, then the routine also
   * maps the particle position to ensure it is in the primary domain. If the
   * container is not periodic, the routine bails out.
   * \param[out] ijk the region index.
   * \param[in,out] (x,y,z) the particle position, remapped into the primary
   *                        domain if necessary.
   * \return True if the particle can be successfully placed into the container,
   * false otherwise. */
  inline bool container_base::put_remap(int &ijk, double &x, double &y, double &z) {
    int l;

    ijk = step_int((x - ax) * xsp);
    if (xperiodic) {
      l = step_mod(ijk, nx);
      x += boxx * (l - ijk);
      ijk = l;
    } else if (ijk < 0 || ijk >= nx) {
      return false;
    }

    int j = step_int((y - ay) * ysp);
    if (yperiodic) {
      l = step_mod(j, ny);
      y += boxy * (l - j);
      j = l;
    } else if (j < 0 || j >= ny) {
      return false;
    }

    int k = step_int((z - az) * zsp);
    if (zperiodic) {
      l = step_mod(k, nz);
      z += boxz * (l - k);
      k = l;
    } else if (k < 0 || k >= nz) {
      return false;
    }

    ijk += nx * j + nxy * k;
    return true;
  }

/** Takes a position vector and attempts to remap it into the primary domain.
   * \param[out] (ai,aj,ak) the periodic image displacement that the vector is in,
   *                       with (0,0,0) corresponding to the primary domain.
   * \param[out] (ci,cj,ck) the index of the block that the position vector is
   *                        within, once it has been remapped.
   * \param[in,out] (x,y,z) the position vector to consider, which is remapped
   *                        into the primary domain during the routine.
   * \param[out] ijk the block index that the vector is within.
   * \return True if the particle is within the container or can be remapped into
   * it, false if it lies outside of the container bounds. */
  inline bool container_base::remap(int &ai, int &aj, int &ak, int &ci, int &cj, int &ck, double &x, double &y, double &z, int &ijk) {
    ci = step_int((x - ax) * xsp);
    if (ci < 0 || ci >= nx) {
      if (xperiodic) {
        ai = step_div(ci, nx);
        x -= ai * (bx - ax);
        ci -= ai * nx;
      } else {
        return false;
      }
    } else {
      ai = 0;
    }

    cj = step_int((y - ay) * ysp);
    if (cj < 0 || cj >= ny) {
      if (yperiodic) {
        aj = step_div(cj, ny);
        y -= aj * (by - ay);
        cj -= aj * ny;
      } else {
        return false;
      }
    } else {
      aj = 0;
    }

    ck = step_int((z - az) * zsp);
    if (ck < 0 || ck >= nz) {
      if (zperiodic) {
        ak = step_div(ck, nz);
        z -= ak * (bz - az);
        ck -= ak * nz;
      } else {
        return false;
      }
    } else {
      ak = 0;
    }

    ijk = ci + nx * cj + nxy * ck;
    return true;
  }

/** Takes a vector and finds the particle whose Voronoi cell contains that
   * vector. This is equivalent to finding the particle which is nearest to the
   * vector. Additional wall classes are not considered by this routine.
   * \param[in] (x,y,z) the vector to test.
   * \param[out] (rx,ry,rz) the position of the particle whose Voronoi cell
   *                        contains the vector. If the container is periodic,
   *                        this may point to a particle in a periodic image of
   *                        the primary domain.
   * \param[out] pid the ID of the particle.
   * \return True if a particle was found. If the container has no particles,
   * then the search will not find a Voronoi cell and false is returned. */
  bool container::find_voronoi_cell(double x, double y, double z, double &rx, double &ry, double &rz, int &pid) {
    int ai;
    int aj;
    int ak;
    int ci;
    int cj;
    int ck;
    int ijk;
    particle_record w;
    double mrs;

    // If the given vector lies outside the domain, but the container
    // is periodic, then remap it back into the domain
    if (!remap(ai, aj, ak, ci, cj, ck, x, y, z, ijk)) {
      return false;
    }
    vc.find_voronoi_cell(x, y, z, ci, cj, ck, ijk, w, mrs);

    if (w.ijk != -1) {
      // Assemble the position vector of the particle to be returned,
      // applying a periodic remapping if necessary
      if (xperiodic) {
        ci += w.di;
        if (ci < 0 || ci >= nx) {
          ai += step_div(ci, nx);
        }
      }
      if (yperiodic) {
        cj += w.dj;
        if (cj < 0 || cj >= ny) {
          aj += step_div(cj, ny);
        }
      }
      if (zperiodic) {
        ck += w.dk;
        if (ck < 0 || ck >= nz) {
          ak += step_div(ck, nz);
        }
      }
      rx = p[w.ijk][3 * w.l] + ai * (bx - ax);
      ry = p[w.ijk][3 * w.l + 1] + aj * (by - ay);
      rz = p[w.ijk][3 * w.l + 2] + ak * (bz - az);
      pid = id[w.ijk][w.l];
      return true;
    }

    // If no particle if found then just return false
    return false;
  }

/** Takes a vector and finds the particle whose Voronoi cell contains that
   * vector. Additional wall classes are not considered by this routine.
   * \param[in] (x,y,z) the vector to test.
   * \param[out] (rx,ry,rz) the position of the particle whose Voronoi cell
   *                        contains the vector. If the container is periodic,
   *                        this may point to a particle in a periodic image of
   *                        the primary domain.
   * \param[out] pid the ID of the particle.
   * \return True if a particle was found. If the container has no particles,
   * then the search will not find a Voronoi cell and false is returned. */
  bool container_poly::find_voronoi_cell(double x, double y, double z, double &rx, double &ry, double &rz, int &pid) {
    int ai;
    int aj;
    int ak;
    int ci;
    int cj;
    int ck;
    int ijk;
    particle_record w;
    double mrs;

    // If the given vector lies outside the domain, but the container
    // is periodic, then remap it back into the domain
    if (!remap(ai, aj, ak, ci, cj, ck, x, y, z, ijk)) {
      return false;
    }
    vc.find_voronoi_cell(x, y, z, ci, cj, ck, ijk, w, mrs);

    if (w.ijk != -1) {
      // Assemble the position vector of the particle to be returned,
      // applying a periodic remapping if necessary
      if (xperiodic) {
        ci += w.di;
        if (ci < 0 || ci >= nx) {
          ai += step_div(ci, nx);
        }
      }
      if (yperiodic) {
        cj += w.dj;
        if (cj < 0 || cj >= ny) {
          aj += step_div(cj, ny);
        }
      }
      if (zperiodic) {
        ck += w.dk;
        if (ck < 0 || ck >= nz) {
          ak += step_div(ck, nz);
        }
      }
      rx = p[w.ijk][4 * w.l] + ai * (bx - ax);
      ry = p[w.ijk][4 * w.l + 1] + aj * (by - ay);
      rz = p[w.ijk][4 * w.l + 2] + ak * (bz - az);
      pid = id[w.ijk][w.l];
      return true;
    }

    // If no particle if found then just return false
    return false;
  }

/** Increase memory for a particular region.
   * \param[in] i the index of the region to reallocate. */
  void container_base::add_particle_memory(int i) const {
    int l;
    const int nmem = mem[i] << 1;

    // Carry out a check on the memory allocation size, and
    // print a status message if requested
    if (nmem > max_particle_memory) {
      voro_fatal_error("Absolute maximum memory allocation exceeded", VOROPP_MEMORY_ERROR);
    }
#if VOROPP_VERBOSE >= 3
    fprintf(stderr, "Particle memory in region %d scaled up to %d\n", i, nmem);
#endif

    // Allocate new memory and copy in the contents of the old arrays
    int *idp = new int[nmem];
    for (l = 0; l < co[i]; l++) {
      idp[l] = id[i][l];
    }
    double *pp = new double[ps * nmem];
    for (l = 0; l < ps * co[i]; l++) {
      pp[l] = p[i][l];
    }

    // Update pointers and delete old arrays
    mem[i] = nmem;
    delete[] id[i];
    id[i] = idp;
    delete[] p[i];
    p[i] = pp;
  }

/** Import a list of particles from an open file stream into the container.
   * Entries of four numbers (Particle ID, x position, y position, z position)
   * are searched for. If the file cannot be successfully read, then the routine
   * causes a fatal error.
   * \param[in] fp the file handle to read from. */
  void container::import(FILE *fp) {
    int i;
    int j;
    double x;
    double y;
    double z;
    while ((j = fscanf(fp, "%d %lg %lg %lg", &i, &x, &y, &z)) == 4) {
      put(i, x, y, z);
    }
    if (j != EOF) {
      voro_fatal_error("File import error", VOROPP_FILE_ERROR);
    }
  }

/** Import a list of particles from an open file stream, also storing the order
   * of that the particles are read. Entries of four numbers (Particle ID, x
   * position, y position, z position) are searched for. If the file cannot be
   * successfully read, then the routine causes a fatal error.
   * \param[in,out] vo a reference to an ordering class to use.
   * \param[in] fp the file handle to read from. */
  void container::import(particle_order &vo, FILE *fp) {
    int i;
    int j;
    double x;
    double y;
    double z;
    while ((j = fscanf(fp, "%d %lg %lg %lg", &i, &x, &y, &z)) == 4) {
      put(vo, i, x, y, z);
    }
    if (j != EOF) {
      voro_fatal_error("File import error", VOROPP_FILE_ERROR);
    }
  }

/** Import a list of particles from an open file stream into the container.
   * Entries of five numbers (Particle ID, x position, y position, z position,
   * radius) are searched for. If the file cannot be successfully read, then the
   * routine causes a fatal error.
   * \param[in] fp the file handle to read from. */
  void container_poly::import(FILE *fp) {
    int i;
    int j;
    double x;
    double y;
    double z;
    double r;
    while ((j = fscanf(fp, "%d %lg %lg %lg %lg", &i, &x, &y, &z, &r)) == 5) {
      put(i, x, y, z, r);
    }
    if (j != EOF) {
      voro_fatal_error("File import error", VOROPP_FILE_ERROR);
    }
  }

/** Import a list of particles from an open file stream, also storing the order
   * of that the particles are read. Entries of four numbers (Particle ID, x
   * position, y position, z position, radius) are searched for. If the file
   * cannot be successfully read, then the routine causes a fatal error.
   * \param[in,out] vo a reference to an ordering class to use.
   * \param[in] fp the file handle to read from. */
  void container_poly::import(particle_order &vo, FILE *fp) {
    int i;
    int j;
    double x;
    double y;
    double z;
    double r;
    while ((j = fscanf(fp, "%d %lg %lg %lg %lg", &i, &x, &y, &z, &r)) == 5) {
      put(vo, i, x, y, z, r);
    }
    if (j != EOF) {
      voro_fatal_error("File import error", VOROPP_FILE_ERROR);
    }
  }

/** Outputs the a list of all the container regions along with the number of
   * particles stored within each. */
  void container_base::region_count() {
    int i;
    int j;
    int k;
    int *cop = co;
    for (k = 0; k < nz; k++) {
      for (j = 0; j < ny; j++) {
        for (i = 0; i < nx; i++) {
          printf("Region (%d,%d,%d): %d particles\n", i, j, k, *(cop++));
        }
      }
    }
  }

/** Clears a container of particles. */
  void container::clear() {
    for (int *cop = co; cop < co + nxyz; cop++) {
      *cop = 0;
    }
  }

/** Clears a container of particles, also clearing resetting the maximum radius
   * to zero. */
  void container_poly::clear() {
    for (int *cop = co; cop < co + nxyz; cop++) {
      *cop = 0;
    }
    max_radius = 0;
  }

/** Computes all the Voronoi cells and saves customized information about them.
   * \param[in] format the custom output string to use.
   * \param[in] fp a file handle to write to. */
  void container::print_custom(const char *format, FILE *fp) {
    c_loop_all vl(*this);
    print_custom(vl, format, fp);
  }

/** Computes all the Voronoi cells and saves customized
   * information about them.
   * \param[in] format the custom output string to use.
   * \param[in] fp a file handle to write to. */
  void container_poly::print_custom(const char *format, FILE *fp) {
    c_loop_all vl(*this);
    print_custom(vl, format, fp);
  }

/** Computes all the Voronoi cells and saves customized information about them.
   * \param[in] format the custom output string to use.
   * \param[in] filename the name of the file to write to. */
  void container::print_custom(const char *format, const char *filename) {
    FILE *fp = safe_fopen(filename, "w");
    print_custom(format, fp);
    fclose(fp);
  }

/** Computes all the Voronoi cells and saves customized
   * information about them
   * \param[in] format the custom output string to use.
   * \param[in] filename the name of the file to write to. */
  void container_poly::print_custom(const char *format, const char *filename) {
    FILE *fp = safe_fopen(filename, "w");
    print_custom(format, fp);
    fclose(fp);
  }

/** Computes all of the Voronoi cells in the container, but does nothing
   * with the output. It is useful for measuring the pure computation time
   * of the Voronoi algorithm, without any additional calculations such as
   * volume evaluation or cell output. */
  void container::compute_all_cells() {
    voronoicell c;
    c_loop_all vl(*this);
    if (vl.start()) {
      do {
        compute_cell(c, vl);
      } while (vl.inc());
    }
  }

/** Computes all of the Voronoi cells in the container, but does nothing
   * with the output. It is useful for measuring the pure computation time
   * of the Voronoi algorithm, without any additional calculations such as
   * volume evaluation or cell output. */
  void container_poly::compute_all_cells() {
    voronoicell c;
    c_loop_all vl(*this);
    if (vl.start()) {
      do {
        compute_cell(c, vl);
      } while (vl.inc());
    }
  }

/** Calculates all of the Voronoi cells and sums their volumes. In most cases
   * without walls, the sum of the Voronoi cell volumes should equal the volume
   * of the container to numerical precision.
   * \return The sum of all of the computed Voronoi volumes. */
  double container::sum_cell_volumes() {
    voronoicell c;
    double vol = 0;
    c_loop_all vl(*this);
    if (vl.start()) {
      do {
        if (compute_cell(c, vl)) {
          vol += c.volume();
        }
      } while (vl.inc());
    }
    return vol;
  }

/** Calculates all of the Voronoi cells and sums their volumes. In most cases
   * without walls, the sum of the Voronoi cell volumes should equal the volume
   * of the container to numerical precision.
   * \return The sum of all of the computed Voronoi volumes. */
  double container_poly::sum_cell_volumes() {
    voronoicell c;
    double vol = 0;
    c_loop_all vl(*this);
    if (vl.start()) {
      do {
        if (compute_cell(c, vl)) {
          vol += c.volume();
        }
      } while (vl.inc());
    }
    return vol;
  }

/** This function tests to see if a given vector lies within the container
   * bounds and any walls.
   * \param[in] (x,y,z) the position vector to be tested.
   * \return True if the point is inside the container, false if the point is
   *         outside. */
  bool container_base::point_inside(double x, double y, double z) {
    if (x < ax || x > bx || y < ay || y > by || z < az || z > bz) {
      return false;
    }
    return point_inside_walls(x, y, z);
  }

/** Draws an outline of the domain in gnuplot format.
   * \param[in] fp the file handle to write to. */
  void container_base::draw_domain_gnuplot(FILE *fp) const {
    fprintf(fp, "%g %g %g\n%g %g %g\n%g %g %g\n%g %g %g\n", ax, ay, az, bx, ay, az, bx, by, az, ax, by, az);
    fprintf(fp, "%g %g %g\n%g %g %g\n%g %g %g\n%g %g %g\n", ax, by, bz, bx, by, bz, bx, ay, bz, ax, ay, bz);
    fprintf(fp, "%g %g %g\n\n%g %g %g\n%g %g %g\n\n", ax, by, bz, ax, ay, az, ax, ay, bz);
    fprintf(fp, "%g %g %g\n%g %g %g\n\n%g %g %g\n%g %g %g\n\n", bx, ay, az, bx, ay, bz, bx, by, az, bx, by, bz);
  }

/** Draws an outline of the domain in POV-Ray format.
   * \param[in] fp the file handle to write to. */
  void container_base::draw_domain_pov(FILE *fp) const {
    fprintf(fp,
            "cylinder{<%g,%g,%g>,<%g,%g,%g>,rr}\n"
            "cylinder{<%g,%g,%g>,<%g,%g,%g>,rr}\n",
            ax, ay, az, bx, ay, az, ax, by, az, bx, by, az);
    fprintf(fp,
            "cylinder{<%g,%g,%g>,<%g,%g,%g>,rr}\n"
            "cylinder{<%g,%g,%g>,<%g,%g,%g>,rr}\n",
            ax, by, bz, bx, by, bz, ax, ay, bz, bx, ay, bz);
    fprintf(fp,
            "cylinder{<%g,%g,%g>,<%g,%g,%g>,rr}\n"
            "cylinder{<%g,%g,%g>,<%g,%g,%g>,rr}\n",
            ax, ay, az, ax, by, az, bx, ay, az, bx, by, az);
    fprintf(fp,
            "cylinder{<%g,%g,%g>,<%g,%g,%g>,rr}\n"
            "cylinder{<%g,%g,%g>,<%g,%g,%g>,rr}\n",
            bx, ay, bz, bx, by, bz, ax, ay, bz, ax, by, bz);
    fprintf(fp,
            "cylinder{<%g,%g,%g>,<%g,%g,%g>,rr}\n"
            "cylinder{<%g,%g,%g>,<%g,%g,%g>,rr}\n",
            ax, ay, az, ax, ay, bz, bx, ay, az, bx, ay, bz);
    fprintf(fp,
            "cylinder{<%g,%g,%g>,<%g,%g,%g>,rr}\n"
            "cylinder{<%g,%g,%g>,<%g,%g,%g>,rr}\n",
            bx, by, az, bx, by, bz, ax, by, az, ax, by, bz);
    fprintf(fp,
            "sphere{<%g,%g,%g>,rr}\nsphere{<%g,%g,%g>,rr}\n"
            "sphere{<%g,%g,%g>,rr}\nsphere{<%g,%g,%g>,rr}\n",
            ax, ay, az, bx, ay, az, ax, by, az, bx, by, az);
    fprintf(fp,
            "sphere{<%g,%g,%g>,rr}\nsphere{<%g,%g,%g>,rr}\n"
            "sphere{<%g,%g,%g>,rr}\nsphere{<%g,%g,%g>,rr}\n",
            ax, ay, bz, bx, ay, bz, ax, by, bz, bx, by, bz);
  }

/** The wall_list constructor sets up an array of pointers to wall classes. */
  wall_list::wall_list() : walls(new wall *[init_wall_size]), wep(walls), wel(walls + init_wall_size), current_wall_size(init_wall_size) {}

/** The wall_list destructor frees the array of pointers to the wall classes.
   */
  wall_list::~wall_list() {
    delete[] walls;
  }

/** Adds all of the walls on another wall_list to this class.
   * \param[in] wl a reference to the wall class. */
  void wall_list::add_wall(wall_list &wl) {
    for (wall **wp = wl.walls; wp < wl.wep; wp++) {
      add_wall(*wp);
    }
  }

/** Deallocates all of the wall classes pointed to by the wall_list. */
  void wall_list::deallocate() const {
    for (wall **wp = walls; wp < wep; wp++) {
      delete *wp;
    }
  }

/** Increases the memory allocation for the walls array. */
  void wall_list::increase_wall_memory() {
    current_wall_size <<= 1;
    if (current_wall_size > max_wall_size) {
      voro_fatal_error("Wall memory allocation exceeded absolute maximum", VOROPP_MEMORY_ERROR);
    }
    wall **nwalls = new wall *[current_wall_size];
    wall **nwp = nwalls;
    wall **wp = walls;
    while (wp < wep) {
      *(nwp++) = *(wp++);
    }
    delete[] walls;
    walls = nwalls;
    wel = walls + current_wall_size;
    wep = nwp;
  }

} // namespace voro
// Voro++, a 3D cell-based Voronoi library
//
// Author   : Chris H. Rycroft (LBL / UC Berkeley)
// Email    : chr@alum.mit.edu
// Date     : August 30th 2011

/** \file unitcell.cc
 * \brief Function implementations for the unitcell class. */

#include <cmath>
#include <queue>

namespace voro {

/** Initializes the unit cell class for a particular non-orthogonal periodic
   * geometry, corresponding to a parallelepiped with sides given by three
   * vectors. The class constructs the unit Voronoi cell corresponding to this
   * geometry.
   * \param[in] (bx_) The x coordinate of the first unit vector.
   * \param[in] (bxy_,by_) The x and y coordinates of the second unit vector.
   * \param[in] (bxz_,byz_,bz_) The x, y, and z coordinates of the third unit
   *                            vector. */
  unitcell::unitcell(double bx_, double bxy_, double by_, double bxz_, double byz_, double bz_) : bx(bx_), bxy(bxy_), by(by_), bxz(bxz_), byz(byz_), bz(bz_) {
    int i;
    int j;
    int l = 1;

    // Initialize the Voronoi cell to be a very large rectangular box
    const double ucx = max_unit_voro_shells * bx;
    const double ucy = max_unit_voro_shells * by;
    const double ucz = max_unit_voro_shells * bz;
    unit_voro.init(-ucx, ucx, -ucy, ucy, -ucz, ucz);

    // Repeatedly cut the cell by shells of periodic image particles
    while (l < 2 * max_unit_voro_shells) {
      // Check to see if any of the planes from the current shell
      // will cut the cell
      if (unit_voro_intersect(l)) {
        // If they do, apply the plane cuts from the current
        // shell
        unit_voro_apply(l, 0, 0);
        for (i = 1; i < l; i++) {
          unit_voro_apply(l, i, 0);
          unit_voro_apply(-l, i, 0);
        }
        for (i = -l; i <= l; i++) {
          unit_voro_apply(i, l, 0);
        }
        for (i = 1; i < l; i++) {
          for (j = -l + 1; j <= l; j++) {
            unit_voro_apply(l, j, i);
            unit_voro_apply(-j, l, i);
            unit_voro_apply(-l, -j, i);
            unit_voro_apply(j, -l, i);
          }
        }
        for (i = -l; i <= l; i++) {
          for (j = -l; j <= l; j++) {
            unit_voro_apply(i, j, l);
          }
        }
      } else {
        // Calculate a bound on the maximum y and z coordinates
        // that could possibly cut the cell. This is based upon
        // a geometric result that particles with z>l can't cut
        // a cell lying within the paraboloid
        // z<=(l*l-x*x-y*y)/(2*l). It is always a tighter bound
        // than the one based on computing the maximum radius
        // of a Voronoi cell vertex.
        max_uv_y = max_uv_z = 0;
        double y;
        double z;
        double q;
        double *pts = unit_voro.pts;
        double *pp = pts;
        while (pp < pts + 3 * unit_voro.p) {
          q = *(pp++);
          y = *(pp++);
          z = *(pp++);
          q = sqrt(q * q + y * y + z * z);
          if (y + q > max_uv_y) {
            max_uv_y = y + q;
          }
          if (z + q > max_uv_z) {
            max_uv_z = z + q;
          }
        }
        max_uv_z *= 0.5;
        max_uv_y *= 0.5;
        return;
      }
      l++;
    }

    // If the routine makes it here, then the unit cell still hasn't been
    // completely bounded by the plane cuts. Give the memory error code,
    // because this is mainly a case of hitting a safe limit, than any
    // inherent problem.
    voro_fatal_error("Periodic cell computation failed", VOROPP_MEMORY_ERROR);
  }

/** Applies a pair of opposing plane cuts from a periodic image point
   * to the unit Voronoi cell.
   * \param[in] (i,j,k) the index of the periodic image to consider. */
  inline void unitcell::unit_voro_apply(int i, int j, int k) {
    const double x = i * bx + j * bxy + k * bxz;
    const double y = j * by + k * byz;
    const double z = k * bz;
    unit_voro.plane(x, y, z);
    unit_voro.plane(-x, -y, -z);
  }

/** Calculates whether the unit Voronoi cell intersects a given periodic image
   * of the domain.
   * \param[in] (dx,dy,dz) the displacement of the periodic image.
   * \param[out] vol the proportion of the unit cell volume within this image,
   *                 only computed in the case that the two intersect.
   * \return True if they intersect, false otherwise. */
  bool unitcell::intersects_image(double dx, double dy, double dz, double &vol) {
    const double bxinv = 1 / bx;
    const double byinv = 1 / by;
    const double bzinv = 1 / bz;
    const double ivol = bxinv * byinv * bzinv;
    voronoicell c;
    c = unit_voro;
    dx *= 2;
    dy *= 2;
    dz *= 2;
    if (!c.plane(0, 0, bzinv, dz + 1)) {
      return false;
    }
    if (!c.plane(0, 0, -bzinv, -dz + 1)) {
      return false;
    }
    if (!c.plane(0, byinv, -byz * byinv * bzinv, dy + 1)) {
      return false;
    }
    if (!c.plane(0, -byinv, byz * byinv * bzinv, -dy + 1)) {
      return false;
    }
    if (!c.plane(bxinv, -bxy * bxinv * byinv, (bxy * byz - by * bxz) * ivol, dx + 1)) {
      return false;
    }
    if (!c.plane(-bxinv, bxy * bxinv * byinv, (-bxy * byz + by * bxz) * ivol, -dx + 1)) {
      return false;
    }
    vol = c.volume() * ivol;
    return true;
  }

/** Computes a list of periodic domain images that intersect the unit Voronoi cell.
   * \param[out] vi a vector containing triplets (i,j,k) corresponding to domain
   *                images that intersect the unit Voronoi cell, when it is
   *                centered in the middle of the primary domain.
   * \param[out] vd a vector containing the fraction of the Voronoi cell volume
   *                within each corresponding image listed in vi. */
  void unitcell::images(std::vector<int> &vi, std::vector<double> &vd) {
    const int ms2 = max_unit_voro_shells * 2 + 1;
    const int mss = ms2 * ms2 * ms2;
    bool *a = new bool[mss];
    bool *ac = a + max_unit_voro_shells * (1 + ms2 * (1 + ms2));
    bool *ap = a;
    int i;
    int j;
    int k;
    double vol;

    // Initialize mask
    while (ap < ac) {
      *(ap++) = true;
    }
    *(ap++) = false;
    while (ap < a + mss) {
      *(ap++) = true;
    }

    // Set up the queue and add (0,0,0) image to it
    std::queue<int> q;
    q.push(0);
    q.push(0);
    q.push(0);

    while (!q.empty()) {
      // Read the next entry on the queue
      i = q.front();
      q.pop();
      j = q.front();
      q.pop();
      k = q.front();
      q.pop();

      // Check intersection of this image
      if (intersects_image(i, j, k, vol)) {
        // Add this entry to the output vectors
        vi.push_back(i);
        vi.push_back(j);
        vi.push_back(k);
        vd.push_back(vol);

        // Add neighbors to the queue if they have not been
        // tested
        ap = ac + i + ms2 * (j + ms2 * k);
        if (k > -max_unit_voro_shells && *(ap - ms2 * ms2)) {
          q.push(i);
          q.push(j);
          q.push(k - 1);
          *(ap - ms2 * ms2) = false;
        }
        if (j > -max_unit_voro_shells && *(ap - ms2)) {
          q.push(i);
          q.push(j - 1);
          q.push(k);
          *(ap - ms2) = false;
        }
        if (i > -max_unit_voro_shells && *(ap - 1)) {
          q.push(i - 1);
          q.push(j);
          q.push(k);
          *(ap - 1) = false;
        }
        if (i < max_unit_voro_shells && *(ap + 1)) {
          q.push(i + 1);
          q.push(j);
          q.push(k);
          *(ap + 1) = false;
        }
        if (j < max_unit_voro_shells && *(ap + ms2)) {
          q.push(i);
          q.push(j + 1);
          q.push(k);
          *(ap + ms2) = false;
        }
        if (k < max_unit_voro_shells && *(ap + ms2 * ms2)) {
          q.push(i);
          q.push(j);
          q.push(k + 1);
          *(ap + ms2 * ms2) = false;
        }
      }
    }

    // Remove mask memory
    delete[] a;
  }

/** Tests to see if a shell of periodic images could possibly cut the periodic
   * unit cell.
   * \param[in] l the index of the shell to consider.
   * \return True if a point in the shell cuts the cell, false otherwise. */
  bool unitcell::unit_voro_intersect(int l) {
    int i;
    int j;
    if (unit_voro_test(l, 0, 0)) {
      return true;
    }
    for (i = 1; i < l; i++) {
      if (unit_voro_test(l, i, 0)) {
        return true;
      }
      if (unit_voro_test(-l, i, 0)) {
        return true;
      }
    }
    for (i = -l; i <= l; i++) {
      if (unit_voro_test(i, l, 0)) {
        return true;
      }
    }
    for (i = 1; i < l; i++) {
      for (j = -l + 1; j <= l; j++) {
        if (unit_voro_test(l, j, i)) {
          return true;
        }
        if (unit_voro_test(-j, l, i)) {
          return true;
        }
        if (unit_voro_test(-l, -j, i)) {
          return true;
        }
        if (unit_voro_test(j, -l, i)) {
          return true;
        }
      }
    }
    for (i = -l; i <= l; i++) {
      for (j = -l; j <= l; j++) {
        if (unit_voro_test(i, j, l)) {
          return true;
        }
      }
    }
    return false;
  }

/** Tests to see if a plane cut from a particular periodic image will cut th
   * unit Voronoi cell.
   * \param[in] (i,j,k) the index of the periodic image to consider.
   * \return True if the image cuts the cell, false otherwise. */
  inline bool unitcell::unit_voro_test(int i, int j, int k) {
    const double x = i * bx + j * bxy + k * bxz;
    const double y = j * by + k * byz;
    const double z = k * bz;
    const double rsq = x * x + y * y + z * z;
    return unit_voro.plane_intersects(x, y, z, rsq);
  }

/** Draws the periodic domain in gnuplot format.
   * \param[in] fp the file handle to write to. */
  void unitcell::draw_domain_gnuplot(FILE *fp) const {
    fprintf(fp, "0 0 0\n%g 0 0\n%g %g 0\n%g %g 0\n", bx, bx + bxy, by, bxy, by);
    fprintf(fp, "%g %g %g\n%g %g %g\n%g %g %g\n%g %g %g\n", bxy + bxz, by + byz, bz, bx + bxy + bxz, by + byz, bz, bx + bxz, byz, bz, bxz, byz, bz);
    fprintf(fp, "0 0 0\n%g %g 0\n\n%g %g %g\n%g %g %g\n\n", bxy, by, bxz, byz, bz, bxy + bxz, by + byz, bz);
    fprintf(fp, "%g 0 0\n%g %g %g\n\n%g %g 0\n%g %g %g\n\n", bx, bx + bxz, byz, bz, bx + bxy, by, bx + bxy + bxz, by + byz, bz);
  }

/** Draws the periodic domain in POV-Ray format.
   * \param[in] fp the file handle to write to. */
  void unitcell::draw_domain_pov(FILE *fp) const {
    fprintf(fp,
            "cylinder{0,0,0>,<%g,0,0>,rr}\n"
            "cylinder{<%g,%g,0>,<%g,%g,0>,rr}\n",
            bx, bxy, by, bx + bxy, by);
    fprintf(fp,
            "cylinder{<%g,%g,%g>,<%g,%g,%g>,rr}\n"
            "cylinder{<%g,%g,%g>,<%g,%g,%g>,rr}\n",
            bxz, byz, bz, bx + bxz, byz, bz, bxy + bxz, by + byz, bz, bx + bxy + bxz, by + byz, bz);
    fprintf(fp,
            "cylinder{<0,0,0>,<%g,%g,0>,rr}\n"
            "cylinder{<%g,0,0>,<%g,%g,0>,rr}\n",
            bxy, by, bx, bx + bxy, by);
    fprintf(fp,
            "cylinder{<%g,%g,%g>,<%g,%g,%g>,rr}\n"
            "cylinder{<%g,%g,%g>,<%g,%g,%g>,rr}\n",
            bxz, byz, bz, bxy + bxz, by + byz, bz, bx + bxz, byz, bz, bx + bxy + bxz, by + byz, bz);
    fprintf(fp,
            "cylinder{<0,0,0>,<%g,%g,%g>,rr}\n"
            "cylinder{<%g,0,0>,<%g,%g,%g>,rr}\n",
            bxz, byz, bz, bx, bx + bxz, byz, bz);
    fprintf(fp,
            "cylinder{<%g,%g,0>,<%g,%g,%g>,rr}\n"
            "cylinder{<%g,%g,0>,<%g,%g,%g>,rr}\n",
            bxy, by, bxy + bxz, by + byz, bz, bx + bxy, by, bx + bxy + bxz, by + byz, bz);
    fprintf(fp,
            "sphere{<0,0,0>,rr}\nsphere{<%g,0,0>,rr}\n"
            "sphere{<%g,%g,0>,rr}\nsphere{<%g,%g,0>,rr}\n",
            bx, bxy, by, bx + bxy, by);
    fprintf(fp,
            "sphere{<%g,%g,%g>,rr}\nsphere{<%g,%g,%g>,rr}\n"
            "sphere{<%g,%g,%g>,rr}\nsphere{<%g,%g,%g>,rr}\n",
            bxz, byz, bz, bx + bxz, byz, bz, bxy + bxz, by + byz, bz, bx + bxy + bxz, by + byz, bz);
  }

} // namespace voro
// Voro++, a 3D cell-based Voronoi library
//
// Author   : Chris H. Rycroft (LBL / UC Berkeley)
// Email    : chr@alum.mit.edu
// Date     : August 30th 2011

/** \file container_prd.cc
 * \brief Function implementations for the container_periodic_base and
 * related classes. */

namespace voro {

/** The class constructor sets up the geometry of container, initializing the
   * minimum and maximum coordinates in each direction, and setting whether each
   * direction is periodic or not. It divides the container into a rectangular
   * grid of blocks, and allocates memory for each of these for storing particle
   * positions and IDs.
   * \param[in] (bx_) The x coordinate of the first unit vector.
   * \param[in] (bxy_,by_) The x and y coordinates of the second unit vector.
   * \param[in] (bxz_,byz_,bz_) The x, y, and z coordinates of the third unit
   *                            vector.
   * \param[in] (nx_,ny_,nz_) the number of grid blocks in each of the three
   *                       coordinate directions.
   * \param[in] init_mem_ the initial memory allocation for each block.
   * \param[in] ps_ the number of floating point entries to store for each
   *                particle. */
  container_periodic_base::container_periodic_base(double bx_, double bxy_, double by_, double bxz_, double byz_, double bz_, int nx_, int ny_, int nz_, int init_mem_, int ps_) :
      unitcell(bx_, bxy_, by_, bxz_, byz_, bz_),
      voro_base(nx_, ny_, nz_, bx_ / nx_, by_ / ny_, bz_ / nz_),
      ey(int(max_uv_y * ysp + 1)),
      ez(int(max_uv_z * zsp + 1)),
      wy(ny + ey),
      wz(nz + ez),
      oy(ny + 2 * ey),
      oz(nz + 2 * ez),
      oxyz(nx * oy * oz),
      id(new int *[oxyz]),
      p(new double *[oxyz]),
      co(new int[oxyz]),
      mem(new int[oxyz]),
      img(new char[oxyz]),
      init_mem(init_mem_),
      ps(ps_) {
    int i;
    int j;
    int k;
    int l;

    // Clear the global arrays
    int *pp = co;
    while (pp < co + oxyz) {
      *(pp++) = 0;
    }
    pp = mem;
    while (pp < mem + oxyz) {
      *(pp++) = 0;
    }
    char *cp = img;
    while (cp < img + oxyz) {
      *(cp++) = 0;
    }

    // Set up memory for the blocks in the primary domain
    for (k = ez; k < wz; k++) {
      for (j = ey; j < wy; j++) {
        for (i = 0; i < nx; i++) {
          l = i + nx * (j + oy * k);
          mem[l] = init_mem;
          id[l] = new int[init_mem];
          p[l] = new double[ps * init_mem];
        }
      }
    }
  }

/** The container destructor frees the dynamically allocated memory. */
  container_periodic_base::~container_periodic_base() {
    for (int l = oxyz - 1; l >= 0; l--) {
      if (mem[l] > 0) {
        delete[] p[l];
        delete[] id[l];
      }
    }
    delete[] img;
    delete[] mem;
    delete[] co;
    delete[] id;
    delete[] p;
  }

/** The class constructor sets up the geometry of container.
   * \param[in] (bx_) The x coordinate of the first unit vector.
   * \param[in] (bxy_,by_) The x and y coordinates of the second unit vector.
   * \param[in] (bxz_,byz_,bz_) The x, y, and z coordinates of the third unit
   *                            vector.
   * \param[in] (nx_,ny_,nz_) the number of grid blocks in each of the three
   *			    coordinate directions.
   * \param[in] init_mem_ the initial memory allocation for each block. */
  container_periodic::container_periodic(double bx_, double bxy_, double by_, double bxz_, double byz_, double bz_, int nx_, int ny_, int nz_, int init_mem_) :
      container_periodic_base(bx_, bxy_, by_, bxz_, byz_, bz_, nx_, ny_, nz_, init_mem_, 3), vc(*this, 2 * nx_ + 1, 2 * ey + 1, 2 * ez + 1) {}

/** The class constructor sets up the geometry of container.
   * \param[in] (bx_) The x coordinate of the first unit vector.
   * \param[in] (bxy_,by_) The x and y coordinates of the second unit vector.
   * \param[in] (bxz_,byz_,bz_) The x, y, and z coordinates of the third unit
   *                            vector.
   * \param[in] (nx_,ny_,nz_) the number of grid blocks in each of the three
   *			    coordinate directions.
   * \param[in] init_mem_ the initial memory allocation for each block. */
  container_periodic_poly::container_periodic_poly(double bx_, double bxy_, double by_, double bxz_, double byz_, double bz_, int nx_, int ny_, int nz_, int init_mem_) :
      container_periodic_base(bx_, bxy_, by_, bxz_, byz_, bz_, nx_, ny_, nz_, init_mem_, 4), vc(*this, 2 * nx_ + 1, 2 * ey + 1, 2 * ez + 1) {
    ppr = p;
  }

/** Put a particle into the correct region of the container.
   * \param[in] n the numerical ID of the inserted particle.
   * \param[in] (x,y,z) the position vector of the inserted particle. */
  void container_periodic::put(int n, double x, double y, double z) {
    int ijk;
    put_locate_block(ijk, x, y, z);
    id[ijk][co[ijk]] = n;
    double *pp = p[ijk] + 3 * co[ijk]++;
    *(pp++) = x;
    *(pp++) = y;
    *pp = z;
  }

/** Put a particle into the correct region of the container.
   * \param[in] n the numerical ID of the inserted particle.
   * \param[in] (x,y,z) the position vector of the inserted particle.
   * \param[in] r the radius of the particle. */
  void container_periodic_poly::put(int n, double x, double y, double z, double r) {
    int ijk;
    put_locate_block(ijk, x, y, z);
    id[ijk][co[ijk]] = n;
    double *pp = p[ijk] + 4 * co[ijk]++;
    *(pp++) = x;
    *(pp++) = y;
    *(pp++) = z;
    *pp = r;
    if (max_radius < r) {
      max_radius = r;
    }
  }

/** Put a particle into the correct region of the container.
   * \param[in] n the numerical ID of the inserted particle.
   * \param[in] (x,y,z) the position vector of the inserted particle.
   * \param[out] (ai,aj,ak) the periodic image displacement that the particle is
   * 			  in, with (0,0,0) corresponding to the primary domain.
   */
  void container_periodic::put(int n, double x, double y, double z, int &ai, int &aj, int &ak) {
    int ijk;
    put_locate_block(ijk, x, y, z, ai, aj, ak);
    id[ijk][co[ijk]] = n;
    double *pp = p[ijk] + 3 * co[ijk]++;
    *(pp++) = x;
    *(pp++) = y;
    *pp = z;
  }

/** Put a particle into the correct region of the container.
   * \param[in] n the numerical ID of the inserted particle.
   * \param[in] (x,y,z) the position vector of the inserted particle.
   * \param[in] r the radius of the particle.
   * \param[out] (ai,aj,ak) the periodic image displacement that the particle is
   * 			  in, with (0,0,0) corresponding to the primary domain.
   */
  void container_periodic_poly::put(int n, double x, double y, double z, double r, int &ai, int &aj, int &ak) {
    int ijk;
    put_locate_block(ijk, x, y, z, ai, aj, ak);
    id[ijk][co[ijk]] = n;
    double *pp = p[ijk] + 4 * co[ijk]++;
    *(pp++) = x;
    *(pp++) = y;
    *(pp++) = z;
    *pp = r;
    if (max_radius < r) {
      max_radius = r;
    }
  }

/** Put a particle into the correct region of the container, also recording
   * into which region it was stored.
   * \param[in] vo the ordering class in which to record the region.
   * \param[in] n the numerical ID of the inserted particle.
   * \param[in] (x,y,z) the position vector of the inserted particle. */
  void container_periodic::put(particle_order &vo, int n, double x, double y, double z) {
    int ijk;
    put_locate_block(ijk, x, y, z);
    id[ijk][co[ijk]] = n;
    vo.add(ijk, co[ijk]);
    double *pp = p[ijk] + 3 * co[ijk]++;
    *(pp++) = x;
    *(pp++) = y;
    *pp = z;
  }

/** Put a particle into the correct region of the container, also recording
   * into which region it was stored.
   * \param[in] vo the ordering class in which to record the region.
   * \param[in] n the numerical ID of the inserted particle.
   * \param[in] (x,y,z) the position vector of the inserted particle.
   * \param[in] r the radius of the particle. */
  void container_periodic_poly::put(particle_order &vo, int n, double x, double y, double z, double r) {
    int ijk;
    put_locate_block(ijk, x, y, z);
    id[ijk][co[ijk]] = n;
    vo.add(ijk, co[ijk]);
    double *pp = p[ijk] + 4 * co[ijk]++;
    *(pp++) = x;
    *(pp++) = y;
    *(pp++) = z;
    *pp = r;
    if (max_radius < r) {
      max_radius = r;
    }
  }

/** Takes a particle position vector and computes the region index into which
   * it should be stored. If the container is periodic, then the routine also
   * maps the particle position to ensure it is in the primary domain. If the
   * container is not periodic, the routine bails out.
   * \param[out] ijk the region index.
   * \param[in,out] (x,y,z) the particle position, remapped into the primary
   *                        domain if necessary.
   * \return True if the particle can be successfully placed into the container,
   * false otherwise. */
  void container_periodic_base::put_locate_block(int &ijk, double &x, double &y, double &z) {
    // Remap particle in the z direction if necessary
    int k = step_int(z * zsp);
    if (k < 0 || k >= nz) {
      const int ak = step_div(k, nz);
      z -= ak * bz;
      y -= ak * byz;
      x -= ak * bxz;
      k -= ak * nz;
    }

    // Remap particle in the y direction if necessary
    int j = step_int(y * ysp);
    if (j < 0 || j >= ny) {
      const int aj = step_div(j, ny);
      y -= aj * by;
      x -= aj * bxy;
      j -= aj * ny;
    }

    // Remap particle in the x direction if necessary
    ijk = step_int(x * xsp);
    if (ijk < 0 || ijk >= nx) {
      const int ai = step_div(ijk, nx);
      x -= ai * bx;
      ijk -= ai * nx;
    }

    // Compute the block index and check memory allocation
    j += ey;
    k += ez;
    ijk += nx * (j + oy * k);
    if (co[ijk] == mem[ijk]) {
      add_particle_memory(ijk);
    }
  }

/** Takes a particle position vector and computes the region index into which
   * it should be stored. If the container is periodic, then the routine also
   * maps the particle position to ensure it is in the primary domain. If the
   * container is not periodic, the routine bails out.
   * \param[out] ijk the region index.
   * \param[in,out] (x,y,z) the particle position, remapped into the primary
   *                        domain if necessary.
   * \param[out] (ai,aj,ak) the periodic image displacement that the particle is
   *                        in, with (0,0,0) corresponding to the primary domain.
   * \return True if the particle can be successfully placed into the container,
   * false otherwise. */
  void container_periodic_base::put_locate_block(int &ijk, double &x, double &y, double &z, int &ai, int &aj, int &ak) {
    // Remap particle in the z direction if necessary
    int k = step_int(z * zsp);
    if (k < 0 || k >= nz) {
      ak = step_div(k, nz);
      z -= ak * bz;
      y -= ak * byz;
      x -= ak * bxz;
      k -= ak * nz;
    } else {
      ak = 0;
    }

    // Remap particle in the y direction if necessary
    int j = step_int(y * ysp);
    if (j < 0 || j >= ny) {
      aj = step_div(j, ny);
      y -= aj * by;
      x -= aj * bxy;
      j -= aj * ny;
    } else {
      aj = 0;
    }

    // Remap particle in the x direction if necessary
    ijk = step_int(x * xsp);
    if (ijk < 0 || ijk >= nx) {
      ai = step_div(ijk, nx);
      x -= ai * bx;
      ijk -= ai * nx;
    } else {
      ai = 0;
    }

    // Compute the block index and check memory allocation
    j += ey;
    k += ez;
    ijk += nx * (j + oy * k);
    if (co[ijk] == mem[ijk]) {
      add_particle_memory(ijk);
    }
  }

/** Takes a position vector and remaps it into the primary domain.
   * \param[out] (ai,aj,ak) the periodic image displacement that the vector is in,
   *                        with (0,0,0) corresponding to the primary domain.
   * \param[out] (ci,cj,ck) the index of the block that the position vector is
   *                        within, once it has been remapped.
   * \param[in,out] (x,y,z) the position vector to consider, which is remapped
   *                        into the primary domain during the routine.
   * \param[out] ijk the block index that the vector is within. */
  inline void container_periodic_base::remap(int &ai, int &aj, int &ak, int &ci, int &cj, int &ck, double &x, double &y, double &z, int &ijk) {
    // Remap particle in the z direction if necessary
    ck = step_int(z * zsp);
    if (ck < 0 || ck >= nz) {
      ak = step_div(ck, nz);
      z -= ak * bz;
      y -= ak * byz;
      x -= ak * bxz;
      ck -= ak * nz;
    } else {
      ak = 0;
    }

    // Remap particle in the y direction if necessary
    cj = step_int(y * ysp);
    if (cj < 0 || cj >= ny) {
      aj = step_div(cj, ny);
      y -= aj * by;
      x -= aj * bxy;
      cj -= aj * ny;
    } else {
      aj = 0;
    }

    // Remap particle in the x direction if necessary
    ci = step_int(x * xsp);
    if (ci < 0 || ci >= nx) {
      ai = step_div(ci, nx);
      x -= ai * bx;
      ci -= ai * nx;
    } else {
      ai = 0;
    }

    cj += ey;
    ck += ez;
    ijk = ci + nx * (cj + oy * ck);
  }

/** Takes a vector and finds the particle whose Voronoi cell contains that
   * vector. This is equivalent to finding the particle which is nearest to the
   * vector.
   * \param[in] (x,y,z) the vector to test.
   * \param[out] (rx,ry,rz) the position of the particle whose Voronoi cell
   *                        contains the vector. This may point to a particle in
   *                        a periodic image of the primary domain.
   * \param[out] pid the ID of the particle.
   * \return True if a particle was found. If the container has no particles,
   * then the search will not find a Voronoi cell and false is returned. */
  bool container_periodic::find_voronoi_cell(double x, double y, double z, double &rx, double &ry, double &rz, int &pid) {
    int ai;
    int aj;
    int ak;
    int ci;
    int cj;
    int ck;
    int ijk;
    particle_record w;
    double mrs;

    // Remap the vector into the primary domain and then search for the
    // Voronoi cell that it is within
    remap(ai, aj, ak, ci, cj, ck, x, y, z, ijk);
    vc.find_voronoi_cell(x, y, z, ci, cj, ck, ijk, w, mrs);

    if (w.ijk != -1) {
      // Assemble the position vector of the particle to be returned,
      // applying a periodic remapping if necessary
      ci += w.di;
      if (ci < 0 || ci >= nx) {
        ai += step_div(ci, nx);
      }
      rx = p[w.ijk][3 * w.l] + ak * bxz + aj * bxy + ai * bx;
      ry = p[w.ijk][3 * w.l + 1] + ak * byz + aj * by;
      rz = p[w.ijk][3 * w.l + 2] + ak * bz;
      pid = id[w.ijk][w.l];
      return true;
    }
    return false;
  }

/** Takes a vector and finds the particle whose Voronoi cell contains that
   * vector. Additional wall classes are not considered by this routine.
   * \param[in] (x,y,z) the vector to test.
   * \param[out] (rx,ry,rz) the position of the particle whose Voronoi cell
   *                        contains the vector. If the container is periodic,
   *                        this may point to a particle in a periodic image of
   *                        the primary domain.
   * \param[out] pid the ID of the particle.
   * \return True if a particle was found. If the container has no particles,
   * then the search will not find a Voronoi cell and false is returned. */
  bool container_periodic_poly::find_voronoi_cell(double x, double y, double z, double &rx, double &ry, double &rz, int &pid) {
    int ai;
    int aj;
    int ak;
    int ci;
    int cj;
    int ck;
    int ijk;
    particle_record w;
    double mrs;

    // Remap the vector into the primary domain and then search for the
    // Voronoi cell that it is within
    remap(ai, aj, ak, ci, cj, ck, x, y, z, ijk);
    vc.find_voronoi_cell(x, y, z, ci, cj, ck, ijk, w, mrs);

    if (w.ijk != -1) {
      // Assemble the position vector of the particle to be returned,
      // applying a periodic remapping if necessary
      ci += w.di;
      if (ci < 0 || ci >= nx) {
        ai += step_div(ci, nx);
      }
      rx = p[w.ijk][4 * w.l] + ak * bxz + aj * bxy + ai * bx;
      ry = p[w.ijk][4 * w.l + 1] + ak * byz + aj * by;
      rz = p[w.ijk][4 * w.l + 2] + ak * bz;
      pid = id[w.ijk][w.l];
      return true;
    }
    return false;
  }

/** Increase memory for a particular region.
   * \param[in] i the index of the region to reallocate. */
  void container_periodic_base::add_particle_memory(int i) const {
    // Handle the case when no memory has been allocated for this block
    if (mem[i] == 0) {
      mem[i] = init_mem;
      id[i] = new int[init_mem];
      p[i] = new double[ps * init_mem];
      return;
    }

    // Otherwise, double the memory allocation for this block. Carry out a
    // check on the memory allocation size, and print a status message if
    // requested.
    int l;
    const int nmem(mem[i] << 1);
    if (nmem > max_particle_memory) {
      voro_fatal_error("Absolute maximum memory allocation exceeded", VOROPP_MEMORY_ERROR);
    }
#if VOROPP_VERBOSE >= 3
    fprintf(stderr, "Particle memory in region %d scaled up to %d\n", i, nmem);
#endif

    // Allocate new memory and copy in the contents of the old arrays
    int *idp = new int[nmem];
    for (l = 0; l < co[i]; l++) {
      idp[l] = id[i][l];
    }
    double *pp = new double[ps * nmem];
    for (l = 0; l < ps * co[i]; l++) {
      pp[l] = p[i][l];
    }

    // Update pointers and delete old arrays
    mem[i] = nmem;
    delete[] id[i];
    id[i] = idp;
    delete[] p[i];
    p[i] = pp;
  }

/** Import a list of particles from an open file stream into the container.
   * Entries of four numbers (Particle ID, x position, y position, z position)
   * are searched for. If the file cannot be successfully read, then the routine
   * causes a fatal error.
   * \param[in] fp the file handle to read from. */
  void container_periodic::import(FILE *fp) {
    int i;
    int j;
    double x;
    double y;
    double z;
    while ((j = fscanf(fp, "%d %lg %lg %lg", &i, &x, &y, &z)) == 4) {
      put(i, x, y, z);
    }
    if (j != EOF) {
      voro_fatal_error("File import error", VOROPP_FILE_ERROR);
    }
  }

/** Import a list of particles from an open file stream, also storing the order
   * of that the particles are read. Entries of four numbers (Particle ID, x
   * position, y position, z position) are searched for. If the file cannot be
   * successfully read, then the routine causes a fatal error.
   * \param[in,out] vo a reference to an ordering class to use.
   * \param[in] fp the file handle to read from. */
  void container_periodic::import(particle_order &vo, FILE *fp) {
    int i;
    int j;
    double x;
    double y;
    double z;
    while ((j = fscanf(fp, "%d %lg %lg %lg", &i, &x, &y, &z)) == 4) {
      put(vo, i, x, y, z);
    }
    if (j != EOF) {
      voro_fatal_error("File import error", VOROPP_FILE_ERROR);
    }
  }

/** Import a list of particles from an open file stream into the container.
   * Entries of five numbers (Particle ID, x position, y position, z position,
   * radius) are searched for. If the file cannot be successfully read, then the
   * routine causes a fatal error.
   * \param[in] fp the file handle to read from. */
  void container_periodic_poly::import(FILE *fp) {
    int i;
    int j;
    double x;
    double y;
    double z;
    double r;
    while ((j = fscanf(fp, "%d %lg %lg %lg %lg", &i, &x, &y, &z, &r)) == 5) {
      put(i, x, y, z, r);
    }
    if (j != EOF) {
      voro_fatal_error("File import error", VOROPP_FILE_ERROR);
    }
  }

/** Import a list of particles from an open file stream, also storing the order
   * of that the particles are read. Entries of four numbers (Particle ID, x
   * position, y position, z position, radius) are searched for. If the file
   * cannot be successfully read, then the routine causes a fatal error.
   * \param[in,out] vo a reference to an ordering class to use.
   * \param[in] fp the file handle to read from. */
  void container_periodic_poly::import(particle_order &vo, FILE *fp) {
    int i;
    int j;
    double x;
    double y;
    double z;
    double r;
    while ((j = fscanf(fp, "%d %lg %lg %lg %lg", &i, &x, &y, &z, &r)) == 5) {
      put(vo, i, x, y, z, r);
    }
    if (j != EOF) {
      voro_fatal_error("File import error", VOROPP_FILE_ERROR);
    }
  }

/** Outputs the a list of all the container regions along with the number of
   * particles stored within each. */
  void container_periodic_base::region_count() {
    int i;
    int j;
    int k;
    int *cop = co;
    for (k = 0; k < nz; k++) {
      for (j = 0; j < ny; j++) {
        for (i = 0; i < nx; i++) {
          printf("Region (%d,%d,%d): %d particles\n", i, j, k, *(cop++));
        }
      }
    }
  }

/** Clears a container of particles. */
  void container_periodic::clear() {
    for (int *cop = co; cop < co + nxyz; cop++) {
      *cop = 0;
    }
  }

/** Clears a container of particles, also clearing resetting the maximum radius
   * to zero. */
  void container_periodic_poly::clear() {
    for (int *cop = co; cop < co + nxyz; cop++) {
      *cop = 0;
    }
    max_radius = 0;
  }

/** Computes all the Voronoi cells and saves customized information about them.
   * \param[in] format the custom output string to use.
   * \param[in] fp a file handle to write to. */
  void container_periodic::print_custom(const char *format, FILE *fp) {
    c_loop_all_periodic vl(*this);
    print_custom(vl, format, fp);
  }

/** Computes all the Voronoi cells and saves customized
   * information about them.
   * \param[in] format the custom output string to use.
   * \param[in] fp a file handle to write to. */
  void container_periodic_poly::print_custom(const char *format, FILE *fp) {
    c_loop_all_periodic vl(*this);
    print_custom(vl, format, fp);
  }

/** Computes all the Voronoi cells and saves customized information about them.
   * \param[in] format the custom output string to use.
   * \param[in] filename the name of the file to write to. */
  void container_periodic::print_custom(const char *format, const char *filename) {
    FILE *fp = safe_fopen(filename, "w");
    print_custom(format, fp);
    fclose(fp);
  }

/** Computes all the Voronoi cells and saves customized
   * information about them
   * \param[in] format the custom output string to use.
   * \param[in] filename the name of the file to write to. */
  void container_periodic_poly::print_custom(const char *format, const char *filename) {
    FILE *fp = safe_fopen(filename, "w");
    print_custom(format, fp);
    fclose(fp);
  }

/** Computes all of the Voronoi cells in the container, but does nothing
   * with the output. It is useful for measuring the pure computation time
   * of the Voronoi algorithm, without any additional calculations such as
   * volume evaluation or cell output. */
  void container_periodic::compute_all_cells() {
    voronoicell c;
    c_loop_all_periodic vl(*this);
    if (vl.start()) {
      do {
        compute_cell(c, vl);
      } while (vl.inc());
    }
  }

/** Computes all of the Voronoi cells in the container, but does nothing
   * with the output. It is useful for measuring the pure computation time
   * of the Voronoi algorithm, without any additional calculations such as
   * volume evaluation or cell output. */
  void container_periodic_poly::compute_all_cells() {
    voronoicell c;
    c_loop_all_periodic vl(*this);
    if (vl.start()) {
      do {
        compute_cell(c, vl);
      } while (vl.inc());
    }
  }

/** Calculates all of the Voronoi cells and sums their volumes. In most cases
   * without walls, the sum of the Voronoi cell volumes should equal the volume
   * of the container to numerical precision.
   * \return The sum of all of the computed Voronoi volumes. */
  double container_periodic::sum_cell_volumes() {
    voronoicell c;
    double vol = 0;
    c_loop_all_periodic vl(*this);
    if (vl.start()) {
      do {
        if (compute_cell(c, vl)) {
          vol += c.volume();
        }
      } while (vl.inc());
    }
    return vol;
  }

/** Calculates all of the Voronoi cells and sums their volumes. In most cases
   * without walls, the sum of the Voronoi cell volumes should equal the volume
   * of the container to numerical precision.
   * \return The sum of all of the computed Voronoi volumes. */
  double container_periodic_poly::sum_cell_volumes() {
    voronoicell c;
    double vol = 0;
    c_loop_all_periodic vl(*this);
    if (vl.start()) {
      do {
        if (compute_cell(c, vl)) {
          vol += c.volume();
        }
      } while (vl.inc());
    }
    return vol;
  }

/** This routine creates all periodic images of the particles. It is meant for
   * diagnostic purposes only, since usually periodic images are dynamically
   * created in when they are referenced. */
  void container_periodic_base::create_all_images() {
    int i;
    int j;
    int k;
    for (k = 0; k < oz; k++) {
      for (j = 0; j < oy; j++) {
        for (i = 0; i < nx; i++) {
          create_periodic_image(i, j, k);
        }
      }
    }
  }

/** Checks that the particles within each block lie within that block's bounds.
   * This is useful for diagnosing problems with periodic image computation. */
  void container_periodic_base::check_compartmentalized() {
    int c;
    int l;
    int i;
    int j;
    int k;
    double mix;
    double miy;
    double miz;
    double max;
    double may;
    double maz;
    double *pp;
    for (k = l = 0; k < oz; k++) {
      for (j = 0; j < oy; j++) {
        for (i = 0; i < nx; i++, l++) {
          if (mem[l] > 0) {
            // Compute the block's bounds, adding in a small tolerance
            mix = i * boxx - tolerance;
            max = mix + boxx + tolerance;
            miy = (j - ey) * boxy - tolerance;
            may = miy + boxy + tolerance;
            miz = (k - ez) * boxz - tolerance;
            maz = miz + boxz + tolerance;

            // Print entries for any particles that lie outside the block's
            // bounds
            for (pp = p[l], c = 0; c < co[l]; c++, pp += ps) {
              if (*pp < mix || *pp > max || pp[1] < miy || pp[1] > may || pp[2] < miz || pp[2] > maz) {
                printf("%d %d %d %d %f %f %f %f %f %f %f %f %f\n", id[l][c], i, j, k, *pp, pp[1], pp[2], mix, max, miy, may, miz, maz);
              }
            }
          }
        }
      }
    }
  }

/** Creates particles within an image block that is aligned with the primary
   * domain in the z axis. In this case, the image block may be comprised of
   * particles from two primary blocks. The routine considers these two primary
   * blocks, and adds the needed particles to the image. The remaining particles
   * from the primary blocks are also filled into the neighboring images.
   * \param[in] (di,dj,dk) the index of the block to consider. The z index must
   *			 satisfy ez<=dk<wz. */
  void container_periodic_base::create_side_image(int di, int dj, int dk) {
    int l;
    const int dijk = di + nx * (dj + oy * dk);
    int odijk;
    const int ima = step_div(dj - ey, ny);
    const int qua = di + step_int(-ima * bxy * xsp);
    const int quadiv = step_div(qua, nx);
    const int fi = qua - quadiv * nx;
    int fijk = fi + nx * (dj - ima * ny + oy * dk);
    double dis = ima * bxy + quadiv * bx;
    double switchx = di * boxx - ima * bxy - quadiv * bx;
    double adis;

    // Left image computation
    if ((img[dijk] & 1) == 0) {
      if (di > 0) {
        odijk = dijk - 1;
        adis = dis;
      } else {
        odijk = dijk + nx - 1;
        adis = dis + bx;
      }
      img[odijk] |= 2;
      for (l = 0; l < co[fijk]; l++) {
        if (p[fijk][ps * l] > switchx) {
          put_image(dijk, fijk, l, dis, by * ima, 0);
        } else {
          put_image(odijk, fijk, l, adis, by * ima, 0);
        }
      }
    }

    // Right image computation
    if ((img[dijk] & 2) == 0) {
      if (fi == nx - 1) {
        fijk += 1 - nx;
        switchx += (1 - nx) * boxx;
        dis += bx;
      } else {
        fijk++;
        switchx += boxx;
      }
      if (di == nx - 1) {
        odijk = dijk - nx + 1;
        adis = dis - bx;
      } else {
        odijk = dijk + 1;
        adis = dis;
      }
      img[odijk] |= 1;
      for (l = 0; l < co[fijk]; l++) {
        if (p[fijk][ps * l] < switchx) {
          put_image(dijk, fijk, l, dis, by * ima, 0);
        } else {
          put_image(odijk, fijk, l, adis, by * ima, 0);
        }
      }
    }

    // All contributions to the block now added, so set both two bits of
    // the image information
    img[dijk] = 3;
  }

/** Creates particles within an image block that is not aligned with the
   * primary domain in the z axis. In this case, the image block may be comprised
   * of particles from four primary blocks. The routine considers these four
   * primary blocks, and adds the needed particles to the image. The remaining
   * particles from the primary blocks are also filled into the neighboring
   * images.
   * \param[in] (di,dj,dk) the index of the block to consider. The z index must
   *			 satisfy dk<ez or dk>=wz. */
  void container_periodic_base::create_vertical_image(int di, int dj, int dk) {
    int l;
    const int dijk = di + nx * (dj + oy * dk);
    int dijkl;
    int dijkr;
    const int ima = step_div(dk - ez, nz);
    const int qj = dj + step_int(-ima * byz * ysp);
    const int qjdiv = step_div(qj - ey, ny);
    int qi = di + step_int((-ima * bxz - qjdiv * bxy) * xsp);
    int qidiv = step_div(qi, nx);
    int fi = qi - qidiv * nx;
    const int fj = qj - qjdiv * ny;
    int fijk = fi + nx * (fj + oy * (dk - ima * nz));
    int fijk2;
    double disy = ima * byz + qjdiv * by;
    double switchy = (dj - ey) * boxy - ima * byz - qjdiv * by;
    double disx = ima * bxz + qjdiv * bxy + qidiv * bx;
    double switchx = di * boxx - ima * bxz - qjdiv * bxy - qidiv * bx;
    double switchx2;
    double disxl;
    double disxr;
    double disx2;
    double disxr2;

    if (di == 0) {
      dijkl = dijk + nx - 1;
      disxl = disx + bx;
    } else {
      dijkl = dijk - 1;
      disxl = disx;
    }

    if (di == nx - 1) {
      dijkr = dijk - nx + 1;
      disxr = disx - bx;
    } else {
      dijkr = dijk + 1;
      disxr = disx;
    }

    // Down-left image computation
    bool y_exist = dj != 0;
    if ((img[dijk] & 1) == 0) {
      img[dijkl] |= 2;
      if (y_exist) {
        img[dijkl - nx] |= 8;
        img[dijk - nx] |= 4;
      }
      for (l = 0; l < co[fijk]; l++) {
        if (p[fijk][ps * l + 1] > switchy) {
          if (p[fijk][ps * l] > switchx) {
            put_image(dijk, fijk, l, disx, disy, bz * ima);
          } else {
            put_image(dijkl, fijk, l, disxl, disy, bz * ima);
          }
        } else {
          if (!y_exist) {
            continue;
          }
          if (p[fijk][ps * l] > switchx) {
            put_image(dijk - nx, fijk, l, disx, disy, bz * ima);
          } else {
            put_image(dijkl - nx, fijk, l, disxl, disy, bz * ima);
          }
        }
      }
    }

    // Down-right image computation
    if ((img[dijk] & 2) == 0) {
      if (fi == nx - 1) {
        fijk2 = fijk + 1 - nx;
        switchx2 = switchx + (1 - nx) * boxx;
        disx2 = disx + bx;
        disxr2 = disxr + bx;
      } else {
        fijk2 = fijk + 1;
        switchx2 = switchx + boxx;
        disx2 = disx;
        disxr2 = disxr;
      }
      img[dijkr] |= 1;
      if (y_exist) {
        img[dijkr - nx] |= 4;
        img[dijk - nx] |= 8;
      }
      for (l = 0; l < co[fijk2]; l++) {
        if (p[fijk2][ps * l + 1] > switchy) {
          if (p[fijk2][ps * l] > switchx2) {
            put_image(dijkr, fijk2, l, disxr2, disy, bz * ima);
          } else {
            put_image(dijk, fijk2, l, disx2, disy, bz * ima);
          }
        } else {
          if (!y_exist) {
            continue;
          }
          if (p[fijk2][ps * l] > switchx2) {
            put_image(dijkr - nx, fijk2, l, disxr2, disy, bz * ima);
          } else {
            put_image(dijk - nx, fijk2, l, disx2, disy, bz * ima);
          }
        }
      }
    }

    // Recomputation of some intermediate quantities for boundary cases
    if (fj == wy - 1) {
      fijk += nx * (1 - ny) - fi;
      switchy += (1 - ny) * boxy;
      disy += by;
      qi = di + step_int(-(ima * bxz + (qjdiv + 1) * bxy) * xsp);
      const int dqidiv = step_div(qi, nx) - qidiv;
      qidiv += dqidiv;
      fi = qi - qidiv * nx;
      fijk += fi;
      disx += bxy + bx * dqidiv;
      disxl += bxy + bx * dqidiv;
      disxr += bxy + bx * dqidiv;
      switchx -= bxy + bx * dqidiv;
    } else {
      fijk += nx;
      switchy += boxy;
    }

    // Up-left image computation
    y_exist = dj != oy - 1;
    if ((img[dijk] & 4) == 0) {
      img[dijkl] |= 8;
      if (y_exist) {
        img[dijkl + nx] |= 2;
        img[dijk + nx] |= 1;
      }
      for (l = 0; l < co[fijk]; l++) {
        if (p[fijk][ps * l + 1] > switchy) {
          if (!y_exist) {
            continue;
          }
          if (p[fijk][ps * l] > switchx) {
            put_image(dijk + nx, fijk, l, disx, disy, bz * ima);
          } else {
            put_image(dijkl + nx, fijk, l, disxl, disy, bz * ima);
          }
        } else {
          if (p[fijk][ps * l] > switchx) {
            put_image(dijk, fijk, l, disx, disy, bz * ima);
          } else {
            put_image(dijkl, fijk, l, disxl, disy, bz * ima);
          }
        }
      }
    }

    // Up-right image computation
    if ((img[dijk] & 8) == 0) {
      if (fi == nx - 1) {
        fijk2 = fijk + 1 - nx;
        switchx2 = switchx + (1 - nx) * boxx;
        disx2 = disx + bx;
        disxr2 = disxr + bx;
      } else {
        fijk2 = fijk + 1;
        switchx2 = switchx + boxx;
        disx2 = disx;
        disxr2 = disxr;
      }
      img[dijkr] |= 4;
      if (y_exist) {
        img[dijkr + nx] |= 1;
        img[dijk + nx] |= 2;
      }
      for (l = 0; l < co[fijk2]; l++) {
        if (p[fijk2][ps * l + 1] > switchy) {
          if (!y_exist) {
            continue;
          }
          if (p[fijk2][ps * l] > switchx2) {
            put_image(dijkr + nx, fijk2, l, disxr2, disy, bz * ima);
          } else {
            put_image(dijk + nx, fijk2, l, disx2, disy, bz * ima);
          }
        } else {
          if (p[fijk2][ps * l] > switchx2) {
            put_image(dijkr, fijk2, l, disxr2, disy, bz * ima);
          } else {
            put_image(dijk, fijk2, l, disx2, disy, bz * ima);
          }
        }
      }
    }

    // All contributions to the block now added, so set all four bits of
    // the image information
    img[dijk] = 15;
  }

/** Copies a particle position from the primary domain into an image block.
   * \param[in] reg the block index within the primary domain that the particle
   *                is within.
   * \param[in] fijk the index of the image block.
   * \param[in] l the index of the particle entry within the primary block.
   * \param[in] (dx,dy,dz) the displacement vector to add to the particle. */
  void container_periodic_base::put_image(int reg, int fijk, int l, double dx, double dy, double dz) {
    if (co[reg] == mem[reg]) {
      add_particle_memory(reg);
    }
    double *p1 = p[reg] + ps * co[reg];
    double *p2 = p[fijk] + ps * l;
    *(p1++) = *(p2++) + dx;
    *(p1++) = *(p2++) + dy;
    *p1 = *p2 + dz;
    if (ps == 4) {
      *(++p1) = *(++p2);
    }
    id[reg][co[reg]++] = id[fijk][l];
  }

} // namespace voro
// Voro++, a 3D cell-based Voronoi library
//
// Author   : Chris H. Rycroft (LBL / UC Berkeley)
// Email    : chr@alum.mit.edu
// Date     : August 30th 2011

/** \file pre_container.cc
 * \brief Function implementations for the pre_container and related classes.
 */

#include <cmath>

namespace voro {

/** The class constructor sets up the geometry of container, initializing the
   * minimum and maximum coordinates in each direction. It allocates an initial
   * chunk into which to store particle information.
   * \param[in] (ax_,bx_) the minimum and maximum x coordinates.
   * \param[in] (ay_,by_) the minimum and maximum y coordinates.
   * \param[in] (az_,bz_) the minimum and maximum z coordinates.
   * \param[in] (xperiodic_,yperiodic_,zperiodic_ ) flags setting whether the
   *                                                container is periodic in each
   *                                                coordinate direction.
   * \param[in] ps_ the number of floating point entries to store for each
   *                particle. */
  pre_container_base::pre_container_base(double ax_, double bx_, double ay_, double by_, double az_, double bz_, bool xperiodic_, bool yperiodic_, bool zperiodic_, int ps_) :
      ax(ax_),
      bx(bx_),
      ay(ay_),
      by(by_),
      az(az_),
      bz(bz_),
      xperiodic(xperiodic_),
      yperiodic(yperiodic_),
      zperiodic(zperiodic_),
      ps(ps_),
      index_sz(init_chunk_size),
      pre_id(new int *[index_sz]),
      end_id(pre_id),
      pre_p(new double *[index_sz]),
      end_p(pre_p) {
    ch_id = *end_id = new int[pre_container_chunk_size];
    l_id = end_id + index_sz;
    e_id = ch_id + pre_container_chunk_size;
    ch_p = *end_p = new double[ps * pre_container_chunk_size];
  }

/** The destructor frees the dynamically allocated memory. */
  pre_container_base::~pre_container_base() {
    delete[] *end_p;
    delete[] *end_id;
    while (end_id != pre_id) {
      end_p--;
      delete[] *end_p;
      end_id--;
      delete[] *end_id;
    }
    delete[] pre_p;
    delete[] pre_id;
  }

/** Makes a guess at the optimal grid of blocks to use, computing in
   * a way that
   * \param[out] (nx,ny,nz) the number of blocks to use. */
  void pre_container_base::guess_optimal(int &nx, int &ny, int &nz) {
    const double dx = bx - ax;
    const double dy = by - ay;
    const double dz = bz - az;
    const double ilscale = pow(total_particles() / (optimal_particles * dx * dy * dz), 1 / 3.0);
    nx = int(dx * ilscale + 1);
    ny = int(dy * ilscale + 1);
    nz = int(dz * ilscale + 1);
  }

/** Stores a particle ID and position, allocating a new memory chunk if
   * necessary. For coordinate directions in which the container is not periodic,
   * the routine checks to make sure that the particle is within the container
   * bounds. If the particle is out of bounds, it is not stored.
   * \param[in] n the numerical ID of the inserted particle.
   * \param[in] (x,y,z) the position vector of the inserted particle. */
  void pre_container::put(int n, double x, double y, double z) {
    if ((xperiodic || (x >= ax && x <= bx)) && (yperiodic || (y >= ay && y <= by)) && (zperiodic || (z >= az && z <= bz))) {
      if (ch_id == e_id) {
        new_chunk();
      }
      *(ch_id++) = n;
      *(ch_p++) = x;
      *(ch_p++) = y;
      *(ch_p++) = z;
    }
#if VOROPP_REPORT_OUT_OF_BOUNDS == 1
    else
      fprintf(stderr, "Out of bounds: (x,y,z)=(%g,%g,%g)\n", x, y, z);
#endif
  }

/** Stores a particle ID and position, allocating a new memory chunk if necessary.
   * \param[in] n the numerical ID of the inserted particle.
   * \param[in] (x,y,z) the position vector of the inserted particle.
   * \param[in] r the radius of the particle. */
  void pre_container_poly::put(int n, double x, double y, double z, double r) {
    if ((xperiodic || (x >= ax && x <= bx)) && (yperiodic || (y >= ay && y <= by)) && (zperiodic || (z >= az && z <= bz))) {
      if (ch_id == e_id) {
        new_chunk();
      }
      *(ch_id++) = n;
      *(ch_p++) = x;
      *(ch_p++) = y;
      *(ch_p++) = z;
      *(ch_p++) = r;
    }
#if VOROPP_REPORT_OUT_OF_BOUNDS == 1
    else
      fprintf(stderr, "Out of bounds: (x,y,z)=(%g,%g,%g)\n", x, y, z);
#endif
  }

/** Transfers the particles stored within the class to a container class.
   * \param[in] con the container class to transfer to. */
  void pre_container::setup(container &con) {
    int **c_id = pre_id;
    int *idp;
    int *ide;
    int n;
    double **c_p = pre_p;
    double *pp;
    double x;
    double y;
    double z;
    while (c_id < end_id) {
      idp = *(c_id++);
      ide = idp + pre_container_chunk_size;
      pp = *(c_p++);
      while (idp < ide) {
        n = *(idp++);
        x = *(pp++);
        y = *(pp++);
        z = *(pp++);
        con.put(n, x, y, z);
      }
    }
    idp = *c_id;
    pp = *c_p;
    while (idp < ch_id) {
      n = *(idp++);
      x = *(pp++);
      y = *(pp++);
      z = *(pp++);
      con.put(n, x, y, z);
    }
  }

/** Transfers the particles stored within the class to a container_poly class.
   * \param[in] con the container_poly class to transfer to. */
  void pre_container_poly::setup(container_poly &con) {
    int **c_id = pre_id;
    int *idp;
    int *ide;
    int n;
    double **c_p = pre_p;
    double *pp;
    double x;
    double y;
    double z;
    double r;
    while (c_id < end_id) {
      idp = *(c_id++);
      ide = idp + pre_container_chunk_size;
      pp = *(c_p++);
      while (idp < ide) {
        n = *(idp++);
        x = *(pp++);
        y = *(pp++);
        z = *(pp++);
        r = *(pp++);
        con.put(n, x, y, z, r);
      }
    }
    idp = *c_id;
    pp = *c_p;
    while (idp < ch_id) {
      n = *(idp++);
      x = *(pp++);
      y = *(pp++);
      z = *(pp++);
      r = *(pp++);
      con.put(n, x, y, z, r);
    }
  }

/** Transfers the particles stored within the class to a container class, also
   * recording the order in which particles were stored.
   * \param[in] vo the ordering class to use.
   * \param[in] con the container class to transfer to. */
  void pre_container::setup(particle_order &vo, container &con) {
    int **c_id = pre_id;
    int *idp;
    int *ide;
    int n;
    double **c_p = pre_p;
    double *pp;
    double x;
    double y;
    double z;
    while (c_id < end_id) {
      idp = *(c_id++);
      ide = idp + pre_container_chunk_size;
      pp = *(c_p++);
      while (idp < ide) {
        n = *(idp++);
        x = *(pp++);
        y = *(pp++);
        z = *(pp++);
        con.put(vo, n, x, y, z);
      }
    }
    idp = *c_id;
    pp = *c_p;
    while (idp < ch_id) {
      n = *(idp++);
      x = *(pp++);
      y = *(pp++);
      z = *(pp++);
      con.put(vo, n, x, y, z);
    }
  }

/** Transfers the particles stored within the class to a container_poly class,
   * also recording the order in which particles were stored.
   * \param[in] vo the ordering class to use.
   * \param[in] con the container_poly class to transfer to. */
  void pre_container_poly::setup(particle_order &vo, container_poly &con) {
    int **c_id = pre_id;
    int *idp;
    int *ide;
    int n;
    double **c_p = pre_p;
    double *pp;
    double x;
    double y;
    double z;
    double r;
    while (c_id < end_id) {
      idp = *(c_id++);
      ide = idp + pre_container_chunk_size;
      pp = *(c_p++);
      while (idp < ide) {
        n = *(idp++);
        x = *(pp++);
        y = *(pp++);
        z = *(pp++);
        r = *(pp++);
        con.put(vo, n, x, y, z, r);
      }
    }
    idp = *c_id;
    pp = *c_p;
    while (idp < ch_id) {
      n = *(idp++);
      x = *(pp++);
      y = *(pp++);
      z = *(pp++);
      r = *(pp++);
      con.put(vo, n, x, y, z, r);
    }
  }

/** Import a list of particles from an open file stream into the container.
   * Entries of four numbers (Particle ID, x position, y position, z position)
   * are searched for. If the file cannot be successfully read, then the routine
   * causes a fatal error.
   * \param[in] fp the file handle to read from. */
  void pre_container::import(FILE *fp) {
    int i;
    int j;
    double x;
    double y;
    double z;
    while ((j = fscanf(fp, "%d %lg %lg %lg", &i, &x, &y, &z)) == 4) {
      put(i, x, y, z);
    }
    if (j != EOF) {
      voro_fatal_error("File import error", VOROPP_FILE_ERROR);
    }
  }

/** Import a list of particles from an open file stream, also storing the order
   * of that the particles are read. Entries of four numbers (Particle ID, x
   * position, y position, z position) are searched for. If the file cannot be
   * successfully read, then the routine causes a fatal error.
   * \param[in] fp the file handle to read from. */
  void pre_container_poly::import(FILE *fp) {
    int i;
    int j;
    double x;
    double y;
    double z;
    double r;
    while ((j = fscanf(fp, "%d %lg %lg %lg %lg", &i, &x, &y, &z, &r)) == 5) {
      put(i, x, y, z, r);
    }
    if (j != EOF) {
      voro_fatal_error("File import error", VOROPP_FILE_ERROR);
    }
  }

/** Allocates a new chunk of memory for storing particles. */
  void pre_container_base::new_chunk() {
    end_id++;
    end_p++;
    if (end_id == l_id) {
      extend_chunk_index();
    }
    ch_id = *end_id = new int[pre_container_chunk_size];
    e_id = ch_id + pre_container_chunk_size;
    ch_p = *end_p = new double[ps * pre_container_chunk_size];
  }

/** Extends the index of chunks. */
  void pre_container_base::extend_chunk_index() {
    index_sz <<= 1;
    if (index_sz > max_chunk_size) {
      voro_fatal_error("Absolute memory limit on chunk index reached", VOROPP_MEMORY_ERROR);
    }
#if VOROPP_VERBOSE >= 2
    fprintf(stderr, "Pre-container chunk index scaled up to %d\n", index_sz);
#endif
    int **n_id = new int *[index_sz];
    int **p_id = n_id;
    int **c_id = pre_id;
    double **n_p = new double *[index_sz];
    double **p_p = n_p;
    double **c_p = pre_p;
    while (c_id < end_id) {
      *(p_id++) = *(c_id++);
      *(p_p++) = *(c_p++);
    }
    delete[] pre_id;
    pre_id = n_id;
    end_id = p_id;
    l_id = pre_id + index_sz;
    delete[] pre_p;
    pre_p = n_p;
    end_p = p_p;
  }

} // namespace voro
// Voro++, a 3D cell-based Voronoi library
//
// Author   : Chris H. Rycroft (LBL / UC Berkeley)
// Email    : chr@alum.mit.edu
// Date     : August 30th 2011

/** \file wall.cc
 * \brief Function implementations for the derived wall classes. */

namespace voro {

/** Tests to see whether a point is inside the sphere wall object.
   * \param[in,out] (x,y,z) the vector to test.
   * \return True if the point is inside, false if the point is outside. */
  bool wall_sphere::point_inside(double x, double y, double z) {
    return (x - xc) * (x - xc) + (y - yc) * (y - yc) + (z - zc) * (z - zc) < rc * rc;
  }

/** Cuts a cell by the sphere wall object. The spherical wall is approximated by
   * a single plane applied at the point on the sphere which is closest to the center
   * of the cell. This works well for particle arrangements that are packed against
   * the wall, but loses accuracy for sparse particle distributions.
   * \param[in,out] c the Voronoi cell to be cut.
   * \param[in] (x,y,z) the location of the Voronoi cell.
   * \return True if the cell still exists, false if the cell is deleted. */
  template <class v_cell> bool wall_sphere::cut_cell_base(v_cell &c, double x, double y, double z) {
    const double xd = x - xc;
    const double yd = y - yc;
    const double zd = z - zc;
    double dq = xd * xd + yd * yd + zd * zd;
    if (dq > 1e-5) {
      dq = 2 * (sqrt(dq) * rc - dq);
      return c.nplane(xd, yd, zd, dq, w_id);
    }
    return true;
  }

/** Tests to see whether a point is inside the plane wall object.
   * \param[in] (x,y,z) the vector to test.
   * \return True if the point is inside, false if the point is outside. */
  bool wall_plane::point_inside(double x, double y, double z) {
    return x * xc + y * yc + z * zc < ac;
  }

/** Cuts a cell by the plane wall object.
   * \param[in,out] c the Voronoi cell to be cut.
   * \param[in] (x,y,z) the location of the Voronoi cell.
   * \return True if the cell still exists, false if the cell is deleted. */
  template <class v_cell> bool wall_plane::cut_cell_base(v_cell &c, double x, double y, double z) {
    const double dq = 2 * (ac - x * xc - y * yc - z * zc);
    return c.nplane(xc, yc, zc, dq, w_id);
  }

/** Tests to see whether a point is inside the cylindrical wall object.
   * \param[in] (x,y,z) the vector to test.
   * \return True if the point is inside, false if the point is outside. */
  bool wall_cylinder::point_inside(double x, double y, double z) {
    double xd = x - xc;
    double yd = y - yc;
    double zd = z - zc;
    const double pa = (xd * xa + yd * ya + zd * za) * asi;
    xd -= xa * pa;
    yd -= ya * pa;
    zd -= za * pa;
    return xd * xd + yd * yd + zd * zd < rc * rc;
  }

/** Cuts a cell by the cylindrical wall object. The cylindrical wall is
   * approximated by a single plane applied at the point on the cylinder which is
   * closest to the center of the cell. This works well for particle arrangements
   * that are packed against the wall, but loses accuracy for sparse particle
   * distributions.
   * \param[in,out] c the Voronoi cell to be cut.
   * \param[in] (x,y,z) the location of the Voronoi cell.
   * \return True if the cell still exists, false if the cell is deleted. */
  template <class v_cell> bool wall_cylinder::cut_cell_base(v_cell &c, double x, double y, double z) {
    double xd = x - xc;
    double yd = y - yc;
    double zd = z - zc;
    double pa = (xd * xa + yd * ya + zd * za) * asi;
    xd -= xa * pa;
    yd -= ya * pa;
    zd -= za * pa;
    pa = xd * xd + yd * yd + zd * zd;
    if (pa > 1e-5) {
      pa = 2 * (sqrt(pa) * rc - pa);
      return c.nplane(xd, yd, zd, pa, w_id);
    }
    return true;
  }

/** Tests to see whether a point is inside the cone wall object.
   * \param[in] (x,y,z) the vector to test.
   * \return True if the point is inside, false if the point is outside. */
  bool wall_cone::point_inside(double x, double y, double z) {
    double xd = x - xc;
    double yd = y - yc;
    double zd = z - zc;
    double pa = (xd * xa + yd * ya + zd * za) * asi;
    xd -= xa * pa;
    yd -= ya * pa;
    zd -= za * pa;
    pa *= gra;
    if (pa < 0) {
      return false;
    }
    pa *= pa;
    return xd * xd + yd * yd + zd * zd < pa;
  }

/** Cuts a cell by the cone wall object. The conical wall is
   * approximated by a single plane applied at the point on the cone which is
   * closest to the center of the cell. This works well for particle arrangements
   * that are packed against the wall, but loses accuracy for sparse particle
   * distributions.
   * \param[in,out] c the Voronoi cell to be cut.
   * \param[in] (x,y,z) the location of the Voronoi cell.
   * \return True if the cell still exists, false if the cell is deleted. */
  template <class v_cell> bool wall_cone::cut_cell_base(v_cell &c, double x, double y, double z) {
    double xd = x - xc;
    double yd = y - yc;
    double zd = z - zc;
    double xf;
    double yf;
    double zf;
    double q;
    double pa = (xd * xa + yd * ya + zd * za) * asi;
    xd -= xa * pa;
    yd -= ya * pa;
    zd -= za * pa;
    pa = xd * xd + yd * yd + zd * zd;
    if (pa > 1e-5) {
      pa = 1 / sqrt(pa);
      q = sqrt(asi);
      xf = -sang * q * xa + cang * pa * xd;
      yf = -sang * q * ya + cang * pa * yd;
      zf = -sang * q * za + cang * pa * zd;
      pa = 2 * (xf * (xc - x) + yf * (yc - y) + zf * (zc - z));
      return c.nplane(xf, yf, zf, pa, w_id);
    }
    return true;
  }

// Explicit instantiation
  template bool wall_sphere::cut_cell_base(voronoicell &, double, double, double);
  template bool wall_sphere::cut_cell_base(voronoicell_neighbor &, double, double, double);
  template bool wall_plane::cut_cell_base(voronoicell &, double, double, double);
  template bool wall_plane::cut_cell_base(voronoicell_neighbor &, double, double, double);
  template bool wall_cylinder::cut_cell_base(voronoicell &, double, double, double);
  template bool wall_cylinder::cut_cell_base(voronoicell_neighbor &, double, double, double);
  template bool wall_cone::cut_cell_base(voronoicell &, double, double, double);
  template bool wall_cone::cut_cell_base(voronoicell_neighbor &, double, double, double);

} // namespace voro
