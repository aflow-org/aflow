// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2021           *
// *                                                                         *
// ***************************************************************************

// aflow_voro.h and aflow_voro.cpp are based on the Voro++ created by Chris H. Rycroft
// under the following original licence

// ********************************************************************************
// Voro++ Copyright (c) 2008, The Regents of the University of California, through
// Lawrence Berkeley National Laboratory (subject to receipt of any required
// approvals from the U.S. Dept. of Energy). All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// (1) Redistributions of source code must retain the above copyright notice, this
// list of conditions and the following disclaimer.
//
// (2) Redistributions in binary form must reproduce the above copyright notice,
// this list of conditions and the following disclaimer in the documentation
// and/or other materials provided with the distribution.
//
// (3) Neither the name of the University of California, Lawrence Berkeley
// National Laboratory, U.S. Dept. of Energy nor the names of its contributors may
// be used to endorse or promote products derived from this software without
// specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
// ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
// ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// You are under no obligation whatsoever to provide any bug fixes, patches, or
// upgrades to the features, functionality or performance of the source code
// ("Enhancements") to anyone; however, if you choose to make your Enhancements
// available either publicly, or directly to Lawrence Berkeley National
// Laboratory, without imposing a separate written license agreement for such
// Enhancements, then you hereby grant the following license: a  non-exclusive,
// royalty-free perpetual license to install, use, modify, prepare derivative
// works, incorporate into other computer software, distribute, and sublicense
// such enhancements or derivative works thereof, in binary and source code form.
// ********************************************************************************

// Main modifications:
//   - reduced scope (user interfaces, cutting planes, ...)
//   - changed to AFLOW style error and output style
//   - reduce the number of compiler warnings

// Voro++, a 3D cell-based Voronoi library
//
// Author   : Chris H. Rycroft (LBL / UC Berkeley)
// Email    : chr@alum.mit.edu
// Date     : August 30th 2011

/** \file voro++.hh
 * \brief A file that loads all of the Voro++ header files. */

#ifndef VOROPP_HH
#define VOROPP_HH

// Voro++, a 3D cell-based Voronoi library
//
// Author   : Chris H. Rycroft (LBL / UC Berkeley)
// Email    : chr@alum.mit.edu
// Date     : August 30th 2011

/** \file config.hh
 * \brief Master configuration file for setting various compile-time options. */

#include "AUROSTD/aurostd_xvector.h"
#ifndef VOROPP_CONFIG_HH
#define VOROPP_CONFIG_HH

#include <cassert>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <vector>

#include <limits.h>
#include <stdint.h>
#include <stdio.h>

#include "AUROSTD/aurostd.h"

namespace voro {

// These constants set the initial memory allocation for the Voronoi cell
/** The initial memory allocation for the number of vertices. */
  const int init_vertices = 256;
/** The initial memory allocation for the maximum vertex order. */
  const int init_vertex_order = 64;
/** The initial memory allocation for the number of regular vertices of order
   * 3. */
  const int init_3_vertices = 256;
/** The initial memory allocation for the number of vertices of higher order.
   */
  const int init_n_vertices = 8;
/** The initial buffer size for marginal cases used by the suretest class. */
  const int init_marginal = 64;
/** The initial size for the delete stack. */
  const int init_delete_size = 256;
/** The initial size for the auxiliary delete stack. */
  const int init_delete2_size = 256;
/** The initial size for the wall pointer array. */
  const int init_wall_size = 32;
/** The default initial size for the ordering class. */
  const int init_ordering_size = 4096;
/** The initial size of the pre_container chunk index. */
  const int init_chunk_size = 256;

// If the initial memory is too small, the program dynamically allocates more.
// However, if the limits below are reached, then the program bails out.
/** The maximum memory allocation for the number of vertices. */
  const int max_vertices = 16777216;
/** The maximum memory allocation for the maximum vertex order. */
  const int max_vertex_order = 2048;
/** The maximum memory allocation for the any particular order of vertex. */
  const int max_n_vertices = 16777216;
/** The maximum buffer size for marginal cases used by the suretest class. */
  const int max_marginal = 16777216;
/** The maximum size for the delete stack. */
  const int max_delete_size = 16777216;
/** The maximum size for the auxiliary delete stack. */
  const int max_delete2_size = 16777216;
/** The maximum amount of particle memory allocated for a single region. */
  const int max_particle_memory = 16777216;
/** The maximum size for the wall pointer array. */
  const int max_wall_size = 2048;
/** The maximum size for the ordering class. */
  const int max_ordering_size = 67108864;
/** The maximum size for the pre_container chunk index. */
  const int max_chunk_size = 65536;

/** The chunk size in the pre_container classes. */
  const int pre_container_chunk_size = 1024;

#ifndef VOROPP_VERBOSE
/** Voro++ can print a number of different status and debugging messages to
 * notify the user of special behavior, and this macro sets the amount which
 * are displayed. At level 0, no messages are printed. At level 1, messages
 * about unusual cases during cell construction are printed, such as when the
 * plane routine bails out due to floating point problems. At level 2, general
 * messages about memory expansion are printed. At level 3, technical details
 * about memory management are printed. */
#define VOROPP_VERBOSE 0
#endif

/** If a point is within this distance of a cutting plane, then the code
   * assumes that point exactly lies on the plane. */
  const double tolerance = 1e-11;

/** If a point is within this distance of a cutting plane, then the code stores
   * whether this point is inside, outside, or exactly on the cutting plane in
   * the marginal cases buffer, to prevent the test giving a different result on
   * a subsequent evaluation due to floating point rounding errors. */
  const double tolerance2 = 2e-11;

/** The square of the tolerance, used when deciding whether some squared
   * quantities are large enough to be used. */
  const double tolerance_sq = tolerance * tolerance;

/** A large number that is used in the computation. */
  const double large_number = 1e30;

/** A radius to use as a placeholder when no other information is available. */
  const double default_radius = 0.5;

/** The maximum number of shells of periodic images to test over. */
  const int max_unit_voro_shells = 10;

/** A guess for the optimal number of particles per block, used to set up the
   * container grid. */
  const double optimal_particles = 5.6;

/** If this is set to 1, then the code reports any instances of particles being
 * put outside of the container geometry. */
#define VOROPP_REPORT_OUT_OF_BOUNDS 0

/** Voro++ returns this status code if there is a file-related error, such as
 * not being able to open file. */
#define VOROPP_FILE_ERROR 1

/** Voro++ returns this status code if there is a memory allocation error, if
 * one of the safe memory limits is exceeded. */
#define VOROPP_MEMORY_ERROR 2

/** Voro++ returns this status code if there is any type of internal error, if
 * it detects that representation of the Voronoi cell is inconsistent. This
 * status code will generally indicate a bug, and the developer should be
 * contacted. */
#define VOROPP_INTERNAL_ERROR 3

/** Voro++ returns this status code if it could not interpret the command line
 * arguments passed to the command line utility. */
#define VOROPP_CMD_LINE_ERROR 4

} // namespace voro

#endif

// Voro++, a 3D cell-based Voronoi library
//
// Author   : Chris H. Rycroft (LBL / UC Berkeley)
// Email    : chr@alum.mit.edu
// Date     : August 30th 2011

/** \file common.hh
 * \brief Header file for the small helper functions. */

#ifndef VOROPP_COMMON_HH
#define VOROPP_COMMON_HH

#include <cstdio>
#include <cstdlib>
#include <vector>

namespace voro {

/** \brief Function for printing fatal error messages and exiting.
   *
   * Function for printing fatal error messages and exiting.
   * \param[in] p a pointer to the message to print.
   * \param[in] status the status code to return with. */
  inline void voro_fatal_error(const char *p, int status) {
    fprintf(stderr, "voro++: %s\n", p);
    exit(status);
  }

/** \brief Prints a vector of positions.
   *
   * Prints a vector of positions as bracketed triplets.
   * \param[in] v the vector to print.
   * \param[in] fp the file stream to print to. */
  inline void voro_print_positions(std::vector<double> &v, FILE *fp = stdout) {
    if (!v.empty()) {
      fprintf(fp, "(%g,%g,%g)", v[0], v[1], v[2]);
      for (int k = 3; (unsigned int) k < v.size(); k += 3) {
        fprintf(fp, " (%g,%g,%g)", v[k], v[k + 1], v[k + 2]);
      }
    }
  }

/** \brief Opens a file and checks the operation was successful.
   *
   * Opens a file, and checks the return value to ensure that the operation
   * was successful.
   * \param[in] filename the file to open.
   * \param[in] mode the cstdio fopen mode to use.
   * \return The file handle. */
  inline FILE *safe_fopen(const char *filename, const char *mode) {
    FILE *fp = fopen(filename, mode);
    if (fp == nullptr) {
      fprintf(stderr, "voro++: Unable to open file '%s'\n", filename);
      exit(VOROPP_FILE_ERROR);
    }
    return fp;
  }

  void voro_print_vector(std::vector<int> &v, FILE *fp = stdout);
  void voro_print_vector(std::vector<double> &v, FILE *fp = stdout);
  void voro_print_face_vertices(std::vector<int> &v, FILE *fp = stdout);

} // namespace voro

#endif

// Voro++, a 3D cell-based Voronoi library
//
// Author   : Chris H. Rycroft (LBL / UC Berkeley)
// Email    : chr@alum.mit.edu
// Date     : August 30th 2011

/** \file cell.hh
 * \brief Header file for the voronoicell and related classes. */

#ifndef VOROPP_CELL_HH
#define VOROPP_CELL_HH

#include <vector>

namespace voro {

/** \brief A class representing a single Voronoi cell.
   *
   * This class represents a single Voronoi cell, as a collection of vertices
   * that are connected by edges. The class contains routines for initializing
   * the Voronoi cell to be simple shapes such as a box, tetrahedron, or octahedron.
   * It the contains routines for recomputing the cell based on cutting it
   * by a plane, which forms the key routine for the Voronoi cell computation.
   * It contains numerous routine for computing statistics about the Voronoi cell,
   * and it can output the cell in several formats.
   *
   * This class is not intended for direct use, but forms the base of the
   * voronoicell and voronoicell_neighbor classes, which extend it based on
   * whether neighboring particle ID information needs to be tracked. */
  class voronoicell_base {
  public:
    /** This holds the current size of the arrays ed and nu, which
     * hold the vertex information. If more vertices are created
     * than can fit in this array, then it is dynamically extended
     * using the add_memory_vertices routine. */
    int current_vertices;
    /** This holds the current maximum allowed order of a vertex,
     * which sets the size of the mem, mep, and mec arrays. If a
     * vertex is created with more vertices than this, the arrays
     * are dynamically extended using the add_memory_vorder routine.
     */
    int current_vertex_order;
    /** This sets the size of the main delete stack. */
    int current_delete_size;
    /** This sets the size of the auxiliary delete stack. */
    int current_delete2_size;
    /** This sets the total number of vertices in the current cell.
     */
    int p;
    /** This is the index of particular point in the cell, which is
     * used to start the tracing routines for plane intersection
     * and cutting. These routines will work starting from any
     * point, but it's often most efficient to start from the last
     * point considered, since in many cases, the cell construction
     * algorithm may consider many planes with similar vectors
     * concurrently. */
    int up;
    /** This is a two dimensional array that holds information
     * about the edge connections of the vertices that make up the
     * cell. The two dimensional array is not allocated in the
     * usual method. To account for the fact the different vertices
     * have different orders, and thus require different amounts of
     * storage, the elements of ed[i] point to one-dimensional
     * arrays in the mep[] array of different sizes.
     *
     * More specifically, if vertex i has order m, then ed[i]
     * points to a one-dimensional array in mep[m] that has 2*m+1
     * entries. The first m elements hold the neighboring edges, so
     * that the jth edge of vertex i is held in ed[i][j]. The next
     * m elements hold a table of relations which is redundant but
     * helps speed up the computation. It satisfies the relation
     * ed[ed[i][j]][ed[i][m+j]]=i. The final entry holds a back
     * pointer, so that ed[i+2*m]=i. The back pointers are used
     * when rearranging the memory. */
    int **ed;
    /** This array holds the order of the vertices in the Voronoi
     * cell. This array is dynamically allocated, with its current
     * size held by current_vertices. */
    int *nu;
    /** This in an array with size 3*current_vertices for holding
     * the positions of the vertices. */
    double *pts;
    voronoicell_base();
    ~voronoicell_base();
    void init_base(double xmin, double xmax, double ymin, double ymax, double zmin, double zmax);
    void init_octahedron_base(double l);
    void init_tetrahedron_base(double x0, double y0, double z0, double x1, double y1, double z1, double x2, double y2, double z2, double x3, double y3, double z3);
    void translate(double x, double y, double z) const;
    void draw_pov(double x, double y, double z, FILE *fp = stdout) const;
    /** Outputs the cell in POV-Ray format, using cylinders for edges
     * and spheres for vertices, to a given file.
     * \param[in] (x,y,z) a displacement to add to the cell's
     *                    position.
     * \param[in] filename the name of the file to write to. */
    inline void draw_pov(double x, double y, double z, const char *filename) const {
      FILE *fp = safe_fopen(filename, "w");
      draw_pov(x, y, z, fp);
      fclose(fp);
    };

    /// AFLOW //HE20220912
    double pbc_distance(const aurostd::xvector<double> &normal, const aurostd::xmatrix<double> &lattice, const aurostd::xvector<double> &p1, const aurostd::xvector<double> &p2);
    double pbc_distance_min(const aurostd::xmatrix<double> &lattice, const aurostd::xvector<double> &p1, const aurostd::xvector<double> &p2);

    /// AFLOW //HE20220908
    void get_vertex_facets(double x, double y, double z, std::vector<aurostd::xvector<double>> &points, std::vector<std::vector<uint>> &facets);

    /// AFLOW //HE20220910
    void fill_laplacian(const int p_id, const std::vector<aurostd::xvector<double>> &vertices, const aurostd::xmatrix<double> &lattice, aurostd::xmatrix<double> &laplacian);

    void draw_pov_mesh(double x, double y, double z, FILE *fp = stdout);
    /** Outputs the cell in POV-Ray format as a mesh2 object to a
     * given file.
     * \param[in] (x,y,z) a displacement to add to the cell's
     *                    position.
     * \param[in] filename the name of the file to write to. */
    inline void draw_pov_mesh(double x, double y, double z, const char *filename) {
      FILE *fp = safe_fopen(filename, "w");
      draw_pov_mesh(x, y, z, fp);
      fclose(fp);
    }
    void draw_gnuplot(double x, double y, double z, FILE *fp = stdout);
    /** Outputs the cell in Gnuplot format a given file.
     * \param[in] (x,y,z) a displacement to add to the cell's
     *                    position.
     * \param[in] filename the name of the file to write to. */
    inline void draw_gnuplot(double x, double y, double z, const char *filename) {
      FILE *fp = safe_fopen(filename, "w");
      draw_gnuplot(x, y, z, fp);
      fclose(fp);
    }
    double volume();
    [[nodiscard]] double max_radius_squared() const;
    [[nodiscard]] double total_edge_distance() const;
    double surface_area();
    void centroid(double &cx, double &cy, double &cz);
    int number_of_faces();
    [[nodiscard]] int number_of_edges() const;
    void vertex_orders(std::vector<int> &v) const;
    void output_vertex_orders(FILE *fp = stdout) const;
    void vertices(std::vector<double> &v) const;
    void output_vertices(FILE *fp = stdout) const;
    void vertices(double x, double y, double z, std::vector<double> &v) const;
    void output_vertices(double x, double y, double z, FILE *fp = stdout) const;
    void face_areas(std::vector<double> &v);
    /** Outputs the areas of the faces.
     * \param[in] fp the file handle to write to. */
    inline void output_face_areas(FILE *fp = stdout) {
      std::vector<double> v;
      face_areas(v);
      voro_print_vector(v, fp);
    }
    void face_orders(std::vector<int> &v);
    /** Outputs a list of the number of sides of each face.
     * \param[in] fp the file handle to write to. */
    inline void output_face_orders(FILE *fp = stdout) {
      std::vector<int> v;
      face_orders(v);
      voro_print_vector(v, fp);
    }
    void face_freq_table(std::vector<int> &v);
    /** Outputs a */
    inline void output_face_freq_table(FILE *fp = stdout) {
      std::vector<int> v;
      face_freq_table(v);
      voro_print_vector(v, fp);
    }
    void face_vertices(std::vector<int> &v);
    /** Outputs the */
    inline void output_face_vertices(FILE *fp = stdout) {
      std::vector<int> v;
      face_vertices(v);
      voro_print_face_vertices(v, fp);
    }
    void face_perimeters(std::vector<double> &v);
    /** Outputs a list of the perimeters of each face.
     * \param[in] fp the file handle to write to. */
    inline void output_face_perimeters(FILE *fp = stdout) {
      std::vector<double> v;
      face_perimeters(v);
      voro_print_vector(v, fp);
    }
    void normals(std::vector<double> &v);
    /** Outputs a list of the perimeters of each face.
     * \param[in] fp the file handle to write to. */
    inline void output_normals(FILE *fp = stdout) {
      std::vector<double> v;
      normals(v);
      voro_print_positions(v, fp);
    }
    /** Outputs a custom string of information about the Voronoi
     * cell to a file. It assumes the cell is at (0,0,0) and has a
     * the default_radius associated with it.
     * \param[in] format the custom format string to use.
     * \param[in] fp the file handle to write to. */
    inline void output_custom(const char *format, FILE *fp = stdout) { output_custom(format, 0, 0, 0, 0, default_radius, fp); }
    void output_custom(const char *format, int i, double x, double y, double z, double r, FILE *fp = stdout);
    template <class vc_class> bool nplane(vc_class &vc, double x, double y, double z, double rsq, int p_id);
    bool plane_intersects(double x, double y, double z, double rsq);
    bool plane_intersects_guess(double x, double y, double z, double rsq);
    void construct_relations() const;
    void check_relations() const;
    void check_duplicates() const;
    void print_edges();
    /** Returns a list of IDs of neighboring particles
     * corresponding to each face.
     * \param[out] v a reference to a vector in which to return the
     *               results. If no neighbor information is
     *               available, a blank vector is returned. */
    virtual void neighbors(std::vector<int> &v) { v.clear(); }
    /** This is a virtual function that is overridden by a routine
     * to print a list of IDs of neighboring particles
     * corresponding to each face. By default, when no neighbor
     * information is available, the routine does nothing.
     * \param[in] fp the file handle to write to. */
    virtual void output_neighbors(FILE *fp = stdout) {}
    /** This a virtual function that is overridden by a routine to
     * print the neighboring particle IDs for a given vertex. By
     * default, when no neighbor information is available, the
     * routine does nothing.
     * \param[in] i the vertex to consider. */
    virtual void print_edges_neighbors(int i) {};
    /** This is a simple inline function for picking out the index
     * of the next edge counterclockwise at the current vertex.
     * \param[in] a the index of an edge of the current vertex.
     * \param[in] p the number of the vertex.
     * \return 0 if a=nu[p]-1, or a+1 otherwise. */
    inline int cycle_up(int a, int p) const { return a == nu[p] - 1 ? 0 : a + 1; }
    /** This is a simple inline function for picking out the index
     * of the next edge clockwise from the current vertex.
     * \param[in] a the index of an edge of the current vertex.
     * \param[in] p the number of the vertex.
     * \return nu[p]-1 if a=0, or a-1 otherwise. */
    inline int cycle_down(int a, int p) const { return a == 0 ? nu[p] - 1 : a - 1; }

  protected:
    /** This a one dimensional array that holds the current sizes
     * of the memory allocations for them mep array.*/
    int *mem;
    /** This is a one dimensional array that holds the current
     * number of vertices of order p that are stored in the mep[p]
     * array. */
    int *mec;
    /** This is a two dimensional array for holding the information
     * about the edges of the Voronoi cell. mep[p] is a
     * one-dimensional array for holding the edge information about
     * all vertices of order p, with each vertex holding 2*p+1
     * integers of information. The total number of vertices held
     * on mep[p] is stored in mem[p]. If the space runs out, the
     * code allocates more using the add_memory() routine. */
    int **mep;
    inline void reset_edges() const;
    template <class vc_class> void check_memory_for_copy(vc_class &vc, voronoicell_base *vb);
    void copy(voronoicell_base *vb);

  private:
    /** This is the delete stack, used to store the vertices which
     * are going to be deleted during the plane cutting procedure.
     */
    int *ds, *stacke;
    /** This is the auxiliary delete stack, which has size set by
     * current_delete2_size. */
    int *ds2, *stacke2;
    /** This stores the current memory allocation for the marginal
     * cases. */
    int current_marginal;
    /** This stores the total number of marginal points which are
     * currently in the buffer. */
    int n_marg;
    /** This array contains a list of the marginal points, and also
     * the outcomes of the marginal tests. */
    int *marg;
    /** The x coordinate of the normal vector to the test plane. */
    double px;
    /** The y coordinate of the normal vector to the test plane. */
    double py;
    /** The z coordinate of the normal vector to the test plane. */
    double pz;
    /** The magnitude of the normal vector to the test plane. */
    double prsq;
    template <class vc_class> void add_memory(vc_class &vc, int i, int *stackp2);
    template <class vc_class> void add_memory_vertices(vc_class &vc);
    template <class vc_class> void add_memory_vorder(vc_class &vc);
    void add_memory_ds(int *&stackp);
    void add_memory_ds2(int *&stackp2);
    template <class vc_class> inline bool collapse_order1(vc_class &vc);
    template <class vc_class> inline bool collapse_order2(vc_class &vc);
    template <class vc_class> inline bool delete_connection(vc_class &vc, int j, int k, bool hand);
    template <class vc_class> inline bool search_for_outside_edge(vc_class &vc, int &up);
    template <class vc_class> inline void add_to_stack(vc_class &vc, int lp, int *&stackp2);
    inline bool plane_intersects_track(double x, double y, double z, double rs, double g);
    inline void normals_search(std::vector<double> &v, int i, int j, int k);
    inline bool search_edge(int l, int &m, int &k) const;
    inline int m_test(int n, double &ans);
    int check_marginal(int n, double &ans);
    friend class voronoicell;
    friend class voronoicell_neighbor;
  };

/** \brief Extension of the voronoicell_base class to represent a Voronoi
   * cell without neighbor information.
   *
   * This class is an extension of the voronoicell_base class, in cases when
   * is not necessary to track the IDs of neighboring particles associated
   * with each face of the Voronoi cell. */
  class voronoicell : public voronoicell_base {
  public:
    using voronoicell_base::nplane;
    /** Copies the information from another voronoicell class into
     * this class, extending memory allocation if necessary.
     * \param[in] c the class to copy. */
    inline void operator=(voronoicell &c) {
      voronoicell_base *vb((voronoicell_base *) &c);
      check_memory_for_copy(*this, vb);
      copy(vb);
    }
    /** Cuts a Voronoi cell using by the plane corresponding to the
     * perpendicular bisector of a particle.
     * \param[in] (x,y,z) the position of the particle.
     * \param[in] rsq the modulus squared of the vector.
     * \param[in] p_id the plane ID, ignored for this case where no
     *                 neighbor tracking is enabled.
     * \return False if the plane cut deleted the cell entirely,
     *         true otherwise. */
    inline bool nplane(double x, double y, double z, double rsq, int p_id) { return nplane(*this, x, y, z, rsq, 0); }
    /** Cuts a Voronoi cell using by the plane corresponding to the
     * perpendicular bisector of a particle.
     * \param[in] (x,y,z) the position of the particle.
     * \param[in] p_id the plane ID, ignored for this case where no
     *                 neighbor tracking is enabled.
     * \return False if the plane cut deleted the cell entirely,
     *         true otherwise. */
    inline bool nplane(double x, double y, double z, int p_id) {
      const double rsq = x * x + y * y + z * z;
      return nplane(*this, x, y, z, rsq, 0);
    }
    /** Cuts a Voronoi cell using by the plane corresponding to the
     * perpendicular bisector of a particle.
     * \param[in] (x,y,z) the position of the particle.
     * \param[in] rsq the modulus squared of the vector.
     * \return False if the plane cut deleted the cell entirely,
     *         true otherwise. */
    inline bool plane(double x, double y, double z, double rsq) { return nplane(*this, x, y, z, rsq, 0); }
    /** Cuts a Voronoi cell using by the plane corresponding to the
     * perpendicular bisector of a particle.
     * \param[in] (x,y,z) the position of the particle.
     * \return False if the plane cut deleted the cell entirely,
     *         true otherwise. */
    inline bool plane(double x, double y, double z) {
      const double rsq = x * x + y * y + z * z;
      return nplane(*this, x, y, z, rsq, 0);
    }
    /** Initializes the Voronoi cell to be rectangular box with the
     * given dimensions.
     * \param[in] (xmin,xmax) the minimum and maximum x coordinates.
     * \param[in] (ymin,ymax) the minimum and maximum y coordinates.
     * \param[in] (zmin,zmax) the minimum and maximum z coordinates. */
    inline void init(double xmin, double xmax, double ymin, double ymax, double zmin, double zmax) { init_base(xmin, xmax, ymin, ymax, zmin, zmax); }
    /** Initializes the cell to be an octahedron with vertices at
     * (l,0,0), (-l,0,0), (0,l,0), (0,-l,0), (0,0,l), and (0,0,-l).
     * \param[in] l a parameter setting the size of the octahedron.
     */
    inline void init_octahedron(double l) { init_octahedron_base(l); }
    /** Initializes the cell to be a tetrahedron.
     * \param[in] (x0,y0,z0) the coordinates of the first vertex.
     * \param[in] (x1,y1,z1) the coordinates of the second vertex.
     * \param[in] (x2,y2,z2) the coordinates of the third vertex.
     * \param[in] (x3,y3,z3) the coordinates of the fourth vertex.
     */
    inline void init_tetrahedron(double x0, double y0, double z0, double x1, double y1, double z1, double x2, double y2, double z2, double x3, double y3, double z3) {
      init_tetrahedron_base(x0, y0, z0, x1, y1, z1, x2, y2, z2, x3, y3, z3);
    }

  private:
    inline void n_allocate(int i, int m) {};
    inline void n_add_memory_vertices(int i) {};
    inline void n_add_memory_vorder(int i) {};
    inline void n_set_pointer(int p, int n) {};
    inline void n_copy(int a, int b, int c, int d) {};
    inline void n_set(int a, int b, int c) {};
    inline void n_set_aux1(int k) {};
    inline void n_copy_aux1(int a, int b) {};
    inline void n_copy_aux1_shift(int a, int b) {};
    inline void n_set_aux2_copy(int a, int b) {};
    inline void n_copy_pointer(int a, int b) {};
    inline void n_set_to_aux1(int j) {};
    inline void n_set_to_aux2(int j) {};
    inline void n_allocate_aux1(int i) {};
    inline void n_switch_to_aux1(int i) {};
    inline void n_copy_to_aux1(int i, int m) {};
    inline void n_set_to_aux1_offset(int k, int m) {};
    inline void n_neighbors(std::vector<int> &v) { v.clear(); };
    friend class voronoicell_base;
  };

/** \brief Extension of the voronoicell_base class to represent a Voronoi cell
   * with neighbor information.
   *
   * This class is an extension of the voronoicell_base class, in cases when the
   * IDs of neighboring particles associated with each face of the Voronoi cell.
   * It contains additional data structures mne and ne for storing this
   * information. */
  class voronoicell_neighbor : public voronoicell_base {
  public:
    using voronoicell_base::nplane;
    /** This two dimensional array holds the neighbor information
     * associated with each vertex. mne[p] is a one dimensional
     * array which holds all of the neighbor information for
     * vertices of order p. */
    int **mne;
    /** This is a two dimensional array that holds the neighbor
     * information associated with each vertex. ne[i] points to a
     * one-dimensional array in mne[nu[i]]. ne[i][j] holds the
     * neighbor information associated with the jth edge of vertex
     * i. It is set to the ID number of the plane that made the
     * face that is clockwise from the jth edge. */
    int **ne;
    voronoicell_neighbor();
    ~voronoicell_neighbor();
    void operator=(voronoicell &c);
    void operator=(voronoicell_neighbor &c);
    /** Cuts the Voronoi cell by a particle whose center is at a
     * separation of (x,y,z) from the cell center. The value of rsq
     * should be initially set to \f$x^2+y^2+z^2\f$.
     * \param[in] (x,y,z) the normal vector to the plane.
     * \param[in] rsq the distance along this vector of the plane.
     * \param[in] p_id the plane ID (for neighbor tracking only).
     * \return False if the plane cut deleted the cell entirely,
     * true otherwise. */
    inline bool nplane(double x, double y, double z, double rsq, int p_id) { return nplane(*this, x, y, z, rsq, p_id); }
    /** This routine calculates the modulus squared of the vector
     * before passing it to the main nplane() routine with full
     * arguments.
     * \param[in] (x,y,z) the vector to cut the cell by.
     * \param[in] p_id the plane ID (for neighbor tracking only).
     * \return False if the plane cut deleted the cell entirely,
     *         true otherwise. */
    inline bool nplane(double x, double y, double z, int p_id) {
      const double rsq = x * x + y * y + z * z;
      return nplane(*this, x, y, z, rsq, p_id);
    }
    /** This version of the plane routine just makes up the plane
     * ID to be zero. It will only be referenced if neighbor
     * tracking is enabled.
     * \param[in] (x,y,z) the vector to cut the cell by.
     * \param[in] rsq the modulus squared of the vector.
     * \return False if the plane cut deleted the cell entirely,
     *         true otherwise. */
    inline bool plane(double x, double y, double z, double rsq) { return nplane(*this, x, y, z, rsq, 0); }
    /** Cuts a Voronoi cell using the influence of a particle at
     * (x,y,z), first calculating the modulus squared of this
     * vector before passing it to the main nplane() routine. Zero
     * is supplied as the plane ID, which will be ignored unless
     * neighbor tracking is enabled.
     * \param[in] (x,y,z) the vector to cut the cell by.
     * \return False if the plane cut deleted the cell entirely,
     *         true otherwise. */
    inline bool plane(double x, double y, double z) {
      const double rsq = x * x + y * y + z * z;
      return nplane(*this, x, y, z, rsq, 0);
    }
    void init(double xmin, double xmax, double ymin, double ymax, double zmin, double zmax);
    void init_octahedron(double l);
    void init_tetrahedron(double x0, double y0, double z0, double x1, double y1, double z1, double x2, double y2, double z2, double x3, double y3, double z3);
    void check_facets();
    virtual void neighbors(std::vector<int> &v);
    virtual void print_edges_neighbors(int i);
    virtual void output_neighbors(FILE *fp = stdout) {
      std::vector<int> v;
      neighbors(v);
      voro_print_vector(v, fp);
    }

  private:
    int *paux1;
    int *paux2;
    inline void n_allocate(int i, int m) const { mne[i] = new int[m * i]; }
    inline void n_add_memory_vertices(int i) {
      int **pp = new int *[i];
      for (int j = 0; j < current_vertices; j++) {
        pp[j] = ne[j];
      }
      delete[] ne;
      ne = pp;
    }
    inline void n_add_memory_vorder(int i) {
      int **p2 = new int *[i];
      for (int j = 0; j < current_vertex_order; j++) {
        p2[j] = mne[j];
      }
      delete[] mne;
      mne = p2;
    }
    inline void n_set_pointer(int p, int n) { ne[p] = mne[n] + n * mec[n]; }
    inline void n_copy(int a, int b, int c, int d) const { ne[a][b] = ne[c][d]; }
    inline void n_set(int a, int b, int c) const { ne[a][b] = c; }
    inline void n_set_aux1(int k) { paux1 = mne[k] + k * mec[k]; }
    inline void n_copy_aux1(int a, int b) { paux1[b] = ne[a][b]; }
    inline void n_copy_aux1_shift(int a, int b) { paux1[b] = ne[a][b + 1]; }
    inline void n_set_aux2_copy(int a, int b) {
      paux2 = mne[b] + b * mec[b];
      for (int i = 0; i < b; i++) {
        ne[a][i] = paux2[i];
      }
    }
    inline void n_copy_pointer(int a, int b) const { ne[a] = ne[b]; }
    inline void n_set_to_aux1(int j) { ne[j] = paux1; }
    inline void n_set_to_aux2(int j) { ne[j] = paux2; }
    inline void n_allocate_aux1(int i) { paux1 = new int[i * mem[i]]; }
    inline void n_switch_to_aux1(int i) {
      delete[] mne[i];
      mne[i] = paux1;
    }
    inline void n_copy_to_aux1(int i, int m) { paux1[m] = mne[i][m]; }
    inline void n_set_to_aux1_offset(int k, int m) { ne[k] = paux1 + m; }
    friend class voronoicell_base;
  };

} // namespace voro

#endif

// Voro++, a 3D cell-based Voronoi library
//
// Author   : Chris H. Rycroft (LBL / UC Berkeley)
// Email    : chr@alum.mit.edu
// Date     : August 30th 2011

/** \file v_base.hh
 * \brief Header file for the base Voronoi container class. */

#ifndef VOROPP_V_BASE_HH
#define VOROPP_V_BASE_HH

// Voro++, a 3D cell-based Voronoi library
//
// Author   : Chris H. Rycroft (LBL / UC Berkeley)
// Email    : chr@alum.mit.edu
// Date     : August 30th 2011

/** \file worklist.hh
 * \brief Header file for setting constants used in the block worklists that are
 * used during cell computation.
 *
 * This file is automatically generated by worklist_gen.pl and it is not
 * intended to be edited by hand. */

#ifndef VOROPP_WORKLIST_HH
#define VOROPP_WORKLIST_HH

namespace voro {

/** Each region is divided into a grid of subregions, and a worklist is
          # constructed for each. This parameter sets is set to half the number of
          # subregions that the block is divided into. */
  const int wl_hgrid = 4;
/** The number of subregions that a block is subdivided into, which is twice
          the value of hgrid. */
  const int wl_fgrid = 8;
/** The total number of worklists, set to the cube of hgrid. */
  const int wl_hgridcu = 64;
/** The number of elements in each worklist. */
  const int wl_seq_length = 64;

} // namespace voro
#endif

namespace voro {

/** \brief Class containing data structures common across all particle container classes.
   *
   * This class contains constants and data structures that are common across all
   * particle container classes. It contains constants setting the size of the
   * underlying subgrid of blocks that forms the basis of the Voronoi cell
   * computations. It also constructs bound tables that are used in the Voronoi
   * cell computation, and contains a number of routines that are common across
   * all container classes. */
  class voro_base {
  public:
    /** The number of blocks in the x direction. */
    const int nx;
    /** The number of blocks in the y direction. */
    const int ny;
    /** The number of blocks in the z direction. */
    const int nz;
    /** A constant, set to the value of nx multiplied by ny, which
     * is used in the routines that step through blocks in
     * sequence. */
    const int nxy;
    /** A constant, set to the value of nx*ny*nz, which is used in
     * the routines that step through blocks in sequence. */
    const int nxyz;
    /** The size of a computational block in the x direction. */
    const double boxx;
    /** The size of a computational block in the y direction. */
    const double boxy;
    /** The size of a computational block in the z direction. */
    const double boxz;
    /** The inverse box length in the x direction. */
    const double xsp;
    /** The inverse box length in the y direction. */
    const double ysp;
    /** The inverse box length in the z direction. */
    const double zsp;
    /** An array to hold the minimum distances associated with the
     * worklists. This array is initialized during container
     * construction, by the initialize_radii() routine. */
    double *mrad;
    /** The pre-computed block worklists. */
    static const unsigned int wl[wl_seq_length * wl_hgridcu];
    bool contains_neighbor(const char *format);
    voro_base(int nx_, int ny_, int nz_, double boxx_, double boxy_, double boxz_);
    ~voro_base() { delete[] mrad; }

  protected:
    /** A custom int function that returns consistent stepping
     * for negative numbers, so that (-1.5, -0.5, 0.5, 1.5) maps
     * to (-2,-1,0,1).
     * \param[in] a the number to consider.
     * \return The value of the custom int operation. */
    inline int step_int(double a) { return a < 0 ? int(a) - 1 : int(a); }
    /** A custom modulo function that returns consistent stepping
     * for negative numbers. For example, (-2,-1,0,1,2) step_mod 2
     * is (0,1,0,1,0).
     * \param[in] (a,b) the input integers.
     * \return The value of a modulo b, consistent for negative
     * numbers. */
    inline int step_mod(int a, int b) { return a >= 0 ? a % b : b - 1 - (b - 1 - a) % b; }
    /** A custom integer division function that returns consistent
     * stepping for negative numbers. For example, (-2,-1,0,1,2)
     * step_div 2 is (-1,-1,0,0,1).
     * \param[in] (a,b) the input integers.
     * \return The value of a div b, consistent for negative
     * numbers. */
    inline int step_div(int a, int b) { return a >= 0 ? a / b : -1 + (a + 1) / b; }

  private:
    void compute_minimum(double &minr, double &xlo, double &xhi, double &ylo, double &yhi, double &zlo, double &zhi, int ti, int tj, int tk) const;
  };

} // namespace voro

#endif

// Voro++, a 3D cell-based Voronoi library
//
// Author   : Chris H. Rycroft (LBL / UC Berkeley)
// Email    : chr@alum.mit.edu
// Date     : August 30th 2011

/** \file rad_option.hh
 * \brief Header file for the classes encapsulating functionality for the
 * regular and radical Voronoi tessellations. */

#ifndef VOROPP_RAD_OPTION_HH
#define VOROPP_RAD_OPTION_HH

#include <cmath>

namespace voro {

/** \brief Class containing all of the routines that are specific to computing
   * the regular Voronoi tessellation.
   *
   * The container and container_periodic classes are derived from this class,
   * and during the Voronoi cell computation, these routines are used to create
   * the regular Voronoi tessellation. */
  class radius_mono {
  protected:
    /** This is called prior to computing a Voronoi cell for a
     * given particle to initialize any required constants.
     * \param[in] ijk the block that the particle is within.
     * \param[in] s the index of the particle within the block. */
    inline void r_init(int ijk, int s) {}
    /** Sets a required constant to be used when carrying out a
     * plane bounds check. */
    inline void r_prime(double rv) {}
    /** Carries out a radius bounds check.
     * \param[in] crs the radius squared to be tested.
     * \param[in] mrs the current maximum distance to a Voronoi
     *                vertex multiplied by two.
     * \return True if particles at this radius could not possibly
     * cut the cell, false otherwise. */
    inline bool r_ctest(double crs, double mrs) { return crs > mrs; }
    /** Scales a plane displacement during a plane bounds check.
     * \param[in] lrs the plane displacement.
     * \return The scaled value. */
    inline double r_cutoff(double lrs) { return lrs; }
    /** Adds the maximum radius squared to a given value.
     * \param[in] rs the value to consider.
     * \return The value with the radius squared added. */
    inline double r_max_add(double rs) { return rs; }
    /** Subtracts the radius squared of a particle from a given
     * value.
     * \param[in] rs the value to consider.
     * \param[in] ijk the block that the particle is within.
     * \param[in] q the index of the particle within the block.
     * \return The value with the radius squared subtracted. */
    inline double r_current_sub(double rs, int ijk, int q) { return rs; }
    /** Scales a plane displacement prior to use in the plane cutting
     * algorithm.
     * \param[in] rs the initial plane displacement.
     * \param[in] ijk the block that the particle is within.
     * \param[in] q the index of the particle within the block.
     * \return The scaled plane displacement. */
    inline double r_scale(double rs, int ijk, int q) { return rs; }
    /** Scales a plane displacement prior to use in the plane
     * cutting algorithm, and also checks if it could possibly cut
     * the cell.
     * \param[in,out] rs the plane displacement to be scaled.
     * \param[in] mrs the current maximum distance to a Voronoi
     *                vertex multiplied by two.
     * \param[in] ijk the block that the particle is within.
     * \param[in] q the index of the particle within the block.
     * \return True if the cell could possibly cut the cell, false
     * otherwise. */
    inline bool r_scale_check(double &rs, double mrs, int ijk, int q) { return rs < mrs; }
  };

/**  \brief Class containing all of the routines that are specific to computing
   * the radical Voronoi tessellation.
   *
   * The container_poly and container_periodic_poly classes are derived from this
   * class, and during the Voronoi cell computation, these routines are used to
   * create the radical Voronoi tessellation. */
  class radius_poly {
  public:
    /** A two-dimensional array holding particle positions and radii. */
    double **ppr;
    /** The current maximum radius of any particle, used to
     * determine when to cut off the radical Voronoi computation.
     * */
    double max_radius;
    /** The class constructor sets the maximum particle radius to
     * be zero. */
    radius_poly() : max_radius(0) {}

  protected:
    /** This is called prior to computing a Voronoi cell for a
     * given particle to initialize any required constants.
     * \param[in] ijk the block that the particle is within.
     * \param[in] s the index of the particle within the block. */
    inline void r_init(int ijk, int s) {
      r_rad = ppr[ijk][4 * s + 3] * ppr[ijk][4 * s + 3];
      r_mul = r_rad - max_radius * max_radius;
    }
    /** Sets a required constant to be used when carrying out a
     * plane bounds check. */
    inline void r_prime(double rv) { r_val = 1 + r_mul / rv; }
    /** Carries out a radius bounds check.
     * \param[in] crs the radius squared to be tested.
     * \param[in] mrs the current maximum distance to a Voronoi
     *                vertex multiplied by two.
     * \return True if particles at this radius could not possibly
     * cut the cell, false otherwise. */
    inline bool r_ctest(double crs, double mrs) const { return crs + r_mul > sqrt(mrs * crs); }
    /** Scales a plane displacement during a plane bounds check.
     * \param[in] lrs the plane displacement.
     * \return The scaled value. */
    inline double r_cutoff(double lrs) const { return lrs * r_val; }
    /** Adds the maximum radius squared to a given value.
     * \param[in] rs the value to consider.
     * \return The value with the radius squared added. */
    inline double r_max_add(double rs) const { return rs + max_radius * max_radius; }
    /** Subtracts the radius squared of a particle from a given
     * value.
     * \param[in] rs the value to consider.
     * \param[in] ijk the block that the particle is within.
     * \param[in] q the index of the particle within the block.
     * \return The value with the radius squared subtracted. */
    inline double r_current_sub(double rs, int ijk, int q) const { return rs - ppr[ijk][4 * q + 3] * ppr[ijk][4 * q + 3]; }
    /** Scales a plane displacement prior to use in the plane cutting
     * algorithm.
     * \param[in] rs the initial plane displacement.
     * \param[in] ijk the block that the particle is within.
     * \param[in] q the index of the particle within the block.
     * \return The scaled plane displacement. */
    inline double r_scale(double rs, int ijk, int q) const { return rs + r_rad - ppr[ijk][4 * q + 3] * ppr[ijk][4 * q + 3]; }
    /** Scales a plane displacement prior to use in the plane
     * cutting algorithm, and also checks if it could possibly cut
     * the cell.
     * \param[in,out] rs the plane displacement to be scaled.
     * \param[in] mrs the current maximum distance to a Voronoi
     *                vertex multiplied by two.
     * \param[in] ijk the block that the particle is within.
     * \param[in] q the index of the particle within the block.
     * \return True if the cell could possibly cut the cell, false
     * otherwise. */
    inline bool r_scale_check(double &rs, double mrs, int ijk, int q) const {
      const double trs = rs;
      rs += r_rad - ppr[ijk][4 * q + 3] * ppr[ijk][4 * q + 3];
      return rs < sqrt(mrs * trs);
    }

  private:
    double r_rad, r_mul, r_val;
  };

} // namespace voro
#endif

// Voro++, a 3D cell-based Voronoi library
//
// Author   : Chris H. Rycroft (LBL / UC Berkeley)
// Email    : chr@alum.mit.edu
// Date     : August 30th 2011

/** \file container.hh
 * \brief Header file for the container_base and related classes. */

#ifndef VOROPP_CONTAINER_HH
#define VOROPP_CONTAINER_HH

#include <cstdio>
#include <vector>

// Voro++, a 3D cell-based Voronoi library
//
// Author   : Chris H. Rycroft (LBL / UC Berkeley)
// Email    : chr@alum.mit.edu
// Date     : August 30th 2011

/** \file c_loops.hh
 * \brief Header file for the loop classes. */

#ifndef VOROPP_C_LOOPS_HH
#define VOROPP_C_LOOPS_HH

namespace voro {

/** A type associated with a c_loop_subset class, determining what type of
   * geometrical region to loop over. */
  enum c_loop_subset_mode { sphere, box, no_check };

/** \brief A class for storing ordering information when particles are added to
   * a container.
   *
   * When particles are added to a container class, they are sorted into an
   * internal computational grid of blocks. The particle_order class provides a
   * mechanism for remembering which block particles were sorted into. The import
   * and put routines in the container class have variants that also take a
   * particle_order class. Each time they are called, they will store the block
   * that the particle was sorted into, plus the position of the particle within
   * the block. The particle_order class can used by the c_loop_order class to
   * specifically loop over the particles that have their information stored
   * within it. */
  class particle_order {
  public:
    /** A pointer to the array holding the ordering. */
    int *o;
    /** A pointer to the next position in the ordering array in
     * which to store an entry. */
    int *op;
    /** The current memory allocation for the class, set to the
     * number of entries which can be stored. */
    int size;
    /** The particle_order constructor allocates memory to store the
     * ordering information.
     * \param[in] init_size the initial amount of memory to
     *                      allocate. */
    particle_order(int init_size = init_ordering_size) : o(new int[init_size << 1]), op(o), size(init_size) {}
    /** The particle_order destructor frees the dynamically allocated
     * memory used to store the ordering information. */
    ~particle_order() { delete[] o; }
    /** Adds a record to the order, corresponding to the memory
     * address of where a particle was placed into the container.
     * \param[in] ijk the block into which the particle was placed.
     * \param[in] q the position within the block where the
     * 		particle was placed. */
    inline void add(int ijk, int q) {
      if (op == o + size) {
        add_ordering_memory();
      }
      *(op++) = ijk;
      *(op++) = q;
    }

  private:
    void add_ordering_memory();
  };

/** \brief Base class for looping over particles in a container.
   *
   * This class forms the base of all classes that can loop over a subset of
   * particles in a contaner in some order. When initialized, it stores constants
   * about the corresponding container geometry. It also contains a number of
   * routines for interrogating which particle currently being considered by the
   * loop, which are common between all of the derived classes. */
  class c_loop_base {
  public:
    /** The number of blocks in the x direction. */
    const int nx;
    /** The number of blocks in the y direction. */
    const int ny;
    /** The number of blocks in the z direction. */
    const int nz;
    /** A constant, set to the value of nx multiplied by ny, which
     * is used in the routines that step through blocks in
     * sequence. */
    const int nxy;
    /** A constant, set to the value of nx*ny*nz, which is used in
     * the routines that step through blocks in sequence. */
    const int nxyz;
    /** The number of floating point numbers per particle in the
     * associated container data structure. */
    const int ps;
    /** A pointer to the particle position information in the
     * associated container data structure. */
    double **p;
    /** A pointer to the particle ID information in the associated
     * container data structure. */
    int **id;
    /** A pointer to the particle counts in the associated
     * container data structure. */
    int *co;
    /** The current x-index of the block under consideration by the
     * loop. */
    int i;
    /** The current y-index of the block under consideration by the
     * loop. */
    int j;
    /** The current z-index of the block under consideration by the
     * loop. */
    int k;
    /** The current index of the block under consideration by the
     * loop. */
    int ijk;
    /** The index of the particle under consideration within the current
     * block. */
    int q;
    /** The constructor copies several necessary constants from the
     * base container class.
     * \param[in] con the container class to use. */
    template <class c_class> c_loop_base(c_class &con) : nx(con.nx), ny(con.ny), nz(con.nz), nxy(con.nxy), nxyz(con.nxyz), ps(con.ps), p(con.p), id(con.id), co(con.co) {}
    /** Returns the position vector of the particle currently being
     * considered by the loop.
     * \param[out] (x,y,z) the position vector of the particle. */
    inline void pos(double &x, double &y, double &z) const {
      double *pp = p[ijk] + ps * q;
      x = *(pp++);
      y = *(pp++);
      z = *pp;
    }
    /** Returns the ID, position vector, and radius of the particle
     * currently being considered by the loop.
     * \param[out] pid the particle ID.
     * \param[out] (x,y,z) the position vector of the particle.
     * \param[out] r the radius of the particle. If no radius
     * 		 information is available the default radius
     * 		 value is returned. */
    inline void pos(int &pid, double &x, double &y, double &z, double &r) const {
      pid = id[ijk][q];
      double *pp = p[ijk] + ps * q;
      x = *(pp++);
      y = *(pp++);
      z = *pp;
      r = ps == 3 ? default_radius : *(++pp);
    }
    /** Returns the x position of the particle currently being
     * considered by the loop. */
    inline double x() const { return p[ijk][ps * q]; }
    /** Returns the y position of the particle currently being
     * considered by the loop. */
    inline double y() const { return p[ijk][ps * q + 1]; }
    /** Returns the z position of the particle currently being
     * considered by the loop. */
    inline double z() const { return p[ijk][ps * q + 2]; }
    /** Returns the ID of the particle currently being considered
     * by the loop. */
    inline int pid() const { return id[ijk][q]; }
  };

/** \brief Class for looping over all of the particles in a container.
   *
   * This is one of the simplest loop classes, that scans the computational
   * blocks in order, and scans all the particles within each block in order. */
  class c_loop_all : public c_loop_base {
  public:
    /** The constructor copies several necessary constants from the
     * base container class.
     * \param[in] con the container class to use. */
    template <class c_class> c_loop_all(c_class &con) : c_loop_base(con) {}
    /** Sets the class to consider the first particle.
     * \return True if there is any particle to consider, false
     * otherwise. */
    inline bool start() {
      i = j = k = ijk = q = 0;
      while (co[ijk] == 0) {
        if (!next_block()) {
          return false;
        }
      }
      return true;
    }
    /** Finds the next particle to test.
     * \return True if there is another particle, false if no more
     * particles are available. */
    inline bool inc() {
      q++;
      if (q >= co[ijk]) {
        q = 0;
        do {
          if (!next_block()) {
            return false;
          }
        } while (co[ijk] == 0);
      }
      return true;
    }

  private:
    /** Updates the internal variables to find the next
     * computational block with any particles.
     * \return True if another block is found, false if there are
     * no more blocks. */
    inline bool next_block() {
      ijk++;
      i++;
      if (i == nx) {
        i = 0;
        j++;
        if (j == ny) {
          j = 0;
          k++;
          if (ijk == nxyz) {
            return false;
          }
        }
      }
      return true;
    }
  };

/** \brief Class for looping over a subset of particles in a container.
   *
   * This class can loop over a subset of particles in a certain geometrical
   * region within the container. The class can be set up to loop over a
   * rectangular box or sphere. It can also rectangular group of internal
   * computational blocks. */
  class c_loop_subset : public c_loop_base {
  public:
    /** The current mode of operation, determining whether tests
     * should be applied to particles to ensure they are within a
     * certain geometrical object. */
    c_loop_subset_mode mode;
    /** The constructor copies several necessary constants from the
     * base container class.
     * \param[in] con the container class to use. */
    template <class c_class>
    c_loop_subset(c_class &con) :
        c_loop_base(con), ax(con.ax), ay(con.ay), az(con.az), sx(con.bx - ax), sy(con.by - ay), sz(con.bz - az), xsp(con.xsp), ysp(con.ysp), zsp(con.zsp), xperiodic(con.xperiodic), yperiodic(con.yperiodic), zperiodic(con.zperiodic) {}
    void setup_sphere(double vx, double vy, double vz, double r, bool bounds_test = true);
    void setup_box(double xmin, double xmax, double ymin, double ymax, double zmin, double zmax, bool bounds_test = true);
    void setup_intbox(int ai_, int bi_, int aj_, int bj_, int ak_, int bk_);
    bool start();
    /** Finds the next particle to test.
     * \return True if there is another particle, false if no more
     * particles are available. */
    inline bool inc() {
      do {
        q++;
        while (q >= co[ijk]) {
          q = 0;
          if (!next_block()) {
            return false;
          }
        }
      } while (mode != no_check && out_of_bounds());
      return true;
    }

  private:
    const double ax, ay, az, sx, sy, sz, xsp, ysp, zsp;
    const bool xperiodic, yperiodic, zperiodic;
    double px, py, pz, apx, apy, apz;
    double v0, v1, v2, v3, v4, v5;
    int ai, bi, aj, bj, ak, bk;
    int ci, cj, ck, di, dj, dk, inc1, inc2;
    inline int step_mod(int a, int b) { return a >= 0 ? a % b : b - 1 - (b - 1 - a) % b; }
    inline int step_div(int a, int b) { return a >= 0 ? a / b : -1 + (a + 1) / b; }
    inline int step_int(double a) { return a < 0 ? int(a) - 1 : int(a); }
    void setup_common();
    bool next_block();
    bool out_of_bounds();
  };

/** \brief Class for looping over all of the particles specified in a
   * pre-assembled particle_order class.
   *
   * The particle_order class can be used to create a specific order of particles
   * within the container. This class can then loop over these particles in this
   * order. The class is particularly useful in cases where the ordering of the
   * output must match the ordering of particles as they were inserted into the
   * container. */
  class c_loop_order : public c_loop_base {
  public:
    /** A reference to the ordering class to use. */
    particle_order &vo;
    /** A pointer to the current position in the ordering class. */
    int *cp;
    /** A pointer to the end position in the ordering class. */
    int *op;
    /** The constructor copies several necessary constants from the
     * base class, and sets up a reference to the ordering class to
     * use.
     * \param[in] con the container class to use.
     * \param[in] vo_ the ordering class to use. */
    template <class c_class> c_loop_order(c_class &con, particle_order &vo_) : c_loop_base(con), vo(vo_), nx(con.nx), nxy(con.nxy) {}
    /** Sets the class to consider the first particle.
     * \return True if there is any particle to consider, false
     * otherwise. */
    inline bool start() {
      cp = vo.o;
      op = vo.op;
      if (cp != op) {
        ijk = *(cp++);
        decode();
        q = *(cp++);
        return true;
      } else {
        return false;
      }
    }
    /** Finds the next particle to test.
     * \return True if there is another particle, false if no more
     * particles are available. */
    inline bool inc() {
      if (cp == op) {
        return false;
      }
      ijk = *(cp++);
      decode();
      q = *(cp++);
      return true;
    }

  private:
    /** The number of computational blocks in the x direction. */
    const int nx;
    /** The number of computational blocks in a z-slice. */
    const int nxy;
    /** Takes the current block index and computes indices in the
     * x, y, and z directions. */
    inline void decode() {
      k = ijk / nxy;
      const int ijkt = ijk - nxy * k;
      j = ijkt / nx;
      i = ijkt - j * nx;
    }
  };

/** \brief A class for looping over all particles in a container_periodic or
   * container_periodic_poly class.
   *
   * Since the container_periodic and container_periodic_poly classes have a
   * fundamentally different memory organization, the regular loop classes cannot
   * be used with them. */
  class c_loop_all_periodic : public c_loop_base {
  public:
    /** The constructor copies several necessary constants from the
     * base periodic container class.
     * \param[in] con the periodic container class to use. */
    template <class c_class> c_loop_all_periodic(c_class &con) : c_loop_base(con), ey(con.ey), ez(con.ez), wy(con.wy), wz(con.wz), ijk0(nx * (ey + con.oy * ez)), inc2(2 * nx * con.ey + 1) {}
    /** Sets the class to consider the first particle.
     * \return True if there is any particle to consider, false
     * otherwise. */
    inline bool start() {
      i = 0;
      j = ey;
      k = ez;
      ijk = ijk0;
      q = 0;
      while (co[ijk] == 0) {
        if (!next_block()) {
          return false;
        }
      }
      return true;
    }
    /** Finds the next particle to test.
     * \return True if there is another particle, false if no more
     * particles are available. */
    inline bool inc() {
      q++;
      if (q >= co[ijk]) {
        q = 0;
        do {
          if (!next_block()) {
            return false;
          }
        } while (co[ijk] == 0);
      }
      return true;
    }

  private:
    /** The lower y index (inclusive) of the primary domain within
     * the block structure. */
    int ey;
    /** The lower y index (inclusive) of the primary domain within
     * the block structure. */
    int ez;
    /** The upper y index (exclusive) of the primary domain within
     * the block structure. */
    int wy;
    /** The upper z index (exclusive) of the primary domain within
     * the block structure. */
    int wz;
    /** The index of the (0,0,0) block within the block structure.
     */
    int ijk0;
    /** A value to increase ijk by when the z index is increased.
     */
    int inc2;
    /** Updates the internal variables to find the next
     * computational block with any particles.
     * \return True if another block is found, false if there are
     * no more blocks. */
    inline bool next_block() {
      i++;
      if (i == nx) {
        i = 0;
        j++;
        if (j == wy) {
          j = ey;
          k++;
          if (k == wz) {
            return false;
          }
          ijk += inc2;
        } else {
          ijk++;
        }
      } else {
        ijk++;
      }
      return true;
    }
  };

/** \brief Class for looping over all of the particles specified in a
   * pre-assembled particle_order class, for use with container_periodic classes.
   *
   * The particle_order class can be used to create a specific order of particles
   * within the container. This class can then loop over these particles in this
   * order. The class is particularly useful in cases where the ordering of the
   * output must match the ordering of particles as they were inserted into the
   * container. */
  class c_loop_order_periodic : public c_loop_base {
  public:
    /** A reference to the ordering class to use. */
    particle_order &vo;
    /** A pointer to the current position in the ordering class. */
    int *cp;
    /** A pointer to the end position in the ordering class. */
    int *op;
    /** The constructor copies several necessary constants from the
     * base class, and sets up a reference to the ordering class to
     * use.
     * \param[in] con the container class to use.
     * \param[in] vo_ the ordering class to use. */
    template <class c_class> c_loop_order_periodic(c_class &con, particle_order &vo_) : c_loop_base(con), vo(vo_), nx(con.nx), oxy(con.nx * con.oy) {}
    /** Sets the class to consider the first particle.
     * \return True if there is any particle to consider, false
     * otherwise. */
    inline bool start() {
      cp = vo.o;
      op = vo.op;
      if (cp != op) {
        ijk = *(cp++);
        decode();
        q = *(cp++);
        return true;
      } else {
        return false;
      }
    }
    /** Finds the next particle to test.
     * \return True if there is another particle, false if no more
     * particles are available. */
    inline bool inc() {
      if (cp == op) {
        return false;
      }
      ijk = *(cp++);
      decode();
      q = *(cp++);
      return true;
    }

  private:
    /** The number of computational blocks in the x direction. */
    const int nx;
    /** The number of computational blocks in a z-slice. */
    const int oxy;
    /** Takes the current block index and computes indices in the
     * x, y, and z directions. */
    inline void decode() {
      k = ijk / oxy;
      const int ijkt = ijk - oxy * k;
      j = ijkt / nx;
      i = ijkt - j * nx;
    }
  };

} // namespace voro

#endif

// Voro++, a 3D cell-based Voronoi library
//
// Author   : Chris H. Rycroft (LBL / UC Berkeley)
// Email    : chr@alum.mit.edu
// Date     : August 30th 2011

/** \file v_compute.hh
 * \brief Header file for the voro_compute template and related classes. */

#ifndef VOROPP_V_COMPUTE_HH
#define VOROPP_V_COMPUTE_HH

namespace voro {

/** \brief Structure for holding information about a particle.
   *
   * This small structure holds information about a single particle, and is used
   * by several of the routines in the voro_compute template for passing
   * information by reference between functions. */
  struct particle_record {
    /** The index of the block that the particle is within. */
    int ijk;
    /** The number of particle within its block. */
    int l;
    /** The x-index of the block. */
    int di;
    /** The y-index of the block. */
    int dj;
    /** The z-index of the block. */
    int dk;
  };

/** \brief Template for carrying out Voronoi cell computations. */
  template <class c_class> class voro_compute {
  public:
    /** A reference to the container class on which to carry out*/
    c_class &con;
    /** The size of an internal computational block in the x
     * direction. */
    const double boxx;
    /** The size of an internal computational block in the y
     * direction. */
    const double boxy;
    /** The size of an internal computational block in the z
     * direction. */
    const double boxz;
    /** The inverse box length in the x direction, set to
     * nx/(bx-ax). */
    const double xsp;
    /** The inverse box length in the y direction, set to
     * ny/(by-ay). */
    const double ysp;
    /** The inverse box length in the z direction, set to
     * nz/(bz-az). */
    const double zsp;
    /** The number of boxes in the x direction for the searching mask. */
    const int hx;
    /** The number of boxes in the y direction for the searching mask. */
    const int hy;
    /** The number of boxes in the z direction for the searching mask. */
    const int hz;
    /** A constant, set to the value of hx multiplied by hy, which
     * is used in the routines which step through mask boxes in
     * sequence. */
    const int hxy;
    /** A constant, set to the value of hx*hy*hz, which is used in
     * the routines which step through mask boxes in sequence. */
    const int hxyz;
    /** The number of floating point entries to store for each
     * particle. */
    const int ps;
    /** This array holds the numerical IDs of each particle in each
     * computational box. */
    int **id;
    /** A two dimensional array holding particle positions. For the
     * derived container_poly class, this also holds particle
     * radii. */
    double **p;
    /** An array holding the number of particles within each
     * computational box of the container. */
    int *co;
    voro_compute(c_class &con_, int hx_, int hy_, int hz_);
    /** The class destructor frees the dynamically allocated memory
     * for the mask and queue. */
    ~voro_compute() {
      delete[] qu;
      delete[] mask;
    }
    template <class v_cell> bool compute_cell(v_cell &c, int ijk, int s, int ci, int cj, int ck);
    void find_voronoi_cell(double x, double y, double z, int ci, int cj, int ck, int ijk, particle_record &w, double &mrs);

  private:
    /** A constant set to boxx*boxx+boxy*boxy+boxz*boxz, which is
     * frequently used in the computation. */
    const double bxsq;
    /** This sets the current value being used to mark tested blocks
     * in the mask. */
    unsigned int mv;
    /** The current size of the search list. */
    int qu_size;
    /** A pointer to the array of worklists. */
    const unsigned int *wl;
    /** An pointer to the array holding the minimum distances
     * associated with the worklists. */
    double *mrad;
    /** This array is used during the cell computation to determine
     * which blocks have been considered. */
    unsigned int *mask;
    /** An array is used to store the queue of blocks to test
     * during the Voronoi cell computation. */
    int *qu;
    /** A pointer to the end of the queue array, used to determine
     * when the queue is full. */
    int *qu_l;
    template <class v_cell> bool corner_test(v_cell &c, double xl, double yl, double zl, double xh, double yh, double zh);
    template <class v_cell> inline bool edge_x_test(v_cell &c, double x0, double yl, double zl, double x1, double yh, double zh);
    template <class v_cell> inline bool edge_y_test(v_cell &c, double xl, double y0, double zl, double xh, double y1, double zh);
    template <class v_cell> inline bool edge_z_test(v_cell &c, double xl, double yl, double z0, double xh, double yh, double z1);
    template <class v_cell> inline bool face_x_test(v_cell &c, double xl, double y0, double z0, double y1, double z1);
    template <class v_cell> inline bool face_y_test(v_cell &c, double x0, double yl, double z0, double x1, double z1);
    template <class v_cell> inline bool face_z_test(v_cell &c, double x0, double y0, double zl, double x1, double y1);
    bool compute_min_max_radius(int di, int dj, int dk, double fx, double fy, double fz, double gx, double gy, double gz, double &crs, double mrs);
    bool compute_min_radius(int di, int dj, int dk, double fx, double fy, double fz, double mrs);
    inline void add_to_mask(int ei, int ej, int ek, int *&qu_e);
    inline void scan_bits_mask_add(unsigned int q, unsigned int *mijk, int ei, int ej, int ek, int *&qu_e);
    inline void scan_all(int ijk, double x, double y, double z, int di, int dj, int dk, particle_record &w, double &mrs);
    void add_list_memory(int *&qu_s, int *&qu_e);
    /** Resets the mask in cases where the mask counter wraps
     * around. */
    inline void reset_mask() {
      for (unsigned int *mp(mask); mp < mask + hxyz; mp++) {
        *mp = 0;
      }
    }
  };

} // namespace voro

#endif

namespace voro {

/** \brief Pure virtual class from which wall objects are derived.
   *
   * This is a pure virtual class for a generic wall object. A wall object
   * can be specified by deriving a new class from this and specifying the
   * functions.*/
  class wall {
  public:
    virtual ~wall() {}
    /** A pure virtual function for testing whether a point is
     * inside the wall object. */
    virtual bool point_inside(double x, double y, double z) = 0;
    /** A pure virtual function for cutting a cell without
     * neighbor-tracking with a wall. */
    virtual bool cut_cell(voronoicell &c, double x, double y, double z) = 0;
    /** A pure virtual function for cutting a cell with
     * neighbor-tracking enabled with a wall. */
    virtual bool cut_cell(voronoicell_neighbor &c, double x, double y, double z) = 0;
  };

/** \brief A class for storing a list of pointers to walls.
   *
   * This class stores a list of pointers to wall classes. It contains several
   * simple routines that make use of the wall classes (such as telling whether a
   * given position is inside all of the walls or not). It can be used by itself,
   * but also forms part of container_base, for associating walls with this
   * class. */
  class wall_list {
  public:
    /** An array holding pointers to wall objects. */
    wall **walls;
    /** A pointer to the next free position to add a wall pointer.
     */
    wall **wep;
    wall_list();
    ~wall_list();
    /** Adds a wall to the list.
     * \param[in] w the wall to add. */
    inline void add_wall(wall *w) {
      if (wep == wel) {
        increase_wall_memory();
      }
      *(wep++) = w;
    }
    /** Adds a wall to the list.
     * \param[in] w a reference to the wall to add. */
    inline void add_wall(wall &w) { add_wall(&w); }
    void add_wall(wall_list &wl);
    /** Determines whether a given position is inside all of the
     * walls on the list.
     * \param[in] (x,y,z) the position to test.
     * \return True if it is inside, false if it is outside. */
    inline bool point_inside_walls(double x, double y, double z) const {
      for (wall **wp = walls; wp < wep; wp++) {
        if (!((*wp)->point_inside(x, y, z))) {
          return false;
        }
      }
      return true;
    }
    /** Cuts a Voronoi cell by all of the walls currently on
     * the list.
     * \param[in] c a reference to the Voronoi cell class.
     * \param[in] (x,y,z) the position of the cell.
     * \return True if the cell still exists, false if the cell is
     * deleted. */
    template <class c_class> bool apply_walls(c_class &c, double x, double y, double z) {
      for (wall **wp = walls; wp < wep; wp++) {
        if (!((*wp)->cut_cell(c, x, y, z))) {
          return false;
        }
      }
      return true;
    }
    void deallocate() const;

  protected:
    void increase_wall_memory();
    /** A pointer to the limit of the walls array, used to
     * determine when array is full. */
    wall **wel;
    /** The current amount of memory allocated for walls. */
    int current_wall_size;
  };

/** \brief Class for representing a particle system in a three-dimensional
   * rectangular box.
   *
   * This class represents a system of particles in a three-dimensional
   * rectangular box. Any combination of non-periodic and periodic coordinates
   * can be used in the three coordinate directions. The class is not intended
   * for direct use, but instead forms the base of the container and
   * container_poly classes that add specialized routines for computing the
   * regular and radical Voronoi tessellations respectively. It contains routines
   * that are commonly between these two classes, such as those for drawing the
   * domain, and placing particles within the internal data structure.
   *
   * The class is derived from the wall_list class, which encapsulates routines
   * for associating walls with the container, and the voro_base class, which
   * encapsulates routines about the underlying computational grid. */
  class container_base : public voro_base, public wall_list {
  public:
    /** The minimum x coordinate of the container. */
    const double ax;
    /** The maximum x coordinate of the container. */
    const double bx;
    /** The minimum y coordinate of the container. */
    const double ay;
    /** The maximum y coordinate of the container. */
    const double by;
    /** The minimum z coordinate of the container. */
    const double az;
    /** The maximum z coordinate of the container. */
    const double bz;
    /** A boolean value that determines if the x coordinate in
     * periodic or not. */
    const bool xperiodic;
    /** A boolean value that determines if the y coordinate in
     * periodic or not. */
    const bool yperiodic;
    /** A boolean value that determines if the z coordinate in
     * periodic or not. */
    const bool zperiodic;
    /** This array holds the numerical IDs of each particle in each
     * computational box. */
    int **id;
    /** A two dimensional array holding particle positions. For the
     * derived container_poly class, this also holds particle
     * radii. */
    double **p;
    /** This array holds the number of particles within each
     * computational box of the container. */
    int *co;
    /** This array holds the maximum amount of particle memory for
     * each computational box of the container. If the number of
     * particles in a particular box ever approaches this limit,
     * more is allocated using the add_particle_memory() function.
     */
    int *mem;
    /** The amount of memory in the array structure for each
     * particle. This is set to 3 when the basic class is
     * initialized, so that the array holds (x,y,z) positions. If
     * the container class is initialized as part of the derived
     * class container_poly, then this is set to 4, to also hold
     * the particle radii. */
    const int ps;
    container_base(double ax_, double bx_, double ay_, double by_, double az_, double bz_, int nx_, int ny_, int nz_, bool xperiodic_, bool yperiodic_, bool zperiodic_, int init_mem, int ps_);
    ~container_base();
    bool point_inside(double x, double y, double z);
    void region_count();
    /** Initializes the Voronoi cell prior to a compute_cell
     * operation for a specific particle being carried out by a
     * voro_compute class. The cell is initialized to fill the
     * entire container. For non-periodic coordinates, this is set
     * by the position of the walls. For periodic coordinates, the
     * space is equally divided in either direction from the
     * particle's initial position. Plane cuts made by any walls
     * that have been added are then applied to the cell.
     * \param[in,out] c a reference to a voronoicell object.
     * \param[in] ijk the block that the particle is within.
     * \param[in] q the index of the particle within its block.
     * \param[in] (ci,cj,ck) the coordinates of the block in the
     * 			 container coordinate system.
     * \param[out] (i,j,k) the coordinates of the test block
     * 		       relative to the voro_compute
     * 		       coordinate system.
     * \param[out] (x,y,z) the position of the particle.
     * \param[out] disp a block displacement used internally by the
     *		    compute_cell routine.
     * \return False if the plane cuts applied by walls completely
     * removed the cell, true otherwise. */
    template <class v_cell> inline bool initialize_voronoicell(v_cell &c, int ijk, int q, int ci, int cj, int ck, int &i, int &j, int &k, double &x, double &y, double &z, int &disp) {
      double x1;
      double x2;
      double y1;
      double y2;
      double z1;
      double z2;
      double *pp = p[ijk] + ps * q;
      x = *(pp++);
      y = *(pp++);
      z = *pp;
      if (xperiodic) {
        x1 = -(x2 = 0.5 * (bx - ax));
        i = nx;
      } else {
        x1 = ax - x;
        x2 = bx - x;
        i = ci;
      }
      if (yperiodic) {
        y1 = -(y2 = 0.5 * (by - ay));
        j = ny;
      } else {
        y1 = ay - y;
        y2 = by - y;
        j = cj;
      }
      if (zperiodic) {
        z1 = -(z2 = 0.5 * (bz - az));
        k = nz;
      } else {
        z1 = az - z;
        z2 = bz - z;
        k = ck;
      }
      c.init(x1, x2, y1, y2, z1, z2);
      if (!apply_walls(c, x, y, z)) {
        return false;
      }
      disp = ijk - i - nx * (j + ny * k);
      return true;
    }
    /** Initializes parameters for a find_voronoi_cell call within
     * the voro_compute template.
     * \param[in] (ci,cj,ck) the coordinates of the test block in
     * 			 the container coordinate system.
     * \param[in] ijk the index of the test block
     * \param[out] (i,j,k) the coordinates of the test block
     * 		       relative to the voro_compute
     * 		       coordinate system.
     * \param[out] disp a block displacement used internally by the
     *		    find_voronoi_cell routine. */
    inline void initialize_search(int ci, int cj, int ck, int ijk, int &i, int &j, int &k, int &disp) {
      i = xperiodic ? nx : ci;
      j = yperiodic ? ny : cj;
      k = zperiodic ? nz : ck;
      disp = ijk - i - nx * (j + ny * k);
    }
    /** Returns the position of a particle currently being computed
     * relative to the computational block that it is within. It is
     * used to select the optimal worklist entry to use.
     * \param[in] (x,y,z) the position of the particle.
     * \param[in] (ci,cj,ck) the block that the particle is within.
     * \param[out] (fx,fy,fz) the position relative to the block.
     */
    inline void frac_pos(double x, double y, double z, double ci, double cj, double ck, double &fx, double &fy, double &fz) {
      fx = x - ax - boxx * ci;
      fy = y - ay - boxy * cj;
      fz = z - az - boxz * ck;
    }
    /** Calculates the index of block in the container structure
     * corresponding to given coordinates.
     * \param[in] (ci,cj,ck) the coordinates of the original block
     * 			 in the current computation, relative
     * 			 to the container coordinate system.
     * \param[in] (ei,ej,ek) the displacement of the current block
     * 			 from the original block.
     * \param[in,out] (qx,qy,qz) the periodic displacement that
     * 			     must be added to the particles
     * 			     within the computed block.
     * \param[in] disp a block displacement used internally by the
     * 		    find_voronoi_cell and compute_cell routines.
     * \return The block index. */
    inline int region_index(int ci, int cj, int ck, int ei, int ej, int ek, double &qx, double &qy, double &qz, int &disp) {
      if (xperiodic) {
        if (ci + ei < nx) {
          ei += nx;
          qx = -(bx - ax);
        } else if (ci + ei >= (nx << 1)) {
          ei -= nx;
          qx = bx - ax;
        } else {
          qx = 0;
        }
      }
      if (yperiodic) {
        if (cj + ej < ny) {
          ej += ny;
          qy = -(by - ay);
        } else if (cj + ej >= (ny << 1)) {
          ej -= ny;
          qy = by - ay;
        } else {
          qy = 0;
        }
      }
      if (zperiodic) {
        if (ck + ek < nz) {
          ek += nz;
          qz = -(bz - az);
        } else if (ck + ek >= (nz << 1)) {
          ek -= nz;
          qz = bz - az;
        } else {
          qz = 0;
        }
      }
      return disp + ei + nx * (ej + ny * ek);
    }
    void draw_domain_gnuplot(FILE *fp = stdout) const;
    /** Draws an outline of the domain in Gnuplot format.
     * \param[in] filename the filename to write to. */
    inline void draw_domain_gnuplot(const char *filename) const {
      FILE *fp = safe_fopen(filename, "w");
      draw_domain_gnuplot(fp);
      fclose(fp);
    }
    void draw_domain_pov(FILE *fp = stdout) const;
    /** Draws an outline of the domain in Gnuplot format.
     * \param[in] filename the filename to write to. */
    inline void draw_domain_pov(const char *filename) const {
      FILE *fp = safe_fopen(filename, "w");
      draw_domain_pov(fp);
      fclose(fp);
    }
    /** Sums up the total number of stored particles.
     * \return The number of particles. */
    inline int total_particles() {
      int tp = *co;
      for (int *cop = co + 1; cop < co + nxyz; cop++) {
        tp += *cop;
      }
      return tp;
    }

  protected:
    void add_particle_memory(int i) const;
    bool put_locate_block(int &ijk, double &x, double &y, double &z);
    inline bool put_remap(int &ijk, double &x, double &y, double &z);
    inline bool remap(int &ai, int &aj, int &ak, int &ci, int &cj, int &ck, double &x, double &y, double &z, int &ijk);
  };

/** \brief Extension of the container_base class for computing regular Voronoi
   * tessellations.
   *
   * This class is an extension of the container_base class that has routines
   * specifically for computing the regular Voronoi tessellation with no
   * dependence on particle radii. */
  class container : public container_base, public radius_mono {
  public:
    container(double ax_, double bx_, double ay_, double by_, double az_, double bz_, int nx_, int ny_, int nz_, bool xperiodic_, bool yperiodic_, bool zperiodic_, int init_mem);
    void clear();
    void put(int n, double x, double y, double z);
    void put(particle_order &vo, int n, double x, double y, double z);
    void import(FILE *fp = stdin);
    void import(particle_order &vo, FILE *fp = stdin);
    /** Imports a list of particles from an open file stream into
     * the container. Entries of four numbers (Particle ID, x
     * position, y position, z position) are searched for. If the
     * file cannot be successfully read, then the routine causes a
     * fatal error.
     * \param[in] filename the name of the file to open and read
     *                     from. */
    inline void import(const char *filename) {
      FILE *fp = safe_fopen(filename, "r");
      import(fp);
      fclose(fp);
    }
    /** Imports a list of particles from an open file stream into
     * the container. Entries of four numbers (Particle ID, x
     * position, y position, z position) are searched for. In
     * addition, the order in which particles are read is saved
     * into an ordering class. If the file cannot be successfully
     * read, then the routine causes a fatal error.
     * \param[in,out] vo the ordering class to use.
     * \param[in] filename the name of the file to open and read
     *                     from. */
    inline void import(particle_order &vo, const char *filename) {
      FILE *fp = safe_fopen(filename, "r");
      import(vo, fp);
      fclose(fp);
    }
    void compute_all_cells();
    double sum_cell_volumes();
    /** Dumps particle IDs and positions to a file.
     * \param[in] vl the loop class to use.
     * \param[in] fp a file handle to write to. */
    template <class c_loop> void draw_particles(c_loop &vl, FILE *fp) {
      double *pp;
      if (vl.start()) {
        do {
          pp = p[vl.ijk] + 3 * vl.q;
          fprintf(fp, "%d %g %g %g\n", id[vl.ijk][vl.q], *pp, pp[1], pp[2]);
        } while (vl.inc());
      }
    }
    /** Dumps all of the particle IDs and positions to a file.
     * \param[in] fp a file handle to write to. */
    inline void draw_particles(FILE *fp = stdout) {
      c_loop_all vl(*this);
      draw_particles(vl, fp);
    }
    /** Dumps all of the particle IDs and positions to a file.
     * \param[in] filename the name of the file to write to. */
    inline void draw_particles(const char *filename) {
      FILE *fp = safe_fopen(filename, "w");
      draw_particles(fp);
      fclose(fp);
    }
    /** Dumps particle positions in POV-Ray format.
     * \param[in] vl the loop class to use.
     * \param[in] fp a file handle to write to. */
    template <class c_loop> void draw_particles_pov(c_loop &vl, FILE *fp) {
      double *pp;
      if (vl.start()) {
        do {
          pp = p[vl.ijk] + 3 * vl.q;
          fprintf(fp, "// id %d\nsphere{<%g,%g,%g>,s}\n", id[vl.ijk][vl.q], *pp, pp[1], pp[2]);
        } while (vl.inc());
      }
    }
    /** Dumps all particle positions in POV-Ray format.
     * \param[in] fp a file handle to write to. */
    inline void draw_particles_pov(FILE *fp = stdout) {
      c_loop_all vl(*this);
      draw_particles_pov(vl, fp);
    }
    /** Dumps all particle positions in POV-Ray format.
     * \param[in] filename the name of the file to write to. */
    inline void draw_particles_pov(const char *filename) {
      FILE *fp = safe_fopen(filename, "w");
      draw_particles_pov(fp);
      fclose(fp);
    }
    /** Computes Voronoi cells and saves the output in gnuplot
     * format.
     * \param[in] vl the loop class to use.
     * \param[in] fp a file handle to write to. */
    template <class c_loop> void draw_cells_gnuplot(c_loop &vl, FILE *fp) {
      voronoicell c;
      double *pp;
      if (vl.start()) {
        do {
          if (compute_cell(c, vl)) {
            pp = p[vl.ijk] + ps * vl.q;
            c.draw_gnuplot(*pp, pp[1], pp[2], fp);
          }
        } while (vl.inc());
      }
    }
    /** Computes all Voronoi cells and saves the output in gnuplot
     * format.
     * \param[in] fp a file handle to write to. */
    inline void draw_cells_gnuplot(FILE *fp = stdout) {
      c_loop_all vl(*this);
      draw_cells_gnuplot(vl, fp);
    }
    /** Computes all Voronoi cells and saves the output in gnuplot
     * format.
     * \param[in] filename the name of the file to write to. */
    inline void draw_cells_gnuplot(const char *filename) {
      FILE *fp = safe_fopen(filename, "w");
      draw_cells_gnuplot(fp);
      fclose(fp);
    }
    /** Computes Voronoi cells and saves the output in POV-Ray
     * format.
     * \param[in] vl the loop class to use.
     * \param[in] fp a file handle to write to. */
    template <class c_loop> void draw_cells_pov(c_loop &vl, FILE *fp) {
      voronoicell c;
      double *pp;
      if (vl.start()) {
        do {
          if (compute_cell(c, vl)) {
            fprintf(fp, "// cell %d\n", id[vl.ijk][vl.q]);
            pp = p[vl.ijk] + ps * vl.q;
            c.draw_pov(*pp, pp[1], pp[2], fp);
          }
        } while (vl.inc());
      }
    }
    /** Computes all Voronoi cells and saves the output in POV-Ray
     * format.
     * \param[in] fp a file handle to write to. */
    inline void draw_cells_pov(FILE *fp = stdout) {
      c_loop_all vl(*this);
      draw_cells_pov(vl, fp);
    }
    /** Computes all Voronoi cells and saves the output in POV-Ray
     * format.
     * \param[in] filename the name of the file to write to. */
    inline void draw_cells_pov(const char *filename) {
      FILE *fp = safe_fopen(filename, "w");
      draw_cells_pov(fp);
      fclose(fp);
    }
    /** Computes the Voronoi cells and saves customized information
     * about them.
     * \param[in] vl the loop class to use.
     * \param[in] format the custom output string to use.
     * \param[in] fp a file handle to write to. */
    template <class c_loop> void print_custom(c_loop &vl, const char *format, FILE *fp) {
      int ijk;
      int q;
      double *pp;
      if (contains_neighbor(format)) {
        voronoicell_neighbor c;
        if (vl.start()) {
          do {
            if (compute_cell(c, vl)) {
              ijk = vl.ijk;
              q = vl.q;
              pp = p[ijk] + ps * q;
              c.output_custom(format, id[ijk][q], *pp, pp[1], pp[2], default_radius, fp);
            }
          } while (vl.inc());
        }
      } else {
        voronoicell c;
        if (vl.start()) {
          do {
            if (compute_cell(c, vl)) {
              ijk = vl.ijk;
              q = vl.q;
              pp = p[ijk] + ps * q;
              c.output_custom(format, id[ijk][q], *pp, pp[1], pp[2], default_radius, fp);
            }
          } while (vl.inc());
        }
      }
    }
    void print_custom(const char *format, FILE *fp = stdout);
    void print_custom(const char *format, const char *filename);
    bool find_voronoi_cell(double x, double y, double z, double &rx, double &ry, double &rz, int &pid);
    /** Computes the Voronoi cell for a particle currently being
     * referenced by a loop class.
     * \param[out] c a Voronoi cell class in which to store the
     * 		 computed cell.
     * \param[in] vl the loop class to use.
     * \return True if the cell was computed. If the cell cannot be
     * computed, if it is removed entirely by a wall or boundary
     * condition, then the routine returns false. */
    template <class v_cell, class c_loop> inline bool compute_cell(v_cell &c, c_loop &vl) { return vc.compute_cell(c, vl.ijk, vl.q, vl.i, vl.j, vl.k); }
    /** Computes the Voronoi cell for given particle.
     * \param[out] c a Voronoi cell class in which to store the
     * 		 computed cell.
     * \param[in] ijk the block that the particle is within.
     * \param[in] q the index of the particle within the block.
     * \return True if the cell was computed. If the cell cannot be
     * computed, if it is removed entirely by a wall or boundary
     * condition, then the routine returns false. */
    template <class v_cell> inline bool compute_cell(v_cell &c, int ijk, int q) {
      int k = ijk / nxy;
      int ijkt = ijk - nxy * k;
      int j = ijkt / nx;
      int i = ijkt - j * nx;
      return vc.compute_cell(c, ijk, q, i, j, k);
    }
    /** Computes the Voronoi cell for a ghost particle at a given
     * location.
     * \param[out] c a Voronoi cell class in which to store the
     * 		 computed cell.
     * \param[in] (x,y,z) the location of the ghost particle.
     * \return True if the cell was computed. If the cell cannot be
     * computed, if it is removed entirely by a wall or boundary
     * condition, then the routine returns false. */
    template <class v_cell> inline bool compute_ghost_cell(v_cell &c, double x, double y, double z) {
      int ijk;
      if (put_locate_block(ijk, x, y, z)) {
        double *pp = p[ijk] + 3 * co[ijk]++;
        *(pp++) = x;
        *(pp++) = y;
        *pp = z;
        bool q = compute_cell(c, ijk, co[ijk] - 1);
        co[ijk]--;
        return q;
      }
      return false;
    }

  private:
    voro_compute<container> vc;
    friend class voro_compute<container>;
  };

/** \brief Extension of the container_base class for computing radical Voronoi
   * tessellations.
   *
   * This class is an extension of container_base class that has routines
   * specifically for computing the radical Voronoi tessellation that depends on
   * the particle radii. */
  class container_poly : public container_base, public radius_poly {
  public:
    container_poly(double ax_, double bx_, double ay_, double by_, double az_, double bz_, int nx_, int ny_, int nz_, bool xperiodic_, bool yperiodic_, bool zperiodic_, int init_mem);
    void clear();
    void put(int n, double x, double y, double z, double r);
    void put(particle_order &vo, int n, double x, double y, double z, double r);
    void import(FILE *fp = stdin);
    void import(particle_order &vo, FILE *fp = stdin);
    /** Imports a list of particles from an open file stream into
     * the container_poly class. Entries of five numbers (Particle
     * ID, x position, y position, z position, radius) are searched
     * for. If the file cannot be successfully read, then the
     * routine causes a fatal error.
     * \param[in] filename the name of the file to open and read
     *                     from. */
    inline void import(const char *filename) {
      FILE *fp = safe_fopen(filename, "r");
      import(fp);
      fclose(fp);
    }
    /** Imports a list of particles from an open file stream into
     * the container_poly class. Entries of five numbers (Particle
     * ID, x position, y position, z position, radius) are searched
     * for. In addition, the order in which particles are read is
     * saved into an ordering class. If the file cannot be
     * successfully read, then the routine causes a fatal error.
     * \param[in,out] vo the ordering class to use.
     * \param[in] filename the name of the file to open and read
     *                     from. */
    inline void import(particle_order &vo, const char *filename) {
      FILE *fp = safe_fopen(filename, "r");
      import(vo, fp);
      fclose(fp);
    }
    void compute_all_cells();
    double sum_cell_volumes();
    /** Dumps particle IDs, positions and radii to a file.
     * \param[in] vl the loop class to use.
     * \param[in] fp a file handle to write to. */
    template <class c_loop> void draw_particles(c_loop &vl, FILE *fp) {
      double *pp;
      if (vl.start()) {
        do {
          pp = p[vl.ijk] + 4 * vl.q;
          fprintf(fp, "%d %g %g %g %g\n", id[vl.ijk][vl.q], *pp, pp[1], pp[2], pp[3]);
        } while (vl.inc());
      }
    }
    /** Dumps all of the particle IDs, positions and radii to a
     * file.
     * \param[in] fp a file handle to write to. */
    inline void draw_particles(FILE *fp = stdout) {
      c_loop_all vl(*this);
      draw_particles(vl, fp);
    }
    /** Dumps all of the particle IDs, positions and radii to a
     * file.
     * \param[in] filename the name of the file to write to. */
    inline void draw_particles(const char *filename) {
      FILE *fp = safe_fopen(filename, "w");
      draw_particles(fp);
      fclose(fp);
    }
    /** Dumps particle positions in POV-Ray format.
     * \param[in] vl the loop class to use.
     * \param[in] fp a file handle to write to. */
    template <class c_loop> void draw_particles_pov(c_loop &vl, FILE *fp) {
      double *pp;
      if (vl.start()) {
        do {
          pp = p[vl.ijk] + 4 * vl.q;
          fprintf(fp, "// id %d\nsphere{<%g,%g,%g>,%g}\n", id[vl.ijk][vl.q], *pp, pp[1], pp[2], pp[3]);
        } while (vl.inc());
      }
    }
    /** Dumps all the particle positions in POV-Ray format.
     * \param[in] fp a file handle to write to. */
    inline void draw_particles_pov(FILE *fp = stdout) {
      c_loop_all vl(*this);
      draw_particles_pov(vl, fp);
    }
    /** Dumps all the particle positions in POV-Ray format.
     * \param[in] filename the name of the file to write to. */
    inline void draw_particles_pov(const char *filename) {
      FILE *fp = safe_fopen(filename, "w");
      draw_particles_pov(fp);
      fclose(fp);
    }
    /** Computes Voronoi cells and saves the output in gnuplot
     * format.
     * \param[in] vl the loop class to use.
     * \param[in] fp a file handle to write to. */
    template <class c_loop> void draw_cells_gnuplot(c_loop &vl, FILE *fp) {
      voronoicell c;
      double *pp;
      if (vl.start()) {
        do {
          if (compute_cell(c, vl)) {
            pp = p[vl.ijk] + ps * vl.q;
            c.draw_gnuplot(*pp, pp[1], pp[2], fp);
          }
        } while (vl.inc());
      }
    }
    /** Compute all Voronoi cells and saves the output in gnuplot
     * format.
     * \param[in] fp a file handle to write to. */
    inline void draw_cells_gnuplot(FILE *fp = stdout) {
      c_loop_all vl(*this);
      draw_cells_gnuplot(vl, fp);
    }
    /** Compute all Voronoi cells and saves the output in gnuplot
     * format.
     * \param[in] filename the name of the file to write to. */
    inline void draw_cells_gnuplot(const char *filename) {
      FILE *fp = safe_fopen(filename, "w");
      draw_cells_gnuplot(fp);
      fclose(fp);
    }
    /** Computes Voronoi cells and saves the output in POV-Ray
     * format.
     * \param[in] vl the loop class to use.
     * \param[in] fp a file handle to write to. */
    template <class c_loop> void draw_cells_pov(c_loop &vl, FILE *fp) {
      voronoicell c;
      double *pp;
      if (vl.start()) {
        do {
          if (compute_cell(c, vl)) {
            fprintf(fp, "// cell %d\n", id[vl.ijk][vl.q]);
            pp = p[vl.ijk] + ps * vl.q;
            c.draw_pov(*pp, pp[1], pp[2], fp);
          }
        } while (vl.inc());
      }
    }
    /** Computes all Voronoi cells and saves the output in POV-Ray
     * format.
     * \param[in] fp a file handle to write to. */
    inline void draw_cells_pov(FILE *fp = stdout) {
      c_loop_all vl(*this);
      draw_cells_pov(vl, fp);
    }
    /** Computes all Voronoi cells and saves the output in POV-Ray
     * format.
     * \param[in] filename the name of the file to write to. */
    inline void draw_cells_pov(const char *filename) {
      FILE *fp = safe_fopen(filename, "w");
      draw_cells_pov(fp);
      fclose(fp);
    }
    /** Computes the Voronoi cells and saves customized information
     * about them.
     * \param[in] vl the loop class to use.
     * \param[in] format the custom output string to use.
     * \param[in] fp a file handle to write to. */
    template <class c_loop> void print_custom(c_loop &vl, const char *format, FILE *fp) {
      int ijk;
      int q;
      double *pp;
      if (contains_neighbor(format)) {
        voronoicell_neighbor c;
        if (vl.start()) {
          do {
            if (compute_cell(c, vl)) {
              ijk = vl.ijk;
              q = vl.q;
              pp = p[ijk] + ps * q;
              c.output_custom(format, id[ijk][q], *pp, pp[1], pp[2], pp[3], fp);
            }
          } while (vl.inc());
        }
      } else {
        voronoicell c;
        if (vl.start()) {
          do {
            if (compute_cell(c, vl)) {
              ijk = vl.ijk;
              q = vl.q;
              pp = p[ijk] + ps * q;
              c.output_custom(format, id[ijk][q], *pp, pp[1], pp[2], pp[3], fp);
            }
          } while (vl.inc());
        }
      }
    }
    /** Computes the Voronoi cell for a particle currently being
     * referenced by a loop class.
     * \param[out] c a Voronoi cell class in which to store the
     * 		 computed cell.
     * \param[in] vl the loop class to use.
     * \return True if the cell was computed. If the cell cannot be
     * computed, if it is removed entirely by a wall or boundary
     * condition, then the routine returns false. */
    template <class v_cell, class c_loop> inline bool compute_cell(v_cell &c, c_loop &vl) { return vc.compute_cell(c, vl.ijk, vl.q, vl.i, vl.j, vl.k); }
    /** Computes the Voronoi cell for given particle.
     * \param[out] c a Voronoi cell class in which to store the
     * 		 computed cell.
     * \param[in] ijk the block that the particle is within.
     * \param[in] q the index of the particle within the block.
     * \return True if the cell was computed. If the cell cannot be
     * computed, if it is removed entirely by a wall or boundary
     * condition, then the routine returns false. */
    template <class v_cell> inline bool compute_cell(v_cell &c, int ijk, int q) {
      int k = ijk / nxy;
      int ijkt = ijk - nxy * k;
      int j = ijkt / nx;
      int i = ijkt - j * nx;
      return vc.compute_cell(c, ijk, q, i, j, k);
    }
    /** Computes the Voronoi cell for a ghost particle at a given
     * location.
     * \param[out] c a Voronoi cell class in which to store the
     * 		 computed cell.
     * \param[in] (x,y,z) the location of the ghost particle.
     * \param[in] r the radius of the ghost particle.
     * \return True if the cell was computed. If the cell cannot be
     * computed, if it is removed entirely by a wall or boundary
     * condition, then the routine returns false. */
    template <class v_cell> inline bool compute_ghost_cell(v_cell &c, double x, double y, double z, double r) {
      int ijk;
      if (put_locate_block(ijk, x, y, z)) {
        double *pp = p[ijk] + 4 * co[ijk]++;
        double tm = max_radius;
        *(pp++) = x;
        *(pp++) = y;
        *(pp++) = z;
        *pp = r;
        if (r > max_radius) {
          max_radius = r;
        }
        bool q = compute_cell(c, ijk, co[ijk] - 1);
        co[ijk]--;
        max_radius = tm;
        return q;
      }
      return false;
    }
    void print_custom(const char *format, FILE *fp = stdout);
    void print_custom(const char *format, const char *filename);
    bool find_voronoi_cell(double x, double y, double z, double &rx, double &ry, double &rz, int &pid);

  private:
    voro_compute<container_poly> vc;
    friend class voro_compute<container_poly>;
  };

} // namespace voro

#endif

// Voro++, a 3D cell-based Voronoi library
//
// Author   : Chris H. Rycroft (LBL / UC Berkeley)
// Email    : chr@alum.mit.edu
// Date     : August 30th 2011

/** \file unitcell.hh
 * \brief Header file for the unitcell class. */

#ifndef VOROPP_UNITCELL_HH
#define VOROPP_UNITCELL_HH

#include <vector>

namespace voro {

/** \brief Class for computation of the unit Voronoi cell associated with
   * a 3D non-rectangular periodic domain. */
  class unitcell {
  public:
    /** The x coordinate of the first vector defining the periodic
     * domain. */
    const double bx;
    /** The x coordinate of the second vector defining the periodic
     * domain. */
    const double bxy;
    /** The y coordinate of the second vector defining the periodic
     * domain. */
    const double by;
    /** The x coordinate of the third vector defining the periodic
     * domain. */
    const double bxz;
    /** The y coordinate of the third vector defining the periodic
     * domain. */
    const double byz;
    /** The z coordinate of the third vector defining the periodic
     * domain. */
    const double bz;
    /** The computed unit Voronoi cell corresponding the given
     * 3D non-rectangular periodic domain geometry. */
    voronoicell unit_voro;
    unitcell(double bx_, double bxy_, double by_, double bxz_, double byz_, double bz_);
    /** Draws an outline of the domain in Gnuplot format.
     * \param[in] filename the filename to write to. */
    inline void draw_domain_gnuplot(const char *filename) const {
      FILE *fp(safe_fopen(filename, "w"));
      draw_domain_gnuplot(fp);
      fclose(fp);
    }
    void draw_domain_gnuplot(FILE *fp = stdout) const;
    /** Draws an outline of the domain in Gnuplot format.
     * \param[in] filename the filename to write to. */
    inline void draw_domain_pov(const char *filename) const {
      FILE *fp(safe_fopen(filename, "w"));
      draw_domain_pov(fp);
      fclose(fp);
    }
    void draw_domain_pov(FILE *fp = stdout) const;
    bool intersects_image(double dx, double dy, double dz, double &vol);
    void images(std::vector<int> &vi, std::vector<double> &vd);

  protected:
    /** The maximum y-coordinate that could possibly cut the
     * computed unit Voronoi cell. */
    double max_uv_y;
    /** The maximum z-coordinate that could possibly cut the
     * computed unit Voronoi cell. */
    double max_uv_z;

  private:
    inline void unit_voro_apply(int i, int j, int k);
    bool unit_voro_intersect(int l);
    inline bool unit_voro_test(int i, int j, int k);
  };

} // namespace voro

#endif

// Voro++, a 3D cell-based Voronoi library
//
// Author   : Chris H. Rycroft (LBL / UC Berkeley)
// Email    : chr@alum.mit.edu
// Date     : August 30th 2011

/** \file container_prd.hh
 * \brief Header file for the container_periodic_base and related classes. */

#ifndef VOROPP_CONTAINER_PRD_HH
#define VOROPP_CONTAINER_PRD_HH

#include <cstdio>
#include <vector>

namespace voro {

/** \brief Class for representing a particle system in a 3D periodic
   * non-orthogonal periodic domain.
   *
   * This class represents a particle system in a three-dimensional
   * non-orthogonal periodic domain. The domain is defined by three periodicity
   * vectors (bx,0,0), (bxy,by,0), and (bxz,byz,bz) that represent a
   * parallelepiped. Internally, the class stores particles in the box 0<x<bx,
   * 0<y<by, 0<z<bz, and constructs periodic images of particles that displaced
   * by the three periodicity vectors when they are necessary for the
   * computation. The internal memory structure for this class is significantly
   * different from the container_base class in order to handle the dynamic
   * construction of these periodic images.
   *
   * The class is derived from the unitcell class, which encapsulates information
   * about the domain geometry, and the voro_base class, which encapsulates
   * information about the underlying computational grid. */
  class container_periodic_base : public unitcell, public voro_base {
  public:
    /** The lower y index (inclusive) of the primary domain within
     * the block structure. */
    int ey;
    /** The lower z index (inclusive) of the primary domain within
     * the block structure. */
    int ez;
    /** The upper y index (exclusive) of the primary domain within
     * the block structure. */
    int wy;
    /** The upper z index (exclusive) of the primary domain within
     * the block structure. */
    int wz;
    /** The total size of the block structure (including images) in
     * the y direction. */
    int oy;
    /** The total size of the block structure (including images) in
     * the z direction. */
    int oz;
    /** The total number of blocks. */
    int oxyz;
    /** This array holds the numerical IDs of each particle in each
     * computational box. */
    int **id;
    /** A two dimensional array holding particle positions. For the
     * derived container_poly class, this also holds particle
     * radii. */
    double **p;
    /** This array holds the number of particles within each
     * computational box of the container. */
    int *co;
    /** This array holds the maximum amount of particle memory for
     * each computational box of the container. If the number of
     * particles in a particular box ever approaches this limit,
     * more is allocated using the add_particle_memory() function.
     */
    int *mem;
    /** An array holding information about periodic image
     * construction at a given location. */
    char *img;
    /** The initial amount of memory to allocate for particles
     * for each block. */
    const int init_mem;
    /** The amount of memory in the array structure for each
     * particle. This is set to 3 when the basic class is
     * initialized, so that the array holds (x,y,z) positions. If
     * the container class is initialized as part of the derived
     * class container_poly, then this is set to 4, to also hold
     * the particle radii. */
    const int ps;
    container_periodic_base(double bx_, double bxy_, double by_, double bxz_, double byz_, double bz_, int nx_, int ny_, int nz_, int init_mem_, int ps);
    ~container_periodic_base();
    /** Prints all particles in the container, including those that
     * have been constructed in image blocks. */
    inline void print_all_particles() const {
      int ijk;
      int q;
      for (ijk = 0; ijk < oxyz; ijk++) {
        for (q = 0; q < co[ijk]; q++) {
          printf("%d %g %g %g\n", id[ijk][q], p[ijk][ps * q], p[ijk][ps * q + 1], p[ijk][ps * q + 2]);
        }
      }
    }
    void region_count();
    /** Initializes the Voronoi cell prior to a compute_cell
     * operation for a specific particle being carried out by a
     * voro_compute class. The cell is initialized to be the
     * pre-computed unit Voronoi cell based on planes formed by
     * periodic images of the particle.
     * \param[in,out] c a reference to a voronoicell object.
     * \param[in] ijk the block that the particle is within.
     * \param[in] q the index of the particle within its block.
     * \param[in] (ci,cj,ck) the coordinates of the block in the
     * 			 container coordinate system.
     * \param[out] (i,j,k) the coordinates of the test block
     * 		       relative to the voro_compute
     * 		       coordinate system.
     * \param[out] (x,y,z) the position of the particle.
     * \param[out] disp a block displacement used internally by the
     *		    compute_cell routine.
     * \return False if the plane cuts applied by walls completely
     * removed the cell, true otherwise. */
    template <class v_cell> inline bool initialize_voronoicell(v_cell &c, int ijk, int q, int ci, int cj, int ck, int &i, int &j, int &k, double &x, double &y, double &z, int &disp) {
      c = unit_voro;
      double *pp = p[ijk] + ps * q;
      x = *(pp++);
      y = *(pp++);
      z = *pp;
      i = nx;
      j = ey;
      k = ez;
      return true;
    }
    /** Initializes parameters for a find_voronoi_cell call within
     * the voro_compute template.
     * \param[in] (ci,cj,ck) the coordinates of the test block in
     * 			 the container coordinate system.
     * \param[in] ijk the index of the test block
     * \param[out] (i,j,k) the coordinates of the test block
     * 		       relative to the voro_compute
     * 		       coordinate system.
     * \param[out] disp a block displacement used internally by the
     *		    find_voronoi_cell routine (but not needed
     *		    in this instance.) */
    inline void initialize_search(int ci, int cj, int ck, int ijk, int &i, int &j, int &k, int &disp) {
      i = nx;
      j = ey;
      k = ez;
    }
    /** Returns the position of a particle currently being computed
     * relative to the computational block that it is within. It is
     * used to select the optimal worklist entry to use.
     * \param[in] (x,y,z) the position of the particle.
     * \param[in] (ci,cj,ck) the block that the particle is within.
     * \param[out] (fx,fy,fz) the position relative to the block.
     */
    inline void frac_pos(double x, double y, double z, double ci, double cj, double ck, double &fx, double &fy, double &fz) {
      fx = x - boxx * ci;
      fy = y - boxy * (cj - ey);
      fz = z - boxz * (ck - ez);
    }
    /** Calculates the index of block in the container structure
     * corresponding to given coordinates.
     * \param[in] (ci,cj,ck) the coordinates of the original block
     * 			 in the current computation, relative
     * 			 to the container coordinate system.
     * \param[in] (ei,ej,ek) the displacement of the current block
     * 			 from the original block.
     * \param[in,out] (qx,qy,qz) the periodic displacement that
     * 			     must be added to the particles
     * 			     within the computed block.
     * \param[in] disp a block displacement used internally by the
     * 		    find_voronoi_cell and compute_cell routines
     * 		    (but not needed in this instance.)
     * \return The block index. */
    inline int region_index(int ci, int cj, int ck, int ei, int ej, int ek, double &qx, double &qy, double &qz, int &disp) {
      int qi = ci + (ei - nx);
      int qj = cj + (ej - ey);
      int qk = ck + (ek - ez);
      const int iv(step_div(qi, nx));
      if (iv != 0) {
        qx = iv * bx;
        qi -= nx * iv;
      } else {
        qx = 0;
      }
      create_periodic_image(qi, qj, qk);
      return qi + nx * (qj + oy * qk);
    }
    void create_all_images();
    void check_compartmentalized();

  protected:
    void add_particle_memory(int i) const;
    void put_locate_block(int &ijk, double &x, double &y, double &z);
    void put_locate_block(int &ijk, double &x, double &y, double &z, int &ai, int &aj, int &ak);
    /** Creates particles within an image block by copying them
     * from the primary domain and shifting them. If the given
     * block is aligned with the primary domain in the z-direction,
     * the routine calls the simpler create_side_image routine
     * where the image block may comprise of particles from up to
     * two primary blocks. Otherwise is calls the more complex
     * create_vertical_image where the image block may comprise of
     * particles from up to four primary blocks.
     * \param[in] (di,dj,dk) the coordinates of the image block to
     *                       create. */
    inline void create_periodic_image(int di, int dj, int dk) {
      if (di < 0 || di >= nx || dj < 0 || dj >= oy || dk < 0 || dk >= oz) {
        voro_fatal_error("Constructing periodic image for nonexistent point", VOROPP_INTERNAL_ERROR);
      }
      if (dk >= ez && dk < wz) {
        if (dj < ey || dj >= wy) {
          create_side_image(di, dj, dk);
        }
      } else {
        create_vertical_image(di, dj, dk);
      }
    }
    void create_side_image(int di, int dj, int dk);
    void create_vertical_image(int di, int dj, int dk);
    void put_image(int reg, int fijk, int l, double dx, double dy, double dz);
    inline void remap(int &ai, int &aj, int &ak, int &ci, int &cj, int &ck, double &x, double &y, double &z, int &ijk);
  };

/** \brief Extension of the container_periodic_base class for computing regular
   * Voronoi tessellations.
   *
   * This class is an extension of the container_periodic_base that has routines
   * specifically for computing the regular Voronoi tessellation with no
   * dependence on particle radii. */
  class container_periodic : public container_periodic_base, public radius_mono {
  public:
    container_periodic(double bx_, double bxy_, double by_, double bxz_, double byz_, double bz_, int nx_, int ny_, int nz_, int init_mem_);
    void clear();
    void put(int n, double x, double y, double z);
    void put(int n, double x, double y, double z, int &ai, int &aj, int &ak);
    void put(particle_order &vo, int n, double x, double y, double z);
    void import(FILE *fp = stdin);
    void import(particle_order &vo, FILE *fp = stdin);
    /** Imports a list of particles from an open file stream into
     * the container. Entries of four numbers (Particle ID, x
     * position, y position, z position) are searched for. If the
     * file cannot be successfully read, then the routine causes a
     * fatal error.
     * \param[in] filename the name of the file to open and read
     *                     from. */
    inline void import(const char *filename) {
      FILE *fp = safe_fopen(filename, "r");
      import(fp);
      fclose(fp);
    }
    /** Imports a list of particles from an open file stream into
     * the container. Entries of four numbers (Particle ID, x
     * position, y position, z position) are searched for. In
     * addition, the order in which particles are read is saved
     * into an ordering class. If the file cannot be successfully
     * read, then the routine causes a fatal error.
     * \param[in,out] vo the ordering class to use.
     * \param[in] filename the name of the file to open and read
     *                     from. */
    inline void import(particle_order &vo, const char *filename) {
      FILE *fp = safe_fopen(filename, "r");
      import(vo, fp);
      fclose(fp);
    }
    void compute_all_cells();
    double sum_cell_volumes();
    /** Dumps particle IDs and positions to a file.
     * \param[in] vl the loop class to use.
     * \param[in] fp a file handle to write to. */
    template <class c_loop> void draw_particles(c_loop &vl, FILE *fp) {
      double *pp;
      if (vl.start()) {
        do {
          pp = p[vl.ijk] + 3 * vl.q;
          fprintf(fp, "%d %g %g %g\n", id[vl.ijk][vl.q], *pp, pp[1], pp[2]);
        } while (vl.inc());
      }
    }
    /** Dumps all of the particle IDs and positions to a file.
     * \param[in] fp a file handle to write to. */
    inline void draw_particles(FILE *fp = stdout) {
      c_loop_all_periodic vl(*this);
      draw_particles(vl, fp);
    }
    /** Dumps all of the particle IDs and positions to a file.
     * \param[in] filename the name of the file to write to. */
    inline void draw_particles(const char *filename) {
      FILE *fp = safe_fopen(filename, "w");
      draw_particles(fp);
      fclose(fp);
    }
    /** Dumps particle positions in POV-Ray format.
     * \param[in] vl the loop class to use.
     * \param[in] fp a file handle to write to. */
    template <class c_loop> void draw_particles_pov(c_loop &vl, FILE *fp) {
      double *pp;
      if (vl.start()) {
        do {
          pp = p[vl.ijk] + 3 * vl.q;
          fprintf(fp, "// id %d\nsphere{<%g,%g,%g>,s}\n", id[vl.ijk][vl.q], *pp, pp[1], pp[2]);
        } while (vl.inc());
      }
    }
    /** Dumps all particle positions in POV-Ray format.
     * \param[in] fp a file handle to write to. */
    inline void draw_particles_pov(FILE *fp = stdout) {
      c_loop_all_periodic vl(*this);
      draw_particles_pov(vl, fp);
    }
    /** Dumps all particle positions in POV-Ray format.
     * \param[in] filename the name of the file to write to. */
    inline void draw_particles_pov(const char *filename) {
      FILE *fp = safe_fopen(filename, "w");
      draw_particles_pov(fp);
      fclose(fp);
    }
    /** Computes Voronoi cells and saves the output in gnuplot
     * format.
     * \param[in] vl the loop class to use.
     * \param[in] fp a file handle to write to. */
    template <class c_loop> void draw_cells_gnuplot(c_loop &vl, FILE *fp) {
      voronoicell c;
      double *pp;
      if (vl.start()) {
        do {
          if (compute_cell(c, vl)) {
            pp = p[vl.ijk] + ps * vl.q;
            c.draw_gnuplot(*pp, pp[1], pp[2], fp);
          }
        } while (vl.inc());
      }
    }
    /** Computes all Voronoi cells and saves the output in gnuplot
     * format.
     * \param[in] fp a file handle to write to. */
    inline void draw_cells_gnuplot(FILE *fp = stdout) {
      c_loop_all_periodic vl(*this);
      draw_cells_gnuplot(vl, fp);
    }
    /** Compute all Voronoi cells and saves the output in gnuplot
     * format.
     * \param[in] filename the name of the file to write to. */
    inline void draw_cells_gnuplot(const char *filename) {
      FILE *fp = safe_fopen(filename, "w");
      draw_cells_gnuplot(fp);
      fclose(fp);
    }
    /** Computes Voronoi cells and saves the output in POV-Ray
     * format.
     * \param[in] vl the loop class to use.
     * \param[in] fp a file handle to write to. */
    template <class c_loop> void draw_cells_pov(c_loop &vl, FILE *fp) {
      voronoicell c;
      double *pp;
      if (vl.start()) {
        do {
          if (compute_cell(c, vl)) {
            fprintf(fp, "// cell %d\n", id[vl.ijk][vl.q]);
            pp = p[vl.ijk] + ps * vl.q;
            c.draw_pov(*pp, pp[1], pp[2], fp);
          }
        } while (vl.inc());
      }
    }
    /** Computes all Voronoi cells and saves the output in POV-Ray
     * format.
     * \param[in] fp a file handle to write to. */
    inline void draw_cells_pov(FILE *fp = stdout) {
      c_loop_all_periodic vl(*this);
      draw_cells_pov(vl, fp);
    }
    /** Computes all Voronoi cells and saves the output in POV-Ray
     * format.
     * \param[in] filename the name of the file to write to. */
    inline void draw_cells_pov(const char *filename) {
      FILE *fp = safe_fopen(filename, "w");
      draw_cells_pov(fp);
      fclose(fp);
    }
    /** Computes the Voronoi cells and saves customized information
     * about them.
     * \param[in] vl the loop class to use.
     * \param[in] format the custom output string to use.
     * \param[in] fp a file handle to write to. */
    template <class c_loop> void print_custom(c_loop &vl, const char *format, FILE *fp) {
      int ijk;
      int q;
      double *pp;
      if (contains_neighbor(format)) {
        voronoicell_neighbor c;
        if (vl.start()) {
          do {
            if (compute_cell(c, vl)) {
              ijk = vl.ijk;
              q = vl.q;
              pp = p[ijk] + ps * q;
              c.output_custom(format, id[ijk][q], *pp, pp[1], pp[2], default_radius, fp);
            }
          } while (vl.inc());
        }
      } else {
        voronoicell c;
        if (vl.start()) {
          do {
            if (compute_cell(c, vl)) {
              ijk = vl.ijk;
              q = vl.q;
              pp = p[ijk] + ps * q;
              c.output_custom(format, id[ijk][q], *pp, pp[1], pp[2], default_radius, fp);
            }
          } while (vl.inc());
        }
      }
    }
    void print_custom(const char *format, FILE *fp = stdout);
    void print_custom(const char *format, const char *filename);
    bool find_voronoi_cell(double x, double y, double z, double &rx, double &ry, double &rz, int &pid);
    /** Computes the Voronoi cell for a particle currently being
     * referenced by a loop class.
     * \param[out] c a Voronoi cell class in which to store the
     * 		 computed cell.
     * \param[in] vl the loop class to use.
     * \return True if the cell was computed. If the cell cannot be
     * computed because it was removed entirely for some reason,
     * then the routine returns false. */
    template <class v_cell, class c_loop> inline bool compute_cell(v_cell &c, c_loop &vl) { return vc.compute_cell(c, vl.ijk, vl.q, vl.i, vl.j, vl.k); }
    /** Computes the Voronoi cell for given particle.
     * \param[out] c a Voronoi cell class in which to store the
     * 		 computed cell.
     * \param[in] ijk the block that the particle is within.
     * \param[in] q the index of the particle within the block.
     * \return True if the cell was computed. If the cell cannot be
     * computed because it was removed entirely for some reason,
     * then the routine returns false. */
    template <class v_cell> inline bool compute_cell(v_cell &c, int ijk, int q) {
      int k(ijk / (nx * oy));
      int ijkt(ijk - (nx * oy) * k);
      int j(ijkt / nx);
      int i(ijkt - j * nx);
      return vc.compute_cell(c, ijk, q, i, j, k);
    }
    /** Computes the Voronoi cell for a ghost particle at a given
     * location.
     * \param[out] c a Voronoi cell class in which to store the
     * 		 computed cell.
     * \param[in] (x,y,z) the location of the ghost particle.
     * \return True if the cell was computed. If the cell cannot be
     * computed, if it is removed entirely by a wall or boundary
     * condition, then the routine returns false. */
    template <class v_cell> inline bool compute_ghost_cell(v_cell &c, double x, double y, double z) {
      int ijk;
      put_locate_block(ijk, x, y, z);
      double *pp = p[ijk] + 3 * co[ijk]++;
      *(pp++) = x;
      *(pp++) = y;
      *(pp++) = z;
      bool q = compute_cell(c, ijk, co[ijk] - 1);
      co[ijk]--;
      return q;
    }

  private:
    voro_compute<container_periodic> vc;
    friend class voro_compute<container_periodic>;
  };

/** \brief Extension of the container_periodic_base class for computing radical
   * Voronoi tessellations.
   *
   * This class is an extension of container_periodic_base that has routines
   * specifically for computing the radical Voronoi tessellation that depends
   * on the particle radii. */
  class container_periodic_poly : public container_periodic_base, public radius_poly {
  public:
    container_periodic_poly(double bx_, double bxy_, double by_, double bxz_, double byz_, double bz_, int nx_, int ny_, int nz_, int init_mem_);
    void clear();
    void put(int n, double x, double y, double z, double r);
    void put(int n, double x, double y, double z, double r, int &ai, int &aj, int &ak);
    void put(particle_order &vo, int n, double x, double y, double z, double r);
    void import(FILE *fp = stdin);
    void import(particle_order &vo, FILE *fp = stdin);
    /** Imports a list of particles from an open file stream into
     * the container_poly class. Entries of five numbers (Particle
     * ID, x position, y position, z position, radius) are searched
     * for. If the file cannot be successfully read, then the
     * routine causes a fatal error.
     * \param[in] filename the name of the file to open and read
     *                     from. */
    inline void import(const char *filename) {
      FILE *fp = safe_fopen(filename, "r");
      import(fp);
      fclose(fp);
    }
    /** Imports a list of particles from an open file stream into
     * the container_poly class. Entries of five numbers (Particle
     * ID, x position, y position, z position, radius) are searched
     * for. In addition, the order in which particles are read is
     * saved into an ordering class. If the file cannot be
     * successfully read, then the routine causes a fatal error.
     * \param[in,out] vo the ordering class to use.
     * \param[in] filename the name of the file to open and read
     *                     from. */
    inline void import(particle_order &vo, const char *filename) {
      FILE *fp = safe_fopen(filename, "r");
      import(vo, fp);
      fclose(fp);
    }
    void compute_all_cells();
    double sum_cell_volumes();
    /** Dumps particle IDs, positions and radii to a file.
     * \param[in] vl the loop class to use.
     * \param[in] fp a file handle to write to. */
    template <class c_loop> void draw_particles(c_loop &vl, FILE *fp) {
      double *pp;
      if (vl.start()) {
        do {
          pp = p[vl.ijk] + 4 * vl.q;
          fprintf(fp, "%d %g %g %g %g\n", id[vl.ijk][vl.q], *pp, pp[1], pp[2], pp[3]);
        } while (vl.inc());
      }
    }
    /** Dumps all of the particle IDs, positions and radii to a
     * file.
     * \param[in] fp a file handle to write to. */
    inline void draw_particles(FILE *fp = stdout) {
      c_loop_all_periodic vl(*this);
      draw_particles(vl, fp);
    }
    /** Dumps all of the particle IDs, positions and radii to a
     * file.
     * \param[in] filename the name of the file to write to. */
    inline void draw_particles(const char *filename) {
      FILE *fp = safe_fopen(filename, "w");
      draw_particles(fp);
      fclose(fp);
    }
    /** Dumps particle positions in POV-Ray format.
     * \param[in] vl the loop class to use.
     * \param[in] fp a file handle to write to. */
    template <class c_loop> void draw_particles_pov(c_loop &vl, FILE *fp) {
      double *pp;
      if (vl.start()) {
        do {
          pp = p[vl.ijk] + 4 * vl.q;
          fprintf(fp, "// id %d\nsphere{<%g,%g,%g>,%g}\n", id[vl.ijk][vl.q], *pp, pp[1], pp[2], pp[3]);
        } while (vl.inc());
      }
    }
    /** Dumps all the particle positions in POV-Ray format.
     * \param[in] fp a file handle to write to. */
    inline void draw_particles_pov(FILE *fp = stdout) {
      c_loop_all_periodic vl(*this);
      draw_particles_pov(vl, fp);
    }
    /** Dumps all the particle positions in POV-Ray format.
     * \param[in] filename the name of the file to write to. */
    inline void draw_particles_pov(const char *filename) {
      FILE *fp(safe_fopen(filename, "w"));
      draw_particles_pov(fp);
      fclose(fp);
    }
    /** Computes Voronoi cells and saves the output in gnuplot
     * format.
     * \param[in] vl the loop class to use.
     * \param[in] fp a file handle to write to. */
    template <class c_loop> void draw_cells_gnuplot(c_loop &vl, FILE *fp) {
      voronoicell c;
      double *pp;
      if (vl.start()) {
        do {
          if (compute_cell(c, vl)) {
            pp = p[vl.ijk] + ps * vl.q;
            c.draw_gnuplot(*pp, pp[1], pp[2], fp);
          }
        } while (vl.inc());
      }
    }
    /** Compute all Voronoi cells and saves the output in gnuplot
     * format.
     * \param[in] fp a file handle to write to. */
    inline void draw_cells_gnuplot(FILE *fp = stdout) {
      c_loop_all_periodic vl(*this);
      draw_cells_gnuplot(vl, fp);
    }
    /** Compute all Voronoi cells and saves the output in gnuplot
     * format.
     * \param[in] filename the name of the file to write to. */
    inline void draw_cells_gnuplot(const char *filename) {
      FILE *fp(safe_fopen(filename, "w"));
      draw_cells_gnuplot(fp);
      fclose(fp);
    }
    /** Computes Voronoi cells and saves the output in POV-Ray
     * format.
     * \param[in] vl the loop class to use.
     * \param[in] fp a file handle to write to. */
    template <class c_loop> void draw_cells_pov(c_loop &vl, FILE *fp) {
      voronoicell c;
      double *pp;
      if (vl.start()) {
        do {
          if (compute_cell(c, vl)) {
            fprintf(fp, "// cell %d\n", id[vl.ijk][vl.q]);
            pp = p[vl.ijk] + ps * vl.q;
            c.draw_pov(*pp, pp[1], pp[2], fp);
          }
        } while (vl.inc());
      }
    }
    /** Computes all Voronoi cells and saves the output in POV-Ray
     * format.
     * \param[in] fp a file handle to write to. */
    inline void draw_cells_pov(FILE *fp = stdout) {
      c_loop_all_periodic vl(*this);
      draw_cells_pov(vl, fp);
    }
    /** Computes all Voronoi cells and saves the output in POV-Ray
     * format.
     * \param[in] filename the name of the file to write to. */
    inline void draw_cells_pov(const char *filename) {
      FILE *fp(safe_fopen(filename, "w"));
      draw_cells_pov(fp);
      fclose(fp);
    }
    /** Computes the Voronoi cells and saves customized information
     * about them.
     * \param[in] vl the loop class to use.
     * \param[in] format the custom output string to use.
     * \param[in] fp a file handle to write to. */
    template <class c_loop> void print_custom(c_loop &vl, const char *format, FILE *fp) {
      int ijk;
      int q;
      double *pp;
      if (contains_neighbor(format)) {
        voronoicell_neighbor c;
        if (vl.start()) {
          do {
            if (compute_cell(c, vl)) {
              ijk = vl.ijk;
              q = vl.q;
              pp = p[ijk] + ps * q;
              c.output_custom(format, id[ijk][q], *pp, pp[1], pp[2], pp[3], fp);
            }
          } while (vl.inc());
        }
      } else {
        voronoicell c;
        if (vl.start()) {
          do {
            if (compute_cell(c, vl)) {
              ijk = vl.ijk;
              q = vl.q;
              pp = p[ijk] + ps * q;
              c.output_custom(format, id[ijk][q], *pp, pp[1], pp[2], pp[3], fp);
            }
          } while (vl.inc());
        }
      }
    }
    /** Computes the Voronoi cell for a particle currently being
     * referenced by a loop class.
     * \param[out] c a Voronoi cell class in which to store the
     * 		 computed cell.
     * \param[in] vl the loop class to use.
     * \return True if the cell was computed. If the cell cannot be
     * computed because it was removed entirely for some reason,
     * then the routine returns false. */
    template <class v_cell, class c_loop> inline bool compute_cell(v_cell &c, c_loop &vl) { return vc.compute_cell(c, vl.ijk, vl.q, vl.i, vl.j, vl.k); }
    /** Computes the Voronoi cell for given particle.
     * \param[out] c a Voronoi cell class in which to store the
     * 		 computed cell.
     * \param[in] ijk the block that the particle is within.
     * \param[in] q the index of the particle within the block.
     * \return True if the cell was computed. If the cell cannot be
     * computed because it was removed entirely for some reason,
     * then the routine returns false. */
    template <class v_cell> inline bool compute_cell(v_cell &c, int ijk, int q) {
      int k(ijk / (nx * oy));
      int ijkt(ijk - (nx * oy) * k);
      int j(ijkt / nx);
      int i(ijkt - j * nx);
      return vc.compute_cell(c, ijk, q, i, j, k);
    }
    /** Computes the Voronoi cell for a ghost particle at a given
     * location.
     * \param[out] c a Voronoi cell class in which to store the
     * 		 computed cell.
     * \param[in] (x,y,z) the location of the ghost particle.
     * \param[in] r the radius of the ghost particle.
     * \return True if the cell was computed. If the cell cannot be
     * computed, if it is removed entirely by a wall or boundary
     * condition, then the routine returns false. */
    template <class v_cell> inline bool compute_ghost_cell(v_cell &c, double x, double y, double z, double r) {
      int ijk;
      put_locate_block(ijk, x, y, z);
      double *pp = p[ijk] + 4 * co[ijk]++;
      double tm = max_radius;
      *(pp++) = x;
      *(pp++) = y;
      *(pp++) = z;
      *pp = r;
      if (r > max_radius) {
        max_radius = r;
      }
      bool q = compute_cell(c, ijk, co[ijk] - 1);
      co[ijk]--;
      max_radius = tm;
      return q;
    }
    void print_custom(const char *format, FILE *fp = stdout);
    void print_custom(const char *format, const char *filename);
    bool find_voronoi_cell(double x, double y, double z, double &rx, double &ry, double &rz, int &pid);

  private:
    voro_compute<container_periodic_poly> vc;
    friend class voro_compute<container_periodic_poly>;
  };

} // namespace voro

#endif

// Voro++, a 3D cell-based Voronoi library
//
// Author   : Chris H. Rycroft (LBL / UC Berkeley)
// Email    : chr@alum.mit.edu
// Date     : August 30th 2011

/** \file pre_container.hh
 * \brief Header file for the pre_container and related classes. */

#ifndef VOROPP_PRE_CONTAINER_HH
#define VOROPP_PRE_CONTAINER_HH

#include <cstdio>

namespace voro {

/** \brief A class for storing an arbitrary number of particles, prior to setting
   * up a container geometry.
   *
   * The pre_container_base class can dynamically import and store an arbitrary
   * number of particles. Once the particles have been read in, an appropriate
   * container class can be set up with the optimal grid size, and the particles
   * can be transferred.
   *
   * The pre_container_base class is not intended for direct use, but forms the
   * base of the pre_container and pre_container_poly classes, that add routines
   * depending on whether particle radii need to be tracked or not. */
  class pre_container_base {
  public:
    /** The minimum x coordinate of the container. */
    const double ax;
    /** The maximum x coordinate of the container. */
    const double bx;
    /** The minimum y coordinate of the container. */
    const double ay;
    /** The maximum y coordinate of the container. */
    const double by;
    /** The minimum z coordinate of the container. */
    const double az;
    /** The maximum z coordinate of the container. */
    const double bz;
    /** A boolean value that determines if the x coordinate in
     * periodic or not. */
    const bool xperiodic;
    /** A boolean value that determines if the y coordinate in
     * periodic or not. */
    const bool yperiodic;
    /** A boolean value that determines if the z coordinate in
     * periodic or not. */
    const bool zperiodic;
    void guess_optimal(int &nx, int &ny, int &nz);
    pre_container_base(double ax_, double bx_, double ay_, double by_, double az_, double bz_, bool xperiodic_, bool yperiodic_, bool zperiodic_, int ps_);
    ~pre_container_base();
    /** Calculates and returns the total number of particles stored
     * within the class.
     * \return The number of particles. */
    inline int total_particles() { return (end_id - pre_id) * pre_container_chunk_size + (ch_id - *end_id); }

  protected:
    /** The number of doubles associated with a single particle
     * (three for the standard container, four when radius
     * information is stored). */
    const int ps;
    void new_chunk();
    void extend_chunk_index();
    /** The size of the chunk index. */
    int index_sz;
    /** A pointer to the chunk index to store the integer particle
     * IDs. */
    int **pre_id;
    /** A pointer to the last allocated integer ID chunk. */
    int **end_id;
    /** A pointer to the end of the integer ID chunk index, used to
     * determine when the chunk index is full. */
    int **l_id;
    /** A pointer to the next available slot on the current
     * particle ID chunk. */
    int *ch_id;
    /** A pointer to the end of the current integer chunk. */
    int *e_id;
    /** A pointer to the chunk index to store the floating point
     * information associated with particles. */
    double **pre_p;
    /** A pointer to the last allocated chunk of floating point
     * information. */
    double **end_p;
    /** A pointer to the next available slot on the current
     * floating point chunk. */
    double *ch_p;
  };

/** \brief A class for storing an arbitrary number of particles without radius
   * information, prior to setting up a container geometry.
   *
   * The pre_container class is an extension of the pre_container_base class for
   * cases when no particle radius information is available. */
  class pre_container : public pre_container_base {
  public:
    /** The class constructor sets up the geometry of container,
     * initializing the minimum and maximum coordinates in each
     * direction.
     * \param[in] (ax_,bx_) the minimum and maximum x coordinates.
     * \param[in] (ay_,by_) the minimum and maximum y coordinates.
     * \param[in] (az_,bz_) the minimum and maximum z coordinates.
     * \param[in] (xperiodic_,yperiodic_,zperiodic_ ) flags setting whether the
     *                                                container is periodic in
     *                                                each coordinate direction. */
    pre_container(double ax_, double bx_, double ay_, double by_, double az_, double bz_, bool xperiodic_, bool yperiodic_, bool zperiodic_) :
        pre_container_base(ax_, bx_, ay_, by_, az_, bz_, xperiodic_, yperiodic_, zperiodic_, 3) {};
    void put(int n, double x, double y, double z);
    void import(FILE *fp = stdin);
    /** Imports particles from a file.
     * \param[in] filename the name of the file to read from. */
    inline void import(const char *filename) {
      FILE *fp = safe_fopen(filename, "r");
      import(fp);
      fclose(fp);
    }
    void setup(container &con);
    void setup(particle_order &vo, container &con);
  };

/** \brief A class for storing an arbitrary number of particles with radius
   * information, prior to setting up a container geometry.
   *
   * The pre_container_poly class is an extension of the pre_container_base class
   * for cases when particle radius information is available. */
  class pre_container_poly : public pre_container_base {
  public:
    /** The class constructor sets up the geometry of container,
     * initializing the minimum and maximum coordinates in each
     * direction.
     * \param[in] (ax_,bx_) the minimum and maximum x coordinates.
     * \param[in] (ay_,by_) the minimum and maximum y coordinates.
     * \param[in] (az_,bz_) the minimum and maximum z coordinates.
     * \param[in] (xperiodic_,yperiodic_,zperiodic_ ) flags setting whether the
     *                                                container is periodic in
     *                                                each coordinate direction. */
    pre_container_poly(double ax_, double bx_, double ay_, double by_, double az_, double bz_, bool xperiodic_, bool yperiodic_, bool zperiodic_) :
        pre_container_base(ax_, bx_, ay_, by_, az_, bz_, xperiodic_, yperiodic_, zperiodic_, 4) {};
    void put(int n, double x, double y, double z, double r);
    void import(FILE *fp = stdin);
    /** Imports particles from a file.
     * \param[in] filename the name of the file to read from. */
    inline void import(const char *filename) {
      FILE *fp = safe_fopen(filename, "r");
      import(fp);
      fclose(fp);
    }
    void setup(container_poly &con);
    void setup(particle_order &vo, container_poly &con);
  };

} // namespace voro

#endif

// Voro++, a 3D cell-based Voronoi library
//
// Author   : Chris H. Rycroft (LBL / UC Berkeley)
// Email    : chr@alum.mit.edu
// Date     : August 30th 2011

/** \file wall.hh
 * \brief Header file for the derived wall classes. */

#ifndef VOROPP_WALL_HH
#define VOROPP_WALL_HH

namespace voro {

/** \brief A class representing a spherical wall object.
   *
   * This class represents a spherical wall object. */
  struct wall_sphere : public wall {
  public:
    /** Constructs a spherical wall object.
     * \param[in] w_id_ an ID number to associate with the wall for
     *		    neighbor tracking.
     * \param[in] (xc_,yc_,zc_) a position vector for the sphere's
     * 			    center.
     * \param[in] rc_ the radius of the sphere. */
    wall_sphere(double xc_, double yc_, double zc_, double rc_, int w_id_ = -99) : w_id(w_id_), xc(xc_), yc(yc_), zc(zc_), rc(rc_) {}
    bool point_inside(double x, double y, double z);
    template <class v_cell> bool cut_cell_base(v_cell &c, double x, double y, double z);
    bool cut_cell(voronoicell &c, double x, double y, double z) { return cut_cell_base(c, x, y, z); }
    bool cut_cell(voronoicell_neighbor &c, double x, double y, double z) { return cut_cell_base(c, x, y, z); }

  private:
    const int w_id;
    const double xc, yc, zc, rc;
  };

/** \brief A class representing a plane wall object.
   *
   * This class represents a single plane wall object. */
  struct wall_plane : public wall {
  public:
    /** Constructs a plane wall object.
     * \param[in] (xc_,yc_,zc_) a normal vector to the plane.
     * \param[in] ac_ a displacement along the normal vector.
     * \param[in] w_id_ an ID number to associate with the wall for
     *		    neighbor tracking. */
    wall_plane(double xc_, double yc_, double zc_, double ac_, int w_id_ = -99) : w_id(w_id_), xc(xc_), yc(yc_), zc(zc_), ac(ac_) {}
    bool point_inside(double x, double y, double z);
    template <class v_cell> bool cut_cell_base(v_cell &c, double x, double y, double z);
    bool cut_cell(voronoicell &c, double x, double y, double z) { return cut_cell_base(c, x, y, z); }
    bool cut_cell(voronoicell_neighbor &c, double x, double y, double z) { return cut_cell_base(c, x, y, z); }

  private:
    const int w_id;
    const double xc, yc, zc, ac;
  };

/** \brief A class representing a cylindrical wall object.
   *
   * This class represents a open cylinder wall object. */
  struct wall_cylinder : public wall {
  public:
    /** Constructs a cylinder wall object.
     * \param[in] (xc_,yc_,zc_) a point on the axis of the
     *			    cylinder.
     * \param[in] (xa_,ya_,za_) a vector pointing along the
     *			    direction of the cylinder.
     * \param[in] rc_ the radius of the cylinder
     * \param[in] w_id_ an ID number to associate with the wall for
     *		    neighbor tracking. */
    wall_cylinder(double xc_, double yc_, double zc_, double xa_, double ya_, double za_, double rc_, int w_id_ = -99) :
        w_id(w_id_), xc(xc_), yc(yc_), zc(zc_), xa(xa_), ya(ya_), za(za_), asi(1 / (xa_ * xa_ + ya_ * ya_ + za_ * za_)), rc(rc_) {}
    bool point_inside(double x, double y, double z);
    template <class v_cell> bool cut_cell_base(v_cell &c, double x, double y, double z);
    bool cut_cell(voronoicell &c, double x, double y, double z) { return cut_cell_base(c, x, y, z); }
    bool cut_cell(voronoicell_neighbor &c, double x, double y, double z) { return cut_cell_base(c, x, y, z); }

  private:
    const int w_id;
    const double xc, yc, zc, xa, ya, za, asi, rc;
  };

/** \brief A class representing a conical wall object.
   *
   * This class represents a cone wall object. */
  struct wall_cone : public wall {
  public:
    /** Constructs a cone wall object.
     * \param[in] (xc_,yc_,zc_) the apex of the cone.
     * \param[in] (xa_,ya_,za_) a vector pointing along the axis of
     *			    the cone.
     * \param[in] ang the angle (in radians) of the cone, measured
     *		  from the axis.
     * \param[in] w_id_ an ID number to associate with the wall for
     *		    neighbor tracking. */
    wall_cone(double xc_, double yc_, double zc_, double xa_, double ya_, double za_, double ang, int w_id_ = -99) :
        w_id(w_id_), xc(xc_), yc(yc_), zc(zc_), xa(xa_), ya(ya_), za(za_), asi(1 / (xa_ * xa_ + ya_ * ya_ + za_ * za_)), gra(tan(ang)), sang(sin(ang)), cang(cos(ang)) {}
    bool point_inside(double x, double y, double z);
    template <class v_cell> bool cut_cell_base(v_cell &c, double x, double y, double z);
    bool cut_cell(voronoicell &c, double x, double y, double z) { return cut_cell_base(c, x, y, z); }
    bool cut_cell(voronoicell_neighbor &c, double x, double y, double z) { return cut_cell_base(c, x, y, z); }

  private:
    const int w_id;
    const double xc, yc, zc, xa, ya, za, asi, gra, sang, cang;
  };

} // namespace voro

#endif

#endif
