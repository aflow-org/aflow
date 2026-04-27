// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *                                                                         *
// ***************************************************************************
// This xplotter class is the evolution of different prior plotting functions in the AFLOW source code.
// hagen.eckert@duke.edu

#ifndef _AUROSTD_xplotter_H_
#define _AUROSTD_xplotter_H_

#include <array>
#include <cstddef>
#include <map>
#include <set>
#include <string>
#include <vector>

#include "AUROSTD/aurostd_xmatrix.h"
#include "aurostd.h"

namespace aurostd {
  class color {
  public:
    color();

    ~color();

    static std::map<std::string, std::array<std::array<float, 3>, 256>> fillColormaps();

    static std::map<std::string, std::array<std::array<float, 3>, 256>> maps;

    static std::map<std::string, std::vector<std::array<float, 3>>> fillCyclics();

    static std::map<std::string, std::vector<std::array<float, 3>>> cyclics;

    static std::string getCyclicTikz(const std::string name, const size_t index);
  };

  // HE20220314 BEGIN X3DWriter
  class x3DWriter {
  private:
    enum class object_types { MATERIAL, SPHERE, FACET, OPEN_CYLINDER };
    void copy(const x3DWriter& x3w);

    ///@brief Helper to store different scene objects (utalizing shared pointers)
    struct storage_object {
      object_types type;
      std::shared_ptr<void> obj;
    };
    std::vector<storage_object> objects; ///< storage for scene objects

    /// @brief Scene object settings for spheres
    struct Sphere {
      xvector<double> center;
      double radius;
      std::string material;
    };

    /// @brief Scene object settings for facets
    struct ConvexFacets {
      std::vector<xvector<double>> vertexes;
      std::vector<std::vector<uint>> facets;
      std::string material;
    };

    /// @brief Scene object settings for cylinders
    struct OpenCylinder {
      xvector<double> base;
      xvector<double> apex;
      xvector<double> center;
      double radius;
      std::string material;
      double height;
      double phi; // rotation around 1 0 0
      double theta; // rotation around 0 0 1
    };

    ///@brief Storage object for material settings
    struct Material {
      std::string name;
      xvector<double> color;
      double ambient;
      double opacity;
      double specular;
      double diffuse;
      /// @brief create a default Material
      Material() : ambient(0.6), opacity(1.0), specular(0.05), diffuse(0.4) {};
    };

    double max_distance = 0.0; ///< maximum distance of scene object from the scene center (for camera position calculation)
    double tachyon_zoom = 2.0; ///< zoom value for rendering the camera rendering with tachyon
    const double tachyon_camera_angle = 45.0; ///< lense angle for the camera rendering with tachyon
    xvector<double> tachyon_camera_position; ///< camera position for tachyon
    std::vector<std::pair<xvector<double>, xvector<double>>> tachyon_lattice_views; ///< collection of camera views based on prepareSceneLattice
    void tachyon_calculate_camera();
    static std::string x3d_material(const std::shared_ptr<x3DWriter::Material>& material);

  public:
    xvector<double> scene_center; ///< center of the scene (should be set before adding objects to the scene)
    double join_threshold = 1E-5; ///< threshould to join vertexes when using joinFacets
    double tachyon_camera_theta = 0.0; ///< tachyon camera angle theta | rotation around 1 0 0
    double tachyon_camera_phi = 0.0; ///<  tachyon camera angle phi | rotation around 0 0 1
    int tachyon_lattice_views_idx = -1; ///< active lattice view (-1 default tachyon viewpoint is used)
    bool tachyon_camera_orthographic = false;
    //@brief available animation formats
    enum class animation_format { GIF, WEBM, MP4 };
    animation_format ani_type = animation_format::GIF; ///< animation type (default GIF)
    std::vector<std::pair<std::string, std::string>> meta; ///<  meta data like title, description, reference (link back to auid)
    x3DWriter();
    x3DWriter(const x3DWriter& x3w);
    ~x3DWriter();
    x3DWriter& operator=(const x3DWriter& x3w);
    void clear();

    void prepareSceneLattice(const xmatrix<double>& lattice);
    void addMaterial(const std::string& name, const xvector<double>& color);
    void addMaterial(const Material& newMaterial);
    std::vector<std::string> addColorSpreadMaterial(uint count, const std::string& cmap = "turbo");
    template <class Container> std::map<std::string, std::string> addNamedColorSpreadMaterial(const Container& names);

    void addSphere(const xvector<double>& center, double radius, const std::string& material = "bm_black");
    void addSpheres(const std::vector<xvector<double>>& center, double radius, const std::string& material = "bm_black");
    void addOpenCylinder(const xvector<double>& base, const xvector<double>& apex, double radius, const std::string& material = "bm_grey");
    void addConvexFacets(const std::vector<xvector<double>>& vertexes, const std::vector<std::vector<uint>>& facets, const std::string& material = "bm_grey_glass", const xvector<double>& shift = {0, 0, 0});
    void addConvexFacets(const std::vector<xvector<double>>& vertexes, const std::vector<std::vector<uint>>& facets, const xvector<double>& shift);
    void addLatticeBox(const xmatrix<double>& lattice, double radius, const std::string& material = "bm_black", bool axis = false);

    void joinFacets(std::vector<xvector<double>>& vertexes, std::vector<std::vector<uint>>& facets, const std::vector<xvector<double>>& new_vertexes, const std::vector<std::vector<uint>>& new_facets) const;

    void animate(float duration, const std::filesystem::path& out_folder, uint fps = 30, bool lr = true);

    std::string toHTML();
    std::string toX3D(bool include_xml = true, bool replace_material = false);
    std::string toTachyon();
  };
// HE20220314 END X3DWriter

  //SD+HE20230217
  // xplotter
  // collection of basic plotting routines
  class xplotter {
  private:
    /// @brief axis type of the figure (set by the first added plot)
    enum class axis_types { TWODIM, BAR, BARSYMBOLIC, THREEDIM, TERNARY, NONE };

    /// @brief available plot types
    enum class object_types { SMOOTH, SCATTER, BAR, LINE, VLINE, HLINE, TERNARY_SCATTER, TERNARY_SURF, TERNARY_CONTOUR, NONE };

    axis_types axisType; ///< current axis type

    /// @brief plot data structure
    struct storage_object {
      object_types type;
      std::shared_ptr<void> obj;
    };

    std::vector<storage_object> objects; ///< save plot data that makes up the figure

    /// @brief data storage for x and y
    struct xy_data {
      std::vector<double> x;
      std::vector<double> y;
      std::map<std::string, std::string> properties;
    };

    /// @brief data storage for x, y, and label list
    struct xyl_data {
      std::vector<double> x;
      std::vector<double> y;
      std::vector<std::string> l;
      std::map<std::string, std::string> properties;
    };

    /// @brief data storage for line position data (two entries for h/v lines and four for custom lines)
    struct line_data {
      std::vector<double> pos;
      std::map<std::string, std::string> properties;
    };

    /// @brief data storage for x, y, z and m
    struct xyzm_data {
      std::vector<double> x;
      std::vector<double> y;
      std::vector<double> z;
      std::vector<double> m;
      std::map<std::string, std::string> properties;
    };

    // property helpers
    static std::set<std::string> excluded_line_properties;
    static std::set<std::string> line_bool_properties;
    static std::set<std::string> excluded_axis_properties;
    static std::set<std::string> string_axis_properties;
    std::vector<std::string> symbolic_coordinates;
    ///< collect symbolic coordinates for all plots - entries should be unique

    std::string prepare_plot_options(std::map<std::string, std::string>& properties, uint& auto_color_count, uint& cyclic_color_count, const object_types plot_types);

    std::string prepare_axis_options(const std::map<std::string, std::string>& properties);

    std::string prepare_colorbar_options(const std::map<std::string, std::string>& properties);

    std::string prepare_meta_csv(const std::set<std::string>& exclude = {});

    // templates and defaults
    static std::string preamble; ///< latex template - preamble
    static std::map<object_types, std::map<std::string, std::string>> object_defaults; ///< defaults for new plots
    static std::map<axis_types, std::map<std::string, std::string>> axis_defaults;
    ///< default for the overall figure
    static std::map<std::string, std::vector<std::pair<std::string, std::string>>> convert_commands;

    // plot helper functions
    void ternary_plot(const std::vector<double>& x, const std::vector<double>& y, const std::vector<double>& z, const std::vector<double>& m, std::map<std::string, std::string> properties, const object_types& ternary_type);

    // basic private class functions
    void copy(const xplotter& xplt);

  public:
    // basic public class functions
    xplotter();

    xplotter(const xplotter& xplt);

    ~xplotter();

    void clear();

    // settings and properties
    std::vector<std::pair<std::string, std::string>> meta; ///< meta information (stored as keywords in .pdf)
    std::map<std::string, std::string> axis_properties; ///< overall figure properties
    std::map<std::string, std::string> colorbar_properties;
    std::string auto_colormap; ///< when {"color", "auto"} - a color is automatically chosen from this colormap
    std::string cyclic_name;
    ///< when {"color", "cyclic"} or {"color", "cyclic5"} - the next or the 5th color in the named cyclic list is used
    std::vector<std::string> nodes; /// compound labels for hull plots

    // plotting functions
    void plot(const std::vector<double>& y, std::map<std::string, std::string> properties = {});

    void plot(const std::vector<double>& x, const std::vector<double>& y, std::map<std::string, std::string> properties = {});

    void scatter(const std::vector<double>& x, const std::vector<double>& y, std::map<std::string, std::string> properties = {});

    void scatter(const std::vector<double>& x, const std::vector<double>& y, const std::vector<double>& z, std::map<std::string, std::string> properties = {});

    void scatter(const std::vector<double>& x, const std::vector<double>& y, const std::vector<double>& z, const std::vector<double>& m, std::map<std::string, std::string> properties = {});

    void bar(const std::vector<double>& values, std::map<std::string, std::string> properties = {});

    void bar(const std::vector<double>& values, const std::vector<double>& edges, std::map<std::string, std::string> properties = {});

    void bar(const std::vector<double>& values, const std::vector<std::string>& labels, std::map<std::string, std::string> properties = {});

    void vline(double x, std::map<std::string, std::string> properties = {});

    void vline(double x, double s, double e, std::map<std::string, std::string> properties = {});

    void hline(double y, std::map<std::string, std::string> properties = {});

    void hline(double y, double s, double e, std::map<std::string, std::string> properties = {});

    void line(const std::vector<double>& pos, std::map<std::string, std::string> properties = {});

    void line(const std::pair<double, double>& start, const std::pair<double, double>& end, std::map<std::string, std::string> properties = {});

    void ternary_scatter(const std::vector<double>& x, const std::vector<double>& y, const std::vector<double>& z, const std::vector<double>& m, std::map<std::string, std::string> properties = {});

    void ternary_contour(const std::vector<double>& x, const std::vector<double>& y, const std::vector<double>& z, const std::vector<double>& m, std::map<std::string, std::string> properties = {});

    void ternary_surf(const std::vector<double>& x, const std::vector<double>& y, const std::vector<double>& z, const std::vector<double>& m, std::map<std::string, std::string> properties = {});

    // export functions
    std::string toString();

    void save(const std::string& file_path);
  };

  // NHA 20250729
  // xtable
  // table creation routines
  class xtable {
  public:
    //member variables
    std::string title;
    std::vector<std::vector<std::string>> data; //simple vector data
    std::map<std::string, std::vector<std::map<std::string, std::vector<std::string>>>> compound_data; //compound data type: each vector has a subtitle in table
    std::vector<std::string> column_titles;
    std::vector<std::string> column_units;
    std::string font;
    int font_size;

    //constructors
    xtable() = default;
    xtable(std::vector<std::vector<std::string>> data, std::vector<std::string> column_titles, std::vector<std::string> column_units, int font_size, std::string title);
    xtable(std::map<std::string, std::vector<std::map<std::string, std::vector<std::string>>>> compound_data, std::vector<std::string> column_titles, std::vector<std::string> column_units, int font_size, std::string title);

    //member functions
    std::string createCompoundTable();
    std::string createTable();
    std::string toString();

    void save(const std::string& file_path); //saves just the table to pdf
    void saveWithPlot(xtable table, xplotter xplt, const std::string& file_path); //appends table to a given plot and saves to pdf
  };
} // namespace aurostd

#endif //_AUROSTD_xplotter_H_
