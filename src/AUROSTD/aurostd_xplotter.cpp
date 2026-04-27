// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2025           *
// *                                                                         *
// ***************************************************************************
// This xplotter class is the evolution of different prior plotting functions in the AFLOW source code.
// hagen.eckert@duke.edu

#include "aurostd_xplotter.h"

#include "config.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <set>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "AUROSTD/aurostd.h"
#include "AUROSTD/aurostd_defs.h"
#include "AUROSTD/aurostd_xerror.h"
#include "AUROSTD/aurostd_xvector.h"
#include "aurostd_time.h"
#include "aurostd_xfile.h"
#include "aurostd_xmatrix.h"

using std::cerr;
using std::endl;
using std::ifstream;
using std::iostream;
using std::istringstream;
using std::map;
using std::ofstream;
using std::ostringstream;
using std::pair;
using std::set;
using std::string;
using std::stringstream;
using std::vector;

//HE20220922  BEGIN x3D
namespace aurostd {
  //***************************************************************************

  /// @class x3DWriter
  /// @brief tools to create 3D scene and export them in different formats
  ///
  /// @authors
  /// @mod{HE,20220922,created x3Dwriter}
  /// @mod{HE,20240530,implemented in AFLOW4}
  /// @mod{HE,20251216,moved to xplotter}
  ///
  /// Basic usage
  /// @code
  /// // load a structure to render
  /// const std::string structure_file = "example.vasp";
  /// xstructure work_structure(structure_file, IOAFLOW_AUTO);
  /// work_structure.ReScale(1.0);
  /// // prepare x3DWriter
  /// aurostd::x3DWriter w;
  /// w.scene_center = (work_structure.lattice.getcol(1) + work_structure.lattice.getcol(2) + work_structure.lattice.getcol(3)) / 3.0;
  /// w.addLatticeBox(work_structure.lattice, 0.1);
  /// for (auto &a: work_structure.atoms)
  ///   w.addSphere({a.cpos[1], a.cpos[2], a.cpos[3]}, 0.3, "bm_blue");
  /// aurostd::string2file(w.toHTML(), "/Users/nathan/Projects/AFLOW4/test_struc/AB_cP2_221_b_a.html");
  /// @endcode
  ///
  /// @see `xstructure::render()`

  /// @brief reconstruct / clear
  void x3DWriter::clear() {
    *this = {};
  } // calling the constructor

  /// @brief create a copy of x3DWriter
  void x3DWriter::copy(const x3DWriter& x3w) {
    objects = x3w.objects;
    max_distance = x3w.max_distance;
    tachyon_zoom = x3w.tachyon_zoom;
    tachyon_camera_position = x3w.tachyon_camera_position;
    tachyon_camera_theta = x3w.tachyon_camera_theta;
    tachyon_camera_phi = x3w.tachyon_camera_phi;
    ani_type = x3w.ani_type;
    meta = x3w.meta;
    scene_center = x3w.scene_center;
    join_threshold = x3w.join_threshold;
  }

  /// @brief class constractor
  x3DWriter::x3DWriter() {
    meta.emplace_back("reference", "https://aflow.org");
    meta.emplace_back("generator", "AFLOW " + static_cast<std::string>(AFLOW_VERSION));
    meta.emplace_back("created", get_datetime_formatted("-", true, " ", ":"));

    // set up basic materials
    addMaterial("bm_black", {0.0, 0.0, 0.0});
    addMaterial("bm_white", {1.0, 1.0, 1.0});
    addMaterial("bm_red", {0.58, 0.07, 0.0});
    addMaterial("bm_blue", {0.01, 0.10, 0.58});
    addMaterial("bm_green", {0.0, 0.56, 0.0});
    addMaterial("bm_grey", {0.3, 0.3, 0.3});
    // set up grey glass
    Material matGreyGlass;
    matGreyGlass.name = "bm_grey_glass";
    matGreyGlass.color = {0.3, 0.3, 0.3};
    matGreyGlass.opacity = 0.3;
    matGreyGlass.specular = 0.0;
    matGreyGlass.ambient = 0.1;
    matGreyGlass.diffuse = 0.7;
    addMaterial(matGreyGlass);

    scene_center = {0.0, 0.0, 0.0};
  }

  /// @brief class copy constructor
  x3DWriter::x3DWriter(const x3DWriter& x3w) {
    copy(x3w);
  }

  /// @brief default class de-constructor
  x3DWriter::~x3DWriter() = default;

  /// @brief assignment operator
  x3DWriter& x3DWriter::operator=(const x3DWriter& x3w) {
    if (this == &x3w) {
      return *this;
    }
    copy(x3w);
    return *this;
  }

  /// @brief prepare a scene to show a lattice by setting the center and a view collection
  void x3DWriter::prepareSceneLattice(const xmatrix<double>& lattice) {
    const xvector<double> endpoint = lattice(1) + lattice(2) + lattice(3);
    scene_center = endpoint / 2.0;
    constexpr std::array<std::pair<int, int>, 3> view_set({
        {{1, 2}, {2, 3}, {3, 1}}
    });
    tachyon_lattice_views.clear();
    for (const auto& [main, secondary] : view_set) {
      const xvector<double> plane_normal = aurostd::vector_product(lattice(main), lattice(secondary));
      const xvector<double> up_dir = std::cos(pi * 0.5) * lattice(main) + std::sin(pi * 0.5) * plane_normal;
      tachyon_lattice_views.emplace_back(lattice(main) * 2, aurostd::normalizeSumToOne(up_dir));
    }
  }

  /// @brief add a Sphere to the scene
  void x3DWriter::addSphere(const xvector<double>& center, double radius, const std::string& material) {
    const std::shared_ptr<x3DWriter::Sphere> newSphere = std::make_shared<x3DWriter::Sphere>();
    newSphere->center = center - scene_center;
    newSphere->radius = radius;
    newSphere->material = material;
    max_distance = std::max(max_distance, aurostd::modulus(center));
    x3DWriter::storage_object so;
    so.obj = newSphere; // cast to a typeless void pointer
    so.type = x3DWriter::object_types::SPHERE;
    objects.emplace_back(so);
  }

  /// @brief add a set of Sphere with a list of center and constant radius and materials
  void x3DWriter::addSpheres(const vector<xvector<double>>& center, const double radius, const std::string& material) {
    for (const xvector<double>& c : center) {
      addSphere(c, radius, material);
    }
  }

  /// @brief add a box forming a unit cell to the scene
  /// @param lattice 3x3 matrix like xstructure.lattice
  /// @param radius radius of the cylinders forming the box
  /// @param material box material
  void x3DWriter::addLatticeBox(const xmatrix<double>& lattice, double radius, const std::string& material, bool axis) {
    vector<xvector<double>> points;
    points.push_back({0, 0, 0});
    for (int i = lattice.lrows; i <= lattice.urows; i++) {
      points.push_back(lattice(i));
    }
    points.push_back(points[1] + points[2]);
    points.push_back(points[1] + points[3]);
    points.push_back(points[2] + points[3]);
    points.push_back(points[1] + points[2] + points[3]);
    addSpheres(points, radius, material);
    vector<std::pair<uint, uint>> connections = {
        {1, 4},
        {1, 5},
        {2, 4},
        {2, 6},
        {3, 5},
        {3, 6},
        {7, 4},
        {7, 5},
        {7, 6}
    };
    ;
    if (!axis) {
      connections.emplace_back(std::pair<uint, uint>({0, 1}));
      connections.emplace_back(std::pair<uint, uint>({0, 2}));
      connections.emplace_back(std::pair<uint, uint>({0, 3}));
    }

    for (auto [base, apex] : connections) {
      addOpenCylinder(points[base], points[apex], radius, material);
    }
    if (axis) {
      addOpenCylinder(points[0], points[1], radius, "bm_red");
      addOpenCylinder(points[0], points[2], radius, "bm_green");
      addOpenCylinder(points[0], points[3], radius, "bm_blue");
    }
  }

  /// @brief add a OpenCylinder to the scene
  /// @param base start point
  /// @param apex end point
  /// @param radius cylinder radius
  /// @param material cylinder material
  void x3DWriter::addOpenCylinder(const xvector<double>& base, const xvector<double>& apex, double radius, const std::string& material) {
    const std::shared_ptr<x3DWriter::OpenCylinder> newCylinder = std::make_shared<x3DWriter::OpenCylinder>();
    const xvector<double> shift_base = base - scene_center;
    const xvector<double> shift_apex = apex - scene_center;
    xvector<double> axis = shift_base - shift_apex;
    newCylinder->base = shift_base;
    newCylinder->apex = shift_apex;
    newCylinder->center = shift_apex + (axis / 2);
    newCylinder->radius = radius;
    newCylinder->material = material;
    newCylinder->height = aurostd::modulus(axis);
    newCylinder->theta = pi / 2 + std::acos(axis(1) / newCylinder->height);
    newCylinder->phi = std::atan2(axis(3), axis(2));
    max_distance = std::max(max_distance, aurostd::modulus(shift_base));
    max_distance = std::max(max_distance, aurostd::modulus(shift_apex));
    x3DWriter::storage_object so;
    so.obj = newCylinder; // cast to a typeless void pointer
    so.type = x3DWriter::object_types::OPEN_CYLINDER;
    objects.emplace_back(so);
  }

  /// @brief add ConvexFacets to the Scene
  /// @param vertexes facet corners
  /// @param facets vertex index that form the facets
  /// @param material facet materials
  /// @param shift
  void x3DWriter::addConvexFacets(const vector<xvector<double>>& vertexes, const vector<vector<uint>>& facets, const std::string& material, const xvector<double>& shift) {
    const std::shared_ptr<x3DWriter::ConvexFacets> newFacet = std::make_shared<x3DWriter::ConvexFacets>();
    newFacet->material = material;

    vector<xvector<double>> shifted_vertexes;
    shifted_vertexes.reserve(vertexes.size());
    for (const xvector<double>& vertex : vertexes) {
      shifted_vertexes.push_back(vertex + shift - scene_center);
    }
    newFacet->vertexes = shifted_vertexes;

    newFacet->facets = facets;
    for (const auto& vertex : vertexes) {
      max_distance = std::max(max_distance, aurostd::modulus(vertex));
    }
    x3DWriter::storage_object so;
    so.obj = newFacet; // cast to a typeless void pointer
    so.type = x3DWriter::object_types::FACET;
    objects.emplace_back(so);
  }

  /// @brief add ConvexFacets with a grey glass material to the Scene
  /// @param vertexes facet corners
  /// @param facets vertex index that form the facets
  /// @param shift
  void x3DWriter::addConvexFacets(const vector<xvector<double>>& vertexes, const vector<vector<uint>>& facets, const xvector<double>& shift) {
    addConvexFacets(vertexes, facets, "bm_grey_glass", shift);
  }

  /// @brief add a new basic Material
  /// @param name name to reference the Material
  /// @param color as R,G,B xvector
  void x3DWriter::addMaterial(const std::string& name, const xvector<double>& color) {
    Material newMat;
    newMat.name = name;
    newMat.color = color;
    addMaterial(newMat);
  }

  /// @brief add a new Material
  /// @param newMaterial
  void x3DWriter::addMaterial(const Material& newMaterial) {
    const std::shared_ptr<x3DWriter::Material> newMat = std::make_shared<x3DWriter::Material>(newMaterial);
    x3DWriter::storage_object so;
    so.obj = newMat; // cast to a typeless void pointer
    so.type = x3DWriter::object_types::MATERIAL;
    objects.emplace_back(so);
  }

  /// @brief create a set of materials using the turbo colormap
  /// @param count number of materials to create
  ///
  /// the materials will be named "auto_cs_{idx}"
  vector<std::string> x3DWriter::addColorSpreadMaterial(const uint count, const std::string& cmap) {
    vector<std::string> created_materials;
    const uint color_step = 256 / count;
    for (uint color_index = 0; color_index < count; color_index++) {
      xvector<double> color;
      color(1) = aurostd::color::maps[cmap][color_index * color_step][0];
      color(2) = aurostd::color::maps[cmap][color_index * color_step][1];
      color(3) = aurostd::color::maps[cmap][color_index * color_step][2];
      created_materials.emplace_back("auto_cs_" + std::to_string(color_index));
      addMaterial(created_materials[color_index], color);
    }
    return created_materials;
  }
  template <typename Container> std::map<std::string, std::string> x3DWriter::addNamedColorSpreadMaterial(const Container& names) {
    std::vector<string> color_spread = addColorSpreadMaterial(names.size());
    std::map<std::string, std::string> mapped_colors;
    std::transform(names.begin(), names.end(), color_spread.begin(), std::inserter(mapped_colors, mapped_colors.end()), std::make_pair<const string&, const string&>);
    return mapped_colors;
  }
  template std::map<std::string, std::string> x3DWriter::addNamedColorSpreadMaterial(const std::vector<std::string>&);
  template std::map<std::string, std::string> x3DWriter::addNamedColorSpreadMaterial(const std::deque<std::string>&);

  /// @brief combine overlapping vertexes and facets
  /// @param vertexes original vertexes positions
  /// @param facets list of vertex indexes forming the facets
  /// @param new_vertexes updated vertexes positions
  /// @param new_facets updated list of vertex indexes forming the facets
  void x3DWriter::joinFacets(vector<xvector<double>>& vertexes, vector<vector<uint>>& facets, const vector<xvector<double>>& new_vertexes, const vector<vector<uint>>& new_facets) const {
    if (vertexes.empty() || facets.empty()) {
      vertexes = new_vertexes;
      facets = new_facets;
      return;
    }
    bool is_set = false;

    vector<uint> vertex_idx_map;
    const size_t original_vertex_count = vertexes.size();
    vector<uint> prepared_facet;
    for (const auto& new_vertexe : new_vertexes) {
      is_set = false;
      for (uint ov_id = 0; ov_id < original_vertex_count; ov_id++) {
        if (distance(vertexes[ov_id], new_vertexe) <= join_threshold) {
          vertex_idx_map.emplace_back(ov_id);
          is_set = true;
          break;
        }
      }
      if (is_set) {
        continue;
      }
      vertex_idx_map.emplace_back(vertexes.size());
      vertexes.emplace_back(new_vertexe);
    }

    for (const auto& new_facet : new_facets) {
      prepared_facet.clear();
      for (const uint v_id : new_facet) {
        prepared_facet.emplace_back(vertex_idx_map[v_id]);
      }
      facets.emplace_back(prepared_facet);
    }

    std::set<uint> facets_to_remove;
    const size_t facet_count = facets.size();
    bool overlap = true;
    // remove double facets to avoid z-fighting
    for (size_t row = 0; row < facet_count; row++) {
      for (size_t col = row + 1; col < facet_count; col++) {
        if (facets_to_remove.find(row) != facets_to_remove.end() || facets_to_remove.find(col) != facets_to_remove.end()) {
          continue;
        }
        overlap = true;
        if (facets[row].size() < facets[col].size()) {
          for (size_t vert_idx = 0; vert_idx < facets[row].size(); vert_idx++) {
            overlap = std::find(facets[col].begin(), facets[col].end(), facets[row][vert_idx]) != facets[col].end();
            if (!overlap) {
              break;
            }
          }
          if (overlap) {
            facets_to_remove.insert(row);
            break;
          }
        } else {
          for (size_t vert_idx = 0; vert_idx < facets[col].size(); vert_idx++) {
            overlap = std::find(facets[row].begin(), facets[row].end(), facets[col][vert_idx]) != facets[row].end();
            if (!overlap) {
              break;
            }
          }
          if (overlap) {
            facets_to_remove.insert(col);
          }
        }
      }
    }
    vector<vector<uint>> cleaned_facets;
    for (uint of_idx = 0; of_idx < facets.size(); of_idx++) {
      if (facets_to_remove.find(of_idx) == facets_to_remove.end()) {
        cleaned_facets.emplace_back(facets[of_idx]);
      }
    }
    facets = cleaned_facets;
  }

  /// @brief prepare an animation
  /// @param duration animation time in second
  /// @param out_folder folder to store the generated files
  /// @param fps frames per second
  /// @param lr move left to right
  ///
  /// In the out_folder a render.sh script that will create a file with the stem _Animation;
  /// for rendering [ffmpeg](https://ffmpeg.org/) and [tachyon](http://jedi.ks.uiuc.edu/~johns/tachyon/) is needed
  void x3DWriter::animate(float duration, const std::filesystem::path& out_folder, uint fps, bool lr) {
    std::filesystem::path file_path;
    std::error_code err;
    if (!std::filesystem::create_directories(out_folder, err)) {
      if (err.value() != 0) // exist is okay
      {
        const string message = "Can not create result folder " + out_folder.string() + ". Error: " + err.message();
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _FILE_ERROR_);
      }
    }
    const uint frame_count = duration * fps;
    std::stringstream command_content;
    const double angle_step = 360.0 / frame_count;
    for (uint frame_index = 0; frame_index < frame_count; frame_index++) {
      std::stringstream filename;
      filename << "frame_" << std::setw(6) << std::setfill('0') << frame_index;
      file_path = out_folder / filename.str();
      if (lr) {
        tachyon_camera_phi = frame_index * angle_step;
      } else {
        tachyon_camera_theta = frame_index * angle_step;
      }
      aurostd::string2file(toTachyon(), file_path.string() + ".dat");
      command_content << "tachyon " << file_path.filename() << ".dat -o " << file_path.filename() << ".ppm -format PPM -auto_skylight 0.7 -fullshade -numthreads 8" << endl;
    }

    switch (ani_type) {
      case (animation_format::MP4): {
        command_content << endl << "ffmpeg -r " << fps << " -f image2 -i frame_%06d.ppm -vcodec libx264 -crf 20 -pix_fmt yuv420p  -an _Animation.mp4" << endl;
        break;
      }
      case (animation_format::GIF): {
        command_content << endl << "ffmpeg -r " << fps << " -f image2 -i frame_%06d.ppm -filter_complex [0:v] fps=" << fps << ", split [a][b];[a] palettegen [p];[b][p] paletteuse _Animation.gif" << endl;
        break;
      }
      case (animation_format::WEBM): {
        command_content << endl << "ffmpeg -r " << fps << " -f image2 -i frame_%06d.ppm -c:v libvpx-vp9 -b:v 0 -crf 20 -pix_fmt yuv420p -pass 1 -an -f null /dev/null" << endl;
        command_content << endl << "ffmpeg -r " << fps << " -f image2 -i frame_%06d.ppm -c:v libvpx-vp9 -b:v 0 -crf 20 -pix_fmt yuv420p -pass 2 -an _Animation.webm" << endl;
        break;
      }
    }
    aurostd::string2file(command_content.str(), (out_folder / "render.sh"));
  }

  /// convert a material to x3d
  /// @return x3d version of a material (xml)
  std::string x3DWriter::x3d_material(const std::shared_ptr<x3DWriter::Material>& material) {
    std::stringstream mat;
    mat << "<Material";
    mat << " diffuseColor='" << material->color << "'";
    mat << " ambientIntensity='" << material->ambient << "'";
    mat << " shininess='" << material->specular << "'";
    mat << " transparency='" << 1.0 - material->opacity << "'/>";
    return mat.str();
  }

  /// @brief calculate the camera position (used fot tachyon)
  void x3DWriter::tachyon_calculate_camera() {
    if (tachyon_camera_orthographic) {
      tachyon_zoom = 1.0 / (max_distance * 1.05);
    } else {
      tachyon_zoom = 2.0;
    }

    const double camera_distance = std::tan((180.0 - (tachyon_camera_angle / tachyon_zoom)) / (180 * 2.0) * pi) * max_distance;
    tachyon_camera_position(1) = camera_distance * std::sin(tachyon_camera_phi / 180.0 * pi) * std::cos(tachyon_camera_theta / 180.0 * pi);
    tachyon_camera_position(2) = camera_distance * std::cos(tachyon_camera_phi / 180.0 * pi) * std::sin(tachyon_camera_theta / 180.0 * pi);
    tachyon_camera_position(3) = camera_distance * std::cos(tachyon_camera_phi / 180.0 * pi);
  }

  /// @brief save scene for the tachyon renderer
  /// use [tachyon](http://jedi.ks.uiuc.edu/~johns/tachyon/) to render the generated file
  std::string x3DWriter::toTachyon() {
    tachyon_calculate_camera();
    stringstream content;
    for (const auto& [name, meta_content] : meta) {
      content << "# " << name << ": " << meta_content << endl;
    }
    content << endl;
    content << "BEGIN_SCENE" << endl;
    content << "  RESOLUTION 720 720" << endl << endl;
    content << "CAMERA" << endl;
    if (tachyon_camera_orthographic) {
      content << "  PROJECTION ORTHOGRAPHIC" << endl;
    }
    content << "  ZOOM " << tachyon_zoom << "  ASPECTRATIO -1.0 ANTIALIASING 2 RAYDEPTH 12" << endl;
    if (tachyon_lattice_views_idx == -1) {
      content << "  CENTER " << tachyon_camera_position << endl;
      content << "  VIEWDIR " << -tachyon_camera_position / max_distance << endl;
      content << "  UPDIR 0 1 0" << endl;
    } else {
      const auto& [center, up] = tachyon_lattice_views[tachyon_lattice_views_idx];
      content << "CENTER " << center << endl;
      content << "VIEWDIR " << -aurostd::normalizeSumToOne(center) << endl;
      content << "UPDIR " << up << endl << endl;
    }
    content << "END_CAMERA" << endl << endl;
    content << "BACKGROUND 1.0 1.0 1.0" << endl;

    for (const auto& [type, obj] : objects) {
      if (type == x3DWriter::object_types::MATERIAL) {
        const std::shared_ptr<x3DWriter::Material> material = std::static_pointer_cast<x3DWriter::Material>(obj);
        content << "TEXDEF " << material->name << endl;
        content << "  AMBIENT " << material->ambient;
        content << "  DIFFUSE " << material->diffuse;
        content << "  SPECULAR " << material->specular;
        content << "  OPACITY " << material->opacity << endl;
        content << "  COLOR " << material->color;
        content << "  TEXFUNC 0" << endl;
      }
    }

    for (const auto& [type, obj] : objects) {
      switch (type) {
        case object_types::MATERIAL: {
          // written to the start of the file
          break;
        }
        case object_types::SPHERE: {
          const std::shared_ptr<x3DWriter::Sphere> sphere = std::static_pointer_cast<x3DWriter::Sphere>(obj);
          content << "SPHERE" << endl;
          content << "  CENTER " << sphere->center << endl;
          content << "  RAD " << sphere->radius << endl;
          content << "  " << sphere->material << endl << endl;
          break;
        }
        case object_types::FACET: {
          const std::shared_ptr<x3DWriter::ConvexFacets> facet_container = std::static_pointer_cast<x3DWriter::ConvexFacets>(obj);
          for (auto facet : facet_container->facets) {
            for (size_t tri_index = 1; tri_index < facet.size() - 1; tri_index++) {
              content << "TRI" << endl;
              content << "  V0 " << facet_container->vertexes[facet[0]] << endl;
              content << "  V1 " << facet_container->vertexes[facet[tri_index]] << endl;
              content << "  V2 " << facet_container->vertexes[facet[tri_index + 1]] << endl;
              content << "  " << facet_container->material << endl << endl;
            }
          }
          break;
        }
        case object_types::OPEN_CYLINDER: {
          const std::shared_ptr<x3DWriter::OpenCylinder> cylinder = std::static_pointer_cast<x3DWriter::OpenCylinder>(obj);
          content << "FCylinder" << endl;
          content << "  BASE " << cylinder->base << endl;
          content << "  APEX " << cylinder->apex << endl;
          content << "  RAD " << cylinder->radius << endl;
          content << "  " << cylinder->material << endl << endl;
          break;
        }
      }
    }

    content << "END_SCENE" << endl;

    return content.str();
  }

  /// @brief save scene as x3d (xml)
  /// @param include_xml create as stand alone xml
  /// @param replace_material switch if Materials are exported
  std::string x3DWriter::toX3D(const bool include_xml, const bool replace_material) {
    stringstream x3d_content;
    if (include_xml) {
      x3d_content << R"(<X3D profile="Immersive" version="3.0">)" << endl << endl;
      x3d_content << "<head>" << endl;
      for (const auto& [key, value] : meta) {
        x3d_content << "  <meta name=\"" << key << "\" content=\"" << value << "\"/>" << endl;
      }
      x3d_content << "</head>" << endl << endl;
    }

    x3d_content << "<Scene>" << endl;
    std::map<std::string, std::string> material_lookup;

    for (const auto& [type, obj] : objects) {
      if (type == x3DWriter::object_types::MATERIAL) {
        const std::shared_ptr<x3DWriter::Material> material = std::static_pointer_cast<x3DWriter::Material>(obj);
        material_lookup[material->name] = x3d_material(material);
      }
    }
    if (!replace_material) {
      for (const auto& [name, content] : material_lookup) {
        x3d_content << "<Appearance DEF='" << name << "'>";
        x3d_content << content;
        x3d_content << "</Appearance>" << endl;
      }
    }

    for (const auto& [type, obj] : objects) {
      switch (type) {
        case object_types::MATERIAL: {
          // written to the start of the file
          break;
        }
        case object_types::SPHERE: {
          const std::shared_ptr<x3DWriter::Sphere> sphere = std::static_pointer_cast<x3DWriter::Sphere>(obj);
          x3d_content << "<Transform translation='" << sphere->center << "'><Shape>" << endl;
          if (replace_material) {
            x3d_content << "  <Appearance>" << material_lookup[sphere->material] << "</Appearance>" << endl;
          } else {
            x3d_content << "  <Appearance USE='" << sphere->material << "'/>" << endl;
          }
          x3d_content << "  <Sphere radius='" << sphere->radius << "'/>" << endl;
          x3d_content << "</Shape></Transform>" << endl << endl;
          break;
        }
        case object_types::FACET: {
          const std::shared_ptr<x3DWriter::ConvexFacets> facet_container = std::static_pointer_cast<x3DWriter::ConvexFacets>(obj);
          x3d_content << "<Shape>" << endl;
          if (replace_material) {
            x3d_content << "  <Appearance>" << material_lookup[facet_container->material] << "</Appearance>" << endl;
          } else {
            x3d_content << "  <Appearance USE='" << facet_container->material << "'/>" << endl;
          }
          x3d_content << R"(  <IndexedFaceSet solid="false" colorPerVertex="false" normalPerVertex="false" coordIndex=")" << endl;
          for (const auto& facet : facet_container->facets) {
            x3d_content << "    ";
            for (const auto point_index : facet) {
              x3d_content << point_index << " ";
            }
            x3d_content << "-1" << endl; // end a facet
          }
          x3d_content << "  \">" << endl;
          x3d_content << "  <Coordinate point=\"" << endl;
          for (const auto& coordinate : facet_container->vertexes) {
            x3d_content << "    " << coordinate << endl;
          }
          x3d_content << "  \"/>" << endl;
          x3d_content << "  </IndexedFaceSet>" << endl;
          x3d_content << "</Shape>" << endl << endl;
          break;
        }
        case object_types::OPEN_CYLINDER: {
          const std::shared_ptr<x3DWriter::OpenCylinder> cylinder = std::static_pointer_cast<x3DWriter::OpenCylinder>(obj);
          x3d_content << "<Transform translation='" << cylinder->center << "'>" << endl;
          x3d_content << "<Transform rotation='1 0 0 " << cylinder->phi << "'><Transform rotation='0 0 1 " << cylinder->theta << "'><Shape>" << endl;
          if (replace_material) {
            x3d_content << "  <Appearance>" << material_lookup[cylinder->material] << "</Appearance>" << endl;
          } else {
            x3d_content << "  <Appearance USE='" << cylinder->material << "'/>" << endl;
          }
          x3d_content << "  <Cylinder radius='" << cylinder->radius << "' height='" << cylinder->height << "' top='FALSE' bottom='FALSE'/>" << endl;
          x3d_content << "</Shape></Transform></Transform></Transform>" << endl << endl;
          break;
        }
      }
    }

    x3d_content << "</Scene>" << endl;
    if (include_xml) {
      x3d_content << "</X3D>" << endl;
    }

    return x3d_content.str();
  }

  /// @brief save scene as x3d and embed it in a html file
  std::string x3DWriter::toHTML() {
    stringstream content;

    content << "<!DOCTYPE html>\n<html lang=\"en\">" << endl << endl;
    content << "<head>" << endl;
    content << "<meta charset=\"UTF-8\">" << endl;
    content << R"(<meta name="viewport" content="width=device-width, initial-scale=1.0">)" << endl;
    for (const auto& [name, meta_content] : meta) {
      content << "  <meta name=\"" << name << "\" content=\"" << meta_content << "\"/>" << endl;
    }
    content << "<title>AFLOW X3D</title>" << endl;
    content << "<script src='https://www.x3dom.org/download/x3dom.js'> </script>" << endl;
    // CDN makes it directly usable
    content << "<link rel='stylesheet' type='text/css' href='https://www.x3dom.org/download/x3dom.css'/>" << endl;
    content << "</head>" << endl << endl;
    content << "<body>" << endl;
    content << "<h1>AFLOW 3D viewer</h1>" << endl << endl;

    content << "<x3d width='1200px' height='1200px'>" << endl;
    content << toX3D(false, true) << endl << endl;
    content << "</x3d>" << endl;
    content << "</body>\n</html>" << endl;

    return content.str();
  }

// Plotter
  /// @class xplotter
  /// @brief collection of basic plotting routines
  ///
  /// @authors
  /// @mod{SD,20230217,creation}
  /// @mod{HE,20230217,creation}
  ///
  /// Basic usage
  /// @code
  /// std::vector<double> x = aurostd::xvector2vector(aurostd::linspace(-2*PI, 2*PI, (int)30));
  /// std::vector<double> y = aurostd::xvector2vector(aurostd::sin(x));
  /// aurostd::xplotter xplt;
  /// xplt.plot(x,y, {{"label", "sin"}, {"line width", "2pt"}, {"mark", "*"}});
  /// xplt.axis_properties["xlabel"] = "$x$-values";
  /// xplt.axis_properties["ylabel"] = "$y$-values";
  /// xplt.axis_properties["title"] = "A test plot";
  /// xplt.save("figure.pdf");
  /// @endcode

  // START Set static members of xplotter

  /// @brief system calls to convert a .pdf to a other format
  /// @note to read the comment string in a picture use "identify -format %c image.png"
  std::map<std::string, std::vector<std::pair<std::string, std::string>>> xplotter::convert_commands = {
      {".pdf",                                                                                                                                                                              {}              },
      {".png",
       {{"pdftocairo", "pdftocairo -q -png -singlefile -r 300 -transp render.pdf render"}, {"convert", "convert -density 300 render.pdf -set comment \"%comment%\" -quality 90 -sharpen 0x1.0  render.png"}}},
      {".tif",
       {
       {"convert", "convert -density 300 render.pdf -set comment \"%comment%\" -quality 90 -sharpen 0x1.0 -compress lzw render.tif"},
       {"pdftocairo", "pdftocairo -q -tiff -singlefile -r 300 -transp render.pdf render"}
       // second choice as not compressed
       }                                                                                                                                                                                                    },
      {".jpg",
       {
       {"pdftocairo", "pdftocairo -q -jpeg -singlefile -r 300 render.pdf render"},
       {"convert", "convert -density 300 render.pdf -set comment \"%comment%\"  -quality 85  -flatten -sharpen 0x1.0  render.jpg"},
       }                                                                                                                                                                                                    },
      {".eps",                                                                                                                                {{"pdftops", "pdftops -q -level3 -eps render.pdf render.eps"}}},
      { ".ps",                                                                                                                                      {{"pdftops", "pdftops -q -level3 render.pdf render.ps"}}},
      {".svg",                                                                                                                              {{"pdftocairo", "pdftocairo -svg -r 300 render.pdf render.svg"}}}
  };

  // property helper sets
  std::set<std::string> xplotter::excluded_line_properties = {"color", "label", "smooth", "line style", "forget plot", "a selected color", "marker scale", "scatter/use mapped color"};
  std::set<std::string> xplotter::line_bool_properties = {"smooth", "forget plot", "scatter", "only marks"};
  std::set<std::string> xplotter::excluded_axis_properties = {"pattern color", "bar type", "symbolic y coords", "symbolic x coords", "symbolic coords"};
  std::set<std::string> xplotter::string_axis_properties = {"xlabel", "ylabel", "title"};

  std::string xplotter::preamble =
      "\\documentclass[english]{revtex4-2}\n"
      "\\usepackage[active,pdftex,tightpage]{preview}\n"
      "\\usepackage[T1]{fontenc}\n"
      "\\usepackage[default]{opensans}\n"
      "\\renewcommand*\\familydefault{\\sfdefault}\n"
      "\\usepackage{tikz}\n"
      "\\usepackage{siunitx}\n"
      "\\usepackage{chemmacros}\n"
      "\\usepackage{pgfplots}\n"
      "\\usetikzlibrary{plotmarks,shapes,positioning,patterns}\n"
      "\\usepgfplotslibrary{units,ternary}\n"
      "\\pgfplotsset{plot coordinates/math parser=false}\n"
      "\\pgfplotsset{compat=newest,height=10cm, width=12cm}\n"
      "\\PreviewEnvironment{tikzpicture}\n"
      "\\pgfdeclarelayer{background}\n"
      "\\pgfdeclarelayer{foreground}\n"
      "\\pgfsetlayers{background,main,foreground}\n\n";

  std::map<xplotter::axis_types, std::map<std::string, std::string>> xplotter::axis_defaults = {
      {     axis_types::TWODIM,
       {{"colormap name", "{viridis}"},
       {"unit markings", "parenthesis"},
       {"xticklabel", R"({$\mathsf{\pgfmathprintnumber[precision=2,fixed,zerofill=false,assume math mode=true]{\tick}}$})"},
       {"yticklabel", R"({$\mathsf{\pgfmathprintnumber[precision=2,fixed,zerofill=false,assume math mode=true]{\tick}}$})"}}},
      {        axis_types::BAR,
       {
       {"bar type", "xbar"},
       {"unit markings", "parenthesis"},
       {"xticklabel", R"({$\mathsf{\pgfmathprintnumber[precision=2,fixed,zerofill=false,assume math mode=true]{\tick}}$})"},
       {"yticklabel", R"({$\mathsf{\pgfmathprintnumber[precision=2,fixed,zerofill=false,assume math mode=true]{\tick}}$})"},
       }                                                                                                                    },
      {axis_types::BARSYMBOLIC,
       {{"symbolic coords", ""},
       {"xticklabel", R"({$\mathsf{\pgfmathprintnumber[precision=2,fixed,zerofill=false,assume math mode=true]{\tick}}$})"},
       {"yticklabel", R"({$\mathsf{\pgfmathprintnumber[precision=2,fixed,zerofill=false,assume math mode=true]{\tick}}$})"}}},
      {    axis_types::TERNARY,
       {{"colorbar", "true"},
       {"colormap name", "{viridis}"},
       {"label style", "{sloped}"},
       {"unit markings", "parenthesis"},
       {"xticklabel", R"({$\mathsf{\pgfmathprintnumber[precision=2,fixed,zerofill=false,assume math mode=true]{\tick}}$})"},
       {"yticklabel", R"({$\mathsf{\pgfmathprintnumber[precision=2,fixed,zerofill=false,assume math mode=true]{\tick}}$})"},
       {"zticklabel", R"({$\mathsf{\pgfmathprintnumber[precision=2,fixed,zerofill=false,assume math mode=true]{\tick}}$})"}}}
  };

  std::map<xplotter::object_types, std::map<std::string, std::string>> xplotter::object_defaults = {
      {         object_types::SMOOTH,                                   {{"label", ""}, {"line width", "1pt"}, {"color", "auto"}, {"mark", ""}, {"mark size", "3pt"}, {"smooth", "true"}}},
      {        object_types::SCATTER, {{"label", ""}, {"scatter", "true"}, {"only marks", "true"}, {"color", "auto"}, {"colormap name", "viridis"}, {"mark", "*"}, {"marker scale", "3"}}},
      {           object_types::LINE,                                                                                                         {{"line width", "1pt"}, {"color", "black"}}},
      {            object_types::BAR,                                                                                                                                 {{"color", "auto"}}},
      {object_types::TERNARY_SCATTER,                                                              {{"point meta", "\\thisrow{m}"}, {"only marks", "true"}, {"nodes near coords*", "{}"}}},
      {object_types::TERNARY_CONTOUR,                                                                        {{"point meta", "\\thisrow{m}"}, {"contour prepared", "{labels over line}"}}},
      {   object_types::TERNARY_SURF,                                                                            {{"point meta", "\\thisrow{m}"}, {"surf", "true"}, {"shader", "interp"}}}
  };

  // END Set static members of xplotter

  /// @brief clean a EntryLoader object
  void xplotter::clear() {
    *this = {};
  } // calling the constructor

  /// @brief create a copy
  void xplotter::copy(const xplotter& xplt) {
    axisType = xplt.axisType;
    objects = xplt.objects;
    meta = xplt.meta;

    axis_properties = xplt.axis_properties;
    colorbar_properties = xplt.colorbar_properties;
    symbolic_coordinates = xplt.symbolic_coordinates;
    auto_colormap = xplt.auto_colormap;
    cyclic_name = xplt.cyclic_name;
  }

  /// @brief initialize the xplotter class
  xplotter::xplotter() {
    axisType = axis_types::NONE;
    auto_colormap = "turbo";
    cyclic_name = "estructure";
    meta.emplace_back("reference", "http://aflow.org");
    meta.emplace_back("generator", "AFLOW " + (std::string) AFLOW_VERSION);
    meta.emplace_back("created", get_datetime_formatted("-", true, " ", ":"));
    colorbar_properties = {
        {"xticklabel", R"({$\mathsf{\pgfmathprintnumber[precision=2,fixed,zerofill=false,assume math mode=true]{\tick}}$})"},
        {"yticklabel", R"({$\mathsf{\pgfmathprintnumber[precision=2,fixed,zerofill=false,assume math mode=true]{\tick}}$})"}
    };
  }

  /// @brief create a copy at initialization
  xplotter::xplotter(const xplotter& xplt) {
    copy(xplt);
  }

  ///@brief class de-constructor
  xplotter::~xplotter() {}

  // plot functions
  /// @brief create a smooth line plot
  /// @note the x values generated are a integer sequence starting at 1
  /// @note for direct line plot add to properties `{"smooth", "false"}`
  void xplotter::plot(const std::vector<double>& y, std::map<std::string, std::string> properties) {
    vector<double> x = aurostd::xvector2vector(aurostd::linspace(1, y.size(), (int) y.size()));
    xplotter::plot(x, y, properties);
  }

  /// @brief create a smooth line plot
  /// @note for direct line plot add to properties `{"smooth", "false"}`
  void xplotter::plot(const std::vector<double>& x, const std::vector<double>& y, std::map<std::string, std::string> properties) {
    if (x.size() != y.size()) {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "x and y vectors must have the same length", _INPUT_ILLEGAL_);
    }
    if (axisType == axis_types::NONE) {
      axisType = axis_types::TWODIM;
    } else if (axisType != axis_types::TWODIM) {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "xplotter.plot not defined for this axis", _RUNTIME_ERROR_);
    }
    std::shared_ptr<xplotter::xy_data> xyo = std::make_shared<xplotter::xy_data>();
    xyo->x = x;
    xyo->y = y;
    xyo->properties = properties;

    xplotter::storage_object so;
    so.obj = xyo; // cast to a typeless void pointer
    so.type = xplotter::object_types::SMOOTH;
    objects.emplace_back(so);
  }

  // scatter functions
  /// @brief create a scatter plot
  void xplotter::scatter(const std::vector<double>& x, const std::vector<double>& y, std::map<std::string, std::string> properties) {
    xplotter::scatter(x, y, {}, {}, properties);
  }

  /// @brief create a scatter plot
  /// \param z mapped to color
  void xplotter::scatter(const std::vector<double>& x, const std::vector<double>& y, const std::vector<double>& z, std::map<std::string, std::string> properties) {
    xplotter::scatter(x, y, z, {}, properties);
  }

  /// @brief create a scatter plot
  /// \param z mapped to color (can be empty '{}' to have just size depends)
  /// \param m mapped to size
  void xplotter::scatter(const std::vector<double>& x, const std::vector<double>& y, const std::vector<double>& z, const std::vector<double>& m, std::map<std::string, std::string> properties) {
    if (axisType == axis_types::NONE) {
      axisType = axis_types::TWODIM;
    }
    if (x.size() != y.size()) {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "x and y vectors must have the same length", _INPUT_ILLEGAL_);
    } else if (axisType != axis_types::TWODIM) {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "xplotter.scatter not defined for this axis", _RUNTIME_ERROR_);
    }
    std::shared_ptr<xplotter::xyzm_data> xyzmo = std::make_shared<xplotter::xyzm_data>();
    xyzmo->x = x;
    xyzmo->y = y;
    xyzmo->z = z;
    xyzmo->m = m;
    xyzmo->properties = properties;

    xplotter::storage_object so;
    so.obj = xyzmo; // cast to a typeless void pointer
    so.type = xplotter::object_types::SCATTER;
    objects.emplace_back(so);
  }

  /// @brief create a bar plot
  /// @param values bar heights
  /// @note the bar positions generated are a integer sequence starting at 1
  void xplotter::bar(const std::vector<double>& values, std::map<std::string, std::string> properties) {
    vector<double> edges = aurostd::xvector2vector(aurostd::linspace(1, values.size(), (int) values.size()));
    xplotter::bar(values, edges, properties);
  }

  /// @brief create a bar plot
  /// @param values bar heights
  /// @param edges bar positions (if values.size()==edges.size() edges are used as center value)
  void xplotter::bar(const std::vector<double>& values, const std::vector<double>& edges, std::map<std::string, std::string> properties) {
    if (values.size() != edges.size()) {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "values and edges must have the same length", _INPUT_ILLEGAL_);
    }
    if (axisType == axis_types::NONE) {
      axisType = axis_types::BAR;
    } else if (axisType != axis_types::BAR) {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "xplotter.bar() not defined for this axis", _RUNTIME_ERROR_);
    }
    std::shared_ptr<xplotter::xyl_data> xylo = std::make_shared<xplotter::xyl_data>();
    xylo->x = edges;
    xylo->y = values;
    xylo->properties = properties;

    xplotter::storage_object so;
    so.obj = xylo; // cast to a typeless void pointer
    so.type = xplotter::object_types::BAR;
    objects.emplace_back(so);
  }

  /// @brief create a symbolic bar plot
  /// @param values bar heights
  /// @param edges symbolic bar labels
  void xplotter::bar(const std::vector<double>& values, const std::vector<std::string>& labels, std::map<std::string, std::string> properties) {
    if (values.size() != labels.size()) {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "values and labels must have the same length", _INPUT_ILLEGAL_);
    }
    if (axisType == axis_types::NONE) {
      axisType = axis_types::BARSYMBOLIC;
    } else if (axisType != axis_types::BARSYMBOLIC) {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "xplotter.bar() not defined for this axis", _RUNTIME_ERROR_);
    }
    std::shared_ptr<xplotter::xyl_data> xylo = std::make_shared<xplotter::xyl_data>();
    xylo->y = values;
    xylo->l = labels;
    xylo->properties = properties;
    for (auto label : labels) {
      if (std::find(symbolic_coordinates.begin(), symbolic_coordinates.end(), label) == symbolic_coordinates.end()) {
        symbolic_coordinates.push_back(label);
      }
    }

    xplotter::storage_object so;
    so.obj = xylo; // cast to a typeless void pointer
    so.type = xplotter::object_types::BAR;
    objects.emplace_back(so);
  }

  /// @brief add a horizontal line
  void xplotter::hline(double y, std::map<std::string, std::string> properties) {
    std::shared_ptr<xplotter::line_data> lineo = std::make_shared<xplotter::line_data>();
    lineo->pos.push_back(y);
    lineo->properties = properties;

    xplotter::storage_object so;
    so.obj = lineo; // cast to a typeless void pointer
    so.type = xplotter::object_types::HLINE;
    objects.emplace_back(so);
  }

  /// @brief add a horizontal line
  /// @param s start
  /// @param e end
  void xplotter::hline(double y, double s, double e, std::map<std::string, std::string> properties) {
    xplotter::line({s, y}, {e, y}, properties);
  }

  /// @brief add a vertical line
  void xplotter::vline(double x, std::map<std::string, std::string> properties) {
    std::shared_ptr<xplotter::line_data> lineo = std::make_shared<xplotter::line_data>();
    lineo->pos.push_back(x);
    lineo->properties = properties;

    xplotter::storage_object so;
    so.obj = lineo; // cast to a typeless void pointer
    so.type = xplotter::object_types::VLINE;
    objects.emplace_back(so);
  }

  /// @brief add a vertical line
  /// @param s start
  /// @param e end
  void xplotter::vline(double x, double s, double e, std::map<std::string, std::string> properties) {
    xplotter::line({x, s}, {x, e}, properties);
  }

  /// @brief add a line
  /// \param pos {start_x, start_y, end_x, end_y}
  void xplotter::line(const vector<double>& pos, std::map<std::string, std::string> properties) {
    std::shared_ptr<xplotter::line_data> lineo = std::make_shared<xplotter::line_data>();
    lineo->pos = pos;
    lineo->properties = properties;

    xplotter::storage_object so;
    so.obj = lineo; // cast to a typeless void pointer
    so.type = xplotter::object_types::LINE;
    objects.emplace_back(so);
  }

  /// @brief add a line
  /// \param start {x, y}
  /// \param end {x, y}
  void xplotter::line(const std::pair<double, double>& start, const std::pair<double, double>& end, std::map<std::string, std::string> properties) {
    xplotter::line({start.first, start.second, end.first, end.second}, properties);
  }

  /// @brief save ternary plot data (private)
  void xplotter::ternary_plot(const std::vector<double>& x, const std::vector<double>& y, const std::vector<double>& z, const std::vector<double>& m, std::map<std::string, std::string> properties, const object_types& ternary_type) {
    if (x.size() != y.size() || x.size() != z.size() || x.size() != m.size()) {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "x, y, z and m vectors must have the same length", _INPUT_ILLEGAL_);
    }
    if (axisType == axis_types::NONE) {
      axisType = axis_types::TERNARY;
    } else if (axisType != axis_types::TERNARY) {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "xplotter.ternary_plot() not defined for this axis", _RUNTIME_ERROR_);
    }
    std::shared_ptr<xplotter::xyzm_data> xyzmo = std::make_shared<xplotter::xyzm_data>();
    xyzmo->x = x;
    xyzmo->y = y;
    xyzmo->z = z;
    xyzmo->m = m;
    xyzmo->properties = properties;

    xplotter::storage_object so;
    so.obj = xyzmo; // cast to a typeless void pointer
    so.type = ternary_type;
    objects.emplace_back(so);
  }

  /// @brief add ternary scatter plot
  /// @param m mapped to color
  /// @note `x+y+z` should be constant
  void xplotter::ternary_scatter(const std::vector<double>& x, const std::vector<double>& y, const std::vector<double>& z, const std::vector<double>& m, std::map<std::string, std::string> properties) {
    xplotter::ternary_plot(x, y, z, m, properties, object_types::TERNARY_SCATTER);
  }

  /// @brief add ternary surface plot
  /// @param m mapped to color
  /// @note `x+y+z` should be constant
  void xplotter::ternary_surf(const std::vector<double>& x, const std::vector<double>& y, const std::vector<double>& z, const std::vector<double>& m, std::map<std::string, std::string> properties) {
    xplotter::ternary_plot(x, y, z, m, properties, object_types::TERNARY_SURF);
  }

  /// @brief add ternary contour plot
  /// @param m mapped to color
  /// @note `x+y+z` should be constant
  void xplotter::ternary_contour(const std::vector<double>& x, const std::vector<double>& y, const std::vector<double>& z, const std::vector<double>& m, std::map<std::string, std::string> properties) {
    xplotter::ternary_plot(x, y, z, m, properties, object_types::TERNARY_CONTOUR);
  }

  /// @brief create the plot option string from given properties
  std::string xplotter::prepare_plot_options(std::map<std::string, std::string>& properties, uint& auto_color_count, uint& cyclic_color_count, const object_types plot_type) {
    stringstream options;
    if (properties["color"] == "auto") {
      properties["a selected color"] = "ac" + std::to_string(auto_color_count);
      // named "a selected color" to ensure it is written out early - avoids a problem in scatter
      auto_color_count += 1;
    } else if (properties["color"].substr(0, 6) == "cyclic") {
      uint cycle_index;
      if (properties["color"].size() > 6) {
        cycle_index = stoi(properties["color"].substr(6));
      } else {
        cycle_index = cyclic_color_count;
        cyclic_color_count += 1;
      }
      properties["a selected color"] = "cyc" + std::to_string(cycle_index);
      // named "a selected color" to ensure it is written out early - avoids a problem in scatter
    } else {
      properties["a selected color"] = properties["color"];
    }
    for (auto entry : properties) {
      if ((!entry.second.empty()) && ((excluded_line_properties.find(entry.first) == excluded_line_properties.end()) && (line_bool_properties.find(entry.first) == line_bool_properties.end()))) {
        options << entry.first << " = " << entry.second << ", ";
      } else if ((entry.first == "a selected color") && (!entry.second.empty())) {
        if (plot_type == object_types::BAR) {
          if (properties.find("pattern") == properties.end()) {
            options << "fill = " << entry.second << ", ";
          } else {
            options << "pattern color = " << entry.second << ", ";
          }
        } else if (plot_type == object_types::SCATTER) {
          if (properties.find("scatter/use mapped color") != properties.end()) {
            options << "scatter/use mapped color = {fill = " << entry.second << "}, ";
          } else {
            options << "mark options = {fill = " << entry.second << "}, ";
          }
        } else {
          options << entry.second << ", ";
        }
      } else if (entry.first == "line style") {
        options << entry.second << ", ";
      } else if ((line_bool_properties.find(entry.first) != line_bool_properties.end()) && (entry.second == "true")) {
        options << entry.first << ", ";
      }
    }
    options.seekp(-2, options.cur);
    return options.str();
  }

  ///@brief create the colorbar option string
  std::string xplotter::prepare_colorbar_options(const std::map<std::string, std::string>& properties) {
    stringstream options;
    options << "{";
    for (auto entry : properties) {
      if (!entry.second.empty()) {
        options << "  " << entry.first << " = " << entry.second << ",\n";
      }
    }
    options.seekp(-2, options.cur);
    options << "} ";
    return options.str();
  }

  ///@brief create the axis option string
  std::string xplotter::prepare_axis_options(const std::map<std::string, std::string>& properties) {
    stringstream options;
    options << "[\n";
    for (auto entry : properties) {
      if (!entry.second.empty() && (excluded_axis_properties.find(entry.first) == excluded_axis_properties.end()) && (string_axis_properties.find(entry.first) == string_axis_properties.end())) {
        options << "  " << entry.first << " = " << entry.second << ",\n";
      } else if (string_axis_properties.find(entry.first) != string_axis_properties.end()) {
        options << "  " << entry.first << " = {" << entry.second << "},\n";
      } else if (entry.first == "bar type") {
        options << entry.second << ",\n";
      } else if (entry.first == "symbolic coords") {
        std::string sym_coord;
        for (auto symbol : symbolic_coordinates) {
          sym_coord += "{" + symbol + "},";
        }
        if (properties.at("bar type") == "xbar") {
          options << "  symbolic y coords = {" << sym_coord << "},\n";
        } else {
          options << "  symbolic x coords = {" << sym_coord << "},\n  xticklabel style={rotate=90},\n";
        }
      }
    }
    std::map<std::string, std::string>::const_iterator cb_iter = properties.find("colorbar");
    if (cb_iter != properties.end()) {
      if (cb_iter->second == "true") {
        options << "colorbar style = " << prepare_colorbar_options(colorbar_properties) << ",\n";
      }
    }
    options.seekp(-2, options.cur);
    options << "]\n";
    return options.str();
  }

  ///@brief convert `meta` into a comma seperated string
  ///@param exclude skip entries with these keys
  std::string xplotter::prepare_meta_csv(const std::set<std::string>& exclude) {
    std::string meta_collection;
    for (auto meta_entry : meta) {
      if (exclude.find(meta_entry.first) == exclude.end()) {
        meta_collection += meta_entry.first + ":" + meta_entry.second + ", ";
      }
    }
    return meta_collection.substr(0, meta_collection.size() - 1);
  }

  /// @brief create figure source tikz + pgfplots
  std::string xplotter::toString() {
    if (objects.empty()) {
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, "No plots added!", _INPUT_MISSING_);
    }
    if (axisType == axis_types::NONE) {
      axisType = axis_types::TWODIM;
    }
    for (auto entry : axis_defaults[axisType]) {
      axis_properties.emplace(entry);
    }

    std::vector<std::string> legend;
    stringstream plot;
    uint auto_color_count = 0;
    uint cyclic_color_count = 0;
    std::string plot_style;
    for (auto so : objects) {
      switch (so.type) {
        case object_types::SMOOTH: {
          std::shared_ptr<xplotter::xy_data> xyo = std::static_pointer_cast<xplotter::xy_data>(so.obj);
          for (auto entry : object_defaults[object_types::SMOOTH]) {
            xyo->properties.emplace(entry);
          }
          plot_style = prepare_plot_options(xyo->properties, auto_color_count, cyclic_color_count, so.type);
          plot << "\\addplot[" << plot_style << "]";
          plot << " table[header=false] {%\n";
          for (size_t value_idx = 0; value_idx < xyo->x.size(); value_idx++) {
            plot << xyo->x[value_idx] << " " << xyo->y[value_idx] << "\n";
          }
          plot << "};\n\n";
          if (!xyo->properties["label"].empty()) {
            legend.push_back("\\addlegendimage{" + plot_style + "}\n\\addlegendentry{" + xyo->properties["label"] + "}\n");
          }
          break;
        }
        case object_types::SCATTER: {
          std::shared_ptr<xplotter::xyzm_data> xyzmo = std::static_pointer_cast<xplotter::xyzm_data>(so.obj);
          for (auto entry : object_defaults[object_types::SCATTER]) {
            xyzmo->properties.emplace(entry);
          }
          if (xyzmo->z.empty()) {
            xyzmo->properties["scatter/use mapped color"] = "";
          } else {
            axis_properties["colorbar"] = "true";
          }
          if (!xyzmo->m.empty()) {
            xyzmo->properties["visualization depends on"] = R"(\thisrow{m}\as\valuem)";
            std::string m_min = std::to_string(aurostd::min(xyzmo->m));
            std::string m_scale = std::to_string(std::sqrt(aurostd::max(xyzmo->m) - aurostd::min(xyzmo->m)));
            xyzmo->properties["scatter/@pre marker code/.append style"] = "{/tikz/mark size=((sqrt(\\valuem-" + m_min + "))/" + m_scale + ")*5+2}";
          }
          if (!xyzmo->z.empty() || !xyzmo->m.empty()) {
            xyzmo->properties["scatter src"] = "explicit";
          }
          plot_style = prepare_plot_options(xyzmo->properties, auto_color_count, cyclic_color_count, so.type);
          plot << "\\addplot[" << plot_style << "]";
          if (xyzmo->z.empty() && xyzmo->m.empty()) {
            plot << " table[header=false] {%\n";
            for (size_t value_idx = 0; value_idx < xyzmo->x.size(); value_idx++) {
              plot << xyzmo->x[value_idx] << " " << xyzmo->y[value_idx] << "\n";
            }
          } else if (xyzmo->m.empty()) {
            plot << " table[x=x,y=y,meta=z] {%\nx y z\n";
            for (size_t value_idx = 0; value_idx < xyzmo->x.size(); value_idx++) {
              plot << xyzmo->x[value_idx] << " " << xyzmo->y[value_idx] << " " << xyzmo->z[value_idx] << "\n";
            }
          } else if (xyzmo->z.empty()) {
            plot << " table[x=x,y=y,meta=m] {%\nx y m\n";
            for (size_t value_idx = 0; value_idx < xyzmo->x.size(); value_idx++) {
              plot << xyzmo->x[value_idx] << " " << xyzmo->y[value_idx] << " " << xyzmo->m[value_idx] << "\n";
            }
          } else {
            plot << " table[x=x,y=y,meta=z] {%\nx y z m\n";

            for (size_t value_idx = 0; value_idx < xyzmo->x.size(); value_idx++) {
              plot << xyzmo->x[value_idx] << " " << xyzmo->y[value_idx] << " " << xyzmo->z[value_idx] << " " << xyzmo->m[value_idx] << "\n";
            }
          }

          plot << "};\n\n";
          if (!xyzmo->properties["label"].empty()) {
            legend.push_back("\\addlegendimage{" + plot_style + "}\n\\addlegendentry{" + xyzmo->properties["label"] + "}\n");
          }
          break;
        }
        case object_types::LINE:
        case object_types::HLINE:
        case object_types::VLINE: {
          std::shared_ptr<xplotter::line_data> lineo = std::static_pointer_cast<xplotter::line_data>(so.obj);
          for (auto entry : object_defaults[object_types::LINE]) {
            lineo->properties.emplace(entry);
          }
          plot_style = prepare_plot_options(lineo->properties, auto_color_count, cyclic_color_count, so.type);
          if (so.type == object_types::VLINE) {
            plot << "\\draw[" << plot_style << "] ({rel axis cs:0,0} -| {axis cs:" << lineo->pos[0] << ",0}) -- ({rel axis cs:0,1} -| {axis cs:" << lineo->pos[0] << ",0});\n\n";
          } else if (so.type == object_types::HLINE) {
            plot << "\\draw[" << plot_style << "] ({rel axis cs:0,0} |- {axis cs:0," << lineo->pos[0] << "}) -- ({rel axis cs:1,0} |- {axis cs:0," << lineo->pos[0] << "});\n\n";
          } else {
            plot << "\\draw[" << plot_style << "] ({axis cs:" << lineo->pos[0] << "," << lineo->pos[1] << "}) -- ({axis cs:" << lineo->pos[2] << "," << lineo->pos[3] << "});\n\n";
          }
          if (!lineo->properties["label"].empty()) {
            legend.push_back("\\addlegendimage{" + plot_style + "}\n\\addlegendentry{" + lineo->properties["label"] + "}\n");
          }
          break;
        }
        case object_types::BAR: {
          std::shared_ptr<xplotter::xyl_data> xylo = std::static_pointer_cast<xplotter::xyl_data>(so.obj);
          for (auto entry : object_defaults[object_types::BAR]) {
            xylo->properties.emplace(entry);
          }

          if ((axis_properties["bar type"] == "ybar") || (axis_properties["bar type"] == "ybar interval")) {
            axis_properties["xticklabel"] = "";
          } else {
            axis_properties["yticklabel"] = "";
          }
          plot_style = prepare_plot_options(xylo->properties, auto_color_count, cyclic_color_count, so.type);
          plot << "\\addplot[" << plot_style << "] table[header=false] {%\n";
          if (!xylo->x.empty()) {
            if ((axis_properties["bar type"] == "ybar") || (axis_properties["bar type"] == "ybar interval")) {
              for (size_t value_idx = 0; value_idx < xylo->y.size(); value_idx++) {
                plot << xylo->x[value_idx] << " " << xylo->y[value_idx] << "\n";
              }
              if (xylo->y.size() + 1 == xylo->x.size()) {
                plot << xylo->x.back() << " " << xylo->y.back() << "\n";
                axis_properties["bar type"] = "ybar interval";
                axis_properties["xticklabel interval boundaries"] = "true";
                axis_properties["x tick label style"] = "{rotate=90, anchor=east}";
              }
            } else {
              for (size_t value_idx = 0; value_idx < xylo->y.size(); value_idx++) {
                plot << xylo->y[value_idx] << " " << xylo->x[value_idx] << "\n";
              }
              if (xylo->y.size() + 1 == xylo->x.size()) {
                plot << xylo->y.back() << " " << xylo->x.back() << "\n";
                axis_properties["bar type"] = "xbar interval";
                axis_properties["yticklabel interval boundaries"] = "true";
              }
            }
          } else if (!xylo->l.empty()) {
            if ((axis_properties["bar type"] == "ybar") || (axis_properties["bar type"] == "ybar interval")) {
              for (size_t value_idx = 0; value_idx < xylo->l.size(); value_idx++) {
                plot << xylo->l[value_idx] << " " << xylo->y[value_idx] << "\n";
              }
            } else {
              for (size_t value_idx = 0; value_idx < xylo->l.size(); value_idx++) {
                plot << xylo->y[value_idx] << " " << xylo->l[value_idx] << "\n";
              }
            }
          }
          plot << "};\n\n";
          if (!xylo->properties["label"].empty()) {
            legend.push_back("\\addlegendimage{" + plot_style + "}\n\\addlegendentry{" + xylo->properties["label"] + "}\n");
          }
          break;
        }
        case object_types::TERNARY_SCATTER:
        case object_types::TERNARY_CONTOUR:
        case object_types::TERNARY_SURF: {
          std::shared_ptr<xplotter::xyzm_data> xyzmo = std::static_pointer_cast<xplotter::xyzm_data>(so.obj);
          for (auto entry : object_defaults[so.type]) {
            xyzmo->properties.emplace(entry);
          }
          plot_style = prepare_plot_options(xyzmo->properties, auto_color_count, cyclic_color_count, so.type);
          plot << "\\addplot3";
          if (so.type == object_types::TERNARY_SCATTER) {
            plot << "+";
          }
          plot << "[" << plot_style << "]";
          plot << " table {%\nx y z m\n";
          for (size_t value_idx = 0; value_idx < xyzmo->x.size(); value_idx++) {
            plot << xyzmo->x[value_idx] << " " << xyzmo->y[value_idx] << " " << xyzmo->z[value_idx] << " " << xyzmo->m[value_idx] << "\n";
          }
          plot << "};\n\n";
          if (!xyzmo->properties["label"].empty()) {
            legend.push_back("\\addlegendimage{" + plot_style + "}\n\\addlegendentry{" + xyzmo->properties["label"] + "}\n");
          }
          break;
        }
        default: {
          break;
        }
      }
    }

    //create TeX file content

    stringstream content;

    //add infile meta data
    for (auto meta_entry : meta) {
      content << "% " << meta_entry.first << " = " << meta_entry.second << endl;
    }
    if (!meta.empty()) {
      content << "% --------------------------------------------- %\n\n";
    }

    //write preamble
    content << preamble;

    //write PDF meta data
    //available fields: /Title /Author /Creator /Producer /CreationDate /ModDate /Subject /Keywords
    if (!meta.empty()) {
      content << "\n\\pdfinfo{\n";
      if (axis_properties.find("title") != axis_properties.end()) {
        content << "  /Title (" << axis_properties["title"] << ")";
      }
      for (auto meta_entry : meta) {
        if (meta_entry.first == "generator") {
          content << "  /Creator (" << meta_entry.second << ")\n";
          break;
        }
      }
      content << "  /CreationDate (D:" << get_datetime_utc_formatted("", true, "", "") << ")\n";
      content << "  /Keywords (" << prepare_meta_csv({"generator"}) << ")\n";
      content << "}\n\n";
    }

    //write colors
    if (auto_color_count > 0) {
      uint color_step = 256 / (auto_color_count * 2);
      double color[3];
      for (uint color_index = 0; color_index < auto_color_count; color_index++) {
        color[0] = aurostd::color::maps[auto_colormap][(color_index * 2 + 1) * color_step][0];
        color[1] = aurostd::color::maps[auto_colormap][(color_index * 2 + 1) * color_step][1];
        color[2] = aurostd::color::maps[auto_colormap][(color_index * 2 + 1) * color_step][2];
        content << "\\definecolor{" << ("ac" + std::to_string(color_index)) << "}{rgb}{" << color[0] << "," << color[1] << "," << color[2] << "}\n";
      }
    }

    uint cyclic_index;
    for (uint color_index = 0; color_index < std::max(aurostd::color::cyclics[cyclic_name].size(), (size_t) cyclic_color_count); color_index++) {
      cyclic_index = color_index % aurostd::color::cyclics[cyclic_name].size();
      content << "\\definecolor{" << ("cyc" + std::to_string(color_index)) << "}" + aurostd::color::getCyclicTikz(cyclic_name, cyclic_index) + "\n";
    }

    // start document
    content << "\n\\begin{document}\n"
               "\\begin{tikzpicture}[]\n";

    if (axisType == axis_types::TERNARY) {
      content << "\n\n\\begin{ternaryaxis}\n";
    } else {
      content << "\n\n\\begin{axis}\n";
    }
    content << prepare_axis_options(axis_properties);
    for (const std::string& entry : legend) {
      content << entry << "\n";
    }
    content << plot.str();

    //add nodes
    if (!nodes.empty()) {
      content << "\\begin{pgfonlayer}{foreground}" << endl;
      ;
      for (const auto& node : nodes) {
        content << node << endl;
      }
      content << "\\end{pgfonlayer}{background}" << endl;
    }
    if (axisType == axis_types::TERNARY) {
      content << "\n\n\\end{ternaryaxis}\n";
    } else {
      content << "\n\n\\end{axis}\n";
    }
    content << "\n\\end{tikzpicture}";
    content << endl << "%LOCATION FOR TABLE:" << endl; //add location for optional table
    content << "\\end{document}";
    return content.str();
  }

  /// @brief save the figure to the filesystem
  /// @note can always create `.tex` file
  /// @note `pdflatex` needs to be installed to render a figure
  /// @note to create image formats (`.png`, `.tif`, `.jpg`) either Imagemagick or Poppler needs to be installed
  /// @note to create vector formats (`.svg`, `.eps`, `.ps`) Poppler needs to be installed
  ///
  /// @see
  /// @xlink{https://imagemagick.org}
  /// @xlink{https://poppler.freedesktop.org}
  /// @xlink{Poppler MAC, https://formulae.brew.sh/formula/poppler}
  /// @xlink{Poppler-utils Linux, https://packages.debian.org/sid/poppler-utils}
  void xplotter::save(const std::string& file_path) {
    // Parse the file path and set defaults
    std::array<std::string, 3> pns = aurostd::splitFilePath(file_path); // parent, name, suffix
    if (pns[2].empty() || pns[2] == ".") {
      pns[2] = ".pdf";
    }
    if (pns[1].empty()) {
      pns[1] = "aflow_plot";
    }
    std::string render_type = pns[2];
    if (render_type == ".tiff") {
      render_type = ".tif";
    } else if (render_type == ".jpeg") {
      render_type = ".jpg";
    }
    // if just .tex is wanted
    if (render_type == ".tex") {
      aurostd::string2file(this->toString(), file_path);
      return;
    }

    // ensure that pdflatex is available
    if (!aurostd::IsCommandAvailable("pdflatex")) {
      cerr << "xplotter: pdflatex not in $PATH - falling back to .tex" << endl;
      aurostd::string2file(this->toString(), pns[0] + pns[1] + ".tex");
      return;
    }

    // ensure that the file type is supported
    if (convert_commands.find(render_type) == convert_commands.end()) {
      cerr << "xplotter: no known conversion to " << pns[2] << " - falling back to .pdf" << endl;
      pns[2] = ".pdf";
      render_type = ".pdf";
    }

    // render PDF in tmp folder
    std::string final_file_path = pns[0] + pns[1] + pns[2];
    std::string active_work_folder = TmpDirectoryCreate("xplotter", "", false);
    aurostd::string2file(this->toString(), active_work_folder + "/render.tex");
    std::string command = "pdflatex -interaction=batchmode -output-dir " + active_work_folder + " " + active_work_folder + "/render.tex";
    aurostd::execute2string(command); // 2string for silence

    // copy the pdf or use convert to generate an alternative format
    if (pns[2] == ".pdf") {
      aurostd::file2file(active_work_folder + "/render.pdf", final_file_path);
    } else {
      bool done = false;
      for (auto convert_option : convert_commands.at(render_type)) {
        if (aurostd::IsCommandAvailable(convert_option.first)) {
          std::string comment = aurostd::StringSubst(prepare_meta_csv(), "\"", "\\\"");
          // command would fail if a " left unescaped
          std::string convert = aurostd::StringSubst(convert_option.second, "%comment%", comment);
          aurostd::execute2string("cd " + active_work_folder + "; " + convert); // 2string for silence
          aurostd::file2file(active_work_folder + "/render" + render_type, final_file_path);
          done = true;
          break;
        } else {
          cerr << "xplotter: " << convert_option.first << " not in $PATH";
        }
      }
      if (!done) {
        cerr << "xplotter: could not convert to " << render_type << " - falling back to .pdf" << endl;
        aurostd::file2file(active_work_folder + "/render.pdf", pns[0] + pns[1] + ".pdf");
      }
    }

    //clean the temp file
    aurostd::RemoveDirectory(active_work_folder);
  }
} // namespace aurostd

//HE20230215 BEGIN color
namespace aurostd {
  std::map<std::string, std::array<std::array<float, 3>, 256>> color::maps = color::fillColormaps();

  std::map<std::string, std::array<std::array<float, 3>, 256>> color::fillColormaps() {
    std::map<std::string, std::array<std::array<float, 3>, 256>> storage;

    // VIRIDIS
    storage.insert({"viridis",
                    {{
                        {0.267004, 0.004874, 0.329415}, {0.26851, 0.009605, 0.335427},  {0.269944, 0.014625, 0.341379}, {0.271305, 0.019942, 0.347269}, {0.272594, 0.025563, 0.353093},
                        {0.273809, 0.031497, 0.358853}, {0.274952, 0.037752, 0.364543}, {0.276022, 0.044167, 0.370164}, {0.277018, 0.050344, 0.375715}, {0.277941, 0.056324, 0.381191},
                        {0.278791, 0.062145, 0.386592}, {0.279566, 0.067836, 0.391917}, {0.280267, 0.073417, 0.397163}, {0.280894, 0.078907, 0.402329}, {0.281446, 0.08432, 0.407414},
                        {0.281924, 0.089666, 0.412415}, {0.282327, 0.094955, 0.417331}, {0.282656, 0.100196, 0.42216},  {0.28291, 0.105393, 0.426902},  {0.283091, 0.110553, 0.431554},
                        {0.283197, 0.11568, 0.436115},  {0.283229, 0.120777, 0.440584}, {0.283187, 0.125848, 0.44496},  {0.283072, 0.130895, 0.449241}, {0.282884, 0.13592, 0.453427},
                        {0.282623, 0.140926, 0.457517}, {0.28229, 0.145912, 0.46151},   {0.281887, 0.150881, 0.465405}, {0.281412, 0.155834, 0.469201}, {0.280868, 0.160771, 0.472899},
                        {0.280255, 0.165693, 0.476498}, {0.279574, 0.170599, 0.479997}, {0.278826, 0.17549, 0.483397},  {0.278012, 0.180367, 0.486697}, {0.277134, 0.185228, 0.489898},
                        {0.276194, 0.190074, 0.493001}, {0.275191, 0.194905, 0.496005}, {0.274128, 0.199721, 0.498911}, {0.273006, 0.20452, 0.501721},  {0.271828, 0.209303, 0.504434},
                        {0.270595, 0.214069, 0.507052}, {0.269308, 0.218818, 0.509577}, {0.267968, 0.223549, 0.512008}, {0.26658, 0.228262, 0.514349},  {0.265145, 0.232956, 0.516599},
                        {0.263663, 0.237631, 0.518762}, {0.262138, 0.242286, 0.520837}, {0.260571, 0.246922, 0.522828}, {0.258965, 0.251537, 0.524736}, {0.257322, 0.25613, 0.526563},
                        {0.255645, 0.260703, 0.528312}, {0.253935, 0.265254, 0.529983}, {0.252194, 0.269783, 0.531579}, {0.250425, 0.27429, 0.533103},  {0.248629, 0.278775, 0.534556},
                        {0.246811, 0.283237, 0.535941}, {0.244972, 0.287675, 0.53726},  {0.243113, 0.292092, 0.538516}, {0.241237, 0.296485, 0.539709}, {0.239346, 0.300855, 0.540844},
                        {0.237441, 0.305202, 0.541921}, {0.235526, 0.309527, 0.542944}, {0.233603, 0.313828, 0.543914}, {0.231674, 0.318106, 0.544834}, {0.229739, 0.322361, 0.545706},
                        {0.227802, 0.326594, 0.546532}, {0.225863, 0.330805, 0.547314}, {0.223925, 0.334994, 0.548053}, {0.221989, 0.339161, 0.548752}, {0.220057, 0.343307, 0.549413},
                        {0.21813, 0.347432, 0.550038},  {0.21621, 0.351535, 0.550627},  {0.214298, 0.355619, 0.551184}, {0.212395, 0.359683, 0.55171},  {0.210503, 0.363727, 0.552206},
                        {0.208623, 0.367752, 0.552675}, {0.206756, 0.371758, 0.553117}, {0.204903, 0.375746, 0.553533}, {0.203063, 0.379716, 0.553925}, {0.201239, 0.38367, 0.554294},
                        {0.19943, 0.387607, 0.554642},  {0.197636, 0.391528, 0.554969}, {0.19586, 0.395433, 0.555276},  {0.1941, 0.399323, 0.555565},   {0.192357, 0.403199, 0.555836},
                        {0.190631, 0.407061, 0.556089}, {0.188923, 0.41091, 0.556326},  {0.187231, 0.414746, 0.556547}, {0.185556, 0.41857, 0.556753},  {0.183898, 0.422383, 0.556944},
                        {0.182256, 0.426184, 0.55712},  {0.180629, 0.429975, 0.557282}, {0.179019, 0.433756, 0.55743},  {0.177423, 0.437527, 0.557565}, {0.175841, 0.44129, 0.557685},
                        {0.174274, 0.445044, 0.557792}, {0.172719, 0.448791, 0.557885}, {0.171176, 0.45253, 0.557965},  {0.169646, 0.456262, 0.55803},  {0.168126, 0.459988, 0.558082},
                        {0.166617, 0.463708, 0.558119}, {0.165117, 0.467423, 0.558141}, {0.163625, 0.471133, 0.558148}, {0.162142, 0.474838, 0.55814},  {0.160665, 0.47854, 0.558115},
                        {0.159194, 0.482237, 0.558073}, {0.157729, 0.485932, 0.558013}, {0.15627, 0.489624, 0.557936},  {0.154815, 0.493313, 0.55784},  {0.153364, 0.497, 0.557724},
                        {0.151918, 0.500685, 0.557587}, {0.150476, 0.504369, 0.55743},  {0.149039, 0.508051, 0.55725},  {0.147607, 0.511733, 0.557049}, {0.14618, 0.515413, 0.556823},
                        {0.144759, 0.519093, 0.556572}, {0.143343, 0.522773, 0.556295}, {0.141935, 0.526453, 0.555991}, {0.140536, 0.530132, 0.555659}, {0.139147, 0.533812, 0.555298},
                        {0.13777, 0.537492, 0.554906},  {0.136408, 0.541173, 0.554483}, {0.135066, 0.544853, 0.554029}, {0.133743, 0.548535, 0.553541}, {0.132444, 0.552216, 0.553018},
                        {0.131172, 0.555899, 0.552459}, {0.129933, 0.559582, 0.551864}, {0.128729, 0.563265, 0.551229}, {0.127568, 0.566949, 0.550556}, {0.126453, 0.570633, 0.549841},
                        {0.125394, 0.574318, 0.549086}, {0.124395, 0.578002, 0.548287}, {0.123463, 0.581687, 0.547445}, {0.122606, 0.585371, 0.546557}, {0.121831, 0.589055, 0.545623},
                        {0.121148, 0.592739, 0.544641}, {0.120565, 0.596422, 0.543611}, {0.120092, 0.600104, 0.54253},  {0.119738, 0.603785, 0.5414},   {0.119512, 0.607464, 0.540218},
                        {0.119423, 0.611141, 0.538982}, {0.119483, 0.614817, 0.537692}, {0.119699, 0.61849, 0.536347},  {0.120081, 0.622161, 0.534946}, {0.120638, 0.625828, 0.533488},
                        {0.12138, 0.629492, 0.531973},  {0.122312, 0.633153, 0.530398}, {0.123444, 0.636809, 0.528763}, {0.12478, 0.640461, 0.527068},  {0.126326, 0.644107, 0.525311},
                        {0.128087, 0.647749, 0.523491}, {0.130067, 0.651384, 0.521608}, {0.132268, 0.655014, 0.519661}, {0.134692, 0.658636, 0.517649}, {0.137339, 0.662252, 0.515571},
                        {0.14021, 0.665859, 0.513427},  {0.143303, 0.669459, 0.511215}, {0.146616, 0.67305, 0.508936},  {0.150148, 0.676631, 0.506589}, {0.153894, 0.680203, 0.504172},
                        {0.157851, 0.683765, 0.501686}, {0.162016, 0.687316, 0.499129}, {0.166383, 0.690856, 0.496502}, {0.170948, 0.694384, 0.493803}, {0.175707, 0.6979, 0.491033},
                        {0.180653, 0.701402, 0.488189}, {0.185783, 0.704891, 0.485273}, {0.19109, 0.708366, 0.482284},  {0.196571, 0.711827, 0.479221}, {0.202219, 0.715272, 0.476084},
                        {0.20803, 0.718701, 0.472873},  {0.214, 0.722114, 0.469588},    {0.220124, 0.725509, 0.466226}, {0.226397, 0.728888, 0.462789}, {0.232815, 0.732247, 0.459277},
                        {0.239374, 0.735588, 0.455688}, {0.24607, 0.73891, 0.452024},   {0.252899, 0.742211, 0.448284}, {0.259857, 0.745492, 0.444467}, {0.266941, 0.748751, 0.440573},
                        {0.274149, 0.751988, 0.436601}, {0.281477, 0.755203, 0.432552}, {0.288921, 0.758394, 0.428426}, {0.296479, 0.761561, 0.424223}, {0.304148, 0.764704, 0.419943},
                        {0.311925, 0.767822, 0.415586}, {0.319809, 0.770914, 0.411152}, {0.327796, 0.77398, 0.40664},   {0.335885, 0.777018, 0.402049}, {0.344074, 0.780029, 0.397381},
                        {0.35236, 0.783011, 0.392636},  {0.360741, 0.785964, 0.387814}, {0.369214, 0.788888, 0.382914}, {0.377779, 0.791781, 0.377939}, {0.386433, 0.794644, 0.372886},
                        {0.395174, 0.797475, 0.367757}, {0.404001, 0.800275, 0.362552}, {0.412913, 0.803041, 0.357269}, {0.421908, 0.805774, 0.35191},  {0.430983, 0.808473, 0.346476},
                        {0.440137, 0.811138, 0.340967}, {0.449368, 0.813768, 0.335384}, {0.458674, 0.816363, 0.329727}, {0.468053, 0.818921, 0.323998}, {0.477504, 0.821444, 0.318195},
                        {0.487026, 0.823929, 0.312321}, {0.496615, 0.826376, 0.306377}, {0.506271, 0.828786, 0.300362}, {0.515992, 0.831158, 0.294279}, {0.525776, 0.833491, 0.288127},
                        {0.535621, 0.835785, 0.281908}, {0.545524, 0.838039, 0.275626}, {0.555484, 0.840254, 0.269281}, {0.565498, 0.84243, 0.262877},  {0.575563, 0.844566, 0.256415},
                        {0.585678, 0.846661, 0.249897}, {0.595839, 0.848717, 0.243329}, {0.606045, 0.850733, 0.236712}, {0.616293, 0.852709, 0.230052}, {0.626579, 0.854645, 0.223353},
                        {0.636902, 0.856542, 0.21662},  {0.647257, 0.8584, 0.209861},   {0.657642, 0.860219, 0.203082}, {0.668054, 0.861999, 0.196293}, {0.678489, 0.863742, 0.189503},
                        {0.688944, 0.865448, 0.182725}, {0.699415, 0.867117, 0.175971}, {0.709898, 0.868751, 0.169257}, {0.720391, 0.87035, 0.162603},  {0.730889, 0.871916, 0.156029},
                        {0.741388, 0.873449, 0.149561}, {0.751884, 0.874951, 0.143228}, {0.762373, 0.876424, 0.137064}, {0.772852, 0.877868, 0.131109}, {0.783315, 0.879285, 0.125405},
                        {0.79376, 0.880678, 0.120005},  {0.804182, 0.882046, 0.114965}, {0.814576, 0.883393, 0.110347}, {0.82494, 0.88472, 0.106217},   {0.83527, 0.886029, 0.102646},
                        {0.845561, 0.887322, 0.099702}, {0.85581, 0.888601, 0.097452},  {0.866013, 0.889868, 0.095953}, {0.876168, 0.891125, 0.09525},  {0.886271, 0.892374, 0.095374},
                        {0.89632, 0.893616, 0.096335},  {0.906311, 0.894855, 0.098125}, {0.916242, 0.896091, 0.100717}, {0.926106, 0.89733, 0.104071},  {0.935904, 0.89857, 0.108131},
                        {0.945636, 0.899815, 0.112838}, {0.9553, 0.901065, 0.118128},   {0.964894, 0.902323, 0.123941}, {0.974417, 0.90359, 0.130215},  {0.983868, 0.904867, 0.136897},
                        {0.993248, 0.906157, 0.143936},
                    }}});
    // PLASMA
    storage.insert({"plasma",
                    {{
                        {0.050383, 0.029803, 0.527975}, {0.063536, 0.028426, 0.533124}, {0.075353, 0.027206, 0.538007}, {0.086222, 0.026125, 0.542658}, {0.096379, 0.025165, 0.547103},
                        {0.10598, 0.024309, 0.551368},  {0.115124, 0.023556, 0.555468}, {0.123903, 0.022878, 0.559423}, {0.132381, 0.022258, 0.56325},  {0.140603, 0.021687, 0.566959},
                        {0.148607, 0.021154, 0.570562}, {0.156421, 0.020651, 0.574065}, {0.16407, 0.020171, 0.577478},  {0.171574, 0.019706, 0.580806}, {0.17895, 0.019252, 0.584054},
                        {0.186213, 0.018803, 0.587228}, {0.193374, 0.018354, 0.59033},  {0.200445, 0.017902, 0.593364}, {0.207435, 0.017442, 0.596333}, {0.21435, 0.016973, 0.599239},
                        {0.221197, 0.016497, 0.602083}, {0.227983, 0.016007, 0.604867}, {0.234715, 0.015502, 0.607592}, {0.241396, 0.014979, 0.610259}, {0.248032, 0.014439, 0.612868},
                        {0.254627, 0.013882, 0.615419}, {0.261183, 0.013308, 0.617911}, {0.267703, 0.012716, 0.620346}, {0.274191, 0.012109, 0.622722}, {0.280648, 0.011488, 0.625038},
                        {0.287076, 0.010855, 0.627295}, {0.293478, 0.010213, 0.62949},  {0.299855, 0.009561, 0.631624}, {0.30621, 0.008902, 0.633694},  {0.312543, 0.008239, 0.6357},
                        {0.318856, 0.007576, 0.63764},  {0.32515, 0.006915, 0.639512},  {0.331426, 0.006261, 0.641316}, {0.337683, 0.005618, 0.643049}, {0.343925, 0.004991, 0.64471},
                        {0.35015, 0.004382, 0.646298},  {0.356359, 0.003798, 0.64781},  {0.362553, 0.003243, 0.649245}, {0.368733, 0.002724, 0.650601}, {0.374897, 0.002245, 0.651876},
                        {0.381047, 0.001814, 0.653068}, {0.387183, 0.001434, 0.654177}, {0.393304, 0.001114, 0.655199}, {0.399411, 0.000859, 0.656133}, {0.405503, 0.000678, 0.656977},
                        {0.41158, 0.000577, 0.65773},   {0.417642, 0.000564, 0.65839},  {0.423689, 0.000646, 0.658956}, {0.429719, 0.000831, 0.659425}, {0.435734, 0.001127, 0.659797},
                        {0.441732, 0.00154, 0.660069},  {0.447714, 0.00208, 0.66024},   {0.453677, 0.002755, 0.66031},  {0.459623, 0.003574, 0.660277}, {0.46555, 0.004545, 0.660139},
                        {0.471457, 0.005678, 0.659897}, {0.477344, 0.00698, 0.659549},  {0.48321, 0.00846, 0.659095},   {0.489055, 0.010127, 0.658534}, {0.494877, 0.01199, 0.657865},
                        {0.500678, 0.014055, 0.657088}, {0.506454, 0.016333, 0.656202}, {0.512206, 0.018833, 0.655209}, {0.517933, 0.021563, 0.654109}, {0.523633, 0.024532, 0.652901},
                        {0.529306, 0.027747, 0.651586}, {0.534952, 0.031217, 0.650165}, {0.54057, 0.03495, 0.64864},    {0.546157, 0.038954, 0.64701},  {0.551715, 0.043136, 0.645277},
                        {0.557243, 0.047331, 0.643443}, {0.562738, 0.051545, 0.641509}, {0.568201, 0.055778, 0.639477}, {0.573632, 0.060028, 0.637349}, {0.579029, 0.064296, 0.635126},
                        {0.584391, 0.068579, 0.632812}, {0.589719, 0.072878, 0.630408}, {0.595011, 0.07719, 0.627917},  {0.600266, 0.081516, 0.625342}, {0.605485, 0.085854, 0.622686},
                        {0.610667, 0.090204, 0.619951}, {0.615812, 0.094564, 0.61714},  {0.620919, 0.098934, 0.614257}, {0.625987, 0.103312, 0.611305}, {0.631017, 0.107699, 0.608287},
                        {0.636008, 0.112092, 0.605205}, {0.640959, 0.116492, 0.602065}, {0.645872, 0.120898, 0.598867}, {0.650746, 0.125309, 0.595617}, {0.65558, 0.129725, 0.592317},
                        {0.660374, 0.134144, 0.588971}, {0.665129, 0.138566, 0.585582}, {0.669845, 0.142992, 0.582154}, {0.674522, 0.147419, 0.578688}, {0.67916, 0.151848, 0.575189},
                        {0.683758, 0.156278, 0.57166},  {0.688318, 0.160709, 0.568103}, {0.69284, 0.165141, 0.564522},  {0.697324, 0.169573, 0.560919}, {0.701769, 0.174005, 0.557296},
                        {0.706178, 0.178437, 0.553657}, {0.710549, 0.182868, 0.550004}, {0.714883, 0.187299, 0.546338}, {0.719181, 0.191729, 0.542663}, {0.723444, 0.196158, 0.538981},
                        {0.72767, 0.200586, 0.535293},  {0.731862, 0.205013, 0.531601}, {0.736019, 0.209439, 0.527908}, {0.740143, 0.213864, 0.524216}, {0.744232, 0.218288, 0.520524},
                        {0.748289, 0.222711, 0.516834}, {0.752312, 0.227133, 0.513149}, {0.756304, 0.231555, 0.509468}, {0.760264, 0.235976, 0.505794}, {0.764193, 0.240396, 0.502126},
                        {0.76809, 0.244817, 0.498465},  {0.771958, 0.249237, 0.494813}, {0.775796, 0.253658, 0.491171}, {0.779604, 0.258078, 0.487539}, {0.783383, 0.2625, 0.483918},
                        {0.787133, 0.266922, 0.480307}, {0.790855, 0.271345, 0.476706}, {0.794549, 0.27577, 0.473117},  {0.798216, 0.280197, 0.469538}, {0.801855, 0.284626, 0.465971},
                        {0.805467, 0.289057, 0.462415}, {0.809052, 0.293491, 0.45887},  {0.812612, 0.297928, 0.455338}, {0.816144, 0.302368, 0.451816}, {0.819651, 0.306812, 0.448306},
                        {0.823132, 0.311261, 0.444806}, {0.826588, 0.315714, 0.441316}, {0.830018, 0.320172, 0.437836}, {0.833422, 0.324635, 0.434366}, {0.836801, 0.329105, 0.430905},
                        {0.840155, 0.33358, 0.427455},  {0.843484, 0.338062, 0.424013}, {0.846788, 0.342551, 0.420579}, {0.850066, 0.347048, 0.417153}, {0.853319, 0.351553, 0.413734},
                        {0.856547, 0.356066, 0.410322}, {0.85975, 0.360588, 0.406917},  {0.862927, 0.365119, 0.403519}, {0.866078, 0.36966, 0.400126},  {0.869203, 0.374212, 0.396738},
                        {0.872303, 0.378774, 0.393355}, {0.875376, 0.383347, 0.389976}, {0.878423, 0.387932, 0.3866},   {0.881443, 0.392529, 0.383229}, {0.884436, 0.397139, 0.37986},
                        {0.887402, 0.401762, 0.376494}, {0.89034, 0.406398, 0.37313},   {0.89325, 0.411048, 0.369768},  {0.896131, 0.415712, 0.366407}, {0.898984, 0.420392, 0.363047},
                        {0.901807, 0.425087, 0.359688}, {0.904601, 0.429797, 0.356329}, {0.907365, 0.434524, 0.35297},  {0.910098, 0.439268, 0.34961},  {0.9128, 0.444029, 0.346251},
                        {0.915471, 0.448807, 0.34289},  {0.918109, 0.453603, 0.339529}, {0.920714, 0.458417, 0.336166}, {0.923287, 0.463251, 0.332801}, {0.925825, 0.468103, 0.329435},
                        {0.928329, 0.472975, 0.326067}, {0.930798, 0.477867, 0.322697}, {0.933232, 0.48278, 0.319325},  {0.93563, 0.487712, 0.315952},  {0.93799, 0.492667, 0.312575},
                        {0.940313, 0.497642, 0.309197}, {0.942598, 0.502639, 0.305816}, {0.944844, 0.507658, 0.302433}, {0.947051, 0.512699, 0.299049}, {0.949217, 0.517763, 0.295662},
                        {0.951344, 0.52285, 0.292275},  {0.953428, 0.52796, 0.288883},  {0.95547, 0.533093, 0.28549},   {0.957469, 0.53825, 0.282096},  {0.959424, 0.543431, 0.278701},
                        {0.961336, 0.548636, 0.275305}, {0.963203, 0.553865, 0.271909}, {0.965024, 0.559118, 0.268513}, {0.966798, 0.564396, 0.265118}, {0.968526, 0.5697, 0.261721},
                        {0.970205, 0.575028, 0.258325}, {0.971835, 0.580382, 0.254931}, {0.973416, 0.585761, 0.25154},  {0.974947, 0.591165, 0.248151}, {0.976428, 0.596595, 0.244767},
                        {0.977856, 0.602051, 0.241387}, {0.979233, 0.607532, 0.238013}, {0.980556, 0.613039, 0.234646}, {0.981826, 0.618572, 0.231287}, {0.983041, 0.624131, 0.227937},
                        {0.984199, 0.629718, 0.224595}, {0.985301, 0.63533, 0.221265},  {0.986345, 0.640969, 0.217948}, {0.987332, 0.646633, 0.214648}, {0.98826, 0.652325, 0.211364},
                        {0.989128, 0.658043, 0.2081},   {0.989935, 0.663787, 0.204859}, {0.990681, 0.669558, 0.201642}, {0.991365, 0.675355, 0.198453}, {0.991985, 0.681179, 0.195295},
                        {0.992541, 0.68703, 0.19217},   {0.993032, 0.692907, 0.189084}, {0.993456, 0.69881, 0.186041},  {0.993814, 0.704741, 0.183043}, {0.994103, 0.710698, 0.180097},
                        {0.994324, 0.716681, 0.177208}, {0.994474, 0.722691, 0.174381}, {0.994553, 0.728728, 0.171622}, {0.994561, 0.734791, 0.168938}, {0.994495, 0.74088, 0.166335},
                        {0.994355, 0.746995, 0.163821}, {0.994141, 0.753137, 0.161404}, {0.993851, 0.759304, 0.159092}, {0.993482, 0.765499, 0.156891}, {0.993033, 0.77172, 0.154808},
                        {0.992505, 0.777967, 0.152855}, {0.991897, 0.784239, 0.151042}, {0.991209, 0.790537, 0.149377}, {0.990439, 0.796859, 0.14787},  {0.989587, 0.803205, 0.146529},
                        {0.988648, 0.809579, 0.145357}, {0.987621, 0.815978, 0.144363}, {0.986509, 0.822401, 0.143557}, {0.985314, 0.828846, 0.142945}, {0.984031, 0.835315, 0.142528},
                        {0.982653, 0.841812, 0.142303}, {0.98119, 0.848329, 0.142279},  {0.979644, 0.854866, 0.142453}, {0.977995, 0.861432, 0.142808}, {0.976265, 0.868016, 0.143351},
                        {0.974443, 0.874622, 0.144061}, {0.97253, 0.88125, 0.144923},   {0.970533, 0.887896, 0.145919}, {0.968443, 0.894564, 0.147014}, {0.966271, 0.901249, 0.14818},
                        {0.964021, 0.90795, 0.14937},   {0.961681, 0.914672, 0.15052},  {0.959276, 0.921407, 0.151566}, {0.956808, 0.928152, 0.152409}, {0.954287, 0.934908, 0.152921},
                        {0.951726, 0.941671, 0.152925}, {0.949151, 0.948435, 0.152178}, {0.946602, 0.95519, 0.150328},  {0.944152, 0.961916, 0.146861}, {0.941896, 0.96859, 0.140956},
                        {0.940015, 0.975158, 0.131326},
                    }}});

    // TURBO
    storage.insert({"turbo",
                    {{{0.18995, 0.07176, 0.23217}, {0.19483, 0.08339, 0.26149}, {0.19956, 0.09498, 0.29024}, {0.20415, 0.10652, 0.31844}, {0.20860, 0.11802, 0.34607}, {0.21291, 0.12947, 0.37314},
                      {0.21708, 0.14087, 0.39964}, {0.22111, 0.15223, 0.42558}, {0.22500, 0.16354, 0.45096}, {0.22875, 0.17481, 0.47578}, {0.23236, 0.18603, 0.50004}, {0.23582, 0.19720, 0.52373},
                      {0.23915, 0.20833, 0.54686}, {0.24234, 0.21941, 0.56942}, {0.24539, 0.23044, 0.59142}, {0.24830, 0.24143, 0.61286}, {0.25107, 0.25237, 0.63374}, {0.25369, 0.26327, 0.65406},
                      {0.25618, 0.27412, 0.67381}, {0.25853, 0.28492, 0.69300}, {0.26074, 0.29568, 0.71162}, {0.26280, 0.30639, 0.72968}, {0.26473, 0.31706, 0.74718}, {0.26652, 0.32768, 0.76412},
                      {0.26816, 0.33825, 0.78050}, {0.26967, 0.34878, 0.79631}, {0.27103, 0.35926, 0.81156}, {0.27226, 0.36970, 0.82624}, {0.27334, 0.38008, 0.84037}, {0.27429, 0.39043, 0.85393},
                      {0.27509, 0.40072, 0.86692}, {0.27576, 0.41097, 0.87936}, {0.27628, 0.42118, 0.89123}, {0.27667, 0.43134, 0.90254}, {0.27691, 0.44145, 0.91328}, {0.27701, 0.45152, 0.92347},
                      {0.27698, 0.46153, 0.93309}, {0.27680, 0.47151, 0.94214}, {0.27648, 0.48144, 0.95064}, {0.27603, 0.49132, 0.95857}, {0.27543, 0.50115, 0.96594}, {0.27469, 0.51094, 0.97275},
                      {0.27381, 0.52069, 0.97899}, {0.27273, 0.53040, 0.98461}, {0.27106, 0.54015, 0.98930}, {0.26878, 0.54995, 0.99303}, {0.26592, 0.55979, 0.99583}, {0.26252, 0.56967, 0.99773},
                      {0.25862, 0.57958, 0.99876}, {0.25425, 0.58950, 0.99896}, {0.24946, 0.59943, 0.99835}, {0.24427, 0.60937, 0.99697}, {0.23874, 0.61931, 0.99485}, {0.23288, 0.62923, 0.99202},
                      {0.22676, 0.63913, 0.98851}, {0.22039, 0.64901, 0.98436}, {0.21382, 0.65886, 0.97959}, {0.20708, 0.66866, 0.97423}, {0.20021, 0.67842, 0.96833}, {0.19326, 0.68812, 0.96190},
                      {0.18625, 0.69775, 0.95498}, {0.17923, 0.70732, 0.94761}, {0.17223, 0.71680, 0.93981}, {0.16529, 0.72620, 0.93161}, {0.15844, 0.73551, 0.92305}, {0.15173, 0.74472, 0.91416},
                      {0.14519, 0.75381, 0.90496}, {0.13886, 0.76279, 0.89550}, {0.13278, 0.77165, 0.88580}, {0.12698, 0.78037, 0.87590}, {0.12151, 0.78896, 0.86581}, {0.11639, 0.79740, 0.85559},
                      {0.11167, 0.80569, 0.84525}, {0.10738, 0.81381, 0.83484}, {0.10357, 0.82177, 0.82437}, {0.10026, 0.82955, 0.81389}, {0.09750, 0.83714, 0.80342}, {0.09532, 0.84455, 0.79299},
                      {0.09377, 0.85175, 0.78264}, {0.09287, 0.85875, 0.77240}, {0.09267, 0.86554, 0.76230}, {0.09320, 0.87211, 0.75237}, {0.09451, 0.87844, 0.74265}, {0.09662, 0.88454, 0.73316},
                      {0.09958, 0.89040, 0.72393}, {0.10342, 0.89600, 0.71500}, {0.10815, 0.90142, 0.70599}, {0.11374, 0.90673, 0.69651}, {0.12014, 0.91193, 0.68660}, {0.12733, 0.91701, 0.67627},
                      {0.13526, 0.92197, 0.66556}, {0.14391, 0.92680, 0.65448}, {0.15323, 0.93151, 0.64308}, {0.16319, 0.93609, 0.63137}, {0.17377, 0.94053, 0.61938}, {0.18491, 0.94484, 0.60713},
                      {0.19659, 0.94901, 0.59466}, {0.20877, 0.95304, 0.58199}, {0.22142, 0.95692, 0.56914}, {0.23449, 0.96065, 0.55614}, {0.24797, 0.96423, 0.54303}, {0.26180, 0.96765, 0.52981},
                      {0.27597, 0.97092, 0.51653}, {0.29042, 0.97403, 0.50321}, {0.30513, 0.97697, 0.48987}, {0.32006, 0.97974, 0.47654}, {0.33517, 0.98234, 0.46325}, {0.35043, 0.98477, 0.45002},
                      {0.36581, 0.98702, 0.43688}, {0.38127, 0.98909, 0.42386}, {0.39678, 0.99098, 0.41098}, {0.41229, 0.99268, 0.39826}, {0.42778, 0.99419, 0.38575}, {0.44321, 0.99551, 0.37345},
                      {0.45854, 0.99663, 0.36140}, {0.47375, 0.99755, 0.34963}, {0.48879, 0.99828, 0.33816}, {0.50362, 0.99879, 0.32701}, {0.51822, 0.99910, 0.31622}, {0.53255, 0.99919, 0.30581},
                      {0.54658, 0.99907, 0.29581}, {0.56026, 0.99873, 0.28623}, {0.57357, 0.99817, 0.27712}, {0.58646, 0.99739, 0.26849}, {0.59891, 0.99638, 0.26038}, {0.61088, 0.99514, 0.25280},
                      {0.62233, 0.99366, 0.24579}, {0.63323, 0.99195, 0.23937}, {0.64362, 0.98999, 0.23356}, {0.65394, 0.98775, 0.22835}, {0.66428, 0.98524, 0.22370}, {0.67462, 0.98246, 0.21960},
                      {0.68494, 0.97941, 0.21602}, {0.69525, 0.97610, 0.21294}, {0.70553, 0.97255, 0.21032}, {0.71577, 0.96875, 0.20815}, {0.72596, 0.96470, 0.20640}, {0.73610, 0.96043, 0.20504},
                      {0.74617, 0.95593, 0.20406}, {0.75617, 0.95121, 0.20343}, {0.76608, 0.94627, 0.20311}, {0.77591, 0.94113, 0.20310}, {0.78563, 0.93579, 0.20336}, {0.79524, 0.93025, 0.20386},
                      {0.80473, 0.92452, 0.20459}, {0.81410, 0.91861, 0.20552}, {0.82333, 0.91253, 0.20663}, {0.83241, 0.90627, 0.20788}, {0.84133, 0.89986, 0.20926}, {0.85010, 0.89328, 0.21074},
                      {0.85868, 0.88655, 0.21230}, {0.86709, 0.87968, 0.21391}, {0.87530, 0.87267, 0.21555}, {0.88331, 0.86553, 0.21719}, {0.89112, 0.85826, 0.21880}, {0.89870, 0.85087, 0.22038},
                      {0.90605, 0.84337, 0.22188}, {0.91317, 0.83576, 0.22328}, {0.92004, 0.82806, 0.22456}, {0.92666, 0.82025, 0.22570}, {0.93301, 0.81236, 0.22667}, {0.93909, 0.80439, 0.22744},
                      {0.94489, 0.79634, 0.22800}, {0.95039, 0.78823, 0.22831}, {0.95560, 0.78005, 0.22836}, {0.96049, 0.77181, 0.22811}, {0.96507, 0.76352, 0.22754}, {0.96931, 0.75519, 0.22663},
                      {0.97323, 0.74682, 0.22536}, {0.97679, 0.73842, 0.22369}, {0.98000, 0.73000, 0.22161}, {0.98289, 0.72140, 0.21918}, {0.98549, 0.71250, 0.21650}, {0.98781, 0.70330, 0.21358},
                      {0.98986, 0.69382, 0.21043}, {0.99163, 0.68408, 0.20706}, {0.99314, 0.67408, 0.20348}, {0.99438, 0.66386, 0.19971}, {0.99535, 0.65341, 0.19577}, {0.99607, 0.64277, 0.19165},
                      {0.99654, 0.63193, 0.18738}, {0.99675, 0.62093, 0.18297}, {0.99672, 0.60977, 0.17842}, {0.99644, 0.59846, 0.17376}, {0.99593, 0.58703, 0.16899}, {0.99517, 0.57549, 0.16412},
                      {0.99419, 0.56386, 0.15918}, {0.99297, 0.55214, 0.15417}, {0.99153, 0.54036, 0.14910}, {0.98987, 0.52854, 0.14398}, {0.98799, 0.51667, 0.13883}, {0.98590, 0.50479, 0.13367},
                      {0.98360, 0.49291, 0.12849}, {0.98108, 0.48104, 0.12332}, {0.97837, 0.46920, 0.11817}, {0.97545, 0.45740, 0.11305}, {0.97234, 0.44565, 0.10797}, {0.96904, 0.43399, 0.10294},
                      {0.96555, 0.42241, 0.09798}, {0.96187, 0.41093, 0.09310}, {0.95801, 0.39958, 0.08831}, {0.95398, 0.38836, 0.08362}, {0.94977, 0.37729, 0.07905}, {0.94538, 0.36638, 0.07461},
                      {0.94084, 0.35566, 0.07031}, {0.93612, 0.34513, 0.06616}, {0.93125, 0.33482, 0.06218}, {0.92623, 0.32473, 0.05837}, {0.92105, 0.31489, 0.05475}, {0.91572, 0.30530, 0.05134},
                      {0.91024, 0.29599, 0.04814}, {0.90463, 0.28696, 0.04516}, {0.89888, 0.27824, 0.04243}, {0.89298, 0.26981, 0.03993}, {0.88691, 0.26152, 0.03753}, {0.88066, 0.25334, 0.03521},
                      {0.87422, 0.24526, 0.03297}, {0.86760, 0.23730, 0.03082}, {0.86079, 0.22945, 0.02875}, {0.85380, 0.22170, 0.02677}, {0.84662, 0.21407, 0.02487}, {0.83926, 0.20654, 0.02305},
                      {0.83172, 0.19912, 0.02131}, {0.82399, 0.19182, 0.01966}, {0.81608, 0.18462, 0.01809}, {0.80799, 0.17753, 0.01660}, {0.79971, 0.17055, 0.01520}, {0.79125, 0.16368, 0.01387},
                      {0.78260, 0.15693, 0.01264}, {0.77377, 0.15028, 0.01148}, {0.76476, 0.14374, 0.01041}, {0.75556, 0.13731, 0.00942}, {0.74617, 0.13098, 0.00851}, {0.73661, 0.12477, 0.00769},
                      {0.72686, 0.11867, 0.00695}, {0.71692, 0.11268, 0.00629}, {0.70680, 0.10680, 0.00571}, {0.69650, 0.10102, 0.00522}, {0.68602, 0.09536, 0.00481}, {0.67535, 0.08980, 0.00449},
                      {0.66449, 0.08436, 0.00424}, {0.65345, 0.07902, 0.00408}, {0.64223, 0.07380, 0.00401}, {0.63082, 0.06868, 0.00401}, {0.61923, 0.06367, 0.00410}, {0.60746, 0.05878, 0.00427},
                      {0.59550, 0.05399, 0.00453}, {0.58336, 0.04931, 0.00486}, {0.57103, 0.04474, 0.00529}, {0.55852, 0.04028, 0.00579}, {0.54583, 0.03593, 0.00638}, {0.53295, 0.03169, 0.00705},
                      {0.51989, 0.02756, 0.00780}, {0.50664, 0.02354, 0.00863}, {0.49321, 0.01963, 0.00955}, {0.47960, 0.01583, 0.01055}}}});

    return storage;
  }

  std::map<std::string, std::vector<std::array<float, 3>>> color::cyclics = color::fillCyclics();

  std::map<std::string, std::vector<std::array<float, 3>>> color::fillCyclics() {
    std::map<std::string, std::vector<std::array<float, 3>>> storage;
    // ESTRUCTURE
    storage.insert({
        "estructure",
        {{0.000000, 0.000000, 0.000000},
          {0.29800, 0.447000, 0.690000},
          {0.333000, 0.658000, 0.407000},
          {0.768000, 0.305000, 0.321000},
          {0.80000, 0.725000, 0.454000},
          {0.505000, 0.447000, 0.698000},
          {0.392000, 0.709000, 0.803000},
          {0.87800, 0.501000, 0.000000},
          {0.000000, 0.376000, 0.376000},
          {0.627000, 0.376000, 0.000000},
          {0.74500, 0.501000, 1.000000}}
    });
    return storage;
  }

  /// @brief get rgb triplet in xcolor format
  /// @param name cyclic list name
  /// @param index
  /// @return xcolor {rgb}{`r`,`g`,`b`}
  std::string color::getCyclicTikz(const std::string name, const size_t index) {
    return "{rgb}{" + std::to_string(cyclics.at(name)[index][0]) + "," + std::to_string(cyclics.at(name)[index][1]) + "," + std::to_string(cyclics.at(name)[index][2]) + "}";
  }
} // namespace aurostd

//NHA 7302025
//xtable implementations

namespace aurostd {
  //constructor
  xtable::xtable(vector<vector<string>> data, vector<string> column_titles, vector<string> column_units, int font_size, string title) :
      data(data), column_titles(column_titles), column_units(column_units), font_size(font_size), title(title) {
    //check to make sure titles and units share size:
    if (column_titles.size() != column_units.size()) {
      cerr << "ERROR IN xtable CONSTRUCTOR: columns titles and column units do not share size!" << endl;
    }
  }
  //constructor
  xtable::xtable(std::map<string, vector<std::map<string, vector<string>>>> compound_data, vector<string> column_titles, vector<string> column_units, int font_size, string title) :
      compound_data(compound_data), column_titles(column_titles), column_units(column_units), font_size(font_size), title(title) {
    //check to make sure titles and units share size:
    if (column_titles.size() != column_units.size()) {
      cerr << "ERROR IN xtable CONSTRUCTOR: columns titles and column units do not share size!" << endl;
    }
  }

  string xtable::createTable() {
    std::stringstream ss_out;
    bool grey_row = true;
    vector<vector<string>> data = this->data;
    vector<string> column_titles = this->column_titles;
    vector<string> column_units = this->column_units;
    int font_size = this->font_size;
    string title = this->title;

    ss_out << "\\begin{longtable}[c]{|";
    for (int i = 0; i < column_titles.size(); i++) {
      ss_out << "c|";
    }
    ss_out << "}";
    ss_out << "\\hline";
    ss_out << "\\multicolumn{" << column_titles.size() << "}{| c |}{" << title << "}\\\\";
    ss_out << "\\hline" << endl;

    for (int i = 0; i < column_titles.size(); i++) {
      ss_out << "{" + column_titles[i];
      if (!column_units[i].empty()) {
        ss_out << " (" + column_units[i] + ")";
      }
      ss_out << "} & ";
    }

    ss_out.seekp(-2, ss_out.cur);
    ss_out << "\\\\" << endl;
    ss_out << "\\hline" << endl;
    ss_out << "\\endfirsthead" << endl;
    ss_out << "\\hline" << endl << endl;

    ss_out << "\\multicolumn{" << column_titles.size() << "}{| c |}{" << title << "}\\\\";
    ss_out << "\\hline" << endl;

    for (int i = 0; i < column_titles.size(); i++) {
      ss_out << "{" + column_titles[i];
      if (!column_units[i].empty()) {
        ss_out << " (" + column_units[i] + ")";
      }
      ss_out << "} & ";
    }
    ss_out.seekp(-2, ss_out.cur);
    ss_out << "\\\\" << endl;
    ss_out << "\\hline" << endl;
    ss_out << "\\endhead" << endl;
    ss_out << "\\hline" << endl;
    ss_out << "\\endfoot" << endl;

    for (auto row : data) {
      ss_out << "\\hline" << endl;
      ss_out << "  {";
      for (auto column : row) {
        ss_out << column << "} & {";
      }
      ss_out.seekp(-3, ss_out.cur);
      ss_out << "\\\\" << endl;
    }

    ss_out << "\\end{longtable}" << endl;
    ss_out << "\\vspace{-20pt}" << endl;

    return ss_out.str();
  }

  //@brief table creation routine for hulls
  string xtable::createCompoundTable() {
    std::stringstream ss_out;
    std::map<string, vector<std::map<string, vector<string>>>> compound_data = this->compound_data; //format: std::map<subtitle, data>
    vector<string> column_titles = this->column_titles;
    vector<string> column_units = this->column_units;
    int font_size = this->font_size;
    string title = this->title;

    //add key:
    ss_out << "\\begin{center}" << endl;
    ss_out << "\\begin{tabular}{ |c| }" << endl;
    ss_out << "\\hline" << endl;
    ss_out << R"(Key: \color{blue} CCE used \\)" << endl;
    ss_out << "\\hline" << endl;
    ss_out << R"($H_f \text{,} \ \Delta H_{hull}$ in (meV/atom) \\)" << endl;
    ss_out << "\\hline" << endl;
    ss_out << "\\end{tabular}" << endl;
    ss_out << "\\end{center}" << endl;
    ss_out << endl;

    //add table:
    ss_out << "\\begin{longtable}[c]{|";
    for (int i = 0; i < column_titles.size(); i++) {
      ss_out << "c|";
    }
    ss_out << "}";
    ss_out << "\\hline";
    ss_out << "\\multicolumn{" << column_titles.size() << "}{| c |}{" << title << "}\\\\";
    ss_out << "\\hline" << endl;

    for (int i = 0; i < column_titles.size(); i++) {
      ss_out << "{" + column_titles[i];
      if (!column_units[i].empty()) {
        ss_out << " (" + column_units[i] + ")";
      }
      ss_out << "} & ";
    }

    ss_out.seekp(-2, ss_out.cur);
    ss_out << "\\\\" << endl;

    for (auto compound : compound_data) {
      ss_out << "\\hline" << endl;
      ss_out << "\\multicolumn{" << column_titles.size() << "}{| c |}{" << compound.first << "}\\\\" << endl;
      for (auto row : compound.second) {
        ss_out << "\\hline" << endl;
        for (auto column : row) {
          if (!column.first.empty()) {
            ss_out << "\\rowcolor{" << column.first << "!25}" << endl;
          }
          ss_out << "  {";
          for (auto c : column.second) {
            ss_out << c << "} & {";
          }
        }
        ss_out.seekp(-3, ss_out.cur);
        ss_out << "\\\\" << endl;
      }
    }
    ss_out << "\\hline" << endl;
    ss_out << "\\end{longtable}" << endl;
    ss_out << "\\vspace{-20pt}" << endl;

    return ss_out.str();
  }

  string xtable::toString() {
    std::stringstream ss_out;
    ss_out << "\\documentclass[12pt]{article}" << endl;
    ss_out << "\\usepackage[a4paper,margin=1.5cm]{geometry}" << endl;
    ss_out << "\\usepackage{longtable}" << endl;
    ss_out << "\\usepackage[table]{xcolor}" << endl;
    ss_out << "\\usepackage{multirow}" << endl;
    ss_out << "\\usepackage{tabularx}" << endl;
    ss_out << "\\usepackage{xcolor}" << endl;
    ss_out << "\\usepackage{amsmath}" << endl;
    ss_out << endl;
    ss_out << "\\begin{document}" << endl;

    if (!this->data.empty() && this->compound_data.empty()) {
      ss_out << this->createTable();
    } else if (!this->compound_data.empty() && this->data.empty()) {
      ss_out << this->createCompoundTable();
    } else {
      cerr << "WARNING IN FUNCTION xtable::toString(): no initialized data!" << endl;
    }

    ss_out << "\\end{document}";

    return ss_out.str();
  }

  void xtable::saveWithPlot(xtable table, xplotter xplt, const std::string& file_path) {
    string plot_store = xplt.toString();
    string table_string = table.createTable();
    string table_locator = "LOCATION FOR TABLE:";

    int pos = plot_store.find(table_locator);
    plot_store.insert(pos + table_locator.size(), "\n" + table_string);

    aurostd::string2file(plot_store, file_path);
  }

  void xtable::save(const std::string& file_path) {
    // Parse the file path and set defaults
    std::array<std::string, 3> pns = aurostd::splitFilePath(file_path); // parent, name, suffix
    if (pns[2].empty() || pns[2] == ".") {
      pns[2] = ".pdf";
    }
    if (pns[1].empty()) {
      pns[1] = "aflow_plot";
    }
    std::string render_type = pns[2];
    if (render_type == ".tiff") {
      render_type = ".tif";
    } else if (render_type == ".jpeg") {
      render_type = ".jpg";
    }
    // if just .tex is wanted
    if (render_type == ".tex") {
      aurostd::string2file(this->toString(), file_path);
      return;
    } else if (render_type == ".pdf") {
      // ensure that pdflatex is available
      if (!aurostd::IsCommandAvailable("pdflatex")) {
        cerr << "xtable: pdflatex not in $PATH - falling back to .tex" << endl;
        aurostd::string2file(this->toString(), pns[0] + pns[1] + ".tex");
        return;
      }
      // render PDF in tmp folder
      std::string final_file_path = pns[0] + pns[1] + pns[2];
      std::string active_work_folder = TmpDirectoryCreate("xplotter", "", false);
      aurostd::string2file(this->toString(), active_work_folder + "/render.tex");
      std::string command = "pdflatex -interaction=batchmode -output-dir " + active_work_folder + " " + active_work_folder + "/render.tex";
      aurostd::execute2string(command); // 2string for silence

      aurostd::file2file(active_work_folder + "/render.pdf", final_file_path);
      //clean the temp file
      aurostd::RemoveDirectory(active_work_folder);
    } else {
      cerr << "ERROR IN xtable::save(): exporting format currently not supported!" << endl;
    }
  }
} // namespace aurostd
