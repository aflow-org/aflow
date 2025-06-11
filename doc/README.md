## AFLOW Documentation

### How to build

- install [doxygen](https://doxygen.nl/download.html)
  - mac: `brew install doxygen`
  - debian/ubuntu: `sudo apt install doxygen`
- change to this `doc` folder
- run `doxygen`
  - `doxygen aflow.doxy`
- open the `index.html`
  - mac: `open ./out/html/index.html`
- create the pdf
  - `make -p ./out/latex`
  
### How to document
```c++
  /// @brief volume of a solid defined by points, facets their normals
  /// @param points collection of points
  /// @param facets collection of ordered point indices describing a facet
  /// @param normals collection of facet normals pointing all either outwards or inwards of the solid
  /// @return volume
  ///
  /// @authors
  /// @mod{HE,20210721,created}
  /// @mod{HE,20220707,changed documentation}
  /// @mod{SD,20220712,added comments}
  ///
  /// A series of pyramids are generated from the solid with the facets as bases and origin as their tips.
  /// Their volumes (1/3 * base area * height) is then summed up.
  /// The height of the pyramids is the scalar product of the normal vector and a point on the facet.
  /// Depending upon the normal direction, the height and, therefore, the volume can be negative.
  /// This ensures that overlapping volumes are handled properly.
  ///
  /// \f[
  /// \frac{1}{3} \left| \sum_F (P0_F \cdot N_F) A_F \right|
  /// \f]
  ///
  /// - \f$ P0_F \f$ first point of a facet (could be any point on facet F)
  /// - \f$ N_F \f$ facet normal vector
  /// - \f$ A_F \f$ facet area
  ///
  /// @warning something dangerous happens when this function is used in the wrong way
  ///
  /// @note some light reminder
  ///
  /// @see
  /// @xlink{Wikipedia Article, https://en.wikipedia.org/wiki/Polyhedron}
  /// @xlink{Wikipedia Article, https://en.wikipedia.org/wiki/Polyhedron#Volume, how to calculate the volume}
  /// @xlink{http://google.com}
  /// @xlink{aurostd::xvector}
  /// @doi{10.1109/ICIP.2001.958278}
  /// @doi{10.1038/s41578-019-0170-8, see the rock-salt structure in Fig. 2}
  ///
  /// @todo a missing piece
  template<class utype>
  double volume(const vector <xvector<utype> > &points, 
                const vector <vector<uint> > &facets,
                const vector <xvector<utype> > &normals) {
    ...
  }
```

- `@brief`
  - short description of the intent usage
- `@param [name]`
  - description of the function parameters
- `@return`
  - description of the data returned by the function (if not `void`)
- `@authors`
  - authors and history of the function
  - for each major change add an `@mod{author,date,change}`
    - `author` is the two letter shorthand used in comment signatures (`//HE20220707`)
    - `date` follows the same format (`YYYYMMDD`)
    - `change` short description of the main change
    - avoid extra spaces between `author` and `date`
- If needed, add a detailed explanation next
  - math symbols in text can be used encapsuled in `\f$`
  - full multiline functions can use `\f[` and `\f]`
- For some function it make sense to add explicit warnings or notes for their usage to avoid misunderstanding
- The last normal section is `@see`
  - `@xlink` can be used for weblinks with a title or as direct link
    - `@xlink{link text, https://link.example.net}` or directly `@xlink{https://link.example.net}`
  - `@xlink` can also be used to point to a different part of the documentation
    - `@xlink{aurostd::xvector}`
  - `@doi` converts the doi to a clickable link 
- During development a `@todo` can help to keep a list of needed changes, but should be removed before merged to staging
