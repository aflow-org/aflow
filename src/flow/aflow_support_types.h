#ifndef AFLOW_TMP_SUPPORT_H
#define AFLOW_TMP_SUPPORT_H

// ME20190628 BEGIN - moved from CHULL for broader access
//  Output formats
enum filetype {   // CO20190629
    // GENERAL FILE TYPES
  txt_ft,         // general plain text
  json_ft,
  aflow_ft,       // ME20210329 - aflowlib.out format
  csv_ft,
  latex_ft,
  gnuplot_ft,
  jupyter2_ft,    // python 2 jupyter
  jupyter3_ft,    // python 3 jupyter
    // CHULL SPECIFIC
  chull_apool_ft,
  chull_web_ft,
};

// Vector reduction types
enum vector_reduction_type {   // CO20190629
  frac_vrt,   // reduce to fractions (normalized to 1)
  gcd_vrt,    // reduce by gcd
  no_vrt,     // no reduction
};
// ME20190628 END

#endif // AFLOW_TMP_SUPPORT_H
