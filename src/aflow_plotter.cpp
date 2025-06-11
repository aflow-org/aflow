// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *            Aflow MARCO ESTERS - Duke University 2019-2021               *
// *                                                                         *
// ***************************************************************************
// 
// A namespace for functions to plot electronic structures, phonon properties,
// and physical properties from AFLOW calculations. The default plotting
// engine is gnuplot with the epslatex terminal.
//
// Each plot type such as a DOS plot has a main function that takes the plot
// options and a stringstringstream to output the plots into the desired
// format (e.g. a gnuplot script). Each format should have its own wrapper
// function (only overload when using gnuplot as the output format).
// 
// Plot options should be handled with xoptions only to allow for as much
// flexibility and customizability as possible.

#include "aflow.h"
#include "aflow_pocc.h"

using std::deque;
using std::string;
using std::stringstream;
using std::vector;
using aurostd::xoption;

#define _DEBUG_PLOTTER_ false  //CO20190116

static const string BANDDOS_SIZE = "8, 4.5";

static const string DEFAULT_IMAGE_FORMAT = "pdf";

//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//                              PLOT FUNCTIONS                              //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////

namespace plotter {

  // Plot options ------------------------------------------------------------

  //getPlotOptions////////////////////////////////////////////////////////////
  // Sets plot options for all plots based on command line arguments.
  xoption getPlotOptions(const aurostd::xoption& xopt, const string& key, bool datasets) {
    xoption plotoptions;

    // Get options
    vector<string> tokens;
    string scheme = xopt.getattachedscheme(key);
    uint ntokens = aurostd::string2tokens(scheme, tokens, ",", true); // Keep empty arguments

    if (ntokens >= 1) {
      if (tokens[0][tokens[0].size() - 1] != '/') tokens[0] += '/';
      plotoptions.push_attached("DIRECTORY", tokens[0]);
    }

    plotoptions.flag("DATATYPE",xopt.flag("ATOMS")); //CO20191010 - which table to look at, "atoms-projected" by default
    plotoptions.flag("PLOT_ALL_ATOMS",xopt.flag("PLOT_ALL_ATOMS")); //CO20191010

    // Need to shift options if partial DOS are plotted
    uint shift = 0;
    if (datasets) { // DATASET = -1: all data sets
      if ((ntokens >= 2) && !tokens[1].empty()) {
        plotoptions.push_attached("DATASET", tokens[1]);
      } else {
        plotoptions.push_attached("DATASET", "-1");
      }
      shift = 1;
    } else {  // DATASET = 0: no data set specified 
      plotoptions.push_attached("DATASET", "0");
    }
    if (ntokens >= 2 + shift) plotoptions.push_attached("XMIN", tokens[1 + shift]);
    if (ntokens >= 3 + shift) plotoptions.push_attached("XMAX", tokens[2 + shift]);
    if (ntokens >= 4 + shift) plotoptions.push_attached("YSCALE", tokens[3 + shift]);

    // Get title if present
    scheme = xopt.getattachedscheme("PLOTTER::TITLE");
    if (!scheme.empty()) plotoptions.push_attached("TITLE", scheme);

    // Get image format
    scheme = xopt.getattachedscheme("PLOTTER::PRINT");
    if (!scheme.empty()) plotoptions.push_attached("IMAGE_FORMAT", aurostd::tolower(scheme));

    //ME20200313 - Get user-defined output file name
    string outfile = xopt.getattachedscheme("PLOTTER::OUTFILE");
    // Remove extension from user-defined file name (extensions will be handled later)
    if (!outfile.empty()) {
      if (scheme.empty()) scheme = "." + DEFAULT_IMAGE_FORMAT;
      else scheme = "." + aurostd::tolower(scheme);
      uint nchar_scheme = scheme.size();
      uint nchar_outfile = outfile.size();
      if (nchar_outfile > nchar_scheme) {
        uint i = 0;
        for (i = 0; i < nchar_scheme; i++) {
          if (scheme[i] != aurostd::tolower(outfile)[nchar_outfile - nchar_scheme + i]) break;
        }
        if (i == nchar_scheme) outfile = outfile.substr(0, nchar_outfile - nchar_scheme);
      }
      plotoptions.push_attached("FILE_NAME_USER", outfile);
    }

    // Set standard background color
    plotoptions.push_attached("BACKGROUND_COLOR", "#FFFFFF");

    // Set standard grid options
    plotoptions.push_attached("GRID_COLOR", "#808080");
    plotoptions.push_attached("GRID_WIDTH", "1");
    plotoptions.push_attached("GRID_LINE_TYPE", "0");

    // Set standard size
    plotoptions.push_attached("PLOT_SIZE", "5.333, 3");

    plotoptions.flag("NOWATERMARK", xopt.flag("PLOTTER::NOWATERMARK"));
    return plotoptions;
  }

  // Electronic structure plots ----------------------------------------------

  //getPlotOptionsEStructure//////////////////////////////////////////////////
  // Sets the plot options that are specific to electronic structure plots.
  xoption getPlotOptionsEStructure(const aurostd::xoption& xopt, const string& key, bool datasets) {
    bool LDEBUG=(FALSE || _DEBUG_PLOTTER_ || XHOST.DEBUG); 

    if(LDEBUG){cerr << __AFLOW_FUNC__ << " BEGIN" << endl;}

    xoption plotoptions = getPlotOptions(xopt, key, datasets);

    // Projection
    string scheme = xopt.getattachedscheme("PLOTTER::PROJECTION");
    if (scheme.empty()) {
      plotoptions.push_attached("PROJECTION", "ORBITALS");
    } else {
      plotoptions.push_attached("PROJECTION", aurostd::toupper(scheme));
    }
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " projection=" << plotoptions.getattachedscheme("PROJECTION") << endl;}

    // No border
    plotoptions.flag("NOBORDER", true);

    // Set to true to not set Efermi to zero
    plotoptions.flag("NOSHIFT", xopt.flag("PLOTTER::NOSHIFT"));

    // Set gray background
    plotoptions.pop_attached("BACKGROUND_COLOR");
    plotoptions.push_attached("BACKGROUND_COLOR", "#E4E4E4");

    // Change grid options
    plotoptions.pop_attached("GRID_COLOR");
    plotoptions.push_attached("GRID_COLOR", "#FFFFFF");
    plotoptions.pop_attached("GRID_WIDTH");
    plotoptions.push_attached("GRID_WIDTH", "2");
    plotoptions.pop_attached("GRID_LINE_TYPE");
    plotoptions.push_attached("GRID_LINE_TYPE", "1");

    return plotoptions;
  }

  //getPlotOptionsPhonons/////////////////////////////////////////////////////
  // Sets the plot options that are specific to phonon dispersions and DOS.
  xoption getPlotOptionsPhonons(const aurostd::xoption& xopt, const string& key) {
    xoption plotoptions = getPlotOptionsEStructure(xopt, key);
    string scheme = xopt.getattachedscheme("PLOTTER::UNIT");
    if (scheme.empty()) {
      plotoptions.push_attached("UNIT", "THZ");
    } else {
      plotoptions.push_attached("UNIT", aurostd::toupper(scheme));
    }
    // There is no Fermi level for phonons, so do not shift
    plotoptions.flag("NOSHIFT", true);

    // Orbital-projections do not exist for phonons
    scheme = plotoptions.getattachedscheme("PROJECTION");
    if ((scheme == "ORBITALS") || (scheme == "LM")) {
      plotoptions.pop_attached("PROJECTION");
      plotoptions.push_attached("PROJECTION", "NONE");
    }
    return plotoptions;
  }

  //AS2020210705 BEGIN
  //getPlotOptionsQHAthermo/////////////////////////////////////////////////////
  // Sets the plot options that are specific to thermal properties obtained by the QHA
  xoption getPlotOptionsQHAthermo(const aurostd::xoption& xopt, const string& key) {
    xoption plotoptions = getPlotOptions(xopt, key);
    string scheme = xopt.getattachedscheme("PLOTTER::EOSMODEL");
    if (scheme.empty()) {
      plotoptions.push_attached("EOSMODEL", "SJ");
    } else {
      plotoptions.push_attached("EOSMODEL", aurostd::toupper(scheme));
    }
    return plotoptions;
  }
  //AS2020210705 END

  // Plot functions ----------------------------------------------------------

  //generateHeader////////////////////////////////////////////////////////////
  // Creates the header in the desired output format.
  void generateHeader(stringstream& out, const aurostd::xoption& plotoptions, bool multiplot) {
    string plottitle = plotoptions.getattachedscheme("PLOT_TITLE");
    string outformat = plotoptions.getattachedscheme("OUTPUT_FORMAT");
    if (outformat == "GNUPLOT") {
      out << "# Generated by AFLOW" << std::endl;
      out << "set terminal epslatex standalone color"
        << " size " << plotoptions.getattachedscheme("PLOT_SIZE") << " linewidth 2" << std::endl;
      out << "set output " << "'" << plotoptions.getattachedscheme("FILE_NAME_LATEX") << ".tex'" << std::endl;
      if (!plottitle.empty())
        out << "set " << (multiplot?"multiplot layout 1,2 ":"")  //ME20200313 - some machines require layout
          << "title '" << plottitle << "' offset 0, -0.5" << std::endl;
      if (plotoptions.flag("NOBORDER")) out << "unset border" << std::endl;
      out << "set object 1 rectangle from graph 0,0 to graph 1,1 fc"
        << " rgb '" << plotoptions.getattachedscheme("BACKGROUND_COLOR") << "' behind fs noborder" << std::endl;
      out << "set grid back lt " << plotoptions.getattachedscheme("GRID_LINE_TYPE")
        << " lc rgb '" << plotoptions.getattachedscheme("GRID_COLOR") << "'"
        << " lw " << plotoptions.getattachedscheme("GRID_WIDTH") << std::endl;
      if (!plotoptions.flag("NOWATERMARK")) {
        out << "set label right '"
          << (plotoptions.flag("BANDDOS")?"":"\\scriptsize") << " "
          << AFLOWLIB_CONSORTIUM_STRING << "' at screen 0.98, 0.025" << std::endl;
      }
      out << std::endl;
    }
  }

  //savePlotGNUPLOT///////////////////////////////////////////////////////////
  // Executes the gnuplot script and converts into the desired image format.
  void savePlotGNUPLOT(const xoption& plotoptions, const stringstream& gpfile) {
    bool LDEBUG=(FALSE || _DEBUG_PLOTTER_ || XHOST.DEBUG); 
    //ME20200327 - Check that all required binaries are available
    // Check that gnuplot is version 5+
    if (XHOST.is_command("gnuplot")) {
      string versionstring = aurostd::execute2string(XHOST.command("gnuplot") + " --version");
      vector<string> tokens;
      aurostd::string2tokens(versionstring, tokens, " ");
      double version = aurostd::string2utype<double>(tokens[1]);
      if (version < 5.0) {
        string message = "Gnuplot needs to be version 5 or newer (found " + tokens[1] + ").";
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_ERROR_);
      }
    }
    // ME20200609
    // Check the pdflatex version - old versions need to use different
    // compilation routes. Since the version string formats can be
    // inconsistent, the copyright year will be used as a proxy.
    uint pdflatex_version = 0;
    if (XHOST.is_command("pdflatex")) {
      string versionstring = aurostd::execute2string(XHOST.command("pdflatex") + " --version");
      vector<string> vstring;
      aurostd::string2vectorstring(versionstring, vstring);
      versionstring = vstring[2];
      aurostd::string2tokens(versionstring, vstring);
      pdflatex_version = aurostd::string2utype<uint>(vstring[1]);
    }

    string binaries = "gnuplot,convert";
    // ME20200609 - old pdfatex versions cannot process eps files
    if (pdflatex_version >= 2010) {
      binaries += ",pdflatex,repstopdf";
    } else {
      binaries += ",latex,dvips,ps2pdf";
    }
    vector<string> missing_binaries, required_binaries;
    aurostd::string2tokens(binaries, required_binaries, ",");
    for (uint i = 0; i < required_binaries.size(); i++) {
      if (!XHOST.is_command(required_binaries[i])) missing_binaries.push_back(required_binaries[i]);
    }

    if (missing_binaries.size() == 0) {
      string directory_work = plotoptions.getattachedscheme("DIRECTORY");
      if(directory_work.empty()){directory_work=aurostd::getPWD();}  //[CO20191112 - OBSOLETE]aurostd::execute2string("pwd")//CO20191004
      if(LDEBUG) {cerr << __AFLOW_FUNC__ << " directory_work=" << directory_work << endl;}
      string filename = plotoptions.getattachedscheme("FILE_NAME");
      if(LDEBUG){cerr << __AFLOW_FUNC__ << " filename=" << filename << endl;}
      string filename_latex = plotoptions.getattachedscheme("FILE_NAME_LATEX");
      // PDF is default since we use pdflatex to compile
      string format = plotoptions.getattachedscheme("IMAGE_FORMAT");
      if (format.empty()) format = "pdf";
      string current_dir = aurostd::getPWD();  //[CO20191112 - OBSOLETE]aurostd::execute2string("pwd")
      // Create temp directory
      string directory_tmp = aurostd::TmpDirectoryCreate("plotLATEX") + "/";
      chdir(directory_tmp.c_str());
      // Execute gnuplot and pdflatex
      string command="",output="";
      aurostd::stringstream2file(gpfile, filename + ".plt");
      command=XHOST.command("gnuplot") + " \"" + filename + ".plt\"";
      if(LDEBUG){cerr << __AFLOW_FUNC__ << " executing command: \"" << command << "\"" << endl;}
      aurostd::execute(command);
      if(LDEBUG) cerr << __AFLOW_FUNC__ << " directory_tmp = " << directory_tmp << endl;
      if(LDEBUG) cerr << __AFLOW_FUNC__ << aurostd::execute("ls -las "+directory_tmp) << endl;
      // ME20200609 - old pdfatex versions cannot process eps files
      if (pdflatex_version >= 2010) {
        command=XHOST.command("pdflatex") + " -interaction=nonstopmode -halt-on-error \"" + filename_latex + ".tex\" 2>&1 > /dev/null";
        if(LDEBUG){cerr << __AFLOW_FUNC__ << " executing command: \"" << command << "\"" << endl;}
        aurostd::execute(command);
      } else {
        command=XHOST.command("latex") + " -interaction=nonstopmode -halt-on-error \"" + filename_latex + ".tex\" 2>&1 > /dev/null";
        if(LDEBUG){cerr << __AFLOW_FUNC__ << " executing command: \"" << command << "\"" << endl;}
        aurostd::execute(command);
        command=XHOST.command("dvips") + " " + filename_latex + ".dvi  > /dev/null 2>&1";
        if(LDEBUG){cerr << __AFLOW_FUNC__ << " executing command: \"" << command << "\"" << endl;}
        aurostd::execute(command);
        command=XHOST.command("ps2pdf") + " " + filename_latex + ".ps";
        if(LDEBUG){cerr << __AFLOW_FUNC__ << " executing command: \"" << command << "\"" << endl;}
        aurostd::execute(command);
      }
      // Convert to the desired format if not pdf
      if (format != "pdf") {
        command=XHOST.command("convert") + " -quiet -density 300 -background white \"" + filename_latex + ".pdf\" convert_output." + format + " 1>/dev/null 2>&1";   // to avoid C: ... Carbon:PBE = C: in window //CO20210701 - io redirection
        if(LDEBUG){cerr << __AFLOW_FUNC__ << " executing command: \"" << command << "\"" << endl;}
        output=aurostd::execute2string(command);
        if(0){  //CO20210701 - does not work because of primitive execute2string functionality, will be fixed with upcoming version
          //start here for a patch:
          //https://www.itechlounge.net/2020/09/web-imagickexception-attempt-to-perform-an-operation-not-allowed-by-the-security-policy-pdf/
          if(output.find("not allowed by the security policy")!=string::npos){
            string message="The ImageMagick policy file disables ghostscript formats. Please see here and update the policy file: https://bugs.archlinux.org/task/60580";
            throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__, message, _RUNTIME_ERROR_);
          }
        }
        if(!aurostd::FileExist("convert_output." + format)){
          string message="The ImageMagick policy file disables ghostscript formats. Please see here and update the policy file: https://bugs.archlinux.org/task/60580";
          throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__, message, _RUNTIME_ERROR_);
        }
        if(LDEBUG) cerr << __AFLOW_FUNC__ << aurostd::execute("ls -las "+directory_tmp) << endl;
        command="mv convert_output." + format + " \"" + filename_latex  + "." + format + "\"";
        if(LDEBUG){cerr << __AFLOW_FUNC__ << " executing command: \"" << command << "\"" << endl;}
        aurostd::execute(command);
        if(LDEBUG) cerr << __AFLOW_FUNC__ << aurostd::execute("ls -las "+directory_tmp) << endl;
      }
      chdir(current_dir.c_str());
      aurostd::CopyFile(directory_tmp + filename_latex + "." + format,directory_work + "/" + filename + "." + format);
      if(LDEBUG) {cerr << __AFLOW_FUNC__ << " moving file to: " << directory_work + "/" + filename + "." + format << endl;}
      // Keep gnuplot file if aflow was called with --keep=gpl
      if (XHOST.vflag_control.flag("KEEP::GPL")) {
        aurostd::CopyFile(directory_tmp + filename + ".plt", directory_work);
      }
      // Clean up
      aurostd::RemoveDirectory(directory_tmp);
      if (!aurostd::FileExist(directory_work + "/" + filename + "." + format)) {
        throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Error while generating plot.", _RUNTIME_ERROR_);
      }
    } else {
      string message = "The following binaries are missing: " + aurostd::joinWDelimiter(missing_binaries, " ") + ".";
      throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message, _RUNTIME_ERROR_);
    }
  }

  // [OBSOLETE] //savePlotGNUPLOT/////////////////////////////////////////////////////////////
  // [OBSOLETE] // Executes the gnuplot script and converts into the desired image format.
  // [OBSOLETE] void savePlotGNUPLOT_OLD(const xoption& plotoptions, const stringstream& gpfile) {
  // [OBSOLETE]   bool LDEBUG=(FALSE || XHOST.DEBUG); 
  // [OBSOLETE]   string directory_work = plotoptions.getattachedscheme("DIRECTORY");
  // [OBSOLETE]   if(directory_work.empty()){directory_work=aurostd::getPWD();}  //[CO20191112 - OBSOLETE]aurostd::execute2string("pwd")//CO20191004
  // [OBSOLETE]   if(LDEBUG) {cerr << __AFLOW_FUNC__ << " directory_work=" << directory_work << endl;}
  // [OBSOLETE]   string filename = plotoptions.getattachedscheme("FILE_NAME");
  // [OBSOLETE]   if(LDEBUG) {cerr << __AFLOW_FUNC__ << " filename=" << filename << endl;}
  // [OBSOLETE]   string filename_latex = plotoptions.getattachedscheme("FILE_NAME_LATEX");
  // [OBSOLETE]   // PDF is default since we use pdflatex to compile
  // [OBSOLETE]   string format = plotoptions.getattachedscheme("IMAGE_FORMAT");
  // [OBSOLETE]   if (format.empty()) format = "pdf";
  // [OBSOLETE]   string current_dir = aurostd::getPWD();  //[CO20191112 - OBSOLETE]aurostd::execute2string("pwd")
  // [OBSOLETE]  // Create temp directory
  // [OBSOLETE]   string directory_tmp = aurostd::TmpDirectoryCreate("plotLATEX") + "/";
  // [OBSOLETE] chdir(directory_tmp.c_str());
  // [OBSOLETE] // Execute gnuplot and pdflatex
  // [OBSOLETE] aurostd::stringstream2file(gpfile, filename + ".plt");
  // [OBSOLETE]   aurostd::execute(XHOST.command("gnuplot") + " " + filename + ".plt");
  // [OBSOLETE]   aurostd::execute(XHOST.command("pdflatex") + " -interaction=nonstopmode -halt-on-error " + filename_latex + ".tex 2>&1 > /dev/null");
  // [OBSOLETE]   // Convert to the desired format if not pdf
  // [OBSOLETE]   if (format != "pdf") {
  // [OBSOLETE]     aurostd::execute(XHOST.command("convert") + " -quiet -density 300 -background white " + filename_latex + ".pdf " + filename_latex  + "." + format);
  // [OBSOLETE]   }
  // [OBSOLETE]   chdir(current_dir.c_str());
  // [OBSOLETE]   aurostd::CopyFile(directory_tmp + filename_latex + "." + format,directory_work + "/" + filename + "." + format);
  // [OBSOLETE]   if(LDEBUG) {cerr << __AFLOW_FUNC__ << " moving file to: " << directory_work + "/" + filename + "." + format << endl;}
  // [OBSOLETE]   // Keep gnuplot file if aflow was called with --keep=gpl
  // [OBSOLETE]   if (XHOST.vflag_control.flag("KEEP::GPL")) {
  // [OBSOLETE]    aurostd::CopyFile(directory_tmp + filename + ".plt", directory_work);
  // [OBSOLETE]   }
  // [OBSOLETE]   // Clean up
  // [OBSOLETE]   aurostd::RemoveDirectory(directory_tmp);
  // [OBSOLETE]   if (!aurostd::FileExist(directory_work + "/" + filename + "." + format)) {
  // [OBSOLETE]      string function = "plotter::savePlotGNUPLOT():";
  // [OBSOLETE]      string message = "Error while generating plot.";
  // [OBSOLETE]     throw aurostd::xerror(__AFLOW_FILE__,function, message, _RUNTIME_ERROR_);
  // [OBSOLETE]   }
  // [OBSOLETE] }

  //setFileName/////////////////////////////////////////////////////////////////
  // Sets the file name of the final plot. FILE_NAME_LATEX is the name of the
  // tex file that is generated by gnuplot, which has different limitations
  // than the output image.
  void setFileName(xoption& plotoptions, string filename) {
    bool LDEBUG=(FALSE || _DEBUG_PLOTTER_ || XHOST.DEBUG);
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " filename_in=" << filename << endl;}
    if (filename.empty()) {
      filename = plotoptions.getattachedscheme("FILE_NAME_USER");  //ME20200313 - user-defined output file
      if (filename.empty()) {
        string default_title = plotoptions.getattachedscheme("DEFAULT_TITLE");
        if(LDEBUG){cerr << __AFLOW_FUNC__ << " default_title=" << default_title << endl;}
        //ME20200228 - Remove ANRL parameters
        string::size_type t = default_title.find(":ANRL=");
        if (t != string::npos) {
          default_title = default_title.substr(0, t);
          if(LDEBUG){std::cerr << __AFLOW_FUNC__ << " default_title (post ANRL)=" << default_title << std::endl;}
        }
        filename = default_title;
        // Get filename
        string ext = plotoptions.getattachedscheme("EXTENSION");
        if (!ext.empty()) {
          if (filename.empty()) filename = ext;
          else filename += "_" + ext;
        }
        filename = aurostd::StringSubst(filename, " ", "_");
        string set = plotoptions.getattachedscheme("DATASET");
        if (aurostd::string2utype<int>(set) > 0) {
          filename += "_" + plotoptions.getattachedscheme("DATALABEL");
          filename += "_" + set;
        }
      }
    }
    plotoptions.push_attached("FILE_NAME", filename);
    //filename_latex is the name of the .tex file
    //the resulting plot file name will be changed (mv) to filename
    string filename_latex=filename;
    aurostd::StringSubst(filename_latex, ".", "_"); // The .tex file created by gnuplot cannot have . or includegraphics will break
    aurostd::StringSubst(filename_latex, ":", "_"); //ME20200409 - Some terminals do not handle : well, which can break convert
    plotoptions.push_attached("FILE_NAME_LATEX", filename_latex);
    if(LDEBUG){
      cerr << __AFLOW_FUNC__ << " filename=" << plotoptions.getattachedscheme("FILE_NAME") << endl;
      cerr << __AFLOW_FUNC__ << " filename_latex=" << plotoptions.getattachedscheme("FILE_NAME_LATEX") << endl;
    }
  }

  //setTitle//////////////////////////////////////////////////////////////////
  // Sets the plot title.
  void setTitle(xoption& plotoptions,ostream& oss) {ofstream FileMESSAGE;return setTitle(plotoptions,FileMESSAGE,oss);} //CO20200404
  void setTitle(xoption& plotoptions,ofstream& FileMESSAGE,ostream& oss) { //CO20200404
    string title = plotoptions.getattachedscheme("TITLE");
    // Format title
    if (title.empty()) title = formatDefaultPlotTitle(plotoptions,FileMESSAGE,oss);
    plotoptions.push_attached("PLOT_TITLE", title);
  }

  //formatDefaultPlotTitle////////////////////////////////////////////////////
  // Checks if the default title is in a known AFLOW format and formats it
  // appropriately.
  string formatDefaultPlotTitle(const xoption& plotoptions,ostream& oss) {ofstream FileMESSAGE;return formatDefaultPlotTitle(plotoptions,FileMESSAGE,oss);} //CO20200404
  string formatDefaultPlotTitle(const xoption& plotoptions,ofstream& FileMESSAGE,ostream& oss) { //CO20200404
    bool LDEBUG=(FALSE || _DEBUG_PLOTTER_ || XHOST.DEBUG);
    string default_title = plotoptions.getattachedscheme("DEFAULT_TITLE");
    if(LDEBUG) {cerr << __AFLOW_FUNC__ << " default_title=" << default_title << endl;}
    if (default_title.empty()) return default_title;
    string title="";
    if (default_title.find("_ICSD_")!=string::npos) {  // Check if AFLOW ICSD format
      vector<string> tokens;
      aurostd::string2tokens(default_title, tokens, "_");
      if (tokens.size() == 3) {
        title = pflow::prettyPrintCompound(tokens[0], no_vrt, true, latex_ft) + " (ICSD \\#" + tokens[2];  //_none_ //_latex_ //CO20190629
        string lattice = plotoptions.getattachedscheme("LATTICE");
        if (lattice.empty()) title += ")";
        else title += ", " + lattice + ")";
      } else { // Title not in ICSD format
        return aurostd::fixStringLatex(default_title, false, false);
      }
    } else if (aurostd::substring2bool(default_title, TAG_TITLE_POCC)) {  // Check if in POCC format
      if(LDEBUG){cerr << __AFLOW_FUNC__ << " found POCC" << endl;}
      title = formatDefaultTitlePOCC(plotoptions,FileMESSAGE,oss); //CO20200404
    } else if (aurostd::substring2bool(default_title, ".")) {  // Check if AFLOW prototype format
      vector<string> tokens;
      aurostd::string2tokens(default_title, tokens, ".");
      //ME20200228 - title may contain ANRL parameters
      if ((tokens.size() > 2) && aurostd::substring2bool(tokens[2], "ANRL")) {
        string::size_type t = tokens[2].find_first_of(":");
        if (t != string::npos) {
          tokens[2] = tokens[2].substr(0, t);
          tokens.erase(tokens.begin() + 3, tokens.end());
        }
      }
      string comp_str=tokens[0];
      string::size_type loc=comp_str.find(":"); //remove LIB1 pp junk
      comp_str=comp_str.substr(0,loc);
      if ((tokens.size() == 2) || (tokens.size() == 3)) {
        string proto = tokens[1];
        if(LDEBUG){cerr << __AFLOW_FUNC__ << " proto=" << proto << endl;}
        vector<string> protos;
        aflowlib::GetAllPrototypeLabels(protos, "anrl");
        if (aurostd::WithinList(protos, proto)) {
          if(LDEBUG){cerr << __AFLOW_FUNC__ << " found proto in ANRL" << endl;}
          if (tokens.size() == 3) proto += "." + tokens[2];
          vector<string> elements = aurostd::getElements(comp_str);
          vector<double> composition = getCompositionFromANRLPrototype(proto);
          if(LDEBUG){
            cerr << __AFLOW_FUNC__ << " elements=" << aurostd::joinWDelimiter(elements,",") << endl;
            cerr << __AFLOW_FUNC__ << " composition=" << aurostd::joinWDelimiter(aurostd::vecDouble2vecString(composition,5),",") << endl;
          }
          proto = aurostd::fixStringLatex(proto, false, false); // Prevent LaTeX errors
          title = pflow::prettyPrintCompound(elements, composition, no_vrt, true, latex_ft) + " (" + proto;  //_none_ //_latex_ //CO20190629
        } else {
          if (tokens.size() == 3) proto += "." + tokens[2];
          vector<string> comp;
          protos.clear();
          aflowlib::GetAllPrototypeLabels(protos, comp, "htqc");
          int index = 0;
          if (aurostd::WithinList(protos, proto, index)) {
            proto = aurostd::fixStringLatex(proto, false, false); // Prevent LaTeX errors
            vector<string> elements = aurostd::getElements(comp_str);
            vector<double> composition = getCompositionFromHTQCPrototype(proto, comp[index]);
            if(LDEBUG){
              cerr << __AFLOW_FUNC__ << " elements=" << aurostd::joinWDelimiter(elements,",") << endl;
              cerr << __AFLOW_FUNC__ << " composition=" << aurostd::joinWDelimiter(aurostd::vecDouble2vecString(composition,5),",") << endl;
            }
            //DX+ME20210729 - the code below is wrapped in a try-catch block because there are unary prototypes in LIB2 (e.g., /common/LIB2/LIB/LaMn_pv/117)
            // in these instances, we will return the default title
            try { title = pflow::prettyPrintCompound(elements, composition, no_vrt, true, latex_ft) + " (" + proto; } //_none_ //_latex_   //CO20190629 //DX+ME20210729 - for unaries in LIB2 errors
            catch(aurostd::xerror& excpt){ return aurostd::fixStringLatex(default_title, false, false); } //DX+ME20210729 - for unaries in LIB2 errors
          } else {  // Title not in prototype format
            return aurostd::fixStringLatex(default_title, false, false);
          }
        }
      } else {
        return aurostd::fixStringLatex(default_title, false, false); //CO20191110
      }
      string lattice = plotoptions.getattachedscheme("LATTICE");
      if (lattice.empty()) title += ")";
      else title += ", " + lattice + ")";
    } else {  // Not an AFLOW-formatted default
      return aurostd::fixStringLatex(default_title, false, false);
    }
    if(LDEBUG) {cerr << __AFLOW_FUNC__ << " title=" << title << endl;}
    // Code only gets here if the title is AFLOW-formatted
    string set = plotoptions.getattachedscheme("DATASET");
    if (aurostd::string2utype<int>(set) > 0) {
      title += " " + aurostd::fixStringLatex(plotoptions.getattachedscheme("SETLABEL"),false,false); //CO20191110
      title += " " + aurostd::fixStringLatex(plotoptions.getattachedscheme("DATALABEL"),false,false);  //CO20191110
      title += " (" + set + ")";
    }
    return title;
  }

  //getCompositionFromHTQCPrototype///////////////////////////////////////////
  // Gets the composition from an HTQC prototype string. The composition
  // string must be retrieved beforehand.
  vector<double> getCompositionFromHTQCPrototype(const string& htqc_prototype,
      const string& composition) {
    string anrl_prototype = composition + "_";
    // Composition already has the correct sequence, so keeping the explicit
    // sequence designator will confuse getCompositonFromANRLPrototype
    string::size_type t = htqc_prototype.find(".");
    anrl_prototype += htqc_prototype.substr(0, t);
    return getCompositionFromANRLPrototype(anrl_prototype);
  }

  //getCompositionFromANRLPrototype///////////////////////////////////////////
  // Gets the composition from an ANRL prototype string.
  vector<double> getCompositionFromANRLPrototype(const string& prototype) {
    // Determine element sequence
    // If there is a . in the prototype string, the element sequence is given explicitly
    string seq;
    string::size_type t = prototype.find(".");
    if (t != string::npos) seq = prototype.substr(t + 1, string::npos);
    t = prototype.find("_");
    string compound = prototype.substr(0, t);
    // If not explicitly given, determine element sequence from the prototype name
    if (seq.empty()) {
      for (uint i = 0; i < compound.size(); i++) {
        if (isalpha(compound[i])) seq += compound[i];
      }
    }
    vector<int> sequence(seq.size());
    for (uint i = 0; i < seq.size(); i++) {
      sequence[i] = (int) seq[i] - 65;
    }

    // Now determine the composition
    vector<double> comp;
    aurostd::getElements(compound, comp);

    // Finally, sort to match sequence
    vector<double> composition(comp.size());
    for (uint i = 0; i < composition.size(); i++) {
      composition[i] = comp[sequence[i]];
    }
    return composition;
  }

  //formatDefaultTitlePOCC//////////////////////////////////////////////////////
  // Converts a POCC-formatted title into a plot title. It currently only works
  // if the POCC string consists only of P-designations.
  string formatDefaultTitlePOCC(const xoption& plotoptions,ostream& oss) {ofstream FileMESSAGE;return formatDefaultTitlePOCC(plotoptions,FileMESSAGE,oss);} //CO20191110  //CO20200404
  string formatDefaultTitlePOCC(const xoption& plotoptions,ofstream& FileMESSAGE,ostream& oss) {return formatDefaultTitlePOCC_20191004(plotoptions,FileMESSAGE,oss);} //CO20191110  //CO20200404
  string formatDefaultTitlePOCC_20191004(const xoption& plotoptions,ostream& oss) {ofstream FileMESSAGE;return formatDefaultTitlePOCC_20191004(plotoptions,FileMESSAGE,oss);}  //CO version //CO20191110  //CO20200404
  string formatDefaultTitlePOCC_20191004(const xoption& plotoptions,ofstream& FileMESSAGE,ostream& oss) {  //CO version //CO20191110  //CO20200404
    bool LDEBUG=(FALSE || _DEBUG_PLOTTER_ || XHOST.DEBUG);
    stringstream message;
    string default_title = plotoptions.getattachedscheme("DEFAULT_TITLE");
    if(LDEBUG) {cerr << __AFLOW_FUNC__ << " default_title=" << default_title << endl;}

    aurostd::xoption pocc_settings;
    xstructure xstr;
    bool found_pocc=pflow::POccInputs2Xstr(default_title,pocc_settings,xstr,FileMESSAGE,oss);
    if(!found_pocc){return aurostd::fixStringLatex(default_title, false, false);}

    string pps=pocc_settings.getattachedscheme("PPS");
    string proto=pocc_settings.getattachedscheme("PROTO");
    string pocc_params=pocc_settings.getattachedscheme("POCC_PARAMS");
    string pocc_tol=pocc_settings.getattachedscheme("POCC_TOL");
    string pocc_arun=pocc_settings.getattachedscheme("POCC_ARUN");
    string module_arun=pocc_settings.getattachedscheme("MODULE_ARUN");

    string new_title="";
    string clean_specie="";
    int comp_prec=(int)ceil(log10(1.0/xstr.partial_occupation_stoich_tol));  //ceil ensures we round up above 1 //CO20181226
    for(uint ispecies=0;ispecies<xstr.species.size();ispecies++){
      clean_specie=KBIN::VASP_PseudoPotential_CleanName(xstr.species[ispecies]);
      if(LDEBUG) {cerr << __AFLOW_FUNC__ << " species[ispecies=" << ispecies << "]=" << clean_specie << endl;}
      new_title+=aurostd::fixStringLatex(clean_specie,false,false);
      new_title+=(aurostd::isequal(xstr.comp_each_type[ispecies],1.0,xstr.partial_occupation_stoich_tol) ? "" : "$_{"+aurostd::utype2string(xstr.comp_each_type[ispecies],comp_prec)+"}$");
    }
    new_title+=" ("+aurostd::fixStringLatex(proto,false,false);

    vector<string> tokens;
    //ARUN.POCC_49
    if(!pocc_arun.empty()){
      aurostd::string2tokens(pocc_arun,tokens,"_");
      string pocc_hash="";
      if(tokens.size()>1){pocc_hash=tokens.back();}
      if(!pocc_hash.empty()){new_title+=":"+aurostd::fixStringLatex(pocc_hash,false,false);}
    }
    //module_arun=ARUN.AGL_6_SF_0.92
    if(!module_arun.empty()){
      aurostd::string2tokens(module_arun,tokens,"_");
      string module_hash="";
      if(tokens.size()>1){
        vector<string> new_tokens;
        for(uint i=2;i<tokens.size();i++){new_tokens.push_back(tokens[i]);}
        module_hash=aurostd::joinWDelimiter(new_tokens,"_");
      }
      if(!module_hash.empty()){new_title+=":"+aurostd::fixStringLatex(module_hash,false,false);}
    }

    new_title+=")";

    if(LDEBUG) {cerr << __AFLOW_FUNC__ << " new_title=" << new_title << endl;}

    return new_title; //aurostd::fixStringLatex(new_title, false, false);  //substs $ for \\$
  }
  string formatDefaultTitlePOCC_20190101(const xoption& plotoptions) {  //ME version
    bool LDEBUG=(FALSE || _DEBUG_PLOTTER_ || XHOST.DEBUG);
    string default_title = plotoptions.getattachedscheme("DEFAULT_TITLE");
    //Get all the pieces of the default title
    string::size_type t = default_title.find(":POCC");
    string proto = default_title.substr(0, t);  // Contains compound and prototype
    string pocc = default_title.substr(t + 5, string::npos);  // POCC string + ARUN
    if (LDEBUG) {
      std::cerr << "proto = " << proto << std::endl;
      std::cerr << "pocc = " << pocc << std::endl;
    }
    bool generic = false;
    // Need _S because S could theoretically also be a decorator
    if (aurostd::substring2bool(pocc, "_S")) generic = true;
    if(LDEBUG) {cerr << __AFLOW_FUNC__ << " found _S tag = " << generic << endl;}
    pocc = pocc.substr(1, pocc.size());  // Remove the leading _

    // Get the HNF matrix string
    vector<string> tokens;
    string hnf = "";
    if (!generic && aurostd::substring2bool(pocc, ":")) {  // Is there an ARUN?
      t = pocc.find(":");
      hnf = pocc.substr(t + 1, string::npos);
      pocc = pocc.substr(0, t);  // Remove ARUN from pocc
      if(LDEBUG) {cerr << __AFLOW_FUNC__ << " hnf=" << hnf << endl;}
      if (!hnf.empty()) {
        aurostd::string2tokens(hnf, tokens, "_");
        hnf = tokens.back();
      }
    }

    // Try to extract the composition from the POCC string and check whether
    // the string is incomplete. While the composition is only supported for
    // P-designations, the algorithm should still run to test for incomplete
    // strings.
    bool broken = false;
    vector<double> composition = getCompositionFromPoccString(pocc, broken);
    if(LDEBUG){  //CO20191110
      cerr << __AFLOW_FUNC__ << " broken=" << broken << endl;
      if(!broken){
        cerr << __AFLOW_FUNC__ << " composition=";
        for(uint i=0;i<composition.size();i++){
          cerr << composition[i] << (i==composition.size()-1?"":",");
        }
        cerr << endl;
      }
    }
    if (broken) {  // If broken, extract the POCC string from the POSCAR title
      try {
        stringstream ss;
        vector<string> vstr;
        string directory = plotoptions.getattachedscheme("DIRECTORY");
        string extension = plotoptions.getattachedscheme("EXTENSION");
        if(LDEBUG){
          cerr << __AFLOW_FUNC__ << " directory=" << directory << endl;
          cerr << __AFLOW_FUNC__ << " extension=" << extension << endl;
        }
        if (aurostd::substring2bool(extension, "phdisp") || (extension == "phdos")) {
          aurostd::efile2vectorstring("PHPOSCAR", vstr);
        } else {
          aurostd::efile2vectorstring(aflowlib::vaspfile2stringstream(directory, "POSCAR", ss), vstr);
        }
        // The POCC string is inside the first item
        t = vstr[0].find(" ");
        string poscartitle = vstr[0].substr(0, t);
        // Only take what is needed by the algorithm
        t = poscartitle.find(":POCC_");
        poscartitle = poscartitle.substr(t + 6, string::npos);
        t = poscartitle.find(":");
        poscartitle = poscartitle.substr(0, t);
        composition = getCompositionFromPoccString(poscartitle, broken);
      } catch (aurostd::xerror& excpt) {
        generic = true;
      }
    }

    // Separate elements and prototype
    string compound, el;
    t = proto.find(".");
    el = proto.substr(0, t);
    proto = proto.substr(t + 1, string::npos);
    if (!generic) {
      vector<string> elements = aurostd::getElements(el);
      if (elements.size() != composition.size()) {
        generic = true;
        broken = true;
      } else {
        compound = pflow::prettyPrintCompound(elements, composition, no_vrt, true, latex_ft);  //_none_ //_latex_ //CO20190629
      }
    }
    if (generic) {  // Broken or unsupported string, so use a very generic title
      return aurostd::fixStringLatex(default_title, false, false); //CO20191110
      //if (!broken) proto += ".POCC:" + pocc;  // Only add POCC string if not broken
    }
    proto = aurostd::fixStringLatex(proto, false, false);

    // Finish title
    string title = compound + " (" + proto;
    if (!hnf.empty()) title += ":" + hnf;
    string lattice = plotoptions.getattachedscheme("LATTICE");
    if (lattice.empty()) title += ")";
    else title += ", " + lattice + ")";
    return title;
  }

  //getCompositionFromPoccString//////////////////////////////////////////////
  // Returns a composition from a POCC string. It also tests whether the POCC
  // string is complete, which is not always the case due to VASP's character
  // limit for titles.
  vector<double> getCompositionFromPoccString(const string& pocc_string, bool& broken) {
    broken = false;
    double SITE_TOL = 0.001;
    // Get each site
    vector<string> tokens;
    aurostd::string2tokens(pocc_string, tokens, "_");
    // Get the composition of each individual site. Do not add yet - it needs
    // to be site-resolved for some checks to work.
    vector<vector<std::pair<int, double> > > sites(tokens.size());
    std::pair<int, double> s;
    vector<string> site;
    int max_index = 0;  // tracks how many decorators can be found in the string
    for (uint i = 0; i < tokens.size(); i++) {
      // Loop over all deocrations on the site
      aurostd::string2tokens(tokens[i], site, "-");
      if (site.size() == 1) { // There must be at least one site (first element is designator)
        broken = true;
        break;
      } else {
        for (uint j = 1; j < site.size(); j++) {
          string::size_type t;
          std::pair<string, string> str_cut;
          t = site[j].find("x");
          if (t != string::npos) {
            str_cut.first = site[j].substr(0, t);
            str_cut.second = site[j].substr(t + 1, string::npos);
          }
          // Format is e.g. Ax0.5, so there must be two elements
          if (str_cut.first.empty() || str_cut.second.empty()) {
            broken = true;
            break;
          } else {
            s.first = (int) str_cut.second[0] - 65;  // [0] is to convert to char
            s.second = aurostd::string2utype<double>(str_cut.first);
            sites[i].push_back(s);
            if (s.first > max_index) max_index = s.first;
          }
        }
      }
    }

    // Determine the composition and check for consistency: no element must
    // can a value of zero and all sites have to add to 1 within a tolerance
    double sum;
    vector<double> composition(max_index + 1);
    for (uint i = 0; i < sites.size(); i++) {
      sum = 0.0;
      for (uint j = 0; j < sites[i].size(); j++) {
        composition[sites[i][j].first] += sites[i][j].second;
        sum += sites[i][j].second;
      }
      if (std::abs(sum - 1.0) > SITE_TOL) broken = true;
    }
    if (!broken) {
      for (uint i = 0; i < composition.size(); i++) {
        if (composition[i] == 0.0) {
          broken = true;
          break;
        }
      }
    }

    return composition;
  }

} // namespace plotter

//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//                           ELECTRONIC STRUCTURE                           //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////

static const string EFERMI_COLOR = "#0000FF";
static const int ESTRUCTURE_NCOLORS = 11;
static const string ESTRUCTURE_COLORS[ESTRUCTURE_NCOLORS] = {
  "#000000",  // black
  "#4C72B0",  // blue
  "#55A868",  // green
  "#C44E52",  // red
  "#CCB974",  // yellow
  "#8172B2",  // purple
  "#64B5CD",  // light blue
  "#E08000",  // orange
  "#006060",  // blue-green
  "#A06000",  // brown
  "#BE80FF"   // light purple
};
static const string ISPIN_COLORS[2] = {"#000000", "#C44E52"};
static const string ORBITALS[4] = {"s", "p", "d", "f"};
static const string LM_ORBITALS_LATEX[16] = {"s", "p_y", "p_z", "p_x",
  "d_{xy}", "d_{yz}", "d_{z^2}", "d_{xz}", "d_{x^2-y^2}",
  "f_{y(3x^2-y^2)}", "f_{xyz}", "f_{yz^2}", "f_{z^3}",
  "f_{xz^2}", "f_{z(x^2-y^2)}", "f_{x(x^2-3y^2)}"};
//AS20201110 as used by E.G. in BANDSDATA_JSON procedure:
static const string LM_ORBITALS[16] = {"s", "py", "pz", "px",
  "dxy", "dyz", "dz2", "dxz", "dx2-y2",
  "f1", "f2", "f3", "f4", "f5", "f6", "f7"};
static const string SPIN_LABEL[2] = {"majority", "minority"};//AS20201102
static const string LS_ORBITALS[16] = {"s_total", "sx", "sy", "sz",
  "p_total", "px", "py", "pz", "d_total", "dx", "dy", "dz",
  "f_total", "fx", "fy", "fz"};//AS20201105

namespace plotter {

  // Plot functions -----------------------------------------------------------

  //PLOT_DOS//////////////////////////////////////////////////////////////////
  // PLots electronic densities of states.
  void PLOT_DOS(xoption& plotoptions,ostream& oss) {ofstream FileMESSAGE;return PLOT_DOS(plotoptions,FileMESSAGE,oss);} //CO20200404
  void PLOT_DOS(xoption& plotoptions,ofstream& FileMESSAGE,ostream& oss) {  //CO20200404
    // Set output format to gnuplot
    plotoptions.push_attached("OUTPUT_FORMAT", "GNUPLOT");

    stringstream out;
    PLOT_DOS(plotoptions, out,FileMESSAGE,oss);
    savePlotGNUPLOT(plotoptions, out);
  }

  void PLOT_DOS(xoption& plotoptions, const xDOSCAR& xdos,ostream& oss) {ofstream FileMESSAGE;return PLOT_DOS(plotoptions,xdos,FileMESSAGE,oss);} //CO20191110  //CO20200404
  void PLOT_DOS(xoption& plotoptions, const xDOSCAR& xdos,ofstream& FileMESSAGE,ostream& oss) { //CO20191110  //CO20200404
    plotoptions.push_attached("OUTPUT_FORMAT", "GNUPLOT");
    stringstream out;
    PLOT_DOS(plotoptions,out,xdos,FileMESSAGE,oss);
    savePlotGNUPLOT(plotoptions, out);
  }

  void PLOT_DOS(xoption& plotoptions, stringstream& out,ostream& oss) {ofstream FileMESSAGE;return PLOT_DOS(plotoptions,out,FileMESSAGE,oss);} //CO20191110 //CO20200404
  void PLOT_DOS(xoption& plotoptions, stringstream& out,ofstream& FileMESSAGE,ostream& oss) { //CO20191110  //CO20200404
    bool LDEBUG=(FALSE || _DEBUG_PLOTTER_ || XHOST.DEBUG);

    // Read files
    string directory = plotoptions.getattachedscheme("DIRECTORY");
    xDOSCAR xdos;
    if(LDEBUG) {cerr << __AFLOW_FUNC__ << " directory=" << directory << endl;}
    //CO20210701 - adding support for POCC
    //check if POCC directory
    vector<string> vfiles,vpocc_doscars;
    aurostd::DirectoryLS(directory,vfiles);
    uint i=0;
    for(i=0;i<vfiles.size();i++){
      if(vfiles[i].find(POCC_DOSCAR_PREFIX)!=string::npos){vpocc_doscars.push_back(aurostd::CleanFileName(directory+"/"+vfiles[i]));}
    }
    if(vpocc_doscars.size()>0){
      std::sort(vpocc_doscars.begin(),vpocc_doscars.end());
      string pscheme=plotoptions.getattachedscheme("PROJECTION");
      if(LDEBUG){cerr << __AFLOW_FUNC__ << " pscheme=" << pscheme << endl;}
      if(pscheme.empty()){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"No projection scheme provided",_INPUT_MISSING_);}
      double temperature=0.0;
      _aflags aflags;aflags.Directory=directory;
      bool see_sub_output=false;//true;
      ostream& oss_empty=cout;if(!see_sub_output){oss_empty.setstate(std::ios_base::badbit);}  //like NULL
      ofstream devnull("/dev/null");  //NULL
      //T0K and T0000K cannot be known from a single file
      //there is logic inside to determine whether we have floats, precision, and how to pad zeros. 
      //it depends on the temperature range specified in the aflow.in (or command line). 
      //however, if this cannot be found (because we only have the DOSCAR's), 
      //then the logic says to understand the files locally as best as possible.
      int temperature_precision=0,zero_padding_temperature=0;
      bool temperatures_int=true;
      //get temperature string parameters, either we grab them from the aflow.in or we guess them from defaults in aflow.rc
      if(aurostd::FileExist(aflags.Directory+"/"+_AFLOWIN_) && pocc::structuresGenerated(aflags.Directory)){
        pocc::POccCalculator pcalc(aflags,devnull,oss_empty);
        vector<double> v_temperatures;
        pcalc.loadDataIntoCalculator();pcalc.setTemperatureStringParameters(v_temperatures); //needed for DOSCAR plots
        pocc::getTemperatureStringParameters(v_temperatures,temperature_precision,temperatures_int,zero_padding_temperature);
        if(pcalc.m_ARUN_directories.size()==0){throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"No ARUN.POCC_* runs found",_FILE_CORRUPT_);}
        plotoptions.push_attached("ARUN_DIRECTORY",pcalc.m_ARUN_directories[0]);
      }else{
        pocc::getTemperatureStringParameters(temperature_precision,temperatures_int,zero_padding_temperature);
        //this is when we don't have the full POCC directory. 
        //maybe someone is trying to plot in a directory only having aflow.in and DOSCAR's (but no ARUN's). 
        //we need the ARUNs if we want to do projection==species, but not if projection==orbital
        if(pscheme=="SPECIES"){ //we need ARUN_DIRECTORY
          throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,"Need ARUN.POCC_* runs for projection==\"species\"",_INPUT_MISSING_);
        }
      }
      if(LDEBUG){cerr << __AFLOW_FUNC__ << " found POCC directory" << endl;}
      for(i=0;i<vpocc_doscars.size();i++){
        if(LDEBUG){cerr << __AFLOW_FUNC__ << " looking at POCC DOSCAR: " << vpocc_doscars[i] << endl;}
        xdos.GetPropertiesFile(vpocc_doscars[i]);
        temperature=pocc::poccDOSCAR2temperature(vpocc_doscars[i]);
        plotoptions.push_attached("NORMALIZATION","ATOM");  //CO20211124
        plotoptions.push_attached("EXTENSION","dos_"+aurostd::tolower(pscheme)+"_T"+pocc::getTemperatureString(temperature,temperature_precision,temperatures_int,zero_padding_temperature)+"K");
        aurostd::StringstreamClean(out);  //CO20211210 - otherwise it writes out of other temperatures
        PLOT_DOS(plotoptions, out, xdos,FileMESSAGE,oss);  //CO20200404
        savePlotGNUPLOT(plotoptions, out);
      }
    }else{
      xdos.GetPropertiesFile(aflowlib::vaspfile2stringstream(directory, "DOSCAR"));
      PLOT_DOS(plotoptions, out, xdos,FileMESSAGE,oss);  //CO20200404
      savePlotGNUPLOT(plotoptions, out);
    }
  }

  void patchDefaultTitleAFLOWIN(xoption& plotoptions) { //CO20191110
    bool LDEBUG=(FALSE || _DEBUG_PLOTTER_ || XHOST.DEBUG);

    const string& directory = plotoptions.getattachedscheme("DIRECTORY");
    if (!aurostd::FileExist(directory + "/" + _AFLOWIN_)) return;  //ME20200922
    if(LDEBUG) {cerr << __AFLOW_FUNC__ << " directory=" << directory << endl;}
    string SYSTEM=KBIN::ExtractSystemName(directory); //CO20200731
    if(!SYSTEM.empty()){
      if(LDEBUG) {cerr << __AFLOW_FUNC__ << " DEFAULT_TITLE(OLD)=" << plotoptions.getattachedscheme("DEFAULT_TITLE") << endl;}
      plotoptions.pop_attached("DEFAULT_TITLE");
      plotoptions.push_attached("DEFAULT_TITLE", SYSTEM);
      if(LDEBUG) {cerr << __AFLOW_FUNC__ << " DEFAULT_TITLE(NEW)=" << plotoptions.getattachedscheme("DEFAULT_TITLE") << endl;}
    }
  }

  void PLOT_DOS(xoption& plotoptions, stringstream& out, const xDOSCAR& xdos,ostream& oss) {ofstream FileMESSAGE;return PLOT_DOS(plotoptions,out,xdos,FileMESSAGE,oss);}  //CO20200404
  void PLOT_DOS(xoption& plotoptions, stringstream& out, const xDOSCAR& xdos,ofstream& FileMESSAGE,ostream& oss) {  //CO20200404
    bool LDEBUG=(FALSE || _DEBUG_PLOTTER_ || XHOST.DEBUG);
    string extension=plotoptions.getattachedscheme("EXTENSION");
    if(extension.empty()) plotoptions.push_attached("EXTENSION", "dos");
    // Make sure the projections are consistent with the DOSCAR file
    if ((plotoptions.getattachedscheme("PROJECTION") == "LM") && !(xdos.lmResolved)) {
      std::cerr << "Found --projection=lm, but DOSCAR is not lm-resolved."
        << " Will choose --projection=orbitals instead." << std::endl;
      plotoptions.pop_attached("PROJECTION");
      plotoptions.push_attached("PROJECTION", "ORBITALS");
    }

    plotoptions.push_attached("DEFAULT_TITLE", xdos.title);
    patchDefaultTitleAFLOWIN(plotoptions);  //CO20191110 - ME, check out and let me know if we should apply everywhere
    setFileName(plotoptions);
    setTitle(plotoptions,FileMESSAGE,oss); //CO20200404

    plotoptions.flag("BANDDOS", false);
    // Set Fermi energy to zero
    if (plotoptions.flag("NOSHIFT")) {
      plotoptions.push_attached("EFERMI", aurostd::utype2string<double>(xdos.Efermi));
    } else {
      plotoptions.push_attached("EFERMI", "0.0");
    }

    if(LDEBUG) {cerr << __AFLOW_FUNC__ << " EFERMI set" << endl;}

    // Set Emin and Emax
    setEMinMax(plotoptions, xdos.energy_min, xdos.energy_max);

    // Plot
    generateHeader(out, plotoptions, false);
    generateDosPlot(out, xdos, plotoptions,FileMESSAGE,oss);  //CO20200404
  }

  //PLOT_PDOS/////////////////////////////////////////////////////////////////
  // Plots projected density of states. If PDOS == -1, the projected DOS of
  // all atoms will be plotted into separate files.
  void PLOT_PDOS(xoption& plotoptions,ostream& oss) {ofstream FileMESSAGE;return PLOT_PDOS(plotoptions,FileMESSAGE,oss);} //CO20200404  //CO20200404
  void PLOT_PDOS(xoption& plotoptions,ofstream& FileMESSAGE,ostream& oss) { //CO20200404  //CO20200404
    // Set output format to gnuplot
    plotoptions.push_attached("OUTPUT_FORMAT", "GNUPLOT");

    stringstream out;
    PLOT_PDOS(plotoptions, out,FileMESSAGE,oss);
  }

  void PLOT_PDOS(xoption& plotoptions, const xDOSCAR& xdos,ostream& oss) {ofstream FileMESSAGE;return PLOT_PDOS(plotoptions,xdos,FileMESSAGE,oss);}  //CO20191110  //CO20200404
  void PLOT_PDOS(xoption& plotoptions, const xDOSCAR& xdos,ofstream& FileMESSAGE,ostream& oss) {  //CO20191110  //CO20200404
    plotoptions.push_attached("OUTPUT_FORMAT", "GNUPLOT");
    stringstream out;
    PLOT_PDOS(plotoptions,out,xdos,FileMESSAGE,oss);
  }

  void PLOT_PDOS(xoption& plotoptions, stringstream& out,ostream& oss) {ofstream FileMESSAGE;return PLOT_PDOS(plotoptions,out,FileMESSAGE,oss);}  //CO20191110  //CO20200404
  void PLOT_PDOS(xoption& plotoptions, stringstream& out,ofstream& FileMESSAGE,ostream& oss) {  //CO20191110  //CO20200404
    // Read files
    string directory = plotoptions.getattachedscheme("DIRECTORY");
    xDOSCAR xdos;
    xdos.GetPropertiesFile(aflowlib::vaspfile2stringstream(directory, "DOSCAR"));
    PLOT_PDOS(plotoptions, out, xdos,FileMESSAGE,oss);
  }

  void PLOT_PDOS(xoption& plotoptions, stringstream& out, const xDOSCAR& xdos,ostream& oss) {ofstream FileMESSAGE;return PLOT_PDOS(plotoptions,out,xdos,FileMESSAGE,oss);} //CO20200404
  void PLOT_PDOS(xoption& plotoptions, stringstream& out, const xDOSCAR& xdos,ofstream& FileMESSAGE,ostream& oss) { //CO20200404
    string set = plotoptions.getattachedscheme("DATASET");
    string datatype=plotoptions.getattachedscheme("DATATYPE");  //CO20191010 - which table to look at, "atoms-projected" by default
    string outformat = plotoptions.getattachedscheme("OUTPUT_FORMAT");
    bool plot_all_atoms = plotoptions.flag("PLOT_ALL_ATOMS"); //CO20191010
    int pdos = -1;
    if (!set.empty()) pdos = aurostd::string2utype<int>(set);
    plotoptions.push_attached("SETLABEL", "PDOS");
    if (pdos == 0) {  // Plot total DOS
      PLOT_DOS(plotoptions,xdos,FileMESSAGE,oss);  //CO20200404
      if (outformat == "GNUPLOT") savePlotGNUPLOT(plotoptions, out);
    } else {
      if (!xdos.partial) {
        std::cerr << "plotter:PLOT_PDOS(): No partial DOS available." << std::endl;
        return;
      }
      xstructure xstr = getStructureWithNames(plotoptions,FileMESSAGE,xdos.carstring,oss); //getStructureWithNames(plotoptions);  //CO20191010 //CO20200404
      if(datatype=="SPECIES"){  //CO20191010 - plot "species-projected"
        if (pdos == -1) {  // Plot partial DOS of all species
          uint nspecies=xstr.num_each_type.size();
          for (uint isp = 0; isp < nspecies; isp++) {
            plotoptions.pop_attached("DATASET");
            plotoptions.pop_attached("DATALABEL");
            plotoptions.push_attached("DATASET", aurostd::utype2string<int>(isp + 1));
            plotoptions.push_attached("DATALABEL", xstr.species[isp]);
            PLOT_DOS(plotoptions, out, xdos,FileMESSAGE,oss);  //CO20200404
            if (outformat == "GNUPLOT") {
              savePlotGNUPLOT(plotoptions, out);
              out.str(string());
              out.clear();
            }
          }
        } else {
          plotoptions.push_attached("DATALABEL", xstr.species[pdos - 1]);
          PLOT_DOS(plotoptions, out, xdos,FileMESSAGE,oss);  //CO20200404
          if (outformat == "GNUPLOT") savePlotGNUPLOT(plotoptions, out);
        }
      }else{  //CO20191010 - plot "atoms-projected"
        if (pdos == -1) {  // Plot partial DOS of all inequivalent atoms
          uint natoms=0;
          if (plot_all_atoms) {natoms = xstr.atoms.size();} //CO20191010 
          else {
            pflow::PerformFullSymmetry(xstr);
            natoms = xstr.iatoms.size();
          }
          int iat=0;
          for (uint i = 0; i < natoms; i++) {
            iat=( plot_all_atoms ? i : xstr.iatoms[i][0] );
            if (!dosDataAvailable(xdos.vDOS, iat + 1)) {
              std::cerr << "plotter:PLOT_PDOS(): No partial DOS available for atom " << iat << "." << std::endl;
              return;
            }
            xstr.atoms[iat].CleanName();
            plotoptions.pop_attached("DATASET");
            plotoptions.pop_attached("DATALABEL");
            plotoptions.push_attached("DATASET", aurostd::utype2string<int>(iat + 1));
            plotoptions.push_attached("DATALABEL", xstr.atoms[iat].cleanname);
            PLOT_DOS(plotoptions, out, xdos,FileMESSAGE,oss);  //CO20200404
            if (outformat == "GNUPLOT") {
              savePlotGNUPLOT(plotoptions, out);
              out.str(string());
              out.clear();
            }
          }
        } else {
          if (!dosDataAvailable(xdos.vDOS, pdos)) {
            std::cerr << "plotter:PLOT_PDOS(): No partial DOS available for index " << pdos << "." << std::endl;
            return;
          }
          xstr.atoms[pdos - 1].CleanName();
          plotoptions.push_attached("DATALABEL", xstr.atoms[pdos - 1].cleanname);
          PLOT_DOS(plotoptions, out, xdos,FileMESSAGE,oss);  //CO20200404
          if (outformat == "GNUPLOT") savePlotGNUPLOT(plotoptions, out);
        }
      }
    }
  }

  //PLOT_BAND/////////////////////////////////////////////////////////////////
  // Plots band structures.
  void PLOT_BAND(xoption& plotoptions,ostream& oss) {ofstream FileMESSAGE;return PLOT_BAND(plotoptions,FileMESSAGE,oss);}  //CO20200404
  void PLOT_BAND(xoption& plotoptions,ofstream& FileMESSAGE,ostream& oss) {  //CO20200404
    // Set k-points format to LaTeX
    plotoptions.push_attached("KPOINT_FORMAT", "LATEX");
    // Set output format to gnuplot
    plotoptions.push_attached("OUTPUT_FORMAT", "GNUPLOT");

    stringstream out;
    PLOT_BAND(plotoptions, out,FileMESSAGE,oss);
    savePlotGNUPLOT(plotoptions, out);
  }

  void PLOT_BAND(xoption& plotoptions, stringstream& out,ostream& oss) {ofstream FileMESSAGE;return PLOT_BAND(plotoptions,out,FileMESSAGE,oss);}  //CO20200404
  void PLOT_BAND(xoption& plotoptions, stringstream& out,ofstream& FileMESSAGE,ostream& oss) {  //CO20200404
    plotoptions.push_attached("EXTENSION", "band");
    plotoptions.push_attached("OUTPUT_FORMAT", "GNUPLOT");
    // Read files
    string directory = plotoptions.getattachedscheme("DIRECTORY");
    xDOSCAR xdos;
    xdos.GetPropertiesFile(aflowlib::vaspfile2stringstream(directory, "DOSCAR"));
    xEIGENVAL xeigen;
    xeigen.GetPropertiesFile(aflowlib::vaspfile2stringstream(directory, "EIGENVAL"));
    xKPOINTS xkpts;
    xkpts.GetPropertiesFile(aflowlib::vaspfile2stringstream(directory, "KPOINTS"));
    std::stringstream poscar;
    aflowlib::vaspfile2stringstream(directory, "POSCAR", poscar);
    xstructure xstr(poscar);

    plotoptions.push_attached("DEFAULT_TITLE", xeigen.title);
    patchDefaultTitleAFLOWIN(plotoptions);  //ME20200217
    plotoptions.push_attached("LATTICE", getLatticeFromKpointsTitle(xkpts.title));
    setFileName(plotoptions);
    setTitle(plotoptions,FileMESSAGE,oss); //CO20200404

    plotoptions.flag("BANDDOS", false);
    // Set Fermi energy to zero
    if (plotoptions.flag("NOSHIFT")) {
      plotoptions.push_attached("EFERMI", aurostd::utype2string<double>(xdos.Efermi));
    } else {
      shiftEfermiToZero(xeigen, xdos.Efermi);
      plotoptions.push_attached("EFERMI", "0.0");
    }

    // Set Emin and Emax
    setEMinMax(plotoptions, xeigen.energy_min, xeigen.energy_max);

    // Plot
    generateHeader(out, plotoptions, false);
    generateBandPlot(out, xeigen, xkpts, xstr, plotoptions);
  }

  //PLOT_BANDDOS//////////////////////////////////////////////////////////////
  // Plots combined band structure + DOS plots.
  void PLOT_BANDDOS(xoption& plotoptions,ostream& oss) {ofstream FileMESSAGE;return PLOT_BANDDOS(plotoptions,FileMESSAGE,oss);} //CO20200404
  void PLOT_BANDDOS(xoption& plotoptions,ofstream& FileMESSAGE,ostream& oss) { //CO20200404
    // Set k-points format to LaTeX
    plotoptions.push_attached("KPOINT_FORMAT", "LATEX");
    // Set output format to gnuplot
    plotoptions.push_attached("OUTPUT_FORMAT", "GNUPLOT");

    stringstream out;
    PLOT_BANDDOS(plotoptions, out,FileMESSAGE,oss);  //CO20200404
    savePlotGNUPLOT(plotoptions, out);
  }

  void PLOT_BANDDOS(xoption& plotoptions, stringstream& out,ostream& oss) {ofstream FileMESSAGE;return PLOT_BANDDOS(plotoptions,out,FileMESSAGE,oss);}  //CO20200404
  void PLOT_BANDDOS(xoption& plotoptions, stringstream& out,ofstream& FileMESSAGE,ostream& oss) {  //CO20200404
    plotoptions.push_attached("EXTENSION", "banddos");
    // Increase plot size
    plotoptions.pop_attached("PLOT_SIZE");
    plotoptions.push_attached("PLOT_SIZE", BANDDOS_SIZE);

    // Read files
    string directory = plotoptions.getattachedscheme("DIRECTORY");
    xDOSCAR xdos;
    xdos.GetPropertiesFile(aflowlib::vaspfile2stringstream(directory, "DOSCAR"));
    xEIGENVAL xeigen;
    xeigen.GetPropertiesFile(aflowlib::vaspfile2stringstream(directory, "EIGENVAL"));
    xKPOINTS xkpts;
    xkpts.GetPropertiesFile(aflowlib::vaspfile2stringstream(directory, "KPOINTS"));
    std::stringstream poscar;
    aflowlib::vaspfile2stringstream(directory, "POSCAR", poscar);
    xstructure xstr(poscar);

    plotoptions.push_attached("DEFAULT_TITLE", xeigen.title);
    patchDefaultTitleAFLOWIN(plotoptions);  //ME20200217
    plotoptions.push_attached("LATTICE", getLatticeFromKpointsTitle(xkpts.title));
    setFileName(plotoptions);
    setTitle(plotoptions,FileMESSAGE,oss); //CO20200404

    plotoptions.flag("BANDDOS", true);
    // Set Fermi energy to zero
    if (plotoptions.flag("NOSHIFT")) {
      plotoptions.push_attached("EFERMI", aurostd::utype2string<double>(xdos.Efermi));
    } else {
      shiftEfermiToZero(xeigen, xdos.Efermi);
      plotoptions.push_attached("EFERMI", "0.0");
    }

    // Set Emin and Emax
    setEMinMax(plotoptions, xdos.energy_min, xdos.energy_max);

    // Plot
    generateHeader(out, plotoptions, true);
    generateBandPlot(out, xeigen, xkpts, xstr, plotoptions);
    generateDosPlot(out, xdos, plotoptions,FileMESSAGE,oss);  //CO20200404
  }

  // Helper functions --------------------------------------------------------

  //getStructureWithNames/////////////////////////////////////////////////////
  // Extracts the structure from VASP input files, including species names.
  xstructure getStructureWithNames(const xoption& plotoptions,const string& carstring,ostream& oss) {ofstream FileMESSAGE;return getStructureWithNames(plotoptions,FileMESSAGE,carstring,oss);} //CO20200404
  xstructure getStructureWithNames(const xoption& plotoptions,ofstream& FileMESSAGE,const string& carstring,ostream& oss) { //CO20200404
    bool LDEBUG=(FALSE || _DEBUG_PLOTTER_ || XHOST.DEBUG);
    string directory = plotoptions.getattachedscheme("DIRECTORY");
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " directory(input)=" << directory << endl;}
    string poscar_file = "POSCAR"; //CO20200404
    std::stringstream poscar;
    //if (plotoptions.getattachedscheme("EXTENSION") == "phdos")
    if (carstring == "PHON")
    { //CO20200106 - patching for auto-indenting
      poscar_file=DEFAULT_APL_PHPOSCAR_FILE; //CO20200404
      aurostd::efile2stringstream(directory+"/"+poscar_file, poscar);  //CO20200404
    } else {
      if(carstring == "POCC") { //CO20191110 //CO20200404
        //[do NOT load in PARTCAR, we need an example ARUN POSCAR, they all have the same num_each_type]aurostd::efile2stringstream(directory+"/PARTCAR", poscar);
        string arun="";
        arun = plotoptions.getattachedscheme("ARUN_DIRECTORY");  //relative path
        if(LDEBUG){cerr << __AFLOW_FUNC__ << " grabbing POSCAR from " << arun << endl;}
        directory+="/"+arun;
        aurostd::StringSubst(directory,"/RAW/","/LIB/");  //read POSCAR from LIB, not RAW (may not have been processed yet)
      }
      if(LDEBUG){cerr << __AFLOW_FUNC__ << " directory(POSCAR)=" << directory << endl;}
      aflowlib::vaspfile2stringstream(directory + "", poscar_file, poscar);  //CO20200404
    }

    if(LDEBUG){cerr << __AFLOW_FUNC__ << " poscar=" << endl << poscar.str() << endl;}
    xstructure xstr(poscar);
    if (xstr.is_vasp4_poscar_format) {  //PARTCAR has species and it is NOT vasp4 format  //CO20191110
      // No special case for phonons needed because PHPOSCAR is always in VASP5 format
      vector<string> atom_names = KBIN::ExtractAtomicSpecies(directory);
      bool not_all_names_given=atom_names.empty(); //CO20200404
      for (uint i = 0; i < atom_names.size() && not_all_names_given==false; i++) { //CO20200404
        if(atom_names[i].empty()){not_all_names_given=true;}
      }
      if(not_all_names_given){ //CO20200404
        stringstream message;
        message << "Species CANNOT be extracted from dir=" << directory << ", no labels can be supplied";pflow::logger(__AFLOW_FILE__, __AFLOW_FUNC__, message, FileMESSAGE, oss, _LOGGER_WARNING_);  //CO20200404
      }else{ //CO20200404
        if(LDEBUG){cerr << __AFLOW_FUNC__ << " patching names of POSCAR with: " << aurostd::joinWDelimiter(atom_names,", ") << endl;}
        //CO20200404 - ALWAYS use AddAtom() for changing atom properties, indices of species/species_pp/etc will be a mess otherwise
        deque<_atom> atoms=xstr.atoms;
        for(uint i=0;i<atoms.size(); i++) {
          atoms[i].name_is_given=true;
          atoms[i].name=atom_names[i];
        }
        xstr.ReplaceAtoms(atoms); //CO20200404 - patches species too, CRITICAL
      }
    }
    if(LDEBUG){cerr << __AFLOW_FUNC__ << " final xstr=" << endl;cerr << xstr << endl;}  //CO20200404
    return xstr;
  }

  //getLatticeFromKpointsTitle////////////////////////////////////////////////
  // If the KPOINTS file is formatted according to the AFLOW standard, return
  // the lattice of the system. Otherwise, return nothing.
  string getLatticeFromKpointsTitle(const string& title) {
    vector<string> tokens;
    aurostd::string2tokens(title, tokens, " ");
    if (tokens.size() >= 3) {  // AFLOW-formatted KPOINTS titles have at least three columns
      if (tokens[0] == "CUB") return tokens[0];
      else if (tokens[0] == "FCC") return tokens[0];
      else if (tokens[0] == "BCC") return tokens[0];
      else if (tokens[0] == "TET") return tokens[0];
      else if (tokens[0] == "BCT") return tokens[0];
      else if (tokens[0] == "BCT1") return tokens[0];
      else if (tokens[0] == "BCT2") return tokens[0];
      else if (tokens[0] == "ORC") return tokens[0];
      else if (tokens[0] == "ORCF") return tokens[0];
      else if (tokens[0] == "ORCF1") return tokens[0];
      else if (tokens[0] == "ORCF2") return tokens[0];
      else if (tokens[0] == "ORCF3") return tokens[0];
      else if (tokens[0] == "ORCI") return tokens[0];
      else if (tokens[0] == "ORCC") return tokens[0];
      else if (tokens[0] == "HEX") return tokens[0];
      else if (tokens[0] == "RHL") return tokens[0];
      else if (tokens[0] == "RHL1") return tokens[0];
      else if (tokens[0] == "RHL2") return tokens[0];
      else if (tokens[0] == "MCL") return tokens[0];
      else if (tokens[0] == "MCLC") return tokens[0];
      else if (tokens[0] == "MCLC1") return tokens[0];
      else if (tokens[0] == "MCLC2") return tokens[0];
      else if (tokens[0] == "MCLC3") return tokens[0];
      else if (tokens[0] == "MCLC4") return tokens[0];
      else if (tokens[0] == "MCLC5") return tokens[0];
      else if (tokens[0] == "TRI") return tokens[0];
      else if (aurostd::toupper(tokens[0]) == "TRI1A") return tokens[0];
      else if (aurostd::toupper(tokens[0]) == "TRI1B") return tokens[0];
      else if (aurostd::toupper(tokens[0]) == "TRI2A") return tokens[0];
      else if (aurostd::toupper(tokens[0]) == "TRI2B") return tokens[0];
      else return "";
    } else {
      return "";
    }
  }


  //shiftEfermiToZero/////////////////////////////////////////////////////////
  // Shift the energies in an xEIGENVAL object so that the Fermi energy is at
  // zero. This is not necessary for xDOSCAR because it has a separate vector
  // for that purpose.
  void shiftEfermiToZero(xEIGENVAL& xeigen, double Efermi) {
    for (uint k = 0; k < xeigen.number_kpoints; k++) {
      for (uint b = 0; b < xeigen.number_bands; b++) {
        for (uint s = 0; s < xeigen.spin + 1; s++) {
          xeigen.venergy[k][b][s] -= Efermi;
        }
      }
    }
  }

  //setEMinMax////////////////////////////////////////////////////////////////
  // Sets the minimum and maximum energy values for electronic structure
  // plots.
  void setEMinMax(xoption& plotoptions, double Emin, double Emax) {
    if (plotoptions.getattachedscheme("XMIN").empty()) {
      if (plotoptions.flag("NOSHIFT")) {
        plotoptions.push_attached("XMIN", aurostd::utype2string<double>(Emin));
      } else {
        plotoptions.push_attached("XMIN", aurostd::utype2string<double>(DEFAULT_DOS_EMIN));
      }
    }
    if (plotoptions.getattachedscheme("XMAX").empty()) {
      if (plotoptions.flag("NOSHIFT")) {
        plotoptions.push_attached("XMAX", aurostd::utype2string<double>(Emax));
      } else {
        plotoptions.push_attached("XMAX", aurostd::utype2string<double>(DEFAULT_DOS_EMAX));
      }
    }
  }

#define BANDS_DOS_JSON_VERSION 1.0

  /// Converts DOS data from xDOSCAR to the file in json format.
  ///
  /// @param xopts controls the following input options:
  /// "DIRECTORY" -- the directory name
  /// "NOSHIFT" -- if true energy is NOT shifted w.r.t Fermi energy
  ///
  /// @param standalone_json_object controls if the output json file is
  /// a standalone object or a part of another json object (i.e. if opening
  /// and closing curly brackets are present or not)
  aurostd::JSONwriter DOS2JSON(xoption &xopt, const xDOSCAR &xdos, ofstream& FileMESSAGE,
      ostream &oss)
  {
    //ME20211015 - Only set directory when none is given
    string directory = xopt.getattachedscheme("DIRECTORY");
    if (directory.empty()) {
      directory = aurostd::getPWD();
      xopt.push_attached("DIRECTORY", directory);
    }

    xstructure xstr = getStructureWithNames(xopt,FileMESSAGE,xdos.carstring,oss);

    string name = KBIN::ExtractSystemName(directory);
    aurostd::JSONwriter dos_json;
    // TDOS header begin
    dos_json.addNumber("version", BANDS_DOS_JSON_VERSION);
    dos_json.addString("name", name);
    dos_json.addVector("species", xstr.species);
    dos_json.addVector("composition", xstr.num_each_type);
    dos_json.addNumber("Emin", xdos.energy_min);
    dos_json.addNumber("Emax", xdos.energy_max);
    dos_json.addNumber("Efermi", xdos.Efermi);
    dos_json.addNumber("DOS_grid", xdos.number_energies);
    // TDOS header end

    aurostd::JSONwriter tdos_data;
    tdos_data.addBool("energies_shifted", !xopt.flag("NOSHIFT"));
    tdos_data.addVector("energy",xopt.flag("NOSHIFT") ? xdos.venergy : xdos.venergyEf);
    if (aurostd::substring2bool(xdos.carstring, "CAR")){
      tdos_data.addString("x_unit", "EV");
      tdos_data.addString("y_unit", "");
    }
    else if (aurostd::substring2bool(xdos.carstring, "PHON")){
      tdos_data.addString("x_unit", "EV");  //ME20211014 - units are always eV
      tdos_data.addString("y_unit", "");
    }

    if (xdos.spin){
      tdos_data.addVector("spin_majority", xdos.vDOS[0][0][0]);
      // negative for minority spin
      deque<double> minority = xdos.vDOS[0][0][1];
      for (uint i=0; i<minority.size(); i++) minority[i] *= -1;
      tdos_data.addVector("spin_minority", minority);

      tdos_data.addVector("sum_spin_majority", xdos.viDOS[0]);
      // negative for minority spin
      minority = xdos.viDOS[1];
      for (uint i=0; i<minority.size(); i++) minority[i] *= -1;
      tdos_data.addVector("sum_spin_minority", minority);
    }
    else{
      tdos_data.addVector("tDOS", xdos.vDOS[0][0][0]);
      tdos_data.addVector("sum", xdos.viDOS[0]);
    }
    dos_json.addJSON("tDOS_data", tdos_data);

    // projected electronic DOS
    if (xdos.partial && aurostd::substring2bool(xdos.carstring, "CAR")){
      // Here we determine what orbital labels are used depending on the type of data
      // present in DOSCAR.
      // For the systems with LS coupling VASP will print at most 16 orbitals:
      // i_total, i_x, i_y, i_z for i in {s,p,d,f}
      // When there is no LS coupling, there are two possibilities:
      // orbitals are classified by their l,m character and, since VASP prints up to
      // f orbital, there would be 16 orbitals at most.
      // Or orbitals are classified only by l character, giving at most 4 orbitals.
      vector<string> orb_labels(xdos.isLSCOUPLING ? 16 : (xdos.lmResolved ? 16 : 4));
      if (xdos.isLSCOUPLING){
        for (uint i=0; i<16; i++) orb_labels[i] = LS_ORBITALS[i];
      }
      else if (xdos.lmResolved){
        for (uint i=0; i<16; i++) orb_labels[i] = LM_ORBITALS[i];
      }
      else{
        for (uint i=0; i<4; i++) orb_labels[i] = ORBITALS[i];
      }

      // for the "orbitals" key, we want to print a list of orbitals up to the
      // highest present in any atom, not the entire list of known/possible orbitals
      int highest_orbital = 0;
      for (uint i=0; i<xdos.vDOS.size(); i++){
        if (xdos.lmResolved){ // there are 1 s, 3 p, 5 d and so on orbitals
          int orbs_total_num = 1, l = 0;
          while (((int)xdos.vDOS[i].size()-1) - orbs_total_num){// 0's index of vDOS stands for total: it should not be counted
            l++;
            orbs_total_num += 2*l+1;
          }
          highest_orbital = std::max(highest_orbital, l+1);
        }
        else{
          highest_orbital = std::max(highest_orbital, (int)xdos.vDOS[i].size()-1);
        }
      }

      // make an array of orbitals labels from s to whichever is highest
      vector<string> orb_labels_out(highest_orbital);
      for (uint i=0; i<orb_labels_out.size(); i++){
        orb_labels_out[i] = ORBITALS[i];
      }

      // header of partial DOS
      aurostd::JSONwriter pdos_data;
      pdos_data.addVector("orbitals", orb_labels_out);
      pdos_data.addBool("spin_polarized", xdos.spin);
      pdos_data.addBool("energies_shifted", !xopt.flag("NOSHIFT"));
      pdos_data.addVector("energy",xopt.flag("NOSHIFT") ? xdos.venergy : xdos.venergyEf);
      pdos_data.addString("x_unit", "EV");
      pdos_data.addString("y_unit", "");

      // create a mapping of species to the ID of the first representative of
      // each group of the symmetry equivalent atoms, i.e. for SG #12 BaBiO3
      // with Ba : {{0,1}}, Bi: {{2},{3}}, O: {{4,5,6,7}, {8,9}} make the
      // following mapping: Ba -> {0}, Bi -> {2,3}, O -> {4,8}
      pflow::PerformFullSymmetry(xstr);
      vector<vector<int> > map_species_to_iatoms(xstr.species.size());
      for (uint species_id=0; species_id<xstr.species.size(); species_id++){
        for (uint iatom=0; iatom<xstr.iatoms.size(); iatom++){
          if (xstr.atoms[xstr.iatoms[iatom][0]].type == (int)species_id)
            map_species_to_iatoms[species_id].push_back(xstr.iatoms[iatom][0]);
        }
      }

      // write atom-projected DOS for each unique atom
      string label = "";
      for (uint species_id=0; species_id<xstr.species.size(); species_id++){
        aurostd::JSONwriter species_json;
        for (uint iatom=0; iatom<map_species_to_iatoms[species_id].size(); iatom++){
          int atom_id = map_species_to_iatoms[species_id][iatom];
          aurostd::JSONwriter atom_json;
          atom_id++; // in vDOS atoms are indexed starting from 1

          if (xdos.isLSCOUPLING){
            for (uint orb=1; orb<xdos.vDOS[atom_id].size(); orb++){
              // there are 4 spin channels
              for (uint spin=0; spin<xdos.vDOS[atom_id][orb].size(); spin++){
                label = orb_labels[4*(orb-1) + spin];
                atom_json.addVector(label, xdos.vDOS[atom_id][orb][spin]);
              }
            }

            atom_json.addVector("total", xdos.vDOS[atom_id][0][0]);
          }
          else{
            for (uint spin=0; spin<=xdos.spin; spin++){
              if (xdos.vDOS[atom_id].size()){
                label = "total" + (xdos.spin ? "_"+ SPIN_LABEL[spin] : "");
                if (spin){
                  // negative for minority spin
                  deque<double> minority = xdos.vDOS[atom_id][0][spin];
                  for (uint i=0; i<minority.size(); i++) minority[i] *= -1;
                  atom_json.addVector(label, minority);
                }
                else{
                  atom_json.addVector(label, xdos.vDOS[atom_id][0][spin]);
                }
              }

              for (uint orb=1; orb<xdos.vDOS[atom_id].size(); orb++){
                label = orb_labels[orb-1] + (xdos.spin ? "_"+ SPIN_LABEL[spin] : "");
                if (spin){
                  // negative for minority spin
                  deque<double> minority = xdos.vDOS[atom_id][orb][spin];
                  for (uint i=0; i<minority.size(); i++) minority[i] *= -1;
                  atom_json.addVector(label, minority);
                }
                else{
                  atom_json.addVector(label, xdos.vDOS[atom_id][orb][spin]);
                }
              }
            }
          }
          species_json.addJSON(aurostd::utype2string<int>(atom_id-1), atom_json);
        }
        pdos_data.addJSON(xstr.species[species_id], species_json); 
      }

      // write the sum of DOS contributions for each orbital (s, p, d and f)
      // for all atoms
      if (xdos.lmResolved){
        deque<deque<double> > orb_dos(xdos.spin+1, deque<double> (xdos.number_energies));
        for (int l=0, norb=1; norb<=(int)xdos.vDOS[0].size()-1; l++, norb += (2*l+1)){
          for (uint spin=0; spin<=xdos.spin; spin++){
            for (uint en=0; en<orb_dos[spin].size(); en++) orb_dos[spin][en] = 0.0;
          }

          for (uint spin=0; spin<=xdos.spin; spin++){
            // orbitals are grouped by 2*l+1 manifolds: loop to sum each group
            for (int i=0; i<2*l+1; i++){
              for (uint en=0; en<orb_dos[spin].size(); en++){
                orb_dos[spin][en] += (xdos.spin ? -1 : 1)*xdos.vDOS[0][norb-i][spin][en]; // negative for minority spin
              }
            }
            label = "sum_" + ORBITALS[l];
            label += xdos.spin ? "_"+SPIN_LABEL[spin] : "";
            pdos_data.addVector(label, orb_dos[spin]);
          }
        }
      }
      else{
        for (uint orb=1; orb<xdos.vDOS[0].size(); orb++){
          for (uint spin=0; spin<=xdos.spin; spin++){
            if (xdos.isLSCOUPLING){
              label = "sum_"+ORBITALS[orb-1];
            }
            else{
              label = "sum_"+orb_labels[orb-1];
              label += xdos.spin ? "_"+SPIN_LABEL[spin] : "";
            }
            if (spin){
              // negative for minority spin
              deque<double> minority = xdos.vDOS[0][orb][spin];
              for (uint i=0; i<minority.size(); i++) minority[i] *= -1;
              pdos_data.addVector(label, minority);
            }
            else{
              pdos_data.addVector(label, xdos.vDOS[0][orb][spin]);
            }
          }
        }
      }
      dos_json.addJSON("pDOS_data", pdos_data);
    }

    // projected phonon DOS
    if (xdos.partial && aurostd::substring2bool(xdos.carstring, "PHON")){
      // header of partial DOS
      aurostd::JSONwriter pdos_data;
      deque<string> projections;
      if (xdos.vDOS.size()>=2){
        for (uint j=1; j<xdos.vDOS[1].size(); j++){
          projections.push_back("projection_" + aurostd::utype2string(j));
        }
      }
      pdos_data.addVector("orbitals", projections);
      pdos_data.addBool("spin_polarized", xdos.spin);
      pdos_data.addBool("energies_shifted", !xopt.flag("NOSHIFT"));
      pdos_data.addVector("energy", xopt.flag("NOSHIFT") ? xdos.venergy : xdos.venergyEf);
      pdos_data.addString("x_unit", "EV");  //ME20211014 - units are always eV
      pdos_data.addString("y_unit", "");

      // create a mapping of species to the id of the first representative of
      // each group of the symmetry equivalent atoms, i.e. for SG #12 BaBiO3
      // with Ba : {{0,1}}, Bi: {{2},{3}}, O: {{4,5,6,7}, {8,9}} make the
      // following mapping: Ba -> {0}, Bi -> {2,3}, O -> {4,8}
      pflow::PerformFullSymmetry(xstr);
      vector<vector<int> > map_species_to_iatoms(xstr.species.size());
      for (uint species_id=0; species_id<xstr.species.size(); species_id++){
        for (uint iatom=0; iatom<xstr.iatoms.size(); iatom++){
          if (xstr.atoms[xstr.iatoms[iatom][0]].type == (int)species_id)
            map_species_to_iatoms[species_id].push_back(xstr.iatoms[iatom][0]);
        }
      }

      // write atom-projected DOS for each unique atom
      string label = "";
      for (uint species_id=0; species_id<xstr.species.size(); species_id++){
        aurostd::JSONwriter species_json;
        for (uint iatom=0; iatom<map_species_to_iatoms[species_id].size(); iatom++){
          int atom_id = map_species_to_iatoms[species_id][iatom];
          aurostd::JSONwriter atom_json;
          atom_id++; // in vDOS atoms are indexed starting from 1

          for (uint p=1; p<xdos.vDOS[atom_id].size(); p++){
            atom_json.addVector(projections[p-1], xdos.vDOS[atom_id][p][0]);
          }
          species_json.addJSON(aurostd::utype2string<int>(atom_id-1), atom_json);
        }
        pdos_data.addJSON(xstr.species[species_id], species_json);
      }

      dos_json.addJSON("pDOS_data", pdos_data);
    }

    return dos_json;
  }

  /// Converts band structure data from xEIGENVAL and xKPOINTS to the file in
  /// json format.
  ///
  /// @plotoptions specifies the following options:
  /// "DIRECTORY" -- the directory name
  /// "NOSHIFT" -- if true energy is NOT shifted w.r.t Fermi energy
  /// "EFERMI" -- the value of Fermi energy
  aurostd::JSONwriter bands2JSON(const xEIGENVAL &xeigen, const xKPOINTS &xkpts,
      const vector<double> &distances, const vector<double> &segment_points,
      const xoption& plotoptions)
  {
    string directory = plotoptions.getattachedscheme("DIRECTORY");
    if (directory.empty()) directory = aurostd::getPWD();
    string name = KBIN::ExtractSystemName(directory);

    // get lattice type as written in KPOINTS.bands file
    vector<string> tokens;
    aurostd::string2tokens(xkpts.title, tokens, " ");
    string LattName = tokens[0];

    // write header
    aurostd::JSONwriter json;
    json.addString("title", name+" ("+LattName+")");
    json.addNumber("n_kpoints", xeigen.number_kpoints);
    json.addNumber("n_bands", xeigen.number_bands);
    if(aurostd::substring2bool(xeigen.carstring, "CAR")){
      json.addString("x_unit", "");
      json.addString("y_unit", "EV");
    }
    else if(aurostd::substring2bool(xeigen.carstring, "PHON")){
      json.addString("x_unit", "");
      json.addString("y_unit", "EV");  //ME20211014 - units are always eV
    }

    static const uint num = 4;
    string tags[num] = {"kpoint_labels", "kpoint_labels_latex", "kpoint_labels_gnuplot", "kpoint_labels_html"};
    string formats[num] = { "", "LATEX", "GNUPLOT", "HTML"};

    // we need clean labels: to be consistent with E.G. format of the JSON file
    // they will be converted into all possible formats
    uint nsegments = xkpts.vpath.size()/2;
    vector<string> labels_formated(nsegments+1);
    for (uint f=0; f<num; f++){
      labels_formated[0] = convertKPointLabel(xkpts.vpath[0], formats[f]);
      for (uint i = 2; i < 2 * nsegments; i += 2) {
        labels_formated[i/2] = convertKPointLabel(xkpts.vpath[i - 1], formats[f]);
        if (xkpts.vpath[i-1] != xkpts.vpath[i]) {
          labels_formated[i/2] += ("GNUPLOT" == formats[f] ? " |" : "|"); //ME, should we add $|$ here too?
          labels_formated[i/2] += convertKPointLabel(xkpts.vpath[i], formats[f]);
        }
      }
      labels_formated.back() = convertKPointLabel(xkpts.vpath.back(), formats[f]);

      // escape backslash symbols in all labels
      for (uint j=0; j<labels_formated.size();j++){
        labels_formated[j] = aurostd::StringSubst(labels_formated[j], "\\", "\\\\");
      }

      json.addVector(tags[f], labels_formated);
    }

    json.addVector("kpoint_positions", segment_points);
    json.addBool("energies_shifted", !plotoptions.flag("NOSHIFT"));

    // write bands data
    string bandslabel = "";
    double Efermi = aurostd::string2utype<double>(plotoptions.getattachedscheme("EFERMI"));
    for (uint s=0; s<=xeigen.spin; s++){
      bandslabel = "bands_data";
      bandslabel += xeigen.spin ? (s ? "_minority" : "_majority") : "";

      vector<vector<double> > bandsdata(distances.size(),
          vector<double> (xeigen.number_bands + 1)); // +1 because the first element is the the distance, the others are energies for each band

      for (uint i=0; i<distances.size(); i++){
        bandsdata[i][0] = distances[i];
        if (plotoptions.flag("NOSHIFT")){
          for (uint band=0; band<xeigen.number_bands; band++){
            bandsdata[i][band+1] = xeigen.venergy[i][band][s];
          }
        }
        else{
          for (uint band=0; band<xeigen.number_bands; band++){
            bandsdata[i][band+1] = xeigen.venergy[i][band][s] - Efermi;
          }
        }
      }

      json.addMatrix(bandslabel, bandsdata);
    }

    return json;
  }

  /// Converts DOS and BANDS data from xDOSCAR, xEIGENVAL and xKPOINTS files
  /// to JSON object
  aurostd::JSONwriter bandsDOS2JSON(const xDOSCAR &xdos, const xEIGENVAL &xeigen,
      const xKPOINTS &xkpts, xoption &xopt, ofstream &FileMESSAGE, ostream &oss)
  {
    // get DOS part of JSON
    aurostd::JSONwriter json = DOS2JSON(xopt, xdos, FileMESSAGE, oss);

    // get BANDS part of JSON
    xstructure xstr = getStructureWithNames(xopt,FileMESSAGE,xdos.carstring,oss);
    xopt.pop_attached("OUTPUT_FORMAT");
    xopt.push_attached("OUTPUT_FORMAT","JSON");
    xopt.pop_attached("EFERMI");
    xopt.push_attached("EFERMI", aurostd::utype2string<double>(xdos.Efermi));
    stringstream json_stream;
    generateBandPlot(json_stream, xeigen, xkpts, xstr, xopt);
    string bands = json_stream.str();
    if (bands.size() > 2){
      bands = bands.substr(1, bands.size() - 2); // remove wrapping curly brackets
    }

    json.mergeRawJSON(bands); //DX20210304 - addRaw -> mergeRawJSON
    return json;
  }

  // DOS ---------------------------------------------------------------------

  //dosDataAvailable//////////////////////////////////////////////////////////
  // ME20200305
  // Checks if the DOS contains data for the requested pdos
  bool dosDataAvailable(const deque<deque<deque<deque<double> > > >& vdos, int pdos) {
    if (pdos < 0) return false;
    if (pdos + 1 > (int) vdos.size()) return false;
    if (vdos[pdos].size() == 0) return false;
    return true;
  }

  //generateDosPlot///////////////////////////////////////////////////////////
  // Generates the data for a DOS plot. 
  // ME20200305 - added DOS data checking
  void generateDosPlot(stringstream& out,const xDOSCAR& xdos,xoption& plotoptions,ostream& oss) {ofstream FileMESSAGE;return generateDosPlot(out,xdos,plotoptions,FileMESSAGE,oss);} //CO20200404
  void generateDosPlot(stringstream& out,const xDOSCAR& xdos,xoption& plotoptions,ofstream& FileMESSAGE,ostream& oss) {  //CO20200404
    bool LDEBUG=(FALSE || _DEBUG_PLOTTER_ || XHOST.DEBUG); 
    deque<deque<deque<double> > > dos;
    string dataset=plotoptions.getattachedscheme("DATASET");
    int pdos = aurostd::string2utype<int>(dataset);
    vector<string> labels;
    labels.push_back("total");  // There is always a total DOS
    if(pdos>0){labels.front()+=" "+plotoptions.getattachedscheme("DATALABEL");} //CO20191010
    string projection = plotoptions.getattachedscheme("PROJECTION");
    string datatype=plotoptions.getattachedscheme("DATATYPE");  //CO20191010 - which table to look at, "atoms-projected" by default
    if(LDEBUG){
      cerr << __AFLOW_FUNC__ << " dataset=" << dataset << endl;
      cerr << __AFLOW_FUNC__ << " pdos=" << pdos << endl;
      cerr << __AFLOW_FUNC__ << " projection=" << projection << endl;
      cerr << __AFLOW_FUNC__ << " datatype=" << datatype << endl;
    }
    if (projection == "ORBITALS") {
      // If the DOSCAR is lm-resolved, the orbital projection is the sum of all individual
      // orbitals with the same quantum number
      int norbitals=0;
      if(datatype=="SPECIES"){  //CO20191010 - plot "species-projected"
        xstructure xstr = getStructureWithNames(plotoptions,FileMESSAGE,xdos.carstring,oss);  //CO20200404
        deque<deque<deque<deque<double> > > > vDOS_species=xdos.GetVDOSSpecies(xstr);
        if (!dosDataAvailable(vDOS_species, pdos)) {
          string message = "DOS data not available for pdos = " + aurostd::utype2string<int>(pdos);
          throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__, message, _RUNTIME_ERROR_);
        }
        norbitals=vDOS_species.front().size();
        dos=vDOS_species[pdos];
      }else{  //CO20191010 - plot "atoms-projected"
        if (!dosDataAvailable(xdos.vDOS, pdos)) {
          string message = "DOS data not available for pdos = " + aurostd::utype2string<int>(pdos);
          throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__, message, _RUNTIME_ERROR_);
        }
        if (xdos.lmResolved){norbitals = (int) std::sqrt(xdos.vDOS[pdos].size());}  // size is either 17 or 10
        else{norbitals = (int) xdos.vDOS[pdos].size() - 1;}
        if(LDEBUG){cerr << __AFLOW_FUNC__ << " norbitals=" << norbitals << endl;}
        if (xdos.lmResolved) {
          // Total DOS and s-orbitals
          dos.push_back(xdos.vDOS[pdos][0]);
          dos.push_back(xdos.vDOS[pdos][1]);
          // p-, d-, and maybe f-orbitals
          deque<deque<deque<double> > > dospart(norbitals - 1, deque<deque<double> >(xdos.spin + 1, deque<double>(xdos.number_energies)));
          for (uint e = 0; e < xdos.number_energies; e++) {
            for (int d = 0; d < norbitals - 1; d++) {
              for (int i = d * (2 + d) + 2; i < d * (d + 4) + 5; i++) {
                for (uint s = 0; s < xdos.spin + 1; s++) {
                  dospart[d][s][e] += xdos.vDOS[pdos][i][s][e];
                }
              }
            }
          }
          for (int d = 0; d < norbitals - 1; d++) dos.push_back(dospart[d]);
        } else {
          dos = xdos.vDOS[pdos];
        }
      }
      if(LDEBUG) {cerr << __AFLOW_FUNC__ << " norbitals=" << norbitals << endl;}
      //CO20191010 - do labels last
      for (int i = 0; i < norbitals; i++) {
        labels.push_back("$" + ORBITALS[i] + "$");
      }
    } else if (projection == "LM") {
      // Safety check
      if (!xdos.lmResolved) {
        string message = "Projection scheme LM chosen, but DOSCAR is not lm-resolved.";
        throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__, message, _RUNTIME_ERROR_);
      }
      if (!dosDataAvailable(xdos.vDOS, pdos)) {
        string message = "DOS data not available for pdos = " + aurostd::utype2string<int>(pdos);
        throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__, message, _RUNTIME_ERROR_);
      }
      for (uint i = 1; i < xdos.vDOS[pdos].size(); i++) {
        labels.push_back("$" + LM_ORBITALS_LATEX[i-1] + "$");
      }
      dos = xdos.vDOS[pdos];
    } else if (projection == "ATOMS") { //CO20191004 - "ATOMS" is really "IATOMS"
      if (!dosDataAvailable(xdos.vDOS, 0)) {
        string message = "Total DOS not available";
        throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__, message, _RUNTIME_ERROR_);
      }
      dos.push_back(xdos.vDOS[0][0]);
      xstructure xstr = getStructureWithNames(plotoptions,FileMESSAGE,xdos.carstring,oss);  //CO20200404
      if (plotoptions.flag("PLOT_ALL_ATOMS")) { //CO20191010 - special mode - ALL atoms
        if (xdos.vDOS.size() > 1) {
          for (uint i = 0; i < xstr.atoms.size(); i++) {
            if (!dosDataAvailable(xdos.vDOS, i + 1)) {
              string message = "DOS data not available for atom " + aurostd::utype2string<uint>(i + 1);
              throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__, message, _RUNTIME_ERROR_);
            }
            dos.push_back(xdos.vDOS[i + 1][0]);
            xstr.atoms[i].CleanName();
            labels.push_back(xstr.atoms[i].cleanname + "(" + aurostd::utype2string<uint>(i + 1) + ")");
          }
        }
      } else {
        if (pdos == 0) {
          if (xdos.vDOS.size() > 1) {
            pflow::PerformFullSymmetry(xstr);
            int iat = 0;
            for (uint i = 0; i < xstr.iatoms.size(); i++) {
              iat = xstr.iatoms[i][0];
              if (!dosDataAvailable(xdos.vDOS, iat + 1)) {
                string message = "DOS data not available for iatom " + aurostd::utype2string<int>(iat + 1);
                throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__, message, _RUNTIME_ERROR_);
              }
              dos.push_back(xdos.vDOS[iat + 1][0]);
              xstr.atoms[iat].CleanName();
              labels.push_back(xstr.atoms[iat].cleanname + "(" + aurostd::utype2string<int>(iat + 1) + ")");
            }
          }
        } else { // In case someone uses pdos option and projection=atoms
          if (!dosDataAvailable(xdos.vDOS, pdos)) {
            string message = "DOS data not available for pdos = " + aurostd::utype2string<int>(pdos);
            throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__, message, _RUNTIME_ERROR_);
          }
          xstr.atoms[pdos - 1].CleanName();
          labels.push_back(xstr.atoms[pdos - 1].cleanname + "(" + aurostd::utype2string<int>(pdos) + ")");
          dos.push_back(xdos.vDOS[pdos][0]);
        }
      }
    } else if (projection == "SPECIES") {  //CO20191110
      xstructure xstr = getStructureWithNames(plotoptions,FileMESSAGE,xdos.carstring,oss);  //CO20200404
      if(LDEBUG){cerr << __AFLOW_FUNC__ << " xstr=" << endl;cerr << xstr << endl;}  //CO20200404
      deque<deque<deque<deque<double> > > > vDOS_species=xdos.GetVDOSSpecies(xstr);
      if (pdos == 0) {
        dos.push_back(vDOS_species[0][0]);
        if (vDOS_species.size() > 1) {
          for (uint i = 0; i < xstr.species.size(); i++) {
            if (!dosDataAvailable(vDOS_species, i + 1)) {
              string message = "DOS data not available for species " + aurostd::utype2string<uint>(i + 1);
              throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__, message, _RUNTIME_ERROR_);
            }
            dos.push_back(vDOS_species[i+1][0]);
            labels.push_back(xstr.species[i]);
          }
        }
      } else { // In case someone uses pdos option and projection=atoms
        //labels.push_back(xstr.species[pdos-1]);
        //dos.push_back(vDOS_species[pdos][0]);
        labels.push_back(xstr.atoms[pdos-1].name);//AS20201028
        int dos_index = xstr.atoms[pdos - 1].type + 1;
        if (!dosDataAvailable(vDOS_species, dos_index)) {
          string message = "DOS data not available for dos_index = " + aurostd::utype2string<int>(dos_index);
          throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__, message, _RUNTIME_ERROR_);
        }
        dos.push_back(vDOS_species[dos_index][0]);//AS20201028
      }
    } else if (projection == "NONE") {  // Total DOS only without projections
      if (!dosDataAvailable(xdos.vDOS, 0)) {
        string message = "Total DOS not available";
        throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__, message, _RUNTIME_ERROR_);
      }
      dos.push_back(xdos.vDOS[0][0]);
    } else {
      string message = "Unknown projection scheme " + projection + ".";
      throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__, message, _INPUT_ILLEGAL_);
    }
    if(plotoptions.flag("LEGEND_HORIZONTAL")==false && labels.size()>4){  //CO20211227, avoid overlap between legend and DOS
      plotoptions.flag("LEGEND_HORIZONTAL",true);
      plotoptions.push_attached("LEGEND_MAXCOLS","5");
    }
    string outformat = plotoptions.getattachedscheme("OUTPUT_FORMAT");
    if (outformat == "GNUPLOT") {
      if (plotoptions.flag("NOSHIFT")) {
        generateDosPlotGNUPLOT(out, xdos, xdos.venergy, dos, labels, plotoptions);
      } else {
        generateDosPlotGNUPLOT(out, xdos, xdos.venergyEf, dos, labels, plotoptions);
      }
    }
  }

  // Bands -------------------------------------------------------------------

  //generateBandPlot//////////////////////////////////////////////////////////
  // Generates the data for a band structure plot.
  void generateBandPlot(stringstream& out, const xEIGENVAL& xeigen, const xKPOINTS& xkpts,
      const xstructure& xstr, const xoption& plotoptions) {
    // Create segments
    uint nsegments = xkpts.vpath.size()/2;
    // Make sure that the number of k-points is consistent with EIGENVAL
    if (xeigen.number_kpoints != nsegments * xkpts.path_grid) {
      string message = "Number of k-points in EIGENVAL and KPOINTS files do not match.";
      throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__, message, _RUNTIME_ERROR_);
    }

    // Labels
    string format = plotoptions.getattachedscheme("KPOINT_FORMAT");
    vector<string> labels(nsegments + 1);
    labels[0] = convertKPointLabel(xkpts.vpath[0], format);
    for (uint i = 2; i < 2 * nsegments; i += 2) {
      labels[i/2] = convertKPointLabel(xkpts.vpath[i - 1], format);
      if (xkpts.vpath[i-1] != xkpts.vpath[i]) {
        labels[i/2] += "$|$" + convertKPointLabel(xkpts.vpath[i], format);  //CO20210701 - | by itself can show up as hyphen, also consider $\\vert$
      }
    }
    labels.back() = convertKPointLabel(xkpts.vpath.back(), format);

    // k-points
    xmatrix<double> f2c = trasp(ReciprocalLattice(xstr.lattice));
    double total_length = 0.0;
    vector<double> segment_points(nsegments + 1, 0.0);
    double dk;
    for (uint i = 0; i < nsegments; i++) {
      dk = aurostd::modulus(f2c * (xkpts.vkpoints[2*i] - xkpts.vkpoints[2*i+1]));
      segment_points[i+1] = segment_points[i] + dk;
      total_length += dk;
    }
    for (uint i = 0; i < nsegments + 1; i++) segment_points[i] /= total_length;

    // Project to path lengths
    vector<double> distances(xeigen.number_kpoints, 0.0);
    int k = 0;
    double ds;
    for (uint i = 0; i < nsegments; i++) {
      ds = (segment_points[i+1] - segment_points[i])/(xkpts.path_grid - 1);
      distances[k] = segment_points[i];
      k++;
      for (int j = 1; j < xkpts.path_grid; j++) {
        distances[k] = distances[k -1] + ds;
        k++;
      }
    }
    string outformat = plotoptions.getattachedscheme("OUTPUT_FORMAT");
    if (outformat == "GNUPLOT") {
      generateBandPlotGNUPLOT(out, xeigen, distances, segment_points, labels, plotoptions);
    }
    else if (outformat == "JSON"){
      out << bands2JSON(xeigen, xkpts, distances, segment_points, plotoptions).toString();
    }
  }

  //convertKPointLabel////////////////////////////////////////////////////////
  // Converts a raw k-point string a k-point label into the desired format.
  string convertKPointLabel(const string& kpoint, const string& format) {
    vector<string> parts;
    string formatted_label;
    aurostd::string2tokens(kpoint, parts, "_");
    if (parts.size() > 2) return kpoint;
    formatted_label = convertKPointLetter(parts[0], format);
    if (parts.size() == 2) {
      if (format == "LATEX") {
        formatted_label += "$_{" + parts[1] + "}$";
        aurostd::RemoveSubStringFirst(kpoint, "$$");
      } else if (format == "HTML") {
        formatted_label += "<sub>" + parts[1] + "</sub>";
      }
      //AS20201103 BEGIN
      else if (format == "GNUPLOT"){
        formatted_label += "_{/"+DEFAULT_GNUPLOT_EPS_FONT+" "+parts[1]+"}";
      }
      else{// return the second part even if no format is specified
        formatted_label += "_" + parts[1];
      }
      //AS20201103 END
    }
    return formatted_label;
  }

  //convertKPointLetter///////////////////////////////////////////////////////
  // Converts a raw k-point letter string into the desired format.
  string convertKPointLetter(string letter, const string& format) {
    if (format == "LATEX") {
      if (aurostd::substring2bool(letter, "\\") && !aurostd::substring2bool(letter, "\\Gamma")) {
        letter = "$\\mathit{" + letter + "}$";
      } else {
        letter = "$" + letter + "$";
      }
    } else if (format == "HTML") {
      if (aurostd::substring2bool(letter, "\\")) {
        aurostd::StringSubst(letter, "\\", "&");
        letter += ";";
      }
    }
    //AS20201103 BEGIN
    else if (format == "GNUPLOT"){
      if (aurostd::substring2bool(letter, "Gamma")){
        letter = "{/Symbol G}";
      }
      else if (aurostd::substring2bool(letter, "Sigma")){
        letter = "{/Symbol S}";
      }
      else{
        letter = "{/"+DEFAULT_GNUPLOT_EPS_FONT_ITALICS+" "+letter+"}";
      }
    }
    //AS20201103 END
    return letter;
  }

  // Gnuplot -----------------------------------------------------------------

  //generateDosPlotGNUPLOT////////////////////////////////////////////////////
  // Generates the gnuplot script for DOS plots.
  void generateDosPlotGNUPLOT(stringstream& out, const xDOSCAR& xdos, const deque<double>& energies,
      const deque<deque<deque<double> > >& dos, const vector<string>& labels,
      const xoption& plotoptions) {
    bool LDEBUG=(FALSE || _DEBUG_PLOTTER_ || XHOST.DEBUG); 
    // Initialize variables
    double Efermi = aurostd::string2utype<double>(plotoptions.getattachedscheme("EFERMI"));
    double Emin = aurostd::string2utype<double>(plotoptions.getattachedscheme("XMIN"));
    double Emax = aurostd::string2utype<double>(plotoptions.getattachedscheme("XMAX"));
    bool banddos = plotoptions.flag("BANDDOS");
    bool swap = (banddos || plotoptions.flag("SWAP_AXES"));
    uint ndos = dos.size();

    double dosmax = getDosLimits(plotoptions, xdos, dos, energies);

    if(LDEBUG){
      cerr << __AFLOW_FUNC__ << " dosmax=" << dosmax << endl;
      cerr << __AFLOW_FUNC__ << " xdos.spin=" << xdos.spin << endl;
    }

    if(aurostd::isequal(dosmax,0.0)){dosmax=1.0;} //CO+ME20210729 - issues with LIB0 pDOS plots, might need to rerun with NEDOS=5000, EMIN=-45, EMAX=30
    string maxdos = aurostd::utype2string<double>(dosmax);
    string mindos="";
    if (xdos.spin == 0) mindos = "0";
    else mindos = "-" + maxdos;

    string unit = plotoptions.getattachedscheme("UNIT");
    if (unit.empty()) unit = "EV";
    string energyLabel;
    if (aurostd::substring2bool(unit, "EV")) energyLabel = "energy";
    else energyLabel = "frequency";
    unit = getFormattedUnit(unit);

    out << std::endl << "# DOS plot" << std::endl;

    // Create data block
    out << std::endl << "$dos_data << EOD" << std::endl;
    for (uint e = 0; e < xdos.number_energies; e++) {
      out << "  " << energies[e];
      for (uint d = 0; d < ndos; d++) {
        out << " " << dos[d][0][e];
      }
      if (xdos.spin == 1) {
        for (uint d = 0; d < ndos; d++) {
          out << " " << -dos[d][1][e];
        }
      }
      out << std::endl;
    }
    out << "EOD" << std::endl << std::endl;

    // Margins
    out << "# Margins" << std::endl;
    if (banddos) {
      out << "set lmargin at screen 0.73" << std::endl;
      out << "set rmargin at screen 0.98" << std::endl;
      out << "set tmargin at screen 0.9" << std::endl;
      out << "set bmargin at screen 0.12" << std::endl;
    } else {
      out << "set tmargin at screen 0.9" << std::endl;
      out << "set bmargin at screen 0.2" << std::endl;
    }

    // Key
    out << std::endl << "# Key" << std::endl;
    if (dos.size() == xdos.spin + 1) { // no need for key when only total DOS is plotted, NO-SPIN: xdos.spin==0, SPIN: xdos.spin==1  //CO20200404
      out << "unset key" << std::endl;
    } else {
      if(plotoptions.flag("LEGEND_HORIZONTAL")){  //CO20200404
        int maxcols=3;
        string maxcols_str=plotoptions.getattachedscheme("LEGEND_MAXCOLS");
        if(!maxcols_str.empty()){maxcols=aurostd::string2utype<int>(maxcols_str);}
        out << "set key horizontal maxcols " << maxcols << std::endl;
      }
      out << "set key samplen 2.5" << std::endl;  // Shorter lines to fit key into image
    }

    // Axes
    out << std::endl << "# Axes" << std::endl;
    if (banddos) {
      out << "unset xtics" << std::endl;
      out << "unset xrange" << std::endl;
    }

    out << "set " << (swap?"x":"y") << "tics " << (dosmax/(2 * (2 - xdos.spin))) << std::endl;
    if(banddos){out << "set ytics format \"\"" << std::endl;}  //CO+ME20210729
    //[CO+ME20210729 - breaks for LIB0 which has swap==false and banddos]out << "set ytics" << (banddos?" format \"\"":"") << std::endl;
    out << "set tic scale 0" << std::endl;
    out << "set " << (swap?"y":"x") << "range [" << Emin << ":" << Emax << "]" << std::endl;
    out << "set " << (swap?"x":"y") << "range [" << mindos << ":" << maxdos << "]" << std::endl;
    string normalization=aurostd::tolower(plotoptions.getattachedscheme("NORMALIZATION"));
    if (banddos) {
      out << "unset ylabel" << std::endl;
      out << "set title 'DOS (states/" << unit << (!normalization.empty()?"/"+normalization:"") << ")' offset 0,-0.7" << std::endl;
    } else {
      out << "set " << (swap?"y":"x") << "label '" << energyLabel << " (" << unit << ")' offset graph 0.00" << std::endl;
      out << "set " << (swap?"x":"y") << "label 'DOS (states/" << unit << (!normalization.empty()?"/"+normalization:"") << ")' offset graph 0.00" << std::endl;
    }

    // Fermi level
    if (Efermi > Emin) {
      out << std::endl << "# Fermi level" << std::endl;
      if (swap) {
        out << "set arrow from " << mindos << "," << Efermi << " to " << maxdos << "," << Efermi;
      } else {
        out << "set arrow from " << Efermi << ", graph 0 to " << Efermi << ", graph 1";
      }
      out << " nohead lt 1 lc rgb '" << EFERMI_COLOR << "' lw 3" << std::endl;
    }

    // Plot data
    out << std::endl << "# Data" << std::endl;
    out << "plot ";
    int xcol, ycol;
    if (swap) ycol = 1;
    else xcol = 1;

    // Majority spin
    for (uint i = 0; i < ndos; i++) {
      if (swap) xcol = i + 2;
      else ycol = i + 2;
      if (i > 0) out << "     ";
      out << "'$dos_data' u " << xcol << ":" << ycol << " w l lt -1 "
        << "lc rgb '" << ESTRUCTURE_COLORS[i % ESTRUCTURE_NCOLORS] << "' lw 2 title '" << labels[i] << "'"
        << (((xdos.spin == 1) || (i < ndos - 1))?",\\":"") << std::endl;
    }
    // Minority spin
    if (xdos.spin == 1) {
      for (uint i = 0; i < ndos; i++) {
        if (swap) xcol = i + ndos + 2;
        else ycol = i + ndos + 2;
        out << "      ";
        out << "'$dos_data' u "  << xcol << ":" << ycol << " w l lt -1 "
          << "lc rgb '" << ESTRUCTURE_COLORS[i % ESTRUCTURE_NCOLORS] << "' lw 2 notitle"  // No title to prevent redundant key entries
          << ((i < ndos - 1)?",\\":"") << std::endl;
      }
    }

    if(LDEBUG){cerr << __AFLOW_FUNC__ << " gnuplot file:" << endl << out.str() << endl;}
  }

  //getDosLimits//////////////////////////////////////////////////////////////
  // Determines the maximum DOS in the plot and sets the limit so that the
  // tics give "nice" numbers. Each plot has four tics, i.e. spin-polarized
  // DOS have two tics per side. This prevents negative numbers from
  // overlapping. This function is fairly primitive but should work for most
  // plots.
  double getDosLimits(const xoption& plotoptions, const xDOSCAR& xdos,
      const deque<deque<deque<double> > >& dos, const deque<double>& energies) {
    bool LDEBUG=(FALSE || _DEBUG_PLOTTER_ || XHOST.DEBUG); 

    double Emin = aurostd::string2utype<double>(plotoptions.getattachedscheme("XMIN"));
    double Emax = aurostd::string2utype<double>(plotoptions.getattachedscheme("XMAX"));
    string dosscale = plotoptions.getattachedscheme("YSCALE");
    uint ndos = dos.size(); //orbital

    if(LDEBUG){
      cerr << __AFLOW_FUNC__ << " Emin=" << Emin << endl;
      cerr << __AFLOW_FUNC__ << " Emax=" << Emax << endl;
    }

    double dosmax = 0.0;
    // First, determine the maximum displayed DOS
    // in the displayed energy range
    // FINDS INDICES OF EMIN AND EMAX
    int e1 = -1, e2 = -1;
    for (uint e = 0; e < xdos.number_energies; e++) {
      if ((e1 < 0) && (energies[e] >= Emin)) e1 = (int) e;
      if ((e2 < 0) && (energies[e] > Emax)) e2 = (int) e;
      if ((e1 > 0) && (e2 > 0)) break;
    }
    if (e1 < 0) e1 = 0;  // Emin not found
    if (e2 < 0) e2 = (int) xdos.number_energies; // Emax not found

    //FINDS DOSMAX BETWEEN EMIN AND EMAX
    for (uint d = 0; d < ndos; d++) { //orbital
      for (uint s = 0; s < xdos.spin + 1; s++) {  //spin
        for (int e = e1; e < e2; e++) { //energy_number
          if (dos[d][s][e] > dosmax) dosmax = dos[d][s][e];
        }
      }
    }

    if(LDEBUG) {cerr << __AFLOW_FUNC__ << " dosmax(FOUND)=" << dosmax << endl;}

    if (!dosscale.empty()) dosmax *= aurostd::string2utype<double>(dosscale);

    if(LDEBUG) {cerr << __AFLOW_FUNC__ << " dosmax(SCALED)=" << dosmax << endl;}

    // l = order of magnitude
    // x = scalar multiple
    // Now round up the numbers to give good tic values.
    int l = (int) log10(dosmax);
    int x = (int) ceil(dosmax/std::pow(10.0, l));

    if(LDEBUG){
      cerr << __AFLOW_FUNC__ << " l=" << l << endl;
      cerr << __AFLOW_FUNC__ << " x=" << x << endl;
    }

    //CO20191110 - add 30% for every orbital bigger than d to account for larger key
    if(0){  //CO20200404 - TOO much padding for LIB6 species-projection dos, better not to include any ad hoc fixes for general plotting
      if(ndos>4){x+=3*(ndos-4);}
      if(LDEBUG){cerr << __AFLOW_FUNC__ << " x(new1)=" << x << endl;}
    }

    //[CO20191110 OBSOLETE]if ((xdos.spin == 1) || (l > 1) || (l < -1) || (2 * x % 4 == 0)) {
    //[CO20191110 OBSOLETE]  dosmax = x * std::pow(10.0, l); 
    //[CO20191110 OBSOLETE]} else {
    //[CO20191110 OBSOLETE]  dosmax = (x + 1) * std::pow(10.0, l);
    //[CO20191110 OBSOLETE]}

    // The DOS axis of the DOS plot should be divided into four tics,
    // which gives a nice grid density. However, this often makes the
    // numbers on the axis look ugly because some will have a decimal
    // point and some won't. To avoid this, the scalar multiple (x) may
    // need to be increased except for the following circumstances:
    //  * The DOS is spin-polarized, in which case we only have two tics per spin.
    //  * The maximum (l) is larger than 10, so there are no decimal points.
    //  * The maximum (l) is smaller than 1, in which case they all have decimal points.
    //  * The number is divisible by 4.
    if (!((xdos.spin == 1) || (l > 1) || (l < -1) || (2 * x % 4 == 0))) {x++;}
    if(LDEBUG) {cerr << __AFLOW_FUNC__ << " x(new2)=" << x << endl;}

    dosmax = x * std::pow(10.0, l); 

    if(LDEBUG){cerr << __AFLOW_FUNC__ << " dosmax=" << dosmax << endl;} //CO20200404

    return dosmax;
  }

  //generateBandPlotGNUPLOT///////////////////////////////////////////////////
  // Generates the gnuplot script for band structure plots.
  void generateBandPlotGNUPLOT(stringstream& out, const xEIGENVAL& xeigen,
      const vector<double>& xvals, const vector<double>& ticvals,
      const vector<string>& ticlabels, const xoption& plotoptions) {
    bool LDEBUG=(FALSE || _DEBUG_PLOTTER_ || XHOST.DEBUG); 

    // Initialize variables
    double Efermi = aurostd::string2utype<double>(plotoptions.getattachedscheme("EFERMI"));
    double Emin = aurostd::string2utype<double>(plotoptions.getattachedscheme("XMIN"));
    double Emax = aurostd::string2utype<double>(plotoptions.getattachedscheme("XMAX"));
    bool banddos = plotoptions.flag("BANDDOS");

    if(LDEBUG){
      cerr << __AFLOW_FUNC__ << " Emin=" << Emin << endl;
      cerr << __AFLOW_FUNC__ << " Emax=" << Emax << endl;
    }

    string unit = plotoptions.getattachedscheme("UNIT");
    if (unit.empty()) unit = "EV";
    string energyLabel;
    //AS20210701 extra check for Grueneisen parameter dispersion plotting
    if (unit=="GRUENEISEN"){
      energyLabel = "$\\gamma$";
      unit = "";
    }
    else if (aurostd::substring2bool(unit, "EV")){
      energyLabel = "energy";
    }
    else{
      energyLabel = "frequency";
    }

    unit = getFormattedUnit(unit);

    out << "# Band structure plot" << std::endl;

    uint kpts_per_segment = xeigen.number_kpoints/(ticvals.size() - 1);
    // Create data block
    out << std::endl << "$band_data << EOD" << std::endl;
    for (uint k = 0; k < xeigen.number_kpoints; k++) {
      out << "  " << xvals[k];
      for (uint s = 0; s < xeigen.spin + 1; s++) {
        for (uint b = 0; b < xeigen.number_bands; b++) {
          out << " " << xeigen.venergy[k][b][s];
        }
      }
      out << std::endl;

      // Put each segment into one block. This prevents discontinuities
      // in the band structure from being connected in the plot
      if (((k + 1) % kpts_per_segment == 0) && (k + 1 < xeigen.number_kpoints)) {
        out << std::endl << std::endl;
      }
    }
    out << "EOD" << std::endl << std::endl;

    // Margins
    out << "# Margins" << std::endl;
    if (banddos) {
      out << "set lmargin at screen 0.08" << std::endl;
      out << "set rmargin at screen 0.70" << std::endl;
      out << "set tmargin at screen 0.9" << std::endl;
      out << "set bmargin at screen 0.12" << std::endl;
    } else {
      out << "set tmargin at screen 0.9" << std::endl;
      out << "set bmargin at screen 0.15" << std::endl;
    }

    // Key
    out << std::endl << "# Key" << std::endl;
    out << "unset key" << std::endl;

    // Axes
    out << std::endl << "# Axes" << std::endl;
    out << "unset xtics" << std::endl;
    out << "set xtics(";
    uint ntics = ticvals.size();
    for (uint i = 0; i < ntics; i++) {
      out << "'" << ticlabels[i] << "' " << ticvals[i];
      if (i < ntics - 1) out << ", ";
    }
    out << ")" << std::endl;
    out << "set tic scale 0" << std::endl;
    out << "set xrange [0:1]" << std::endl;
    out << "set yrange [" << Emin << ":" << Emax << "]" << std::endl;
    if (unit.empty()){//AS20210701 this might be a Grueneisen parameter dispersion plot, unitless
      out << "set ylabel '" << energyLabel << std::endl;
    }
    else{
      out << "set ylabel '" << energyLabel << " (" << unit << ")'" << std::endl;
    }

    // Fermi level
    if (Efermi > Emin) {
      out << std::endl << "# Fermi level" << std::endl;
      out << "set arrow from 0, " << Efermi << " to graph 1, first " << Efermi
        << " nohead lt 1 lc rgb '" << EFERMI_COLOR << "' lw 3" << std::endl;
    }

    // Plot data
    out << std::endl << "# Data" << std::endl;
    out << "plot ";
    // Majority spin
    for (uint b = 0; b < xeigen.number_bands; b++) {
      if (b > 0) out << "     ";
      out << "'$band_data' u 1:" << (b + 2)
        << " w l lt -1 lc rgb '" <<  ISPIN_COLORS[0] << "' lw 2"
        << (((xeigen.spin == 1) || (b < xeigen.number_bands - 1))?",\\":"") << std::endl;
    }

    // Minority spin
    if (xeigen.spin == 1) {
      for (uint b = 0; b < xeigen.number_bands; b++) {
        out << "      ";
        out << "'$band_data' u 1:" << (b + xeigen.number_bands + 2)
          << " w l lt 1 lc rgb '" << ISPIN_COLORS[1] << "' lw 2"
          << ((b < xeigen.number_bands - 1)?",\\":"") << std::endl;
      }
    }
  }

  //getFormattedUnit//////////////////////////////////////////////////////////
  // Formats the energy/frequency unit for band structures. This is
  // especially useful for phonons.
  string getFormattedUnit(const string& unit) {
    if (unit == "EV") return "eV";
    if (unit == "MEV") return "meV";
    if (unit == "THZ") return "THz";
    if (unit == "HZ") return "Hz";
    if ((unit == "CM-1") || (unit == "RCM")) return "cm$^{-1}$";
    return unit;
  }

}  // namespace plotter

//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//                                 PHONONS                                  //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////

namespace plotter {

  //PLOT_PHDOS////////////////////////////////////////////////////////////////
  // Plots phonon DOS.
  void PLOT_PHDOS(xoption& plotoptions,ostream& oss) {ofstream FileMESSAGE;return PLOT_PHDOS(plotoptions,FileMESSAGE,oss);} //CO20200404
  void PLOT_PHDOS(xoption& plotoptions,ofstream& FileMESSAGE,ostream& oss) {  //CO20200404
    // Set k-points format to LaTeX
    plotoptions.push_attached("KPOINT_FORMAT", "LATEX");
    // Set output format to gnuplot
    plotoptions.push_attached("OUTPUT_FORMAT", "GNUPLOT");

    stringstream out;
    PLOT_PHDOS(plotoptions, out,FileMESSAGE,oss);  //CO20200404
  }

  void PLOT_PHDOS(xoption& plotoptions, stringstream& out,ostream& oss) {ofstream FileMESSAGE;return PLOT_PHDOS(plotoptions,out,FileMESSAGE,oss);}  //CO20200404
  void PLOT_PHDOS(xoption& plotoptions, stringstream& out,ofstream& FileMESSAGE,ostream& oss) {  //CO20200404
    plotoptions.push_attached("EXTENSION", "phdos");
    plotoptions.push_attached("OUTPUT_FORMAT", "GNUPLOT");
    // Read files
    string directory = plotoptions.getattachedscheme("DIRECTORY");
    //ME20211008 - adding support for POCC
    //Check if directory contains PHDOSCAR.pocc files
    vector<string> vfiles, vpocc_doscars;
    aurostd::DirectoryLS(directory, vfiles);
    for (uint i = 0; i < vfiles.size(); i++) {
      // File name must start with prefix
      if (vfiles[i].find(POCC_PHDOSCAR_PREFIX) == 0) vpocc_doscars.push_back(vfiles[i]);
    }
    xDOSCAR xdos;
    if (vpocc_doscars.size() > 0) {
      string T = "";
      plotoptions.flag("PLOT_ALL_ATOMS", true);  // PARTCAR has no iatoms
      for (uint i = 0; i < vpocc_doscars.size(); i++) {
        // Grab the temperature string. Format is always prefix + T + "K" + extension
        T = vpocc_doscars[i].substr(POCC_PHDOSCAR_PREFIX.length());
        T = T.substr(0, T.find("K"));
        plotoptions.push_attached("EXTENSION", "phdos_T" + T + "K");
        xdos.GetPropertiesFile(directory + "/" + vpocc_doscars[i]);
        PLOT_PHDOS(plotoptions, out, xdos, FileMESSAGE, oss);
      }
    } else {
      xdos.GetPropertiesFile(directory+"/"+DEFAULT_APL_PHDOSCAR_FILE); //CO20191110
      PLOT_PHDOS(plotoptions, out, xdos, FileMESSAGE, oss);
    }
  }

  // ME20210927
  void PLOT_PHDOS(xoption& plotoptions, const xDOSCAR& xdos, ostream& oss) {ofstream FileMESSAGE; return PLOT_PHDOS(plotoptions,xdos,FileMESSAGE,oss);}
  void PLOT_PHDOS(xoption& plotoptions, const xDOSCAR& xdos, ofstream& FileMESSAGE, ostream& oss) {
    // Set k-points format to LaTeX
    plotoptions.push_attached("KPOINT_FORMAT", "LATEX");
    // Set output format to gnuplot
    plotoptions.push_attached("OUTPUT_FORMAT", "GNUPLOT");

    stringstream out;
    PLOT_PHDOS(plotoptions, out, xdos, FileMESSAGE, oss);
  }

  // Create a copy of xDOSCAR here because the energy units may be changed, which
  // should not propagate outside
  void PLOT_PHDOS(xoption& plotoptions, stringstream& out, xDOSCAR xdos, ofstream& FileMESSAGE, ostream& oss) {
    plotoptions.push_attached("DEFAULT_TITLE", xdos.title);
    setFileName(plotoptions);
    setTitle(plotoptions,FileMESSAGE,oss); //CO20200404

    plotoptions.flag("BANDDOS", false);

    // Convert energies if necessary
    string unit = plotoptions.getattachedscheme("UNIT");
    if (!unit.empty() && (unit != "EV")) {
      convertEnergies(xdos, unit);
    }

    // Set Emin and Emax
    setEMinMax(plotoptions, xdos.energy_min, xdos.energy_max);

    generateHeader(out, plotoptions, false);
    generateDosPlot(out, xdos, plotoptions,FileMESSAGE,oss);  //CO20200404
    savePlotGNUPLOT(plotoptions, out);
  }

  //PLOT_PHDISP///////////////////////////////////////////////////////////////
  // Plots phonon dispersion curves.
  void PLOT_PHDISP(xoption& plotoptions,ostream& oss) {ofstream FileMESSAGE;return PLOT_PHDISP(plotoptions,FileMESSAGE,oss);}  //CO20200404
  void PLOT_PHDISP(xoption& plotoptions,ofstream& FileMESSAGE,ostream& oss) {  //CO20200404
    // Set k-points format to LaTeX
    plotoptions.push_attached("KPOINT_FORMAT", "LATEX");
    // Set output format to gnuplot
    plotoptions.push_attached("OUTPUT_FORMAT", "GNUPLOT");

    stringstream out;
    PLOT_PHDISP(plotoptions, out,FileMESSAGE,oss); //CO20200404
  }

  void PLOT_PHDISP(xoption& plotoptions, stringstream& out,ostream& oss) {ofstream FileMESSAGE;return PLOT_PHDISP(plotoptions,out,FileMESSAGE,oss);} //CO20200404
  void PLOT_PHDISP(xoption& plotoptions, stringstream& out,ofstream& FileMESSAGE,ostream& oss) { //CO20200404
    plotoptions.push_attached("EXTENSION", "phdisp");
    plotoptions.push_attached("OUTPUT_FORMAT", "GNUPLOT");
    // Read files
    string directory = plotoptions.getattachedscheme("DIRECTORY");
    xEIGENVAL xeigen;
    xeigen.GetPropertiesFile(directory+"/"+DEFAULT_APL_PHEIGENVAL_FILE); //CO20191110
    xKPOINTS xkpts;
    xkpts.GetPropertiesFile(directory+"/"+DEFAULT_APL_PHKPOINTS_FILE); //CO20191110
    stringstream poscar;
    aurostd::efile2stringstream(directory+"/"+DEFAULT_APL_PHPOSCAR_FILE, poscar);  //CO20191110
    xstructure xstr(poscar);

    plotoptions.push_attached("DEFAULT_TITLE", xeigen.title);
    plotoptions.push_attached("LATTICE", getLatticeFromKpointsTitle(xkpts.title));
    setFileName(plotoptions);
    setTitle(plotoptions,FileMESSAGE,oss); //CO20200404

    plotoptions.flag("BANDDOS", false);

    // Convert energies if necessary
    string unit = plotoptions.getattachedscheme("UNIT");
    if (!unit.empty() && (unit != "EV")) {
      convertEnergies(xeigen, unit);
    }

    // Set Emin and Emax
    setEMinMax(plotoptions, xeigen.energy_min, xeigen.energy_max);

    generateHeader(out, plotoptions, false);
    generateBandPlot(out, xeigen, xkpts, xstr, plotoptions);
    savePlotGNUPLOT(plotoptions, out);
  }

  //PLOT_PHDISPDOS////////////////////////////////////////////////////////////
  // Plots combined phonon band structure + DOS plots.
  void PLOT_PHDISPDOS(xoption& plotoptions,ostream& oss) {ofstream FileMESSAGE;return PLOT_PHDISPDOS(plotoptions,FileMESSAGE,oss);} //CO20200404
  void PLOT_PHDISPDOS(xoption& plotoptions,ofstream& FileMESSAGE,ostream& oss) { //CO20200404
    // Set k-points format to LaTeX
    plotoptions.push_attached("KPOINT_FORMAT", "LATEX");
    // Set output format to gnuplot
    plotoptions.push_attached("OUTPUT_FORMAT", "GNUPLOT");

    stringstream out;
    PLOT_PHDISPDOS(plotoptions, out,FileMESSAGE,oss);  //CO20200404
  }

  void PLOT_PHDISPDOS(xoption& plotoptions, stringstream& out,ostream& oss) {ofstream FileMESSAGE;return PLOT_PHDISPDOS(plotoptions,out,FileMESSAGE,oss);}  //CO20200404
  void PLOT_PHDISPDOS(xoption& plotoptions, stringstream& out,ofstream& FileMESSAGE,ostream& oss) {  //CO20200404
    plotoptions.push_attached("EXTENSION", "phdispdos");
    plotoptions.push_attached("OUTPUT_FORMAT", "GNUPLOT");
    plotoptions.push_attached("PLOT_SIZE", BANDDOS_SIZE);
    // Read files
    string directory = plotoptions.getattachedscheme("DIRECTORY");
    xDOSCAR xdos;
    xdos.GetPropertiesFile(directory+"/"+DEFAULT_APL_PHDOSCAR_FILE); //CO20191110
    xEIGENVAL xeigen;
    xeigen.GetPropertiesFile(directory+"/"+DEFAULT_APL_PHEIGENVAL_FILE); //CO20191110
    xKPOINTS xkpts;
    xkpts.GetPropertiesFile(directory+"/"+DEFAULT_APL_PHKPOINTS_FILE); //CO20191110
    stringstream poscar;
    aurostd::efile2stringstream(directory+"/"+DEFAULT_APL_PHPOSCAR_FILE, poscar);  //CO20191110
    xstructure xstr(poscar);

    plotoptions.push_attached("DEFAULT_TITLE", xeigen.title);
    plotoptions.push_attached("LATTICE", getLatticeFromKpointsTitle(xkpts.title));
    setFileName(plotoptions);
    setTitle(plotoptions,FileMESSAGE,oss); //CO20200404

    plotoptions.flag("BANDDOS", true);

    // Convert energies if necessary
    string unit = plotoptions.getattachedscheme("UNIT");
    if (!unit.empty() && (unit != "EV")) {
      convertEnergies(xdos, unit);
      convertEnergies(xeigen, unit);
    }

    // Set Emin and Emax
    setEMinMax(plotoptions, xdos.energy_min, xdos.energy_max);

    generateHeader(out, plotoptions, true);
    generateBandPlot(out, xeigen, xkpts, xstr, plotoptions);
    generateDosPlot(out, xdos, plotoptions,FileMESSAGE,oss);  //CO20200404
    savePlotGNUPLOT(plotoptions, out);
  }

  //convertEnergies///////////////////////////////////////////////////////////
  // Converts the energies in an electronic structure
  // object (xDOSCAR/xEIGENVAL) into the desired energy/frequency unit.
  void convertEnergies(xEIGENVAL& xeigen, const string& unit) {
    double conversion_factor = getEnergyConversionFactor(unit);
    for (uint k = 0; k < xeigen.number_kpoints; k++) {
      for (uint b = 0; b < xeigen.number_bands; b++) {
        for (uint s = 0; s < xeigen.spin + 1; s++) {
          xeigen.venergy[k][b][s] *= conversion_factor;
        }
      }
    }
    xeigen.energy_min *= conversion_factor;
    xeigen.energy_max *= conversion_factor;
  }

  void convertEnergies(xDOSCAR& xdos, const string& unit) {
    double conversion_factor = getEnergyConversionFactor(unit);
    for (uint i = 0; i < xdos.number_energies; i++) {
      xdos.venergy[i] *= conversion_factor;
      xdos.venergyEf[i] *= conversion_factor;
    }
    xdos.energy_min *= conversion_factor;
    xdos.energy_max *= conversion_factor;
  }

  //getEnergyConversionFactor/////////////////////////////////////////////////
  // Returns the factor to convert eV into the desired energy/frequency unit.
  // Supported units are meV, THz, Hz, and reciprocal cm (CM-1/RCM).
  //ME20200121 - Replaced with constants from xscalar.
  double getEnergyConversionFactor(const string& unit) {
    if (unit == "MEV") return 1000.0;
    if (unit == "THZ") return (eV2Hz * Hz2THz);
    if (unit == "HZ") return eV2Hz;
    if ((unit == "CM-1") || (unit == "RCM")) return eV2rcm;
    return 1.0;
  }

}  // namespace plotter

//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//                           PROPERTIES PLOTTERS                            //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////

namespace plotter {

  //PLOT_THERMO///////////////////////////////////////////////////////////////
  // Plots APL thermal properties.
  void PLOT_THERMO(xoption& plotoptions,ostream& oss) {ofstream FileMESSAGE;return PLOT_THERMO(plotoptions,FileMESSAGE,oss);}  //CO20200404
  void PLOT_THERMO(xoption& plotoptions,ofstream& FileMESSAGE,ostream& oss) {  //CO20200404
    stringstream out;
    plotoptions.push_attached("OUTPUT_FORMAT", "GNUPLOT");
    plotoptions.push_attached("COLOR", "#000000");
    plotoptions.push_attached("LINETYPES", "-1");
    PLOT_THERMO(plotoptions, out,FileMESSAGE,oss); //CO20200404
  }

  void PLOT_THERMO(xoption& plotoptions, stringstream& out,ostream& oss) {ofstream FileMESSAGE;return PLOT_THERMO(plotoptions,out,FileMESSAGE,oss);} //CO20200404
  void PLOT_THERMO(xoption& plotoptions, stringstream& out,ofstream& FileMESSAGE,ostream& oss) { //CO20200404
    // Set labels
    static const int nprops = 4;
    string ylabels[nprops] = {"U", "F", "S", "c_V"};
    string extensions[nprops] = {"vib_internal_energy", "vib_free_energy", "vib_entropy", "apl_cV"};
    string yunits[nprops] = {"meV/cell", "meV/cell", "$k_B$/cell", "$k_B$/cell"};
    string ymin[nprops] = {"", "", "0", "0"};

    // Get data
    string directory = plotoptions.getattachedscheme("DIRECTORY");
    string thermo_file = directory + "/" + DEFAULT_APL_FILE_PREFIX + DEFAULT_APL_THERMO_FILE;
    //ME20200413 - Since multiple data files are plotted, the user file
    // name functions as base file name.
    string user_file_name = plotoptions.getattachedscheme("FILE_NAME_USER");
    plotoptions.pop_attached("FILE_NAME_USER");
    if (aurostd::EFileExist(thermo_file)) {
      string outformat = plotoptions.getattachedscheme("OUTPUT_FORMAT");
      plotoptions.push_attached("DATA_FILE", thermo_file);
      plotoptions.push_attached("KEYWORD", "APL_THERMO");
      vector<vector<double> > data = readAflowDataFile(plotoptions);
      if (!user_file_name.empty()) plotoptions.push_attached("DEFAULT_TITLE", user_file_name);  //ME20200413
      for (int i = 0; i < nprops; i++) {
        plotoptions.pop_attached("YMIN");
        if (!ymin[i].empty()) plotoptions.push_attached("YMIN", ymin[i]);
        plotoptions.push_attached("EXTENSION", extensions[i]);
        setPlotLabels(plotoptions, "T", "K", ylabels[i], yunits[i]);
        plotSingleFromSet(plotoptions, out, data, i + 2,FileMESSAGE,oss);  //CO20200404
        if (outformat == "GNUPLOT") {
          savePlotGNUPLOT(plotoptions, out);
        }
        out.str("");  //ME20200513 - reset stringstream
      }
    } else {
      string message = "Could not find file " + thermo_file + ".";
      throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__, message, _FILE_NOT_FOUND_);
    }
  }

  //AS20200909 BEGIN
  //PLOT_THERMO_QHA///////////////////////////////////////////////////////////////
  // Plots QHA thermal properties.
  void PLOT_THERMO_QHA(xoption& plotoptions,ostream& oss) {ofstream FileMESSAGE; PLOT_THERMO_QHA(plotoptions,FileMESSAGE,oss);}  //CO20200404
  void PLOT_THERMO_QHA(xoption& plotoptions,ofstream& FileMESSAGE,ostream& oss) {  //CO20200404
    stringstream out;
    plotoptions.push_attached("OUTPUT_FORMAT", "GNUPLOT");
    plotoptions.push_attached("COLOR", "#000000");
    plotoptions.push_attached("LINETYPES", "-1");
    PLOT_THERMO_QHA(plotoptions, out,FileMESSAGE,oss); //CO20200404
  }

  void PLOT_THERMO_QHA(xoption& plotoptions, stringstream& out,ostream& oss) {ofstream FileMESSAGE; PLOT_THERMO_QHA(plotoptions,out,FileMESSAGE,oss);} //CO20200404
  void PLOT_THERMO_QHA(xoption& plotoptions, stringstream& out,ofstream& FileMESSAGE,ostream& oss) 
  {
    string msg = "";

    // Set labels
    static const int nprops = 7;
    string ylabels[nprops] = {"V", "F", "B", "\\beta", "c_V", "c_P", "\\gamma"};
    string extensions[nprops] = {"volume_equilibrium_qha", "energy_free_qha",
      "modulus_bulk_qha", "thermal_expansion_qha", "cV_qha", "cP_qha",
      "gruneisen_parameter_qha"};
    string yunits[nprops] = {"\\AA$^{3}$/atom", "eV/atom", "GPa", "$10^{-5}K^{-1}$",
      "$k_B$/atom", "$k_B$/atom", ""};
    string ymin[nprops] = {"", "", "", "", "0", "0", ""};

    string eos_model = plotoptions.getattachedscheme("EOSMODEL");
    if (eos_model != "SJ"  && eos_model != "BM2" && eos_model != "BM3" &&
        eos_model != "BM4" && eos_model != "M"){
      msg = "Wrong name of the EOS model was specified. ";
      msg += "Only SJ, BM2, BM3, BM4 or M labels are allowed.";
      throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__, msg, _INPUT_ILLEGAL_);
    }
    string keyword = "QHA_" + eos_model + "_THERMO";

    // Get data
    string directory = plotoptions.getattachedscheme("DIRECTORY");
    string thermo_file = directory+"/"+DEFAULT_QHA_FILE_PREFIX+DEFAULT_QHA_THERMO_FILE;
    //ME20200413 - Since multiple data files are plotted, the user file
    // name functions as base file name.
    string user_file_name = plotoptions.getattachedscheme("FILE_NAME_USER");
    plotoptions.pop_attached("FILE_NAME_USER");
    if (aurostd::EFileExist(thermo_file)) {
      string outformat = plotoptions.getattachedscheme("OUTPUT_FORMAT");
      plotoptions.push_attached("DATA_FILE", thermo_file);
      plotoptions.push_attached("KEYWORD", keyword);
      vector<vector<double> > data = readAflowDataFile(plotoptions);
      if (!user_file_name.empty()) plotoptions.push_attached("DEFAULT_TITLE", user_file_name);  //ME20200413
      for (int i = 0; i < nprops; i++) {
        plotoptions.pop_attached("YMIN");
        if (!ymin[i].empty()) plotoptions.push_attached("YMIN", ymin[i]);
        plotoptions.push_attached("EXTENSION", extensions[i] + '_' + eos_model);
        setPlotLabels(plotoptions, "T", "K", ylabels[i], yunits[i]);
        plotSingleFromSet(plotoptions, out, data, i + 1,FileMESSAGE,oss);  //CO20200404
        if (outformat == "GNUPLOT") {
          savePlotGNUPLOT(plotoptions, out);
        }
        out.str("");  //ME20200513 - reset stringstream
      }
    } else {
      msg = "Could not find file " + thermo_file + ".";
      throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__, msg, _FILE_NOT_FOUND_);
    }
  }
  //AS20200909 END

  //AS20210701 BEGIN
  //PLOT_GRUENEISEN_DISPERSION///////////////////////////////////////////////////////////////
  /// Plots Grueneisen parameter dispersion curves.
  /// Follows PLOT_GRUENEISEN_DISPERSION function.
  void PLOT_GRUENEISEN_DISPERSION(xoption& plotoptions,ostream& oss) {ofstream FileMESSAGE;return PLOT_GRUENEISEN_DISPERSION(plotoptions,FileMESSAGE,oss);}  //CO20200404
  void PLOT_GRUENEISEN_DISPERSION(xoption& plotoptions,ofstream& FileMESSAGE,ostream& oss) {  //CO20200404
    // Set k-points format to LaTeX
    plotoptions.push_attached("KPOINT_FORMAT", "LATEX");
    // Set output format to gnuplot
    plotoptions.push_attached("OUTPUT_FORMAT", "GNUPLOT");

    stringstream out;
    PLOT_GRUENEISEN_DISPERSION(plotoptions, out,FileMESSAGE,oss); //CO20200404
    savePlotGNUPLOT(plotoptions, out);
  }

  void PLOT_GRUENEISEN_DISPERSION(xoption& plotoptions, stringstream& out,ostream& oss) {ofstream FileMESSAGE;return PLOT_GRUENEISEN_DISPERSION(plotoptions,out,FileMESSAGE,oss);} //CO20200404
  void PLOT_GRUENEISEN_DISPERSION(xoption& plotoptions, stringstream& out,ofstream& FileMESSAGE,ostream& oss) { //CO20200404
    // Grueneisen parameters for acoustic modes at the Gamma point are ill-defined,
    // so for the plot to be pretty, one needs to substitute it with NaN
    static double nan = std::numeric_limits<double>::quiet_NaN();

    plotoptions.push_attached("EXTENSION", "grdisp");
    plotoptions.push_attached("OUTPUT_FORMAT", "GNUPLOT");
    // Read files
    string directory = plotoptions.getattachedscheme("DIRECTORY");
    xEIGENVAL xeigen;
    xeigen.GetPropertiesFile(directory+"/"+ DEFAULT_QHA_FILE_PREFIX
        + DEFAULT_QHA_GP_PATH_FILE);
    xKPOINTS xkpts;
    xkpts.GetPropertiesFile(directory+"/"+ DEFAULT_QHA_FILE_PREFIX +
        DEFAULT_QHA_KPOINTS_FILE);
    stringstream poscar;
    aurostd::efile2stringstream(directory+"/"+DEFAULT_APL_PHPOSCAR_FILE, poscar);  //CO20191110
    xstructure xstr(poscar);

    // substituting the values of Grueneisen parameters for acoustic modes at
    // Gamma point with NaN.
    // Grueneisen parameters at Gamma may be large due to numerical noise, so
    // recalculate energy_min and energy_max
    xeigen.energy_max = -AUROSTD_MAX_DOUBLE;
    xeigen.energy_min =  AUROSTD_MAX_DOUBLE;
    double eigval = 0.0;
    for (uint k=0; k<xeigen.number_kpoints; k++){
      for (uint b=0; b<xeigen.number_bands; b++){
        if (b<=3){// acoustic modes should be the first three ones
          // if we are close enough to Gamma, substitute with NaN
          if (aurostd::modulus(xeigen.vkpoint[k]) < AUROSTD_IDENTITY_TOL){
            xeigen.venergy[k][b][0] = nan;
          }
        }

        eigval = xeigen.venergy[k][b][0];
        if (eigval < xeigen.energy_min) xeigen.energy_min = eigval;
        if (eigval > xeigen.energy_max) xeigen.energy_max = eigval;
      }
    }

    // proceed with plot setup
    plotoptions.push_attached("DEFAULT_TITLE", xeigen.title);
    plotoptions.push_attached("LATTICE", getLatticeFromKpointsTitle(xkpts.title));
    setFileName(plotoptions);
    setTitle(plotoptions,FileMESSAGE,oss); //CO20200404

    plotoptions.flag("BANDDOS", false);

    plotoptions.push_attached("UNIT", "GRUENEISEN");

    // Set Emin and Emax
    setEMinMax(plotoptions, xeigen.energy_min, xeigen.energy_max);

    generateHeader(out, plotoptions, false);
    generateBandPlot(out, xeigen, xkpts, xstr, plotoptions);
  }
  //AS20210701 END

  //PLOT_TCOND////////////////////////////////////////////////////////////////
  // Plots AAPL thermal conductivity tensors
  void PLOT_TCOND(xoption& plotoptions,ostream& oss) {ofstream FileMESSAGE;return PLOT_TCOND(plotoptions,FileMESSAGE,oss);} //CO20200404
  void PLOT_TCOND(xoption& plotoptions,ofstream& FileMESSAGE,ostream& oss) { //CO20200404
    stringstream out;
    plotoptions.push_attached("OUTPUT_FORMAT", "GNUPLOT");
    plotoptions.flag("LINESPOINTS", true);
    PLOT_TCOND(plotoptions, out,FileMESSAGE,oss);  //CO20200404
    savePlotGNUPLOT(plotoptions, out);
  }

  void PLOT_TCOND(xoption& plotoptions, stringstream& out,ostream& oss) {ofstream FileMESSAGE;return PLOT_TCOND(plotoptions,out,FileMESSAGE,oss);}  //CO20200404
  void PLOT_TCOND(xoption& plotoptions, stringstream& out,ofstream& FileMESSAGE,ostream& oss) {  //CO20200404
    plotoptions.push_attached("EXTENSION", "thermal_conductivity");
    string directory = plotoptions.getattachedscheme("DIRECTORY");
    string tcond_file = directory + "/" + DEFAULT_AAPL_FILE_PREFIX + DEFAULT_AAPL_TCOND_FILE;
    if (aurostd::EFileExist(tcond_file)) {
      plotoptions.push_attached("DATA_FILE", tcond_file);
      plotoptions.push_attached("KEYWORD", "AAPL_THERMAL_CONDUCTIVITY");
      plotoptions.flag("CONTRAVARIANT", true);
      plotoptions.push_attached("YMIN", "0");
      plotoptions.flag("LEGEND_HORIZONTAL", true); //CO20200404
      plotoptions.push_attached("LEGEND_MAXCOLS","3");  //CO20200404
      setPlotLabels(plotoptions, "T", "K", "\\kappa", "W/m K");
      plotMatrix(plotoptions, out,FileMESSAGE,oss);  //CO20200404
    } else {
      string message = "Could not find file " + tcond_file + ".";
      throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__, message, _FILE_NOT_FOUND_);
    }
  }

}  // namespace plotter

//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//                              GENERAL PLOTS                               //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////

namespace plotter {

  // Color palettes ----------------------------------------------------------

  // The color palette was designed to be accessible for people with color
  // vision deficiencies. When changing the colors, please make sure that they
  // are still distinguishable for everyone!
  static const string MATRIX_COLORS = "#000000,#004949,#009292,#490092,#B66DFF,#6DB6FF,#924900,#D55E00,#EDB120";

  // Point styles ------------------------------------------------------------

  static const string MATRIX_POINT_STYLES = "17,35,51,44,18,9,60,11,20";

  // Line types --------------------------------------------------------------

  static const string MATRIX_LINE_TYPES = "-1";

  // Pre-set labels ---------------------------------------------------------

  static const string MATRIX_LABELS[9] = {"xx", "yx", "zx",
    "xy", "yy", "zy",
    "xz", "yz", "zz"};


  //plotSingleFromSet/////////////////////////////////////////////////////////
  // Plots a single column from a dataset.
  void plotSingleFromSet(xoption& plotoptions, stringstream& out,const vector<vector<double> >& data_set, int col,ostream& oss) {ofstream FileMESSAGE;return plotSingleFromSet(plotoptions,out,data_set,col,FileMESSAGE,oss);} //CO20200404
  void plotSingleFromSet(xoption& plotoptions, stringstream& out,const vector<vector<double> >& data_set, int col,ofstream& FileMESSAGE,ostream& oss) { //CO20200404
    setFileName(plotoptions);
    setTitle(plotoptions,FileMESSAGE,oss); //CO20200404

    vector<vector<double> > data(data_set.size(), vector<double>(2));
    for (uint i = 0; i < data.size(); i++) {
      data[i][0] = data_set[i][0];
      data[i][1] = data_set[i][col];
    }

    generateHeader(out, plotoptions);
    generatePlotGNUPLOT(out, plotoptions, data);
  }

  //plotMatrix////////////////////////////////////////////////////////////////
  // Plots a 3 x 3 matrix.
  void plotMatrix(xoption& plotoptions, stringstream& out,ostream& oss) {ofstream FileMESSAGE;return plotMatrix(plotoptions,out,FileMESSAGE,oss);}  //CO20200404
  void plotMatrix(xoption& plotoptions, stringstream& out,ofstream& FileMESSAGE,ostream& oss) {  //CO20200404
    vector<vector<double> > data = readAflowDataFile(plotoptions);
    setFileName(plotoptions);
    setTitle(plotoptions,FileMESSAGE,oss); //CO20200404

    // Set additional plot options
    bool contravariant = plotoptions.flag("CONTRAVARIANT");
    string title_base = plotoptions.getattachedscheme("YLABEL");
    vector<string> titles(9);
    vector<string> colors(9);
    vector<int> point_styles(9);
    for (int i = 0; i < 9; i++) {
      if (contravariant) titles[i] = "$" + title_base + "^{" + MATRIX_LABELS[i] + "}$";
      else titles[i] = "$" + title_base + "_{" + MATRIX_LABELS[i] + "}$";
    }
    plotoptions.push_attached("TITLES", aurostd::joinWDelimiter(titles, ","));
    plotoptions.push_attached("COLORS", MATRIX_COLORS);
    plotoptions.push_attached("POINTSTYLES", MATRIX_POINT_STYLES);
    plotoptions.push_attached("LINETYPES", MATRIX_LINE_TYPES);

    generateHeader(out, plotoptions);
    generatePlotGNUPLOT(out, plotoptions, data);
  }

  //setPlotLabels/////////////////////////////////////////////////////////////
  // Stores the labels and units for the plots.
  void setPlotLabels(xoption& plotoptions,
      const string& xlabel, const string& xunit,
      const string& ylabel, const string& yunit) {
    plotoptions.pop_attached("XLABEL");
    plotoptions.push_attached("XLABEL", xlabel);
    plotoptions.pop_attached("XUNIT");
    plotoptions.push_attached("XUNIT", xunit);
    plotoptions.pop_attached("YLABEL");
    plotoptions.push_attached("YLABEL", ylabel);
    plotoptions.pop_attached("YUNIT");
    plotoptions.push_attached("YUNIT", yunit);
  }

  //readDataFile//////////////////////////////////////////////////////////////
  // Reads data from an AFLOW data file. Requires a START and STOP string to
  // be present so that it can skip headers and other data sets.
  vector<vector<double> > readAflowDataFile(xoption& plotoptions) {
    string message = "";
    vector<vector<double> > data;
    vector<double> row;
    vector<string> vcontent;
    string keyword = plotoptions.getattachedscheme("KEYWORD");
    string path_to_file = plotoptions.getattachedscheme("DATA_FILE");
    string startstring = "[" + keyword + "]START";
    string stopstring = "[" + keyword + "]STOP";
    string systemstring = "[" + keyword + "]SYSTEM=";
    aurostd::efile2vectorstring(path_to_file, vcontent);
    string line;
    uint nlines = vcontent.size();
    for (uint iline = 0; iline < nlines; iline++) {
      line = aurostd::RemoveWhiteSpacesFromTheFrontAndBack(vcontent[iline]);
      if (aurostd::substring2bool(line, systemstring)) {
        vector<string> tokens;
        aurostd::string2tokens(line, tokens, "=");
        if (tokens.size() == 2) plotoptions.push_attached("DEFAULT_TITLE", tokens.back());
      }
      if (vcontent[iline] == startstring) {
        iline++;
        while ((vcontent[iline] != stopstring) && (iline < nlines)) {
          line = aurostd::RemoveWhiteSpacesFromTheFront(vcontent[iline]);
          if (line[0] != '#') {
            aurostd::string2tokens(vcontent[iline], row, " ");
            data.push_back(row);
          }
          iline++;
        }
      }
      if (vcontent[iline] == stopstring) break;
      if (iline == nlines) {
        message = "Wrong file format. No STOP tag found.";
        throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__, message, _FILE_WRONG_FORMAT_);
      }
    }
    if (data.size() == 0) {
      message = "No data extracted from file " + path_to_file + ".";
      message += "File is either empty or has the wrong format.";
      throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__, message, _FILE_WRONG_FORMAT_);
    }
    return data;
  }

  //generatePlotGNUPLOT///////////////////////////////////////////////////////
  // Generate the gnuplot scripts for general plots.
  void generatePlotGNUPLOT(stringstream& out, const xoption& plotoptions,
      const vector<vector<double> >& data) {
    bool LDEBUG=(FALSE || _DEBUG_PLOTTER_ || XHOST.DEBUG); 
    uint ndata = data[0].size() - 1;

    // Axes settings
    string xmax = plotoptions.getattachedscheme("XMAX");
    string xmin = plotoptions.getattachedscheme("XMIN");
    string ymax = plotoptions.getattachedscheme("YMAX");
    string ymin = plotoptions.getattachedscheme("YMIN");
    string xlabel = plotoptions.getattachedscheme("XLABEL");
    string xunit = plotoptions.getattachedscheme("XUNIT");
    if (!xunit.empty()) xunit = " (" + xunit + ")";
    string ylabel = plotoptions.getattachedscheme("YLABEL");
    string yunit = plotoptions.getattachedscheme("YUNIT");
    if (!yunit.empty()) yunit = " (" + yunit + ")";

    if(LDEBUG){
      cerr << __AFLOW_FUNC__ << " ymin=" << ymin << endl;
      cerr << __AFLOW_FUNC__ << " ymax=" << ymax << endl;
    }

    // Plot types: lines (default), linespoints, or points only (nolines)
    bool linespoints = plotoptions.flag("LINESPOINTS");
    bool points = linespoints || plotoptions.flag("NOLINES");
    bool lines = (linespoints || !points);
    string plotstyle;
    if (lines) plotstyle += "l";
    if (points) plotstyle += "p";

    // Colors, point styles, line types
    vector<string> colors;
    vector<int> point_styles, line_types;
    aurostd::string2tokens(plotoptions.getattachedscheme("COLORS"), colors, ", ");
    aurostd::string2tokens(plotoptions.getattachedscheme("POINTSTYLES"), point_styles, ", ");
    aurostd::string2tokens(plotoptions.getattachedscheme("LINETYPES"), line_types, ", ");
    // Assume the same color, etc. for all plots if only one is given
    if ((ndata > 1) && (colors.size() == 1)) colors.assign(ndata, colors[0]);
    if ((ndata > 1) && (point_styles.size() == 1)) point_styles.assign(ndata, point_styles[0]);
    if ((ndata > 1) && (line_types.size() == 1)) line_types.assign(ndata, line_types[0]);
    bool colors_set = (colors.size() == ndata);
    bool points_set = ((linespoints || points) && (point_styles.size() == ndata));
    bool lines_set = (lines && (line_types.size() == ndata));

    // Titles
    vector<string> titles;
    aurostd::string2tokens(plotoptions.getattachedscheme("TITLES"), titles, ", ");
    if (titles.size() < ndata) {  // Fill up with empty titles if there aren't enough
      for (uint i = titles.size(); i < ndata; i++) titles.push_back("");
    }

    // Data block
    out << "$matrix_data << EOD" << std::endl;
    for (uint i = 0; i < data.size(); i++) {
      for (uint j = 0; j < ndata + 1; j++) {  // ndata + 1 to include x-values
        out << " " << data[i][j];
      }
      out << std::endl;
    }
    out << "EOD" << std::endl << std::endl;

    // Margins
    out << "# Margins" << std::endl;
    out << "set tmargin at screen 0.9" << std::endl;
    out << "set bmargin at screen 0.22" << std::endl;

    // Axes
    out << "# Axes" << std::endl;
    out << "set xrange [" << xmin << ":" << xmax << "]" << std::endl;
    out << "set yrange [" << ymin << ":" << ymax << "]" << std::endl;
    out << "set xlabel '$" << xlabel << "$" << xunit << "'" << std::endl;
    out << "set ylabel '$" << ylabel << "$" << yunit << "'" << std::endl;
    out << "set tics nomirror out" << std::endl;

    // Key
    out << std::endl << "# Key" << std::endl;
    if (ndata == 1) {  // No need for legend if only one set of data to plot
      out << "unset key" << std::endl;
    } else {
      if(plotoptions.flag("LEGEND_HORIZONTAL")){  //CO20200404
        int maxcols=3;
        string maxcols_str=plotoptions.getattachedscheme("LEGEND_MAXCOLS");
        if(!maxcols_str.empty()){maxcols=aurostd::string2utype<int>(maxcols_str);}
        out << "set key horizontal maxcols " << maxcols << std::endl;
      }
    }

    // Plot
    out << std::endl << "# Plot" << std::endl;
    out << "plot ";
    for (uint i = 0; i < ndata; i++) {
      if (i > 0) out << "     ";
      out << "'$matrix_data' u 1:" << (i + 2) << " w " << plotstyle;
      if (lines) out << " lw 2";
      if (lines_set) out << " lt " << line_types[i];
      if (points) out << " ps 1.5";
      if (colors_set) out << " lc rgb '" << colors[i] << "'";
      if (points_set) out << " pt " << point_styles[i];
      if (titles[i].empty()) {
        out << " notitle";
      } else {
        out << " title '" << titles[i] << "'";
      }
      out << ((i < ndata - 1)?",\\":"") << std::endl;
    }
  }

}  // namespace plotter

// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *            Aflow MARCO ESTERS - Duke University 2019-2021               *
// *                                                                         *
// ***************************************************************************
