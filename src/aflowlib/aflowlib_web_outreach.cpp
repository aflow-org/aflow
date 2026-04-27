// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2024           *
// *                                                                         *
// ***************************************************************************
// Stefano Curtarolo - Duke

#include "aflowlib/aflowlib_web_outreach.h"

#include "config.h"

#include <algorithm>
#include <ctime>
#include <iostream>
#include <ostream>
#include <sstream>
#include <string>
#include <vector>

#include "AUROSTD/aurostd.h"
#include "AUROSTD/aurostd_xerror.h"
#include "AUROSTD/aurostd_xfile.h"
#include "AUROSTD/aurostd_xscalar.h"

#include "aflow_xhost.h"
#include "aflowlib/aflowlib.h"
#include "flow/aflow_ivasp.h"

#define TOC 0

#define WEB_PDF string("http://" + XHOST.AFLOW_MATERIALS_SERVER + "/auro/AUROARTICULA/")
#define WEB_BIBTEX string("http://" + XHOST.AFLOW_MATERIALS_SERVER + "/auro/AUROARTICULA/BIBITEMS/")
#define WEB_BIBTEX_FILE string("/www/auro/AUROARTICULA/BIBITEMS/")
#define WEB_DOI string("http://dx.doi.org/")

#define THRUST_RECENT_ARTICLES 20
#define THRUST_RECENT_YEARS 5
#define AUTHOR_RECENT_ARTICLES 300
#define AUTHOR_RECENT_YEARS 30

// use MAX_YEAR_PRESENTATIONS if you don't want the current year as cutoff for presentations
// #define MAX_YEAR_PRESENTATIONS 2025

using std::cerr;
using std::cout;
using std::endl;
using std::istream;
using std::ostream;
using std::string;
using std::stringstream;
using std::vector;

// ******************************************************************************************************************************************************
// _OUTREACH CLASS
// namespace web

_outreach::_outreach() {
  newflag = false;
  wnumber = 0;
  year = 0;
  vauthor.clear();
  vcorrespondingauthor.clear();
  title = "";
  journal = "";
  arxiv = "";
  supplementary = "";
  supplementary_url = "";
  place = "";
  date = "";
  link = "";
  type = "";
  _isinvited = false;
  _isonline = false;
  host = "";
  abstract = "";
  pdf = "";
  doi = "";
  bibtex = "";
  bibtex_journal = "";
  bibtex_volume = "";
  bibtex_issue = "";
  bibtex_pages = "";
  bibtex_year = "";
  vextra_html.clear();
  vextra_latex.clear();
  vkeyword.clear();
  vsponsor.clear();
  valloy.clear();
}

// destructor
_outreach::~_outreach() {
  free();
}

void _outreach::free() {}

void _outreach::copy(const _outreach& b) {
  // const _outreach& _outreach::operator=(const _outreach& b)       // operator=
  free();
  newflag = b.newflag;
  wnumber = b.wnumber;
  year = b.year;
  vauthor.clear();
  for (size_t i = 0; i < b.vauthor.size(); i++) {
    vauthor.push_back(b.vauthor[i]);
  }
  vcorrespondingauthor.clear();
  for (size_t i = 0; i < b.vcorrespondingauthor.size(); i++) {
    vcorrespondingauthor.push_back(b.vcorrespondingauthor[i]);
  }
  title = b.title;
  journal = b.journal;
  arxiv = b.arxiv;
  supplementary = b.supplementary;
  supplementary_url = b.supplementary_url;
  place = b.place;
  date = b.date;
  link = b.link;
  type = b.type;
  _isinvited = b._isinvited;
  _isonline = b._isonline;
  host = b.host;
  abstract = b.abstract;
  pdf = b.pdf;
  doi = b.doi;
  bibtex = b.bibtex;
  bibtex_journal = b.bibtex_journal;
  bibtex_volume = b.bibtex_volume;
  bibtex_issue = b.bibtex_issue;
  bibtex_pages = b.bibtex_pages;
  bibtex_year = b.bibtex_year;
  vextra_html.clear();
  for (size_t i = 0; i < b.vextra_html.size(); i++) {
    vextra_html.push_back(b.vextra_html[i]);
  }
  vextra_latex.clear();
  for (size_t i = 0; i < b.vextra_latex.size(); i++) {
    vextra_latex.push_back(b.vextra_latex[i]);
  }
  vkeyword.clear();
  for (size_t i = 0; i < b.vkeyword.size(); i++) {
    vkeyword.push_back(b.vkeyword[i]);
  }
  vsponsor.clear();
  for (size_t i = 0; i < b.vsponsor.size(); i++) {
    vsponsor.push_back(b.vsponsor[i]);
  }
  valloy.clear();
  for (size_t i = 0; i < b.valloy.size(); i++) {
    valloy.push_back(b.valloy[i]);
  }
}

// copy
_outreach::_outreach(const _outreach& b) {
  copy(b);
}

// copies xtructures: b=a
const _outreach& _outreach::operator=(const _outreach& b) {  // operator=
  if (this != &b) {
    free();
    copy(b);
  }
  return *this;
}

void _outreach::clear() {
  const _outreach outreach_temp;
  copy(outreach_temp);
}

string _outreach::print_string() const {
  stringstream z;
  z << *this;
  return z.str();
}

// ******************************************************************************************************************************************************
// ******************************************************************************************************************************************************
// ARTICLES/PRESENTATIONS
// ******************************************************************************************************************************************************
// ******************************************************************************************************************************************************

uint fixlabel(const vector<string>& vlabel, string& labelfixed) {
  for (size_t j = 0; j < vlabel.size(); j += 2) {
    //    if(labelfixed==vlabel.at(j))
    aurostd::StringSubstInPlace(labelfixed, vlabel.at(j), vlabel.at(j + 1));
  }
  aurostd::StringSubstInPlace(labelfixed, "_WEB_PDF_", WEB_PDF);
  aurostd::StringSubstInPlace(labelfixed, "_WEB_DOI_", WEB_DOI);
  return true;
}

uint fixlabel(const vector<string>& vlabel, vector<string>& vlabelfixed) {
  for (size_t i = 0; i < vlabelfixed.size(); i++) {
    for (size_t j = 0; j < vlabel.size(); j += 2) {
      //      if(vlabelfixed[i]==vlabel.at(j))
      aurostd::StringSubstInPlace(vlabelfixed[i], vlabel.at(j), vlabel.at(j + 1));
    }
    aurostd::StringSubstInPlace(vlabelfixed[i], "_WEB_PDF_", WEB_PDF);
    aurostd::StringSubstInPlace(vlabelfixed[i], "_WEB_DOI_", WEB_DOI);
  }
  return vlabelfixed.size();
}

uint bibtex2file(string bibtex, string _authors, string _title, string journal, string volume, string issue, string pages, string year, string abstract, string doi, string bibfile) { // SC20201229
  stringstream bibcontent;
  // FIX THE AUTHORS
  string authors = _authors;
  const vector<string> vfix{"Buongiorno Nardelli", "van Roekeghem", "Aspuru-Guzik", "Hattrick-Simpers", "DeCost",      "de Coss",   "De Santo", "De Gennaro", "Al Rahal Al Orabi", "de Jong", "D'Amico",
                            "van der Zwaag",       "van de Walle",  "Di Stefano",   "Ojeda Mota",       "Simmons Jr.", "Mattos Jr."};
  for (size_t i = 0; i < vfix.size(); i++) {
    aurostd::StringSubstInPlace(authors, vfix[i], string("{" + vfix[i] + "}")); // FIX
  }
  aurostd::StringSubstInPlace(authors, ".", ".~");
  aurostd::StringSubstInPlace(authors, "~ ", " ");
  aurostd::StringSubstInPlace(authors, "~ ", " ");
  aurostd::StringSubstInPlace(authors, "~-", "-");
  aurostd::StringSubstInPlace(authors, ",", " and");
  aurostd::StringSubstInPlace(authors, "and and", "and");
  authors = aurostd::html2latex(authors);
  // FIX THE TITLES
  string title = "{" + _title + "}";
  // vfix = {"SnSe", "SnTe", "GeTe", "CoSb", "FePt", "AgPt", "LSO:Ce", "Eu:SrI", "Mo-Ru-Ta-W", "BaSnO", "Heusler", "CsI(Tl)", "NaI(Tl)", "LaBr", "L1", "IrV", "RhV", "Pt$_{8}$Ti", "Mg-B", "Fe:Mo:C", "Fe-C", "LiB", "MgB", "AFLUX", "LUX", "ACBN0", "PAOFLOW", "Ag", "Au", "Cu", "Mg", "Bi", "In", "Sb", "Fe", "Mo", "Kr", "Ar", "Ne", "Cs", "Xe"};
  for (size_t i = 0; i < vfix.size(); i++) {
    aurostd::StringSubstInPlace(title, vfix[i], string("{" + vfix[i] + "}")); // FIX
  }
  // NOW build
  bibcontent << "@article{" << bibtex << "," << endl;
  if (!authors.empty()) {
    bibcontent << " author={" << authors << "}"; // AUTHORS
  }
  if (!title.empty()) {
    bibcontent << "," << endl << " title={" << aurostd::html2latex(title) << "}"; // TITLE
  }
  if (!journal.empty()) {
    bibcontent << "," << endl << " journal={" << journal << "}"; // JOURNAL
  }
  if (!volume.empty()) {
    bibcontent << "," << "volume={" << volume << "}"; // VOLUME
  }
  if (!issue.empty()) {
    bibcontent << "," << "issue={" << issue << "}"; // ISSUE
  }
  if (!pages.empty()) {
    bibcontent << "," << "pages={" << pages << "}"; // PAGES
  }
  if (!year.empty()) {
    bibcontent << "," << "year={" << year << "}"; // YEAR
  }
  if (!abstract.empty()) {
    bibcontent << "," << endl << " abstract={" << abstract << "}"; // ABSTRACT
  }
  if (!doi.empty()) {
    bibcontent << "," << "doi={" << doi << "}"; // DOI
  }
  bibcontent << endl << "}" << endl;
  bibcontent << "% Automatically generated - AFLOW " << AFLOW_VERSION << endl;

  // save
  aurostd::stringstream2file(bibcontent, bibfile);
  return bibcontent.str().size();
}

ostream& operator<<(ostream& oss, const _outreach& outreach) {
  // ******************************************************************************************************************************************************
  const bool compact = true;
  stringstream newline;
  newline << endl;

  // ***************************************************************************
  // ARTICLE
  if (aurostd::substring2bool(outreach.type, "ARTICLE")) {
    // generate authors_json
    string authors_json;
    authors_json = "[";
    for (size_t iauth = 0; iauth < outreach.vauthor.size(); iauth++) {
      authors_json += "\"" + outreach.vauthor[iauth] + "\"";
      if (iauth != outreach.vauthor.size() - 1) {
        authors_json += ",";
      }
    }
    authors_json += "]";

    // generate authors_txt
    string authors_txt;
    authors_txt = "";
    for (size_t iauth = 0; iauth < outreach.vauthor.size(); iauth++) {
      authors_txt += outreach.vauthor[iauth];
      if (outreach.vauthor.size() == 2 && iauth == outreach.vauthor.size() - 2) {
        authors_txt += " and ";
      }
      if (outreach.vauthor.size() != 2 && iauth == outreach.vauthor.size() - 2) {
        authors_txt += ", and ";
      }
      if (iauth != outreach.vauthor.size() - 2 && iauth != outreach.vauthor.size() - 1) {
        authors_txt += ", ";
      }
    }

    string wnumber = aurostd::utype2string(outreach.wnumber);
    if (wnumber.length() < 2) {
      wnumber = "0" + wnumber;
    }

    //  generate journal HTML style
    string journal;
    if (outreach.bibtex_journal.empty()) { // with journal
      journal = outreach.journal; // with journal
    } else {
      journal = outreach.bibtex_journal;
      if (!outreach.bibtex_volume.empty()) {
        journal += " <b>" + outreach.bibtex_volume + "</b>";
      }
      if (!outreach.bibtex_issue.empty()) {
        journal += "(" + outreach.bibtex_issue + ")";
      }
      if (!outreach.bibtex_volume.empty()) {
        journal += ",";
      }
      if (!outreach.bibtex_pages.empty()) {
        journal += " " + outreach.bibtex_pages;
      }
      if (!outreach.bibtex_year.empty()) {
        journal += " (" + outreach.bibtex_year + ")";
      }
    }

    // ***************************************************************************
    // TXT
    if (XHOST.vflag_control.flag("PRINT_MODE::TXT")) {
      // 3rd line AUTHORS
      aurostd::StringSubstInPlace(authors_txt, "Csányi", "Csanyi");
      oss << aurostd::html2txt(authors_txt) << ", ";
      // 4th line TITLE
      oss << "\"" << aurostd::html2txt(outreach.title) << "\", ";
      // 5th line JOURNAL with year
      oss << "" << aurostd::html2txt(journal) << ". ";

      //   oss  << endl;
      // EXTRA STUFF FOR PHP WEB
      // 7th line DOI
      if (!outreach.doi.empty()) {
        oss << "doi: " << outreach.doi;
      }
      //     // 6th line PDF
      //     if(outreach.pdf.size())
      //       oss << "[<a href="+WEB_PDF+outreach.pdf << "><b>pdf</b></a>] " << endl;
      //     // 6th line SUPPLEMENTARY
      //     if(outreach.supplementary.size())
      //       oss << "[<a href="+WEB_PDF+outreach.supplementary << "><b>suppl</b></a>] " << endl;
      //     // 6th line ARXIV
      //     if(outreach.arxiv.size())
      //       oss << "[<a href="+outreach.arxiv << "><b>arxiv</b></a>] " << endl;
      //     // 6th line LINK
      //     if(outreach.link.size())
      //       oss << "[<a href="+outreach.link << "><b>link</b></a>] " << endl;
      //     // Xth line EXTRA
      //     for(size_t ix=0;ix<outreach.vextra_html.size();ix++)
      //       oss << outreach.vextra_html.at(ix) << endl;
      //     oss << "<br>";
      //     oss << endl; not needed in new SC20210616
    }
    // ***************************************************************************
    // JSON
    if (XHOST.vflag_control.flag("PRINT_MODE::JSON")) {
      oss << "{";
      if (outreach.newflag) {
        if (XHOST.vflag_control.flag("PRINT_MODE::NEW") && outreach.newflag) {
          oss << "\"new\":true,";
        }
        if (XHOST.vflag_control.flag("PRINT_MODE::DATE") && outreach.newflag && !outreach.date.empty()) {
          oss << R"("date":")" << outreach.date << "\",";
        }
      }
      if (XHOST.vflag_control.flag("PRINT_MODE::NUMBER")) {
        if (outreach.wnumber == 0) { // SHORTCUT for CHAPTER
          oss << "\"number\":" << outreach.wnumber << ",";
          oss << "\"chapter\":true,";
        } else {
          oss << "\"number\":" << outreach.wnumber << ",";
          oss << "\"article\":true,";
        }
      }
      // 3rd line AUTHORS
      aurostd::StringSubstInPlace(authors_json, "Csányi", "Csanyi");
      oss << "\"authors\":" << authors_json << ",";
      // 4th line TITLE
      oss << R"("title":")" << outreach.title << "\",";
      // 5th line JOURNAL with year
      oss << R"("journal":")" << journal << "\",";
      // EXTRA STUFF FOR PHP WEB
      // 7th line DOI
      if (XHOST.vflag_control.flag("PRINT_MODE::DOI") && !outreach.doi.empty()) {
        oss << R"("doi":")" << outreach.doi << "\",";
      }
      // 6th line PDF
      if (XHOST.vflag_control.flag("PRINT_MODE::PDF") && !outreach.pdf.empty()) {
        oss << R"("pdf":")" << WEB_PDF << outreach.pdf << "\",";
      }
      // 6th line SUPPLEMENTARY
      if (XHOST.vflag_control.flag("PRINT_MODE::PDF") && !outreach.supplementary.empty()) {
        oss << R"("supplementary":")" << WEB_PDF << outreach.supplementary << "\",";
      }
      // 6th line SUPPLEMENTARY_URL
      if (XHOST.vflag_control.flag("PRINT_MODE::DOI") && !outreach.supplementary_url.empty()) {
        oss << R"("supplementary_url":")" << outreach.supplementary_url << "\",";
      }
      // 6th line ARXIV
      if (!outreach.arxiv.empty()) {
        oss << R"("arxiv":")" << outreach.arxiv << "\",";
      }
      // 6th line LINK
      if (!outreach.link.empty()) {
        oss << R"("link":")" << outreach.link << "\",";
      }
      // 7th line BIBTEX
      if (XHOST.vflag_control.flag("PRINT_MODE::BIBTEX") && !outreach.bibtex.empty() && !outreach.bibtex_journal.empty() && !outreach.doi.empty() && !outreach.title.empty() &&
          !outreach.vauthor.empty()) { // needs to have bibtex pdf doi title and authors
        bibtex2file(outreach.bibtex, authors_txt, outreach.title, outreach.bibtex_journal, outreach.bibtex_volume, outreach.bibtex_issue, outreach.bibtex_pages, outreach.bibtex_year, outreach.abstract,
                    outreach.doi, string(WEB_BIBTEX_FILE + outreach.bibtex + ".txt"));
        oss << R"("bibtex":")" << WEB_BIBTEX << outreach.bibtex << ".txt" << "\",";
      }
      // Xth line EXTRA
      if (!outreach.vextra_html.empty()) {
        oss << R"("vextra_html":")";
        for (size_t ix = 0; ix < outreach.vextra_html.size(); ix++) {
          string extrahtml = outreach.vextra_html[ix];
          aurostd::StringSubstInPlace(extrahtml, "\"", "\\\"");
          if (!aurostd::substring2bool(outreach.vextra_html[ix], WEB_PDF) && !aurostd::substring2bool(outreach.vextra_html[ix], WEB_DOI)) {
            oss << extrahtml << (compact ? " " : "<br>");
          } else {
            if (aurostd::substring2bool(outreach.vextra_html[ix], WEB_DOI) && XHOST.vflag_control.flag("PRINT_MODE::DOI")) {
              oss << extrahtml << (compact ? " " : "<br>");
            } else {
              if (aurostd::substring2bool(outreach.vextra_html[ix], WEB_PDF) && XHOST.vflag_control.flag("PRINT_MODE::PDF")) {
                oss << extrahtml << (compact ? " " : "<br>");
              }
            }
          }
        }
        oss << "\",";
      }
      // Print year
      oss << "\"year\":" << outreach.year;
      oss << "}";
    }
    // ***************************************************************************
    //  HTML
    if (XHOST.vflag_control.flag("PRINT_MODE::HTML")) {
      oss << "<li>";
      if (outreach.newflag) {
        if (XHOST.vflag_control.flag("PRINT_MODE::NEW") && outreach.newflag) {
          oss << "<blink><b>NEW </b></blink>";
        }
        if (XHOST.vflag_control.flag("PRINT_MODE::DATE") && outreach.newflag && !outreach.date.empty()) {
          oss << "<b>" << outreach.date << "</b>";
        }
        oss << " "; // endl;
      }
      if (XHOST.vflag_control.flag("PRINT_MODE::NUMBER")) {
        if (outreach.wnumber == 0) {
          oss << "<span class=\"pubNumber\">" << "Chapter" << ". </span>";
        } else {
          oss << "<span class=\"pubNumber\">" << outreach.wnumber << ". </span>";
        }
      }
      // 3rd line AUTHORS
      aurostd::StringSubstInPlace(authors_txt, "Csányi", "Csanyi");
      oss << "<span class=\"pubAuthors\">" << authors_txt << "</span>, " << (compact ? " " : newline.str());
      // 4th line TITLE
      oss << "<span class=\"pubTitle\">" << outreach.title << "</span>, " << (compact ? " " : newline.str());
      // 5th line JOURNAL with year
      oss << "<span class=\"pubJournal\">" << journal << "</span>. " << (compact ? " " : newline.str());
      // EXTRA STUFF FOR PHP WEB
      // 7th line DOI
      if (XHOST.vflag_control.flag("PRINT_MODE::DOI") && !outreach.doi.empty()) {
        oss << "[<a href=" + WEB_DOI + outreach.doi << "><b>doi" << "=" << outreach.doi << "</b></a>] " << (compact ? " " : newline.str());
      }
      // 6th line PDF
      if (XHOST.vflag_control.flag("PRINT_MODE::PDF") && !outreach.pdf.empty()) {
        oss << "[<a href=" + WEB_PDF + outreach.pdf << "><b>pdf</b></a>] " << (compact ? " " : newline.str());
      }
      // 6th line SUPPLEMENTARY
      if (XHOST.vflag_control.flag("PRINT_MODE::PDF") && !outreach.supplementary.empty()) {
        oss << "[<a href=" + WEB_PDF + outreach.supplementary << "><b>suppl</b></a>] " << (compact ? " " : newline.str());
      }
      // 6th line SUPPLEMENTARY_URL
      if (XHOST.vflag_control.flag("PRINT_MODE::DOI") && !outreach.supplementary_url.empty()) {
        oss << "[<a href=" + outreach.supplementary_url << "><b>suppl</b></a>] " << (compact ? " " : newline.str());
      }
      // 6th line ARXIV
      if (!outreach.arxiv.empty()) {
        oss << "[<a href=" + outreach.arxiv << "><b>arxiv</b></a>] " << (compact ? " " : newline.str());
      }
      // 6th line LINK
      if (!outreach.link.empty()) {
        oss << "[<a href=" + outreach.link << "><b>link</b></a>] " << (compact ? " " : newline.str());
      }
      // 7th line BIBTEX
      if (XHOST.vflag_control.flag("PRINT_MODE::BIBTEX") && !outreach.bibtex.empty() && !outreach.bibtex_journal.empty() && !outreach.doi.empty() && !outreach.title.empty() &&
          !outreach.vauthor.empty()) { // needs to have bibtex pdf doi title and authors
        bibtex2file(outreach.bibtex, authors_txt, outreach.title, outreach.bibtex_journal, outreach.bibtex_volume, outreach.bibtex_issue, outreach.bibtex_pages, outreach.bibtex_year, outreach.abstract,
                    outreach.doi, string(WEB_BIBTEX_FILE + outreach.bibtex + ".txt"));
        oss << "[<a href=" << WEB_BIBTEX << outreach.bibtex << ".txt" << "><b>bibtex</b></a>] " << (compact ? " " : newline.str());
      }
      // Xth line EXTRA
      for (size_t ix = 0; ix < outreach.vextra_html.size(); ix++) {
        if (!aurostd::substring2bool(outreach.vextra_html[ix], WEB_PDF) && !aurostd::substring2bool(outreach.vextra_html[ix], WEB_DOI)) {
          oss << outreach.vextra_html[ix] << (compact ? " " : newline.str());
        } else {
          if (aurostd::substring2bool(outreach.vextra_html[ix], WEB_DOI) && XHOST.vflag_control.flag("PRINT_MODE::DOI")) {
            oss << outreach.vextra_html[ix] << (compact ? " " : newline.str());
          } else {
            if (aurostd::substring2bool(outreach.vextra_html[ix], WEB_PDF) && XHOST.vflag_control.flag("PRINT_MODE::PDF")) {
              oss << outreach.vextra_html[ix] << (compact ? " " : newline.str());
            }
          }
        }
      }
      oss << "</li>";
    }
    // ***************************************************************************
    // LATEX
    if (XHOST.vflag_control.flag("PRINT_MODE::LATEX")) {
      // 3rd line AUTHORS
      authors_txt = aurostd::html2latex(authors_txt) + ", ";
      oss << " " << authors_txt;
      // 4th line TITLE
      const string title = "\\textit{" + aurostd::html2latex(outreach.title) + "}, ";
      oss << " " << title;
      // 5th line JOURNAL with year
      oss << " " << aurostd::html2latex(journal) << ". ";
      // EXTRA STUFF FOR LATEX
      bool link = false;
      // 7th line DOI
      if (XHOST.vflag_control.flag("PRINT_MODE::DOI") && !outreach.doi.empty()) {
        string doi;
        doi = R"(\ifthenelse {\equal{\hyperlinks}{true}}{)";
        doi += R"({\newline \textsf{\href{http://dx.doi.org/)" + outreach.doi + "}{DOI: " + aurostd::html2latex(outreach.doi) + "}}}";
        doi += "}{{\\newline \\textsf{DOI: " + aurostd::html2latex(outreach.doi) + "}}}";
        oss << " " << doi;
        link = true;
      } else {
        if (outreach.wnumber != 0 && outreach.wnumber != 27 && outreach.wnumber != 14 && outreach.wnumber != 11 && outreach.wnumber != 8 && outreach.wnumber != 2 && outreach.wnumber != 1) {
          const string doi; // ="{\\newline \\textsf{\\href{http://dx.doi.org/}{DOI: N/A}}}";
          oss << " " << doi;
          link = true;
        }
      }
      // 7th line PDF
      if ((XHOST.vflag_control.flag("PRINT_MODE::PDF") || (!XHOST.vflag_control.flag("PRINT_MODE::PDF") && XHOST.vflag_control.flag("PRINT_MODE::DOI"))) && outreach.doi.empty() && !outreach.pdf.empty()) {
        string pdf;
        pdf = R"(\ifthenelse {\equal{\hyperlinks}{true}}{\textsf{[\href{)" + WEB_PDF + outreach.pdf + "}{pdf}]}}{}";
        oss << " " << pdf;
        link = true;
      }
      if (link == false) {
        ;
      }; // if(!XHOST.QUIET_CERR) cerr  << wnumber << endl;
      // LATEX EXTRA
      for (size_t ix = 0; ix < outreach.vextra_latex.size(); ix++) {
        if (!aurostd::substring2bool(outreach.vextra_latex[ix], WEB_PDF) && !aurostd::substring2bool(outreach.vextra_latex[ix], WEB_DOI)) {
          oss << " " << outreach.vextra_latex[ix];
        } else {
          if (aurostd::substring2bool(outreach.vextra_latex[ix], WEB_DOI) && XHOST.vflag_control.flag("PRINT_MODE::DOI")) {
            oss << " " << outreach.vextra_latex[ix];
          } else {
            if (aurostd::substring2bool(outreach.vextra_latex[ix], WEB_PDF) && XHOST.vflag_control.flag("PRINT_MODE::PDF")) {
              oss << " " << outreach.vextra_latex[ix];
            }
          }
        }
      }
      // oss << " % ART" << wnumber << " - aflow " << string(AFLOW_VERSION) << endl; OLD NEEDS endls SC20210616
      oss << " % ART" << wnumber << " - aflow " << string(AFLOW_VERSION); // << endl; OLD NEEDS endls SC20210616
    }
  }
  // ***************************************************************************
  // PRESENTATION
  if (aurostd::substring2bool(outreach.type, "PRESENTATION")) {
    if (XHOST.vflag_control.flag("PRINT_MODE::TXT")) {
      if (outreach._isinvited && outreach.type == "PRESENTATION_TALK") {
        oss << "Invited talk";
      }
      if (outreach._isinvited && outreach.type == "PRESENTATION_SEMINAR") {
        oss << "Invited seminar";
      }
      if (outreach._isinvited && outreach.type == "PRESENTATION_COLLOQUIUM") {
        oss << "Invited colloquium";
      }
      if (outreach._isinvited && outreach.type == "PRESENTATION_PLENARY") {
        oss << "Plenary Speaker";
      }
      if (outreach._isinvited && outreach.type == "PRESENTATION_KEYNOTE") {
        oss << "Keynote Speaker";
      }
      if (outreach._isinvited && outreach.type == "PRESENTATION_TUTORIAL") {
        oss << "Tutorial";
      }
      if (outreach._isinvited && outreach.type == "PRESENTATION_PANELIST") {
        oss << "Invited panelist";
      }
      if (outreach._isonline) {
        oss << " (online):";
      } else {
        oss << ":";
      }
      oss << " " << aurostd::latex2txt(outreach.title) << "; ";
      oss << "" << aurostd::latex2txt(outreach.place) << ", ";
      oss << "" << aurostd::latex2txt(outreach.date) << ". ";
      if (!outreach.link.empty()) {
        oss << "link: " << outreach.link;
      }
    }
    if (XHOST.vflag_control.flag("PRINT_MODE::HTML")) {
      if (outreach._isinvited && outreach.type == "PRESENTATION_TALK") {
        oss << "<b> Invited talk";
      }
      if (outreach._isinvited && outreach.type == "PRESENTATION_SEMINAR") {
        oss << "<b> Invited seminar";
      }
      if (outreach._isinvited && outreach.type == "PRESENTATION_COLLOQUIUM") {
        oss << "<b> Invited colloquium";
      }
      if (outreach._isinvited && outreach.type == "PRESENTATION_PLENARY") {
        oss << "<b> Plenary Speaker";
      }
      if (outreach._isinvited && outreach.type == "PRESENTATION_KEYNOTE") {
        oss << "<b> Keynote Speaker";
      }
      if (outreach._isinvited && outreach.type == "PRESENTATION_TUTORIAL") {
        oss << "<b> Tutorial";
      }
      if (outreach._isinvited && outreach.type == "PRESENTATION_PANELIST") {
        oss << "<b> Invited panelist";
      }
      if (outreach._isonline) {
        oss << " (online):</b>";
      } else {
        oss << ":</b>";
      }
      //    oss << "<br>" << endl;
      oss << "<i> " << aurostd::latex2html(outreach.title) << "</i>; " << endl;
      oss << "" << aurostd::latex2html(outreach.place) << ", " << endl;
      oss << "" << aurostd::latex2html(outreach.date) << ". "; // << endl;
      oss << "<br>" << endl;
      if (!outreach.link.empty()) {
        oss << "[<a href=" + outreach.link << "><b>link" << "=" << outreach.link << "</b></a>] " << (compact ? " " : newline.str());
        oss << "<br>" << endl;
      }
    }
    if (XHOST.vflag_control.flag("PRINT_MODE::LATEX")) {
      if (outreach._isinvited && outreach.type == "PRESENTATION_TALK") {
        oss << "\\textbf{Invited talk";
      }
      if (outreach._isinvited && outreach.type == "PRESENTATION_SEMINAR") {
        oss << "\\textbf{Invited seminar";
      }
      if (outreach._isinvited && outreach.type == "PRESENTATION_COLLOQUIUM") {
        oss << "\\textbf{Invited colloquium";
      }
      if (outreach._isinvited && outreach.type == "PRESENTATION_PLENARY") {
        oss << "\\textbf{Plenary Speaker";
      }
      if (outreach._isinvited && outreach.type == "PRESENTATION_KEYNOTE") {
        oss << "\\textbf{Keynote Speaker";
      }
      if (outreach._isinvited && outreach.type == "PRESENTATION_TUTORIAL") {
        oss << "\\textbf{Tutorial";
      }
      if (outreach._isinvited && outreach.type == "PRESENTATION_PANELIST") {
        oss << "\\textbf{Invited panelist";
      }
      if (outreach._isonline) {
        oss << " (online):}";
      } else {
        oss << ":}";
      }
      oss << endl;
      oss << "\\textit{" << outreach.title << "}; " << endl;
      oss << "" << outreach.place << ", " << endl;
      oss << "" << outreach.date << ". "; // << endl;
      if (!outreach.link.empty()) {
        oss << endl;
        oss << R"(\ifthenelse {\equal{\hyperlinks}{true}}{)";
        oss << R"({\newline \textsf{\href{)" << outreach.link << "}{LINK: " << aurostd::html2latex(outreach.link) << "}}}";
        oss << "}{{\\newline \\textsf{LINK: " << aurostd::html2latex(outreach.link) << "}}}";
      }
    }
  }
  return oss;
}

// ******************************************************************************************************************************************************
void voutreach_print_publications_EXTRA(ostream& oss);

vector<_outreach> voutreach_global_list;
uint voutreach_global_max_year = 0, voutreach_global_min_year = 9999;

// ******************************************************************************************************************************************************
// ******************************************************************************************************************************************************
// ARTICLES
// ******************************************************************************************************************************************************
// ******************************************************************************************************************************************************

uint voutreach_remove_duplicate(vector<_outreach>& voutreach) {
  const bool LOCAL_VERBOSE = true; // false;
  if (LOCAL_VERBOSE && !XHOST.QUIET_CERR) {
    cerr << "voutreach_remove_duplicate: voutreach.size()=" << voutreach.size() << endl;
  }
  for (size_t i = 0; i < voutreach.size(); i++) {
    for (size_t j = i + 1; j < voutreach.size(); j++) {
      if (i != j && voutreach[j].year == voutreach[i].year) { // same year
        // if(LOCAL_VERBOSE && !XHOST.QUIET_CERR) cerr << "voutreach_remove_duplicate: found same year: (i,j)=(" << i << "," << j << ")  " << voutreach.size() << endl;
        if (voutreach[j].vauthor.size() == voutreach[i].vauthor.size()) { // same number of authors
          // if(LOCAL_VERBOSE && !XHOST.QUIET_CERR) cerr << "voutreach_remove_duplicate: found same vauthor.size: (i,j)=(" << i << "," << j << ")  " << voutreach.size() << endl;
          if (voutreach[j].vauthor.at(0) == voutreach[i].vauthor.at(0)) { // same first author
            // if(LOCAL_VERBOSE && !XHOST.QUIET_CERR) cerr << "voutreach_remove_duplicate: found same vauthor.at(0): (i,j)=(" << i << "," << j << ")  " << voutreach.size() << endl;
            if (voutreach[j].vauthor.at(voutreach[j].vauthor.size() - 1) == voutreach[i].vauthor.at(voutreach[i].vauthor.size() - 1)) { // same last author
              // if(LOCAL_VERBOSE && !XHOST.QUIET_CERR) cerr << "voutreach_remove_duplicate: found same vauthor.at(N): (i,j)=(" << i << "," << j << ")  " << voutreach.size() << endl;
              vector<string> vwordi;
              aurostd::string2tokens(voutreach[i].title, vwordi, " ");
              vector<string> vwordj;
              aurostd::string2tokens(voutreach[j].title, vwordj, " ");
              if (vwordj.size() == vwordi.size()) { // same spaces
                // if(LOCAL_VERBOSE && !XHOST.QUIET_CERR) cerr << "voutreach_remove_duplicate: found same number of title words: (i,j)=(" << i << "," << j << ")  " << vwordj.size() << endl;
                if (aurostd::abs(voutreach[i].title.length() - voutreach[j].title.length()) < 3) {
                  if (LOCAL_VERBOSE && !XHOST.QUIET_CERR) {
                    cerr << "voutreach_remove_duplicate: found similar title length: (i,j)=(" << i << "," << j << ")  " << voutreach[i].title.length() << "," << voutreach[j].title.length() << endl;
                  }
                  { // ACT
                    if (LOCAL_VERBOSE && !XHOST.QUIET_CERR) {
                      cerr << "voutreach_remove_duplicate: found same year: (i,j)=(" << i << "," << j << ")  " << voutreach.size() << endl;
                    }
                    if (LOCAL_VERBOSE && !XHOST.QUIET_CERR) {
                      cerr << "voutreach_remove_duplicate: found same vauthor.size: (i,j)=(" << i << "," << j << ")  " << voutreach.size() << endl;
                    }
                    if (LOCAL_VERBOSE && !XHOST.QUIET_CERR) {
                      cerr << "voutreach_remove_duplicate: found same vauthor.at(0): (i,j)=(" << i << "," << j << ")  " << voutreach.size() << endl;
                    }
                    if (LOCAL_VERBOSE && !XHOST.QUIET_CERR) {
                      cerr << "voutreach_remove_duplicate: found same vauthor.at(N): (i,j)=(" << i << "," << j << ")  " << voutreach.size() << endl;
                    }
                    if (LOCAL_VERBOSE && !XHOST.QUIET_CERR) {
                      cerr << "voutreach_remove_duplicate: found same number of title words: (i,j)=(" << i << "," << j << ")  " << vwordj.size() << endl;
                    }
                    if (LOCAL_VERBOSE && !XHOST.QUIET_CERR) {
                      cerr << "voutreach_remove_duplicate: found similar title length: (i,j)=(" << i << "," << j << ")  " << voutreach[i].title.length() << "," << voutreach[j].title.length() << endl;
                    }
                    if (LOCAL_VERBOSE && !XHOST.QUIET_CERR) {
                      cerr << "voutreach_remove_duplicate: removing (i,j)=(" << i << "," << j << ")  title1=" << voutreach[j].title << " | " << voutreach[i].title << endl;
                    }
                    if (i != j) {
                      voutreach.erase(voutreach.begin() + j);
                      i = 0;
                      j = 0;
                    }
                    if (LOCAL_VERBOSE && !XHOST.QUIET_CERR) {
                      cerr << "voutreach_remove_duplicate: Article: found duplicate (i,j)=(" << i << "," << j << ")  " << voutreach.size() << endl;
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  return voutreach.size();
}

// ******************************************************************************************************************************************************
void voutreach_print_publications_EXTRA(ostream& oss) {
  oss << "<! *************************************************************************************************************>" << endl;
  oss << "<! *************************************************************************************************************>" << endl;
  oss << "" << endl;
  oss << "<img border=\"0\" width=100% height=3 src=http://" << XHOST.AFLOW_MATERIALS_SERVER << "/auro/images/line.gif>" << endl;
  oss << "" << endl;
  // **********************************************************************************************
  oss << "<li><span class=\"pubYears\">Patents</span></li>" << endl;
  oss << "<li>" << endl;
  oss << "<span class=\"pubAuthors\">G. Ceder, C. Fischer, K. Tibbetts, D. Morgan, and S. Curtarolo</span>, " << endl;
  oss << "<span class=\"pubTitle\">Systems and Methods for predicting materials properties</span>," << endl;
  oss << "<span class=\"pubJournal\">US Patent #7292958 (2007)</span>. [<a HREF=" << WEB_PDF << "PAT7292958.pdf>pdf</a>]" << endl;
  oss << "</li>" << endl;

  // **********************************************************************************************
  oss << "<li><span class=\"pubYears\">Theses</span></li>" << endl;
  oss << "<li>" << endl;
  oss << R"(<span class="pubNumber">4. </span><span class="pubAuthors">S. Curtarolo</span>,)" << endl;
  oss << "<span class=\"pubTitle\">Coarse-Graining and Data Mining Approaches to the Prediction of Structures and their Dynamics</span>, " << endl;
  oss << "<span class=\"pubJournal\">Ph.D., Massachusetts Institute of Technology (2003)</span>." << endl;
  oss << "[<a href=" << WEB_PDF << "scurtarolo4.pdf>pdf</a>]<br>" << endl;
  oss << "</li>" << endl;
  oss << "<li>" << endl;
  oss << R"(<span class="pubNumber">3. </span><span class="pubAuthors">S. Curtarolo</span>,)" << endl;
  oss << "<span class=\"pubTitle\">Adsorption Problems investigated with Computer Simulation</span>, " << endl;
  oss << "<span class=\"pubJournal\">M.S., Pennsylvania State University (1999)</span>." << endl;
  oss << "[<a href=" << WEB_PDF << "scurtarolo3.pdf>pdf</a>]<br>" << endl;
  oss << "</li>" << endl;
  oss << "<li>" << endl;
  oss << R"(<span class="pubNumber">2. </span><span class="pubAuthors">S. Curtarolo</span>,)" << endl;
  oss << "<span class=\"pubTitle\">Influenza della rugosita' sul prewetting di Neon su Magnesio</span>, " << endl;
  oss << "<span class=\"pubJournal\">Laurea in Fisica, Universita` di Padova (1998)</span>." << endl;
  oss << "[<a href=" << WEB_PDF << "scurtarolo2.pdf>pdf</a>] (in Italian). <br>" << endl;
  oss << "</li>" << endl;
  oss << "<li>" << endl;
  oss << R"(<span class="pubNumber">1. </span><span class="pubAuthors">S. Curtarolo</span>,)" << endl;
  oss << "<span class=\"pubTitle\">Approccio analitico e numerico allo studio degli adattatori dielettrici</span>, " << endl;
  oss << "<span class=\"pubJournal\">Laurea in Ingegneria Elettronica, Universita` di Padova (1995)</span>." << endl;
  oss << "[<!A href=" << WEB_PDF << "scurtarolo1.pdf>pdf, not available<!/a>] (in Italian). <br>" << endl;
  oss << "</li>" << endl;
}

// ***************************************************************************
void SystemReferences(const string& system_in, vector<uint>& voutreach_wnumber) { // ADD REFERENCES
  voutreach_wnumber.clear();

  const string system = aurostd::VASP_PseudoPotential_CleanName(system_in);
  vector<_outreach> voutreach;
  voutreach_load(voutreach, "PUBLICATIONS");
  // if(!XHOST.QUIET_CERR) cerr << "LOADED " << voutreach.size() << " articles" << endl;
  for (size_t iart = 0; iart < voutreach.size(); iart++) {
    if (!voutreach[iart].valloy.empty()) {
      for (size_t itoken = 0; itoken < voutreach[iart].valloy.size(); itoken++) {
        if (system == voutreach[iart].valloy[itoken]) {
          voutreach_wnumber.push_back(voutreach[iart].wnumber);
        }
      }
    }
  }
}

// ***************************************************************************

void SystemReferences(const string& system_in, vector<string>& vrefs, bool AUTHOR_ETAL) { // ADD REFERENCES
  vrefs.clear();
  vector<uint> voutreach_wnumber;
  SystemReferences(system_in, voutreach_wnumber);
  vector<_outreach> voutreach;
  voutreach_load(voutreach, "PUBLICATIONS");

  stringstream oss;

  for (size_t iarticle = 0; iarticle < voutreach_wnumber.size(); iarticle++) {
    for (size_t i = 0; i < voutreach.size(); i++) {
      if (voutreach_wnumber[iarticle] == voutreach[i].wnumber) {
        oss.clear();
        oss.str(std::string());
        // GOT THE ARTICLE
        //	if(iarticle+1<10)  oss << "<sup>&nbsp;&nbsp;" << iarticle+1 << "</sup> "; // LABEL
        //	if(iarticle+1>=10) oss << "<sup>" << iarticle+1 << "</sup> "; // LABEL
        // AUTHORS
        string authors;
        if (!AUTHOR_ETAL || voutreach[i].vauthor.size() <= 4) { // all authors OR <=4
          for (size_t iauth = 0; iauth < voutreach[i].vauthor.size(); iauth++) {
            authors += voutreach[i].vauthor[iauth];
            if (voutreach[i].vauthor.size() == 2 && iauth == voutreach[i].vauthor.size() - 2) {
              authors += " and ";
            }
            if (voutreach[i].vauthor.size() != 2 && iauth == voutreach[i].vauthor.size() - 2) {
              authors += ", and ";
            }
            if (iauth != voutreach[i].vauthor.size() - 2) {
              authors += ", ";
            }
          }
          //	aurostd::StringSubst(authors,"Csányi","Csanyi");
        } else { // first et al..
          authors += voutreach[i].vauthor.at(0) + " et al., ";
        }
        oss << authors; // << endl;
        // TITLE
        string title = voutreach[i].title;
        title = "" + title + ", ";
        oss << title;
        // JOURNAL with year
        const string journal = voutreach[i].journal;
        oss << "" << journal << ". ";
        vrefs.push_back(aurostd::html2txt(oss.str()));
        oss.clear();
        oss.str(std::string());
      }
    }
  }
}

// ******************************************************************************************************************************************************

bool SystemInSystems(const string& system, const string& systems) {
  vector<string> tokens;
  aurostd::string2tokens(systems, tokens, ",");
  for (size_t i = 0; i < tokens.size(); i++) {
    if (system == tokens[i]) {
      return true;
    }
  }
  return false;
}

bool SystemInSystems(const string& system, const vector<string>& vsystems) {
  for (size_t i = 0; i < vsystems.size(); i++) {
    if (system == vsystems[i]) {
      return true;
    }
  }
  return false;
}

// ******************************************************************************************************************************************************

bool AlloyAlphabeticLIBRARY(const string& s) {
  // never anymore this problem
  if (s == "MoMg" || s == "NaMg" || s == "NbMg" || s == "OsMg" || s == "PbMg" || s == "RbMg" || s == "RbPd" || s == "RbPt" || s == "ReMg" || s == "RhMg" || s == "RuMg" || s == "ScMg" || s == "SiPd") {
    return false;
  }
  if (s == "SiPt" || s == "SnMg" || s == "SnPd" || s == "SnPt" || s == "SrMg" || s == "SrPd" || s == "SrPt" || s == "TaMg" || s == "TiMg" || s == "VMg" || s == "WMg" || s == "YMg" || s == "ZnMg" || s == "ZrMg") {
    return false;
  }
  return true;
}

// ******************************************************************************************************************************************************
// ******************************************************************************************************************************************************
// ARTICLES/PRESENTATIONS
// ******************************************************************************************************************************************************
// ******************************************************************************************************************************************************

void voutreach_print(uint _mode, ostream& oss, string what2print) {
#ifdef MAX_YEAR_PRESENTATIONS
  const uint max_year_presentations = MAX_YEAR_PRESENTATIONS;
#else
  const std::time_t t = std::time(nullptr);
  const std::tm* now = std::localtime(&t);
  const uint max_year_presentations = 1900 + now->tm_year;
#endif
  string message;
  const aurostd::xoption vflag = XHOST.vflag_outreach;
  const bool LDEBUG = (false || XHOST.DEBUG);
  uint mode = _mode;

  // ******************************************************************************************************************************************************
  // ARTICLES
  if (aurostd::substring2bool(what2print, "ARTICLE") || aurostd::substring2bool(what2print, "PUBLICATION")) {
    vector<_outreach> voutreach;
    voutreach_load(voutreach, "PUBLICATIONS");
    const bool LOCAL_VERBOSE = true; // false;
    if (!XHOST.QUIET_CERR) {
      cerr << "LOADED " << voutreach.size() << " articles" << endl;
    }
    const bool HTRESOURCE_MODE_PHP_PREAMBLE = false;

    vector<string> voss;

    // if(LDEBUG && !XHOST.QUIET_CERR) cerr << voutreach.at(0);// << endl;
    // if(LDEBUG && !XHOST.QUIET_CERR) cerr << voutreach.at(1);// << endl;
    // if(LDEBUG && !XHOST.QUIET_CERR) cerr << voutreach.at(2);// << endl;
    // if(LDEBUG && !XHOST.QUIET_CERR) cerr << voutreach.at(3);// << endl;
    // if(LDEBUG && !XHOST.QUIET_CERR) cerr << voutreach.at(4);// << endl;
    // if(LDEBUG && !XHOST.QUIET_CERR) cerr << voutreach.at(5);// << endl;
    // if(LDEBUG && !XHOST.QUIET_CERR) cerr << voutreach.at(6);// << endl;
    // if(LDEBUG && !XHOST.QUIET_CERR) cerr << voutreach.at(7);// << endl;
    // if(LDEBUG && !XHOST.QUIET_CERR) cerr << voutreach.at(8);// << endl;

    vector<_outreach> voutreach_local;
    vector<string> vauthor;
    vector<string> vkeyword;
    vector<string> valloy;
    //  vector<string> vargument;
    const bool flag_simple = false;

    if (LOCAL_VERBOSE && !XHOST.QUIET_CERR) {
      cerr << "voutreach_print [begin]" << endl;
    }
    //   LDEBUG=true;

    if (LDEBUG) {
      oss << "XHOST.vflag_control.flags" << endl;
    }
    if (LDEBUG) {
      oss << "XHOST.vflag_control.flag(\"PRINT_MODE::JSON\")=" << XHOST.vflag_control.flag("PRINT_MODE::JSON") << endl;
    }
    if (LDEBUG) {
      oss << "XHOST.vflag_control.flag(\"PRINT_MODE::HTML\")=" << XHOST.vflag_control.flag("PRINT_MODE::HTML") << endl;
    }
    if (LDEBUG) {
      oss << "XHOST.vflag_control.flag(\"PRINT_MODE::TXT\")=" << XHOST.vflag_control.flag("PRINT_MODE::TXT") << endl;
    }
    if (LDEBUG) {
      oss << "XHOST.vflag_control.flag(\"PRINT_MODE::LATEX\")=" << XHOST.vflag_control.flag("PRINT_MODE::LATEX") << endl;
    }
    if (LDEBUG) {
      oss << "XHOST.vflag_control.flag(\"PRINT_MODE::YEAR\")=" << XHOST.vflag_control.flag("PRINT_MODE::YEAR") << endl;
    }
    if (LDEBUG) {
      oss << "XHOST.vflag_control.flag(\"PRINT_MODE::DOI\")=" << XHOST.vflag_control.flag("PRINT_MODE::DOI") << endl;
    }
    if (LDEBUG) {
      oss << "XHOST.vflag_control.flag(\"PRINT_MODE::BIBTEX\")=" << XHOST.vflag_control.flag("PRINT_MODE::BIBTEX") << endl;
    }
    if (LDEBUG) {
      oss << "XHOST.vflag_control.flag(\"PRINT_MODE::EXTRA\")=" << XHOST.vflag_control.flag("PRINT_MODE::EXTRA") << endl;
    }
    if (LDEBUG) {
      oss << "XHOST.vflag_control.flag(\"PRINT_MODE::NUMBER\")=" << XHOST.vflag_control.flag("PRINT_MODE::NUMBER") << endl;
    }
    if (LDEBUG) {
      oss << "XHOST.vflag_control.flag(\"PRINT_MODE::HYPERLINKS\")=" << XHOST.vflag_control.flag("PRINT_MODE::HYPERLINKS") << endl;
    }
    if (LDEBUG) {
      oss << "XHOST.vflag_control.flag(\"PRINT_MODE::NOTE\")=" << XHOST.vflag_control.flag("PRINT_MODE::NOTE") << endl;
    }
    if (LDEBUG) {
      oss << "XHOST.vflag_control.flag(\"PRINT_MODE::NEW\")=" << XHOST.vflag_control.flag("PRINT_MODE::NEW") << endl;
    }
    if (LDEBUG) {
      oss << "XHOST.vflag_control.flag(\"PRINT_MODE::DATE\")=" << XHOST.vflag_control.flag("PRINT_MODE::DATE") << endl;
    }

    // big hack in chinatown ! if Frisco gives me nothing...
    if (XHOST.vflag_control.flag("PHP::PUBS_ALLOY")) {
      if (XHOST.vflag_control.getattachedscheme("PHP::PUBS_ALLOY") == "--print=html") {
        mode = HTRESOURCE_MODE_PHP_ALLOY;
        XHOST.vflag_control.flag("PRINT_MODE::LATEX", false);
        XHOST.vflag_control.flag("PRINT_MODE::TXT", false);
        XHOST.vflag_control.flag("PRINT_MODE::JSON", false);
        XHOST.vflag_control.flag("PRINT_MODE::HTML", true); // DEFAULT
        valloy.emplace_back("AFLOW");
        valloy.emplace_back("AFLOWLIB");
        valloy.emplace_back("nmatHT");
      }
    }

    if (mode == HTRESOURCE_MODE_NONE) {
      mode = HTRESOURCE_MODE_PHP_AUTHOR;
    } // by default

    if (XHOST.vflag_control.flag("CV::AUTHOR") || mode == HTRESOURCE_MODE_PHP_AUTHOR) {
      mode = HTRESOURCE_MODE_PHP_AUTHOR;
      if (LOCAL_VERBOSE && !XHOST.QUIET_CERR) {
        cerr << "voutreach_print: [" << XHOST.vflag_control.flag("CV::AUTHOR") << "]" << endl;
      }
      if (XHOST.vflag_control.flag("CV::AUTHOR")) {
        if (LOCAL_VERBOSE && !XHOST.QUIET_CERR) {
          cerr << "voutreach_print: [" << XHOST.vflag_control.getattachedscheme("CV::AUTHOR") << "]" << endl;
        }
        aurostd::string2tokens(XHOST.vflag_control.getattachedscheme("CV::AUTHOR"), vauthor, ",");
        if (vauthor.empty()) {
          message = "No authors specified.";
          throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message);
        }
      }
    }

    if (LOCAL_VERBOSE && !XHOST.QUIET_CERR) {
      cerr << "voutreach_print: vauthor.size()=" << vauthor.size() << endl;
    }

    if (XHOST.vflag_control.flag("PHP::PUBS_ALLOY") || mode == HTRESOURCE_MODE_PHP_ALLOY) {
      if (XHOST.vflag_control.getattachedscheme("PHP::PUBS_ALLOY") != "--print=html") {
        mode = HTRESOURCE_MODE_PHP_ALLOY;
        aurostd::string2tokens(XHOST.vflag_control.getattachedscheme("PHP::PUBS_ALLOY"), valloy, ",");
        if (valloy.empty()) {
          message = "valloy.size() == 0";
          throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message);
        }
      }
    }
    if (XHOST.vflag_control.flag("PHP::PUBS_KEYWORD") || mode == HTRESOURCE_MODE_PHP_THRUST) {
      mode = HTRESOURCE_MODE_PHP_THRUST;
      aurostd::string2tokens(XHOST.vflag_control.getattachedscheme("PHP::PUBS_KEYWORD"), vkeyword, ",");
      if (vkeyword.empty()) {
        message = "No keywords specified.";
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message);
      }
      for (size_t i = 0; i < vkeyword.size(); i++) {
        if (!XHOST.QUIET_CERR) {
          cerr << "vkeyword.at(" << i << ")=" << vkeyword[i] << endl;
        }
      }
    }

    if (!XHOST.QUIET_CERR) {
      cerr << "voutreach_print: vauthor.size()=" << vauthor.size() << endl;
    }
    if (!XHOST.QUIET_CERR) {
      cerr << "voutreach_print: valloy.size()=" << valloy.size() << endl;
    }
    if (!XHOST.QUIET_CERR) {
      cerr << "voutreach_print: vkeyword.size()=" << vkeyword.size() << endl;
    }
    if (!XHOST.QUIET_CERR) {
      cerr << "voutreach_print: mode=" << mode << endl;
    }
    // ************************************************************************************************************************
    // ************************************************************************************************************************
    // OPERATION ON MULTIPLE VAUTHOR
    if (!vauthor.empty()) {
      // if(!XHOST.QUIET_CERR) cerr << "OPERATION ON MULTIPLE VAUTHOR" << endl;
      uint recent_articles = 0;
      for (size_t iauthor = 0; iauthor < vauthor.size(); iauthor++) {
        voutreach_local.clear();
        if (!XHOST.vflag_control.flag("PRINT_MODE::JSON")) {
          oss << endl; // just some  leeway but not in json SC20210616
        }
        if (HTRESOURCE_MODE_PHP_PREAMBLE == true && mode == HTRESOURCE_MODE_PHP_AUTHOR) {
          oss << "<?php" << endl;
          if (flag_simple == false) {
            oss << "if($author==\"" << vauthor[iauthor] << "\")" << endl;
          }
          oss << "{" << endl;
          oss << "?>" << endl;
        }
        if (XHOST.vflag_control.flag("PRINT_MODE::LATEX")) {
          oss << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl;
          oss << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl;
          // oss << "% PUBLICATIONS" << endl;
          oss << "" << endl;
          oss << "\\section{Publications} \\label{publications}" << endl;
          oss << "[Articles can be accessed through the embedded link. Only published, submitted and ``in press'' articles are listed. The list might be slightly out of order.]" << endl;

          // oss << "\\begin{list}{\\labelitemi}{\\leftmargin=1em}" << endl;
        }

        if (LOCAL_VERBOSE && !XHOST.QUIET_CERR) {
          cerr << "voutreach_print: LOADED " << voutreach.size() << " articles" << endl;
        }

        for (size_t i = 0; i < voutreach.size(); i++) {
          if (aurostd::substring2bool(voutreach[i].vauthor, vauthor[iauthor])) {
            voutreach_local.push_back(voutreach[i]); // push authors
          }
        }

        if (LOCAL_VERBOSE && !XHOST.QUIET_CERR) {
          cerr << "voutreach_print: LOADED_LOCAL " << voutreach_local.size() << " articles" << endl;
        }
        voutreach_remove_duplicate(voutreach_local);
        if (LOCAL_VERBOSE && !XHOST.QUIET_CERR) {
          cerr << "voutreach_print: LOADED_DUPLICATE " << voutreach_local.size() << " articles" << endl;
        }
        voutreach_rsort_wnumber(voutreach_local); // sort TOP numbers first
        if (LOCAL_VERBOSE && !XHOST.QUIET_CERR) {
          cerr << "voutreach_print: LOADED_SORTED " << voutreach_local.size() << " articles" << endl;
        }

        if (mode == HTRESOURCE_MODE_PHP_AUTHOR && !XHOST.vflag_control.flag("PRINT_MODE::LATEX") && !XHOST.vflag_control.flag("PRINT_MODE::TXT") && !XHOST.vflag_control.flag("PRINT_MODE::JSON")) {
          oss << "<br><br>" << endl;
        }
        if (XHOST.vflag_control.flag("PRINT_MODE::LATEX")) {
          oss << endl;
        }

        // ************************************************************************************************************************
        // HTML/JSON MODE
        if (XHOST.vflag_control.flag("PRINT_MODE::HTML") || XHOST.vflag_control.flag("PRINT_MODE::JSON")) {
          //	if(mode==HTRESOURCE_MODE_PHP_AUTHOR && !XHOST.vflag_control.flag("PRINT_MODE::LATEX")&& !XHOST.vflag_control.flag("PRINT_MODE::TXT"))
          if (XHOST.vflag_control.flag("PRINT_MODE::HTML")) {
            voss.push_back("<li><span class=\"aflow\"> AFLOW V" + string(AFLOW_VERSION) + " </span></li>");
          }
          if (XHOST.vflag_control.flag("PRINT_MODE::JSON")) {
            voss.emplace_back("[");
          }
          for (uint year = voutreach_global_max_year; year >= voutreach_global_min_year; year--) {
            bool flag_year = false;
            for (size_t i = 0; i < voutreach_local.size() && !flag_year; i++) {
              if (year == voutreach_local[i].year) {
                flag_year = true;
              }
            }
            if (flag_year && recent_articles < AUTHOR_RECENT_ARTICLES && year >= voutreach_global_max_year - AUTHOR_RECENT_YEARS && year <= voutreach_global_max_year) {
              if (XHOST.vflag_control.flag("PRINT_MODE::YEAR")) {
                if (XHOST.vflag_control.flag("PRINT_MODE::HTML")) {
                  voss.push_back("<li><span class=\"pubYears\">" + aurostd::utype2string(year) + "</span></li>");
                }
              }
              for (size_t i = 0; i < voutreach_local.size(); i++) {
                string wnumber = aurostd::utype2string(voutreach_local[i].wnumber);
                if (wnumber.length() < 3) {
                  wnumber = "0" + wnumber;
                }
                if (wnumber.length() < 2) {
                  wnumber = "0" + wnumber;
                }
                if (voutreach_local[i].year == year && recent_articles < AUTHOR_RECENT_ARTICLES && year >= voutreach_global_max_year - AUTHOR_RECENT_YEARS && year <= voutreach_global_max_year) {
                  voss.push_back(voutreach_local[i].print_string());
                  recent_articles++;
                }
              }
            }
          }
          if (XHOST.vflag_control.flag("PRINT_MODE::HTML")) { // this is only for HTML, json has all preables elsewhere
            if (HTRESOURCE_MODE_PHP_PREAMBLE == true && mode == HTRESOURCE_MODE_PHP_AUTHOR) {
              voss.emplace_back("\n <?php");
              voss.emplace_back("}");
              voss.emplace_back("?>");
            }
            for (size_t j = 0; j < voss.size(); j++) {
              oss << voss[j] << endl;
            }
          }
          if (XHOST.vflag_control.flag("PRINT_MODE::JSON")) {
            voss.emplace_back("]");
            for (size_t j = 0; j < voss.size(); j++) {
              oss << voss[j];
              if (j > 0 && j < voss.size() - 2) {
                oss << ",";
              }
              // oss << endl;
            }
            oss << endl;
          }
          if (XHOST.vflag_control.flag("PRINT_MODE::EXTRA")) {
            voutreach_print_publications_EXTRA(oss);
          }
        } // HTML/JSON MODE
        // ************************************************************************************************************************

        // ************************************************************************************************************************
        /* OLD MODE
        // LATEX MODE
        if(XHOST.vflag_control.flag("PRINT_MODE::LATEX")) {
        if(!XHOST.vflag_control.flag("PRINT_MODE::YEAR")) oss << "\\begin{itemize}" << endl;

        for(uint year=voutreach_global_max_year;year>=voutreach_global_min_year;year--) {
        bool flag_year=false;
        for(size_t i=0;i<voutreach_local.size()&&!flag_year;i++)
        if(year==voutreach_local.at(i).year)
        flag_year=true;
        if(flag_year && year<=voutreach_global_max_year) {
        oss << endl;
        if(XHOST.vflag_control.flag("PRINT_MODE::YEAR")) oss << "\\ \\\\  \\textbf{ \\color{blue}{[" << year << "]}}" << endl;
        if(XHOST.vflag_control.flag("PRINT_MODE::YEAR")) oss << "\\begin{itemize}" << endl;
        for(size_t i=0;i<voutreach_local.size();i++) {
        string wnumber=aurostd::utype2string(voutreach_local.at(i).wnumber);
        if(wnumber.length()<3) wnumber="0"+wnumber;
        if(wnumber.length()<2) wnumber="0"+wnumber;
        if(voutreach_local.at(i).year==year && year>=voutreach_global_max_year-AUTHOR_RECENT_YEARS && year<=voutreach_global_max_year) {
        // print 1 article
        if(!XHOST.vflag_control.flag("PRINT_MODE::NUMBER")) {
        oss << "" << "\\item[$\\bullet$]{}";// << "% ART" << wnumber << " - aflow " << string(AFLOW_VERSION);
        } else {
        if(voutreach_local.at(i).wnumber==0) {
        oss << "" << "\\item[\\textbf{Chapter.\\,}]{}";// << "% ART" << wnumber << " - aflow " << string(AFLOW_VERSION);
        } else {
        oss << "" << "\\item[\\textbf{"+aurostd::utype2string(voutreach_local.at(i).wnumber)+".\\,}]{}";// << "% ART" << wnumber << " - aflow " << string(AFLOW_VERSION);
        }
        }
        // oss << endl;
        oss << voutreach_local.at(i);// << endl;
        }
        }
        if(XHOST.vflag_control.flag("PRINT_MODE::YEAR")) oss << "\\end{itemize}" << endl;
        }
        }
        oss << "" << endl;
        if(!XHOST.vflag_control.flag("PRINT_MODE::YEAR")) oss << "\\end{itemize}" << endl;
        //    oss << "\\end{list}" << endl;
        oss << "" << endl;
        oss << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl;
        oss << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl;
        } // LATEX MODE
        */
        // ************************************************************************************************************************
        // LATEX MODE
        if (XHOST.vflag_control.flag("PRINT_MODE::LATEX")) {
          if (!XHOST.vflag_control.flag("PRINT_MODE::YEAR")) {
            voss.emplace_back("\\begin{itemize}");
          }
          for (uint year = voutreach_global_max_year; year >= voutreach_global_min_year; year--) {
            bool flag_year = false;
            for (size_t i = 0; i < voutreach_local.size() && !flag_year; i++) {
              if (year == voutreach_local[i].year) {
                flag_year = true;
              }
            }
            if (flag_year && year <= voutreach_global_max_year) {
              voss.emplace_back("");
              if (XHOST.vflag_control.flag("PRINT_MODE::YEAR")) {
                voss.push_back(R"(\ \\  \textbf{ \color{blue}{[)" + aurostd::utype2string(year) + "]}}");
              }
              if (XHOST.vflag_control.flag("PRINT_MODE::YEAR")) {
                voss.emplace_back("\\begin{itemize}");
              }
              for (size_t i = 0; i < voutreach_local.size(); i++) {
                string wnumber = aurostd::utype2string(voutreach_local[i].wnumber);
                if (wnumber.length() < 3) {
                  wnumber = "0" + wnumber;
                }
                if (wnumber.length() < 2) {
                  wnumber = "0" + wnumber;
                }
                if (voutreach_local[i].year == year && year >= voutreach_global_max_year - AUTHOR_RECENT_YEARS && year <= voutreach_global_max_year) {
                  // print 1 article
                  if (!XHOST.vflag_control.flag("PRINT_MODE::NUMBER")) {
                    voss.push_back("\\item[$\\bullet$]{}" + voutreach_local[i].print_string());
                  } else {
                    if (voutreach_local[i].wnumber == 0) {
                      voss.push_back(R"(\item[\textbf{Chapter.\,}]{})" + voutreach_local[i].print_string());
                    } else {
                      voss.push_back("\\item[\\textbf{" + aurostd::utype2string(voutreach_local[i].wnumber) + ".\\,}]{}" + voutreach_local[i].print_string());
                    }
                  }
                }
              }
              if (XHOST.vflag_control.flag("PRINT_MODE::YEAR")) {
                voss.emplace_back("\\end{itemize}");
              }
            }
          }
          voss.emplace_back("");
          if (!XHOST.vflag_control.flag("PRINT_MODE::YEAR")) {
            voss.emplace_back("\\end{itemize}");
          }
          //    voss.push_back("\\end{list}");
          voss.emplace_back("");
          voss.emplace_back("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%");
          voss.emplace_back("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%");
          for (size_t j = 0; j < voss.size(); j++) {
            oss << voss[j] << endl;
          }
        } // LATEX MODE
        // ************************************************************************************************************************

        // ************************************************************************************************************************
        // OLD MODE
        /* TXT MODE
           if(XHOST.vflag_control.flag("PRINT_MODE::TXT")) {
           for(uint year=voutreach_global_max_year;year>=voutreach_global_min_year;year--) {
           bool flag_year=false;
           for(size_t i=0;i<voutreach_local.size()&&!flag_year;i++)
           if(year==voutreach_local.at(i).year)
           flag_year=true;
           if(flag_year && year<=voutreach_global_max_year) {
           oss << endl;
           if(XHOST.vflag_control.flag("PRINT_MODE::YEAR")) oss << "[" << year << "]" << endl;
           for(size_t i=0;i<voutreach_local.size();i++) {
           string wnumber=aurostd::utype2string(voutreach_local.at(i).wnumber);
           if(wnumber.length()<3) wnumber="0"+wnumber;
           if(wnumber.length()<2) wnumber="0"+wnumber;
           if(voutreach_local.at(i).year==year && year>=voutreach_global_max_year-AUTHOR_RECENT_YEARS && year<=voutreach_global_max_year) {
        // print 1 article
        if(XHOST.vflag_control.flag("PRINT_MODE::NUMBER")) {
        oss << "" << "["+aurostd::utype2string(voutreach_local.at(i).wnumber)+".] ";
        }
        oss << voutreach_local.at(i);
        }
        }
        }
        }
        oss << "" << endl;
        oss << "" << endl;
        oss << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl;
        oss << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl;
        } // TXT MODE
        */
        // ************************************************************************************************************************

        // ************************************************************************************************************************
        // TXT MODE
        if (XHOST.vflag_control.flag("PRINT_MODE::TXT")) {
          for (uint year = voutreach_global_max_year; year >= voutreach_global_min_year; year--) {
            bool flag_year = false;
            for (size_t i = 0; i < voutreach_local.size() && !flag_year; i++) {
              if (year == voutreach_local[i].year) {
                flag_year = true;
              }
            }
            if (flag_year && year <= voutreach_global_max_year) {
              if (XHOST.vflag_control.flag("PRINT_MODE::YEAR")) {
                voss.emplace_back("");
                voss.push_back("[" + aurostd::utype2string(year) + "]");
              }
              for (size_t i = 0; i < voutreach_local.size(); i++) {
                string wnumber = aurostd::utype2string(voutreach_local[i].wnumber);
                if (wnumber.length() < 3) {
                  wnumber = "0" + wnumber;
                }
                if (wnumber.length() < 2) {
                  wnumber = "0" + wnumber;
                }
                if (voutreach_local[i].year == year && year >= voutreach_global_max_year - AUTHOR_RECENT_YEARS && year <= voutreach_global_max_year) {
                  // print 1 article
                  if (XHOST.vflag_control.flag("PRINT_MODE::NUMBER")) {
                    voss.push_back("[" + aurostd::utype2string(voutreach_local[i].wnumber) + ".] " + voutreach_local[i].print_string());
                  } else {
                    voss.push_back(voutreach_local[i].print_string());
                  }
                }
              }
            }
          }
          voss.emplace_back("");
          voss.emplace_back("");
          voss.emplace_back("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%");
          voss.emplace_back("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%");
          for (size_t j = 0; j < voss.size(); j++) {
            oss << voss[j] << endl;
          }
        } // TXT MODE
        // ************************************************************************************************************************
      } // iauthor
    }
    // ************************************************************************************************************************
    // ************************************************************************************************************************
    // OPERATION ON MULTIPLE VTHRUST
    if (!vkeyword.empty()) {
      // if(!XHOST.QUIET_CERR) cerr << "OPERATION ON MULTIPLE VTHRUST" << endl;
      uint recent_articles = 0;
      voutreach_local.clear();
      for (size_t ithrust = 0; ithrust < vkeyword.size(); ithrust++) {
        oss << endl;
        if (HTRESOURCE_MODE_PHP_PREAMBLE == true && mode == HTRESOURCE_MODE_PHP_THRUST) {
          oss << "<?php" << endl;
          if (flag_simple == false) {
            oss << "if($thrust==\"" << vkeyword[ithrust] << "\")" << endl;
          }
          oss << "{" << endl;
          oss << "?>" << endl;
        }
        for (size_t i = 0; i < voutreach.size(); i++) {
          // if(!XHOST.QUIET_CERR) cerr << "voutreach.at(" << i << ").vkeyword.size()=" << voutreach.at(i).vkeyword.size() << endl;
          for (size_t j = 0; j < voutreach[i].vkeyword.size(); j++) {
            if (voutreach[i].vkeyword[j] == vkeyword[ithrust]) {
              voutreach_local.push_back(voutreach[i]); // push thrust
            }
          }
        }
      } // ithrust
      voutreach_remove_duplicate(voutreach_local);
      voutreach_sort_year(voutreach_local); // sort TOP numbers first

      //    if(mode==HTRESOURCE_MODE_PHP_THRUST) {oss << "<!br><!br>" << endl;}
      // ************************************************************************************************************************
      // HTML MODE
      if (XHOST.vflag_control.flag("PRINT_MODE::HTML")) {
        if (mode == HTRESOURCE_MODE_PHP_THRUST) {
          if (LOCAL_VERBOSE && !XHOST.QUIET_CERR) {
            oss << "<li><span class=\"aflow\">" << "AFLOW V" << string(AFLOW_VERSION) << " </span></li>" << endl;
          }
          for (uint year = voutreach_global_max_year; year >= voutreach_global_min_year; year--) {
            bool flag_year = false;
            for (size_t i = 0; i < voutreach_local.size() && !flag_year; i++) {
              if (year == voutreach_local[i].year) {
                flag_year = true;
              }
            }
            if (flag_year && recent_articles < THRUST_RECENT_ARTICLES && year >= voutreach_global_max_year - THRUST_RECENT_YEARS && year <= voutreach_global_max_year) {
              for (size_t i = 0; i < voutreach_local.size(); i++) {
                string wnumber = aurostd::utype2string(voutreach_local[i].wnumber);
                if (wnumber.length() < 3) {
                  wnumber = "0" + wnumber;
                }
                if (wnumber.length() < 2) {
                  wnumber = "0" + wnumber;
                }
                if (voutreach_local[i].year == year && recent_articles < THRUST_RECENT_ARTICLES && year >= voutreach_global_max_year - THRUST_RECENT_YEARS && year <= voutreach_global_max_year) {
                  XHOST.vflag_control.flag("PRINT_MODE::NUMBER", false);
                  oss << voutreach_local[i] << endl;
                  recent_articles++;
                }
              }
            }
          }
          if (HTRESOURCE_MODE_PHP_PREAMBLE == true && mode == HTRESOURCE_MODE_PHP_THRUST) {
            oss << "<?php" << endl;
            oss << "}" << endl;
            oss << "?>" << endl;
          }
        }
      } // HTML MODE
      // ************************************************************************************************************************
    }
    // ************************************************************************************************************************
    // ************************************************************************************************************************
    // OPERATION ON MULTIPLE VALLOY
    if (!valloy.empty()) {
      // if(!XHOST.QUIET_CERR) cerr << "OPERATION ON MULTIPLE VALLOY" << endl;
      //   uint recent_articles=0;
      voutreach_local.clear();
      for (size_t ialloy = 0; ialloy < valloy.size(); ialloy++) {
        vector<uint> voutreach_wnumber;
        SystemReferences(aurostd::VASP_PseudoPotential_CleanName(valloy[ialloy]), voutreach_wnumber);
        for (size_t iwnumber = 0; iwnumber < voutreach_wnumber.size(); iwnumber++) {
          for (size_t iarticle = 0; iarticle < voutreach.size(); iarticle++) {
            if (voutreach[iarticle].wnumber == voutreach_wnumber[iwnumber]) {
              //	  if(voutreach.at(iart)!=65 && voutreach.at(iart)!=75)     // remove aflow and aflowlib papers
              voutreach_local.push_back(voutreach[iarticle]); // remove aflow and aflowlib papers
            }
          }
        }
        // if(!XHOST.QUIET_CERR) cerr << voutreach_local.size() << endl;

      } // ialloy
      voutreach_remove_duplicate(voutreach_local);
      voutreach_sort_year(voutreach_local); // sort TOP numbers first

      //    if(mode==HTRESOURCE_MODE_PHP_ALLOY) {oss << "<!br><!br>" << endl;}
      // ************************************************************************************************************************
      // PHP AND WEB MODE
      if (mode == HTRESOURCE_MODE_PHP_ALLOY) {
        if (LOCAL_VERBOSE && !XHOST.QUIET_CERR) {
          if (XHOST.vflag_control.flag("PRINT_MODE::HTML") || mode == HTRESOURCE_MODE_PHP_AUTHOR) {
            oss << "<li><span class=\"aflow\">" << "AFLOW V" << string(AFLOW_VERSION) << " </span></li>" << endl;
          }
        }
        oss << "<ol class=\"reference\">" << endl;
        for (size_t i = 0; i < voutreach_local.size(); i++) {
          //	string wnumber=voutreach_local[i].wnumber;
          //	if(wnumber.length()<2) wnumber="0"+wnumber;
          XHOST.vflag_control.flag("PRINT_MODE::NUMBER", false);
          voutreach_local[i].vextra_html.clear(); // delete it
          oss << voutreach_local[i] << endl;
          // recent_articles++;
        }
        oss << "</ol>" << endl;
        if (HTRESOURCE_MODE_PHP_PREAMBLE == true && mode == HTRESOURCE_MODE_PHP_ALLOY) {
          oss << "<?php" << endl;
          oss << "}" << endl;
          oss << "?>" << endl;
        }
      } // PHP WEB
      // ************************************************************************************************************************
    }

    //   for(size_t i=0;i<voutreach_local.size();i++)
    //   oss << voutreach_local.at(i) << endl;
  }

  // ******************************************************************************************************************************************************
  // PRESENTATIOS
  if (aurostd::substring2bool(what2print, "PRESENTATION")) {
    vector<_outreach> voutreach;
    voutreach_load(voutreach, "PRESENTATIONS");
    if (!XHOST.QUIET_CERR) {
      cerr << "LOADED " << voutreach.size() << " presentations" << endl;
    }

    if (mode == HTRESOURCE_MODE_NONE) {
      mode = HTRESOURCE_MODE_PHP_AUTHOR;
    }

    vector<_outreach> voutreach_local;
    vector<string> vauthors;
    vauthors.emplace_back("Curtarolo");

    uint noutreach = 0;

    // ************************************************************************************************************************
    // TXT MODE
    if (XHOST.vflag_control.flag("PRINT_MODE::TXT")) {
      for (size_t iauthor = 0; iauthor < vauthors.size(); iauthor++) {
        voutreach_local.clear();

        oss << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl;
        for (size_t i = 0; i < voutreach.size(); i++) {
          if (aurostd::substring2bool(voutreach[i].vauthor, vauthors[iauthor])) {
            voutreach_local.push_back(voutreach[i]);
          }
        }
        for (uint year = max_year_presentations; year >= 1998; year--) {
          bool flag_year = false;
          for (size_t i = 0; i < voutreach_local.size() && !flag_year; i++) {
            if (year == voutreach_local[i].year) {
              flag_year = true;
            }
          }
          if (flag_year) {
            for (size_t i = 0; i < voutreach_local.size(); i++) {
              if (voutreach_local[i].year == year) {
                // print 1 outreach
                oss << endl;
                oss << voutreach_local[i] << endl;
              }
            }
          }
        }
        oss << "" << endl;
      } // iauthor
    } // TXT MODE
    // ************************************************************************************************************************
    // LATEX MODE
    if (XHOST.vflag_control.flag("PRINT_MODE::LATEX")) {
      for (size_t iauthor = 0; iauthor < vauthors.size(); iauthor++) {
        voutreach_local.clear();

        oss << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl;
        oss << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl;
        oss << "% TALKS" << endl;
        oss << "" << endl;
        oss << "\\section{Invited Talks, Seminars, Plenaries}" << endl;
        oss << "\\label{italks}" << endl;

        oss << "" << endl;
        oss << "\\begin{itemize}" << endl;
        // oss << "\\begin{list}{\\labelitemi}{\\leftmargin=1em}" << endl;

        for (size_t i = 0; i < voutreach.size(); i++) {
          if (aurostd::substring2bool(voutreach[i].vauthor, vauthors[iauthor])) {
            voutreach_local.push_back(voutreach[i]);
          }
        }

        for (uint year = max_year_presentations; year >= 1998; year--) {
          bool flag_year = false;
          for (size_t i = 0; i < voutreach_local.size() && !flag_year; i++) {
            if (year == voutreach_local[i].year) {
              flag_year = true;
            }
          }
          if (flag_year) {
            //	oss << endl << "<b><font size=\"3\"><font COLOR=blue><i>" << year << ".</i> </font></b><br><br>" << endl;
            for (size_t i = 0; i < voutreach_local.size(); i++) {
              if (voutreach_local[i].year == year) {
                // print 1 outreach
                oss << endl;
                if (vauthors[iauthor] == "Curtarolo") {
                  oss << "\\item[\\textbf{" << voutreach_local.size() - noutreach++ << ".\\,}]{}" << endl;
                }
                oss << voutreach_local[i] << endl;
              }
            }
          }
        }
        oss << "\\end{itemize}" << endl;
        // 	oss << "\\end{list}" << endl;

        oss << "" << endl;
        oss << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl;
        oss << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl;
      } // iauthor
    } // TXT MODE
    // ************************************************************************************************************************
    // HTML MODE
    if (XHOST.vflag_control.flag("PRINT_MODE::HTML")) {
      if (mode == HTRESOURCE_MODE_PHP_AUTHOR) {
        for (size_t iauthor = 0; iauthor < vauthors.size(); iauthor++) {
          voutreach_local.clear();

          for (size_t i = 0; i < voutreach.size(); i++) {
            if (aurostd::substring2bool(voutreach[i].vauthor, vauthors[iauthor])) {
              voutreach_local.push_back(voutreach[i]);
            }
          }

          for (uint year = max_year_presentations; year >= 1998; year--) {
            bool flag_year = false;
            for (size_t i = 0; i < voutreach_local.size() && !flag_year; i++) {
              if (year == voutreach_local[i].year) {
                flag_year = true;
              }
            }
            if (flag_year) {
              //	oss << endl << "<b><font size=\"3\"><font COLOR=blue><i>" << year << ".</i> </font></b><br><br>" << endl;
              for (size_t i = 0; i < voutreach_local.size(); i++) {
                if (voutreach_local[i].year == year) {
                  // print 1 outreach
                  oss << endl;
                  if (vauthors[iauthor] == "Curtarolo") {
                    oss << "<font COLOR=red><b>" << voutreach_local.size() - noutreach++ << ". </b></font>";
                  }
                  XHOST.vflag_control.flag("PRINT_MODE::LATEX", false);
                  XHOST.vflag_control.flag("PRINT_MODE::TXT", false);
                  XHOST.vflag_control.flag("PRINT_MODE::JSON", false);
                  XHOST.vflag_control.flag("PRINT_MODE::HTML", true);
                  oss << voutreach_local[i] << endl;
                }
              }
            }
          }
        }
      }
    } // HTML MODE
  } // PRESENTATIONS
}

// ******************************************************************************************************************************************************
void voutreach_print_everything(ostream& oss, const vector<string>& vitems, string msg1, string msg2, string sectionlabel) {
  const bool flag_ITMZ = false;
  const bool flag_LIST = true;
  if (!XHOST.QUIET_CERR) {
    cerr << "LOADED " << vitems.size() << " " << msg1 << endl;
  }

  if (XHOST.vflag_control.flag("PRINT_MODE::LATEX")) {
    oss << "% " << msg2 << endl;
    oss << sectionlabel << endl;
    if (flag_ITMZ) {
      oss << "\\begin{itemize}" << endl;
    }
    if (flag_LIST) {
      oss << R"(\begin{list}{\labelitemi}{\leftmargin=1em})" << endl;
    }
    for (size_t i = 0; i < vitems.size(); i++) {
      oss << "\\item{}" << " " << vitems[i] << "  % " << msg2 << endl;
    }
    if (flag_ITMZ) {
      oss << "\\end{itemize}" << endl;
    }
    if (flag_LIST) {
      oss << "\\end{list}" << endl;
    }
    //    oss << "" << endl;
  }
}

// ******************************************************************************************************************************************************
vector<_outreach> voutreach_presentations;
vector<_outreach> voutreach_publications;
int voutreach_call = 0;
uint voutreach_load(vector<_outreach>& voutreach, const string& what2print) {
  // if(!XHOST.QUIET_CERR) cerr << voutreach_call++ << endl;
  // check cache
  const bool LDEBUG = false;
  voutreach.clear();
  vector<string> valabel;
  vector<string> vjlabel;
  if (aurostd::substring2bool(what2print, "ARTICLE") || aurostd::substring2bool(what2print, "PUBLICATION")) {
    if (!voutreach_publications.empty()) {
      for (const _outreach& voutreach_publication : voutreach_publications) {
        voutreach.push_back(voutreach_publication);
      }
      return voutreach.size();
    }
  }
  if (aurostd::substring2bool(what2print, "PRESENTATION")) {
    if (!voutreach_presentations.empty()) {
      for (const _outreach& voutreach_presentation : voutreach_presentations) {
        voutreach.push_back(voutreach_presentation);
      }
      return voutreach.size();
    }
  }
  // no cache
  // <b><font size="3"><font COLOR=blue><i>2010-current.</i> </font></b><br><br>
  _outreach ptmp;

  if (LDEBUG) {
    cerr << "voutreach_load [0]" << endl;
  }
  // collect all txt files in the CV folder and concatenate them
  const fs::path cv_dir(XHOST.vflag_control.getattachedscheme("CV::DIR"));
  std::vector<fs::path> cv_files;
  for (const auto& entry : fs::directory_iterator(cv_dir)) {
    if (aurostd::tolower(entry.path().extension()) == ".txt") {
      cv_files.push_back(entry.path());
    }
  }
  std::string raw_cv_data;

  for (const auto& file : cv_files) {
    raw_cv_data += aurostd::file2string(file) + "\n";
  }
  raw_cv_data = aurostd::RemoveComments(raw_cv_data);

  if (LDEBUG) {
    cerr << "voutreach_load [1]" << endl;
  }
  vector<string> vpres;
  vector<string> ktokens;
  vector<string> tokens;
  std::vector<std::string> allowed_what2print = {"ARTICLE",        "ARTICLES",      "PUBLICATION", "PUBLICATIONS", "PRESENTATION", "PRESENTATIONS", "EDUCATION", "RESEARCH", "ACADEMIC",
                                                 "SERVICEOUTSIDE", "SERVICEINSIDE", "TEACHING",    "ADVISING",     "PATENT",       "PATENTS",       "PRESS",     "AWARD",    "AWARDS"};

  if (std::find(allowed_what2print.begin(), allowed_what2print.end(), what2print) != allowed_what2print.end()) {
    aurostd::string2vectorstring(raw_cv_data, vpres);
  } else {
    return 0;
  }

  if (LDEBUG) {
    cerr << "voutreach_load [2]" << endl;
  }

  string iline;
  string jline;
  string kline;
  size_t i;
  size_t j;
  size_t k;
  size_t l;

  for (i = 0; i < vpres.size(); i++) {
    //  cout << vpres[i] << endl;

    iline = vpres[i];

    if (LDEBUG) {
      cerr << "voutreach_load [2.line]=" << iline << endl;
    }

    aurostd::StringSubstInPlace(iline, " ", "");
    if (iline == "OBJECT={") { // found an object
      for (j = i; i < vpres.size(); j++) {
        jline = vpres.at(j);
        aurostd::StringSubstInPlace(jline, " ", "");
        if (jline.empty()) {
          continue;
        }
        if (jline.at(0) == '}') {
          break;
        }
      }
      // now I have an object between i+1 and j-1
      ptmp.clear();
      for (k = i + 1; k < j; k++) {
        kline = vpres.at(k);
        aurostd::StringSubstInPlace(kline, " ", "");
        aurostd::string2tokens(kline, ktokens, "=");
        aurostd::string2tokens(vpres.at(k), tokens, "=");
        // if(!XHOST.QUIET_CERR) cerr << tokens.at(1) << endl;
        if (!tokens.empty() && !ktokens.empty()) {
          string object;
          string _ktoken0_aus;
          object = tokens.at(1);
          for (l = 2; l < tokens.size(); l++) {
            object += "=" + tokens[l];
          }

          // clean the start and end of the string
          size_t clean_start = 0;
          size_t clean_stop = object.size();
          for (size_t pos = 0; pos < object.size(); pos++) {
            if (object[pos] == ' ' || object[pos] == '\"') {
              clean_start = pos + 1;
            } else {
              break;
            }
          }
          for (size_t pos = object.size() - 1; pos > clean_start; pos--) {
            if (object[pos] == ' ' || object[pos] == ';' || object[pos] == '\"') {
              clean_stop = pos;
            } else {
              break;
            }
          }
          object = object.substr(clean_start, clean_stop - clean_start);

          _ktoken0_aus = ktokens.at(0); // SC20240117
          if (_ktoken0_aus == "YEAR" || _ktoken0_aus == "year") {
            ptmp.year = aurostd::string2utype<uint>(ktokens.at(1)); // check year // SC20240117
          }
          if (_ktoken0_aus == "TYPE" || _ktoken0_aus == "type") { // check type // SC20240117
            if (object == "TALK" || object == "PRESENTATION_TALK") {
              ptmp.type = "PRESENTATION_TALK";
              ptmp._isinvited = true;
            }
            if (object == "SEMINAR" || object == "PRESENTATION_SEMINAR") {
              ptmp.type = "PRESENTATION_SEMINAR";
              ptmp._isinvited = true;
            }
            if (object == "COLLOQUIUM" || object == "PRESENTATION_COLLOQUIUM") {
              ptmp.type = "PRESENTATION_COLLOQUIUM";
              ptmp._isinvited = true;
            }
            if (object == "KEYNOTE" || object == "PRESENTATION_KEYNOTE") {
              ptmp.type = "PRESENTATION_KEYNOTE";
              ptmp._isinvited = true;
            }
            if (object == "PLENARY" || object == "PRESENTATION_PLENARY") {
              ptmp.type = "PRESENTATION_PLENARY";
              ptmp._isinvited = true;
            }
            if (object == "TUTORIAL" || object == "PRESENTATION_TUTORIAL") {
              ptmp.type = "PRESENTATION_TUTORIAL";
              ptmp._isinvited = true;
            }
            if (object == "PANELIST" || object == "PRESENTATION_PANELIST") {
              ptmp.type = "PRESENTATION_PANELIST";
              ptmp._isinvited = true;
            }
            if (object == "CONTRIBUTED" || object == "PRESENTATION_CONTRIBUTED") {
              ptmp.type = "PRESENTATION_CONTRIBUTED";
            }
            if (object == "POSTER" || object == "PRESENTATION_POSTER") {
              ptmp.type = "PRESENTATION_POSTER";
            }
            if (object == "ARTICLE" || object == "article") {
              ptmp.type = "ARTICLE";
            }
            if (object == "ALABEL" || object == "alabel") {
              ptmp.type = "ALABEL";
            }
            if (object == "JLABEL" || object == "jlabel") {
              ptmp.type = "JLABEL";
            }
            if (object == "PUBLICATION" || object == "publication") {
              ptmp.type = "ARTICLE";
            }
            if (object == "EDUCATION" || object == "education") {
              ptmp.type = "EDUCATION";
            }
            if (object == "RESEARCH" || object == "research") {
              ptmp.type = "RESEARCH";
            }
            if (object == "ACADEMIC" || object == "academic") {
              ptmp.type = "ACADEMIC";
            }
            if (object == "SERVICEOUTSIDE" || object == "serviceoutside") {
              ptmp.type = "SERVICEOUTSIDE";
            }
            if (object == "SERVICEINSIDE" || object == "serviceinside") {
              ptmp.type = "SERVICEINSIDE";
            }
            if (object == "TEACHING" || object == "teaching") {
              ptmp.type = "TEACHING";
            }
            if (object == "ADVISING" || object == "advising") {
              ptmp.type = "ADVISING";
            }
            if (object == "PATENTS" || object == "patents") {
              ptmp.type = "PATENTS";
            }
            if (object == "PRESS" || object == "press") {
              ptmp.type = "PRESS";
            }
            if (object == "AWARDS" || object == "awards") {
              ptmp.type = "AWARDS";
            }
          }
          if (_ktoken0_aus == "ONLINE" || _ktoken0_aus == "online") { // check ONLINESS  // SC20240117
            if (object == "YES" || object == "yes" || object == "true" || object == "true" || object == "1") {
              ptmp._isonline = true;
            }
            if (object == "NO" || object == "no" || object == "false" || object == "false" || object == "0") {
              ptmp._isonline = false;
            }
          }

          if (_ktoken0_aus == "TITLE" || _ktoken0_aus == "title") {
            ptmp.title = object; // check title // SC20240117
          }
          if (_ktoken0_aus == "JOURNAL" || _ktoken0_aus == "journal") {
            ptmp.journal = object; // check journal // SC20240117
          }
          if (_ktoken0_aus == "LINK" || _ktoken0_aus == "link") {
            ptmp.link = object; // check link // SC20240117
          }
          if (_ktoken0_aus == "ARXIV" || _ktoken0_aus == "arxiv") {
            ptmp.arxiv = object; // check arxiv // SC20240117
          }
          if (_ktoken0_aus == "SUPPLEMENTARY" || _ktoken0_aus == "supplementary") {
            ptmp.supplementary = object; // check supplementary // SC20240117
          }
          if (_ktoken0_aus == "SUPPLEMENTARY_URL" || _ktoken0_aus == "supplementary_url") {
            ptmp.supplementary_url = object; // check supplementary_url // SC20240117
          }
          if (_ktoken0_aus == "BIBTEX" || _ktoken0_aus == "bibtex") {
            ptmp.bibtex = object; // check bibtex // SC20240117
          }
          if (_ktoken0_aus == "BIBTEX_JOURNAL" || _ktoken0_aus == "bibtex_journal") {
            ptmp.bibtex_journal = object; // check bibtex_journal // SC20240117
          }
          if (_ktoken0_aus == "BIBTEX_VOLUME" || _ktoken0_aus == "bibtex_volume") {
            ptmp.bibtex_volume = object; // check bibtex_volume // SC20240117
          }
          if (_ktoken0_aus == "BIBTEX_ISSUE" || _ktoken0_aus == "bibtex_issue") {
            ptmp.bibtex_issue = object; // check bibtex_issue // SC20240117
          }
          if (_ktoken0_aus == "BIBTEX_PAGES" || _ktoken0_aus == "bibtex_pages") {
            ptmp.bibtex_pages = object; // check bibtex_pages // SC20240117
          }
          if (_ktoken0_aus == "BIBTEX_YEAR" || _ktoken0_aus == "bibtex_year") {
            ptmp.bibtex_year = object; // check bibtex_year // SC20240117
          }
          if (_ktoken0_aus == "PLACE" || _ktoken0_aus == "place") {
            ptmp.place = object; // check place // SC20240117
          }
          if (_ktoken0_aus == "DATE" || _ktoken0_aus == "date") {
            ptmp.date = object; // check date // SC20240117
          }
          if (_ktoken0_aus == "HOST" || _ktoken0_aus == "host") {
            ptmp.host = object; // check host // SC20240117
          }
          if (_ktoken0_aus == "ABSTRACT" || _ktoken0_aus == "abstract") {
            ptmp.abstract = object; // check abstract // SC20240117
          }
          if (_ktoken0_aus == "PDF" || _ktoken0_aus == "pdf") {
            ptmp.pdf = object; // check pdf // SC20240117
          }
          if (_ktoken0_aus == "DOI" || _ktoken0_aus == "doi") {
            ptmp.doi = object; // check doi // SC20240117
          }
          if (_ktoken0_aus == "NEW" || _ktoken0_aus == "new") {
            ptmp.newflag = true; // check doi // SC20240117
          }
          if (_ktoken0_aus == "NEWFLAG" || _ktoken0_aus == "newflag") {
            ptmp.newflag = true; // check doi // SC20240117
          }
          if (_ktoken0_aus == "WNUMBER" || _ktoken0_aus == "wnumber") {
            ptmp.wnumber = aurostd::string2utype<uint>(object); // check wnumber   // SC20240117
          }
          if (_ktoken0_aus == "AUTHOR" || _ktoken0_aus == "author") {
            aurostd::string2tokensAdd(object, ptmp.vauthor, ","); // SC20240117
          }
          if (_ktoken0_aus == "EXTRA_HTML" || _ktoken0_aus == "extra_html") {
            ptmp.vextra_html.push_back(object); // SC20240117
          }
          if (_ktoken0_aus == "EXTRA_LATEX" || _ktoken0_aus == "extra_latex") {
            ptmp.vextra_latex.push_back(object); // SC20240117
          }
          if (_ktoken0_aus == "KEYWORD" || _ktoken0_aus == "keyword") {
            aurostd::string2tokensAdd(object, ptmp.vkeyword, ","); // SC20240117
          }
          if (_ktoken0_aus == "SPONSOR" || _ktoken0_aus == "sponsor") {
            aurostd::string2tokensAdd(object, ptmp.vsponsor, ","); // SC20240117
          }
          if (_ktoken0_aus == "ALLOY" || _ktoken0_aus == "alloy") {
            aurostd::string2tokensAdd(object, ptmp.valloy, ","); // SC20240117
          }
          if (_ktoken0_aus == "ALABEL" || _ktoken0_aus == "alabel") {
            aurostd::string2tokensAdd(object, valabel, ","); // SC20240117
          }
          if (_ktoken0_aus == "JLABEL" || _ktoken0_aus == "jlabel") {
            aurostd::string2tokensAdd(object, vjlabel, ","); // SC20240117
          }
        }
      }

      if (ptmp.bibtex.empty() && !ptmp.bibtex_volume.empty()) {
        stringstream oss;
        oss << "curtarolo:art" << ptmp.wnumber;
        ptmp.bibtex = oss.str();
      }

      if (ptmp.type.empty()) {
        for (k = i; k <= j; k++) {
          if (!XHOST.QUIET_CERR) {
            cerr << "entry=[" << vpres.at(k) << "]" << endl;
          }
        }
      } else {
        if ((aurostd::substring2bool(ptmp.type, "ARTICLE") || aurostd::substring2bool(ptmp.type, "PUBLICATION")) &&
            (aurostd::substring2bool(what2print, "ARTICLE") || aurostd::substring2bool(what2print, "PUBLICATION"))) {
          voutreach.push_back(ptmp);
        }
        if (aurostd::substring2bool(ptmp.type, "PRESENTATION") && aurostd::substring2bool(what2print, "PRESENTATION")) {
          if (ptmp.vauthor.empty()) {
            ptmp.vauthor.emplace_back("S. Curtarolo"); // for SC CV
          }
          voutreach.push_back(ptmp);
        }
        if (aurostd::substring2bool(ptmp.type, "EDUCATION") && aurostd::substring2bool(what2print, "EDUCATION")) {
          voutreach.push_back(ptmp);
        }
        if (aurostd::substring2bool(ptmp.type, "RESEARCH") && aurostd::substring2bool(what2print, "RESEARCH")) {
          voutreach.push_back(ptmp);
        }
        if (aurostd::substring2bool(ptmp.type, "ACADEMIC") && aurostd::substring2bool(what2print, "ACADEMIC")) {
          voutreach.push_back(ptmp);
        }
        if (aurostd::substring2bool(ptmp.type, "SERVICEOUTSIDE") && aurostd::substring2bool(what2print, "SERVICEOUTSIDE")) {
          voutreach.push_back(ptmp);
        }
        if (aurostd::substring2bool(ptmp.type, "SERVICEINSIDE") && aurostd::substring2bool(what2print, "SERVICEINSIDE")) {
          voutreach.push_back(ptmp);
        }
        if (aurostd::substring2bool(ptmp.type, "TEACHING") && aurostd::substring2bool(what2print, "TEACHING")) {
          voutreach.push_back(ptmp);
        }
        if (aurostd::substring2bool(ptmp.type, "ADVISING") && aurostd::substring2bool(what2print, "ADVISING")) {
          voutreach.push_back(ptmp);
        }
        if (aurostd::substring2bool(ptmp.type, "PATENTS") && aurostd::substring2bool(what2print, "PATENTS")) {
          voutreach.push_back(ptmp);
        }
        if (aurostd::substring2bool(ptmp.type, "PRESS") && aurostd::substring2bool(what2print, "PRESS")) {
          voutreach.push_back(ptmp);
        }
        if (aurostd::substring2bool(ptmp.type, "AWARDS") && aurostd::substring2bool(what2print, "AWARDS")) {
          voutreach.push_back(ptmp);
        }
      }
    }
  }

  if (LDEBUG) {
    cerr << "voutreach_load [3]" << endl;
  }

  for (_outreach& entry : voutreach) {
    fixlabel(valabel, entry.vauthor);
    //   cout << voutreach[i].vauthor.size() << endl;
    for (string& author : entry.vauthor) { // fix corresponding  authors
      aurostd::StringSubstInPlace(author, "*", ""); // fix corresponding authors
    }
  }

  if (LDEBUG) {
    cerr << "voutreach_load [4]" << endl;
  }

  for (size_t i = 0; i < voutreach.size(); i++) {
    fixlabel(valabel, voutreach[i].vextra_html);
  }
  for (size_t i = 0; i < voutreach.size(); i++) {
    fixlabel(valabel, voutreach[i].vextra_latex);
  }
  for (size_t i = 0; i < voutreach.size(); i++) {
    fixlabel(vjlabel, voutreach[i].journal);
  }
  for (size_t i = 0; i < voutreach.size(); i++) {
    fixlabel(vjlabel, voutreach[i].bibtex_journal);
  }
  for (size_t i = 0; i < voutreach.size(); i++) {
    fixlabel(vjlabel, voutreach[i].vextra_html);
  }
  for (size_t i = 0; i < voutreach.size(); i++) {
    fixlabel(vjlabel, voutreach[i].vextra_latex);
  }
  if (aurostd::substring2bool(ptmp.type, "ARTICLE") || aurostd::substring2bool(ptmp.type, "PUBLICATION")) {
    voutreach_remove_duplicate(voutreach); // talks can be duplicated
  }
  // SAVE the global one
  if (voutreach_global_list.empty()) {
    voutreach_global_list = voutreach;
  }
  // CHECK MAX YEAR
  for (size_t i = 0; i < voutreach.size(); i++) {
    if (voutreach_global_max_year < voutreach[i].year) {
      voutreach_global_max_year = voutreach[i].year;
    }
  }
  // CHECK MIN YEAR
  for (size_t i = 0; i < voutreach.size(); i++) {
    if (voutreach_global_min_year > voutreach[i].year) {
      voutreach_global_min_year = voutreach[i].year;
    }
  }
  // DONE

  if (LDEBUG) {
    cerr << "voutreach_load [5]" << endl;
  }

  // save cache
  if (aurostd::substring2bool(ptmp.type, "ARTICLE") || aurostd::substring2bool(ptmp.type, "PUBLICATION")) {
    voutreach_publications.clear();
    for (const _outreach& entry : voutreach) {
      voutreach_publications.push_back(entry);
      // if(!XHOST.QUIET_CERR) cerr << voutreach[i] << endl;
    }
  }
  if (LDEBUG) {
    cerr << "voutreach_load [6]" << endl;
  }

  if (aurostd::substring2bool(ptmp.type, "PRESENTATION")) {
    voutreach_presentations.clear();
    for (const _outreach& entry : voutreach) {
      voutreach_presentations.push_back(entry);
    }
  }

  if (LDEBUG) {
    cerr << "voutreach_load [9]" << endl;
  }

  return voutreach.size();
}

// ******************************************************************************************************************************************************
void HT_CHECK_GRANTS(ostream& oss) { //,const vector<string>& vitems,string msg1,string msg2,string sectionlabel)
  const aurostd::xoption vflag = XHOST.vflag_outreach;
  vector<_outreach> voutreach;
  voutreach_load(voutreach, "PUBLICATIONS");
  if (!vflag.flag("GRANTS")) {
    const string message = "No grants.";
    throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message);
  }
  if (!XHOST.QUIET_CERR) {
    cerr << "LOADED " << voutreach.size() << " " << endl;
  }
  const string grant = vflag.getattachedscheme("GRANTS");
  //  if(mode==HTRESOURCE_MODE_NONE) {mode=HTRESOURCE_MODE_PHP_AUTHOR;} // by default

  if (XHOST.vflag_control.flag("PRINT_MODE::LATEX")) {
    for (size_t i = 0; i < voutreach.size(); i++) {
      bool found = false;
      for (size_t j = 0; j < voutreach[i].vsponsor.size() && found == false; j++) {
        if (aurostd::substring2bool(voutreach[i].vsponsor[j], grant)) {
          found = true;
        }
      }
      if (found) {
        XHOST.vflag_control.flag("PRINT_MODE::TXT", true);
        oss << voutreach[i] << endl; // << endl;
      }
    }
  }
}

// ******************************************************************************************************************************************************
bool ProcessPhpLatexCv() {
  const aurostd::xoption vflag = XHOST.vflag_outreach;

  bool Arun = false;
  vector<_outreach> voutreach;
  if (!XHOST.vflag_control.flag("PRINT_MODE::JSON") && !XHOST.vflag_control.flag("PRINT_MODE::HTML") && !XHOST.vflag_control.flag("PRINT_MODE::LATEX") && !XHOST.vflag_control.flag("PRINT_MODE::TXT")) {
    XHOST.vflag_control.flag("PRINT_MODE::HTML", true);
  }
  if (XHOST.vflag_control.flag("PRINT_MODE::JSON")) {
    XHOST.vflag_control.flag("PRINT_MODE::LATEX", false);
    XHOST.vflag_control.flag("PRINT_MODE::TXT", false);
    XHOST.vflag_control.flag("PRINT_MODE::HTML", false);
  }
  if (XHOST.vflag_control.flag("PRINT_MODE::HTML")) {
    XHOST.vflag_control.flag("PRINT_MODE::LATEX", false);
    XHOST.vflag_control.flag("PRINT_MODE::TXT", false);
    XHOST.vflag_control.flag("PRINT_MODE::JSON", false);
  }
  if (XHOST.vflag_control.flag("PRINT_MODE::LATEX")) {
    XHOST.vflag_control.flag("PRINT_MODE::TXT", false);
    XHOST.vflag_control.flag("PRINT_MODE::HTML", false);
    XHOST.vflag_control.flag("PRINT_MODE::JSON", false);
  }
  if (XHOST.vflag_control.flag("PRINT_MODE::TXT")) {
    XHOST.vflag_control.flag("PRINT_MODE::HTML", false);
    XHOST.vflag_control.flag("PRINT_MODE::LATEX", false);
    XHOST.vflag_control.flag("PRINT_MODE::JSON", false);
  }

  if (!Arun && XHOST.vflag_control.flag("PHP::CENTER_MISSION")) {
    center_print(HTRESOURCE_MODE_PHP_AUTHOR, cout);
    Arun = true;
  }
  if (!Arun && XHOST.vflag_control.flag("PHP::PUBS_ALLOY")) {
    voutreach_print(HTRESOURCE_MODE_PHP_ALLOY, cout, "ARTICLES");
    Arun = true;
  }
  if (!Arun && XHOST.vflag_control.flag("PHP::PUBS_KEYWORD")) {
    voutreach_print(HTRESOURCE_MODE_PHP_THRUST, cout, "ARTICLES");
    Arun = true;
  }
  if (!Arun && XHOST.vflag_control.flag("PHP::ITALKS")) {
    voutreach_print(HTRESOURCE_MODE_PHP_AUTHOR, cout, "PRESENTATIONS");
    Arun = true;
  }
  if (!Arun && XHOST.vflag_control.flag("CV::PUBS")) {
    voutreach_print(HTRESOURCE_MODE_PHP_AUTHOR, cout, "ARTICLES");
    Arun = true;
  }
  if (!Arun && vflag.flag("GRANTS")) {
    HT_CHECK_GRANTS(cout);
    Arun = true;
  }
  if (!Arun && XHOST.vflag_control.flag("CV::ITALKS")) {
    voutreach_print(HTRESOURCE_MODE_PHP_AUTHOR, cout, "PRESENTATIONS");
    Arun = true;
  }
  if (!Arun && XHOST.vflag_control.flag("CV::ACADEMIC")) {
    voutreach_load(voutreach, "ACADEMIC");
    voutreach_print_everything(cout, voutreach.at(0).vextra_latex, "Academic", "ACADEMIC", "\\section{Academic Positions}\\label{academic_positions}");
    Arun = true;
  }
  if (!Arun && XHOST.vflag_control.flag("CV::RESEARCH")) {
    voutreach_load(voutreach, "RESEARCH");
    voutreach_print_everything(cout, voutreach.at(0).vextra_latex, "Research", "RESEARCH", "\\section{Research Experience}\\label{research_experience}");
    Arun = true;
  }
  if (!Arun && XHOST.vflag_control.flag("CV::EDUCATION")) {
    voutreach_load(voutreach, "EDUCATION");
    voutreach_print_everything(cout, voutreach.at(0).vextra_latex, "Education", "EDUCATION", "\\section{Education}\\label{education}");
    Arun = true;
  }
  if (!Arun && XHOST.vflag_control.flag("CV::TEACHING")) {
    voutreach_load(voutreach, "TEACHING");
    voutreach_print_everything(cout, voutreach.at(0).vextra_latex, "Teaching", "TEACHING", "\\section{Teaching Experience}\\label{teaching_experience}");
    Arun = true;
  }
  if (!Arun && XHOST.vflag_control.flag("CV::ADVISING")) {
    voutreach_load(voutreach, "ADVISING");
    voutreach_print_everything(cout, voutreach.at(0).vextra_latex, "Advising", "ADVISING", "\\section{Advising Experience}\\label{advising_experience}");
    Arun = true;
  }
  if (!Arun && XHOST.vflag_control.flag("CV::AWARDS")) {
    voutreach_load(voutreach, "AWARDS");
    voutreach_print_everything(cout, voutreach.at(0).vextra_latex, "Awards", "AWARDS", "\\section{Awards and Honors}\\label{awards}");
    Arun = true;
  }
  if (!Arun && XHOST.vflag_control.flag("CV::PRESS")) {
    voutreach_load(voutreach, "PRESS");
    voutreach_print_everything(cout, voutreach.at(0).vextra_latex, "Press", "PRESS", "\\section{Press and news releases}\\label{pressreleases}");
    Arun = true;
  }
  if (!Arun && XHOST.vflag_control.flag("CV::PATENTS")) {
    voutreach_load(voutreach, "PATENTS");
    voutreach_print_everything(cout, voutreach.at(0).vextra_latex, "Patents", "PATENTS", "\\section{Patents}\\label{patents}");
    Arun = true;
  }
  if (!Arun && XHOST.vflag_control.flag("CV::SERVICE_OUTSIDE")) {
    voutreach_load(voutreach, "SERVICEOUTSIDE");
    voutreach_print_everything(cout, voutreach.at(0).vextra_latex, "ServiceOutside", "SERVICE OUTSIDE", "\\section{Outreach and Professional Activities}\\label{outreach}");
    Arun = true;
  }
  if (!Arun && XHOST.vflag_control.flag("CV::SERVICE_INSIDE")) {
    voutreach_load(voutreach, "SERVICEINSIDE");
    voutreach_print_everything(cout, voutreach.at(0).vextra_latex, "ServiceInside", "SERVICE INSIDE", "\\section{Duke University - Academic Service Activities}\\label{duke_service}");
    Arun = true;
  }
  if (!Arun && XHOST.vflag_control.flag("CV::AUTHOR")) {
    // something need to be specified
    voutreach_print(HTRESOURCE_MODE_PHP_AUTHOR, cout, "ARTICLES");
    Arun = true;
  }

  return Arun;
}

// ******************************************************************************************************************************************************

uint voutreach_sort_year(vector<_outreach>& voutreach) {
  if (!voutreach.empty()) {
    if (voutreach.at(0).year != 0) {
      sort(voutreach.begin(), voutreach.end(), _sort_outreach_outreach_year_());
    }
  }
  return voutreach.size();
}
uint voutreach_rsort_year(vector<_outreach>& voutreach) {
  if (!voutreach.empty()) {
    if (voutreach.at(0).year != 0) {
      sort(voutreach.begin(), voutreach.end(), _rsort_outreach_outreach_year_());
    }
  }
  return voutreach.size();
}
uint voutreach_sort_wnumber(vector<_outreach>& voutreach) {
  if (!voutreach.empty()) {
    if (voutreach.at(0).wnumber != 0) {
      sort(voutreach.begin(), voutreach.end(), _sort_outreach_outreach_wnumber_());
    }
  }
  return voutreach.size();
}
uint voutreach_rsort_wnumber(vector<_outreach>& voutreach) {
  if (!voutreach.empty()) {
    if (voutreach.at(0).wnumber != 0) {
      sort(voutreach.begin(), voutreach.end(), _rsort_outreach_outreach_wnumber_());
    }
  }
  return voutreach.size();
}

// ******************************************************************************************************************************************************
void center_print(uint mode, ostream& oss) {
  const aurostd::xoption vflag = XHOST.vflag_outreach;

  if (mode) {
    ;
  } // dummy load

  if (XHOST.vflag_control.flag("PHP::CENTER_MISSION")) {
    //  if(mode=print_mode && mode=HTRESOURCE_MODE_LATEX)
    oss << "  <br>" << endl;
  }
}

// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2024           *
// *                                                                         *
// ***************************************************************************
