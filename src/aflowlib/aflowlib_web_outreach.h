
#ifndef AFLOWLIB_WEB_OUTREACH_H
#define AFLOWLIB_WEB_OUTREACH_H

#include <ostream>
#include <string>
#include <vector>

#include <sys/types.h>

// ***************************************************************************
// _OUTREACH CLASS
class _outreach {
public:
    // constructors/destructors
  _outreach();    // do nothing
  _outreach(const _outreach& b);    // do nothing
  ~_outreach();        // do nothing
    // OPERATORS                                              // --------------------------------------
  const _outreach& operator=(const _outreach& b);             // some operators
  void clear();                                         // clear
    // CONTENT
  uint wnumber;
  bool newflag;
  uint year;
  std::vector<std::string> vauthor;
  std::vector<std::string> vcorrespondingauthor;
  std::string title;
  std::string journal, link, arxiv, supplementary, supplementary_url;
  std::string place, date;
  std::string type;   // ARTICLE PRESENTATION_TALK PRESENTATION_SEMINAR PRESENTATION_COLLOQUIUM PRESENTATION_KEYNOTE PRESENTATION_PLENARY PRESENTATION_TUTORIAL PRESENTATION_CONTRIBUTED PRESENTATION_POSTER
  bool _isinvited;       // YES
  bool _isonline;       // YES
  std::string host;           // who invited
  std::string abstract;       // if available and in LaTeX
  std::string pdf;
  std::string doi;
  std::string bibtex, bibtex_journal, bibtex_volume, bibtex_issue, bibtex_pages, bibtex_year;  // SC20201228
  std::vector<std::string> vextra_html, vextra_latex;
  std::vector<std::string> vkeyword, vsponsor, valloy;
    // operators/functions                                    // operator/functions
  friend std::ostream& operator<<(std::ostream&, const _outreach&);   // print
  [[nodiscard]] std::string print_string() const;                                // print string from <<
private:                                                    // ---------------------------------------
  void free();                                              // to free everything
  void copy(const _outreach& b);                            // the flag is necessary because sometimes you need to allocate the space.
};

uint voutreach_load(std::vector<_outreach>& voutreach, const std::string& what2print);
void voutreach_print(uint _mode, std::ostream& oss, std::string what2print);
void voutreach_print_everything(std::ostream& oss, const std::vector<std::string>& vitems, std::string msg1, std::string msg2, std::string sectionlabel);

// sort

class _sort_outreach_outreach_year_ {
public:
  bool operator()(const _outreach& a1, const _outreach& a2) const { return (bool) (a1.year < a2.year); }
}; // sorting through reference
class _rsort_outreach_outreach_year_ {
public:
  bool operator()(const _outreach& a1, const _outreach& a2) const { return (bool) (a1.year > a2.year); }
}; // sorting through reference
class _sort_outreach_outreach_wnumber_ {
public:
  bool operator()(const _outreach& a1, const _outreach& a2) const { return (bool) (a1.wnumber < a2.wnumber); }
}; // sorting through reference
class _rsort_outreach_outreach_wnumber_ {
public:
  bool operator()(const _outreach& a1, const _outreach& a2) const { return (bool) (a1.wnumber > a2.wnumber); }
}; // sorting through reference

uint voutreach_sort_year(std::vector<_outreach>& voutreach);
uint voutreach_rsort_year(std::vector<_outreach>& voutreach);
uint voutreach_sort_wnumber(std::vector<_outreach>& voutreach);
uint voutreach_rsort_wnumber(std::vector<_outreach>& voutreach);

uint voutreach_remove_duplicate(std::vector<_outreach>& voutreach);

// automatic load up
bool ProcessPhpLatexCv();
// ***************************************************************************
void center_print(uint mode, std::ostream& oss);

// ***************************************************************************
// references search
void SystemReferences(const std::string& system_in, std::vector<std::string>& vrefs, bool = false);  // if true then put et_al
void SystemReferences(const std::string& system_in, std::vector<uint>& vwnumber);
bool SystemInSystems(const std::string& system, const std::string& systems);
bool SystemInSystems(const std::string& system, const std::vector<std::string>& vsystems);
bool AlloyAlphabeticLIBRARY(const std::string& system);

#endif // AFLOWLIB_WEB_OUTREACH_H
