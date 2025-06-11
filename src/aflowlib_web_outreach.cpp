// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *                                                                         *
// ***************************************************************************
// Stefano Curtarolo - Duke

#ifndef _AFLOWLIB_WEB_OUTREACH_CPP_
#define _AFLOWLIB_WEB_OUTREACH_CPP_
#include "aflow.h"

#define TOC 0

#define WEB_PDF  string("http://"+XHOST.AFLOW_MATERIALS_SERVER+"/auro/AUROARTICULA/")
#define WEB_BIBTEX  string("http://"+XHOST.AFLOW_MATERIALS_SERVER+"/auro/AUROARTICULA/BIBITEMS/")
#define WEB_BIBTEX_FILE  string("/www/auro/AUROARTICULA/BIBITEMS/")
#define WEB_DOI  string("http://dx.doi.org/")

#define THRUST_RECENT_ARTICLES  20
#define THRUST_RECENT_YEARS     5
#define AUTHOR_RECENT_ARTICLES  300
#define AUTHOR_RECENT_YEARS     30

#define MAX_YEAR_PRESENTATIONS  2023

// ******************************************************************************************************************************************************
// _OUTREACH CLASS
// namespace web

_outreach::_outreach() {
  newflag=FALSE;
  wnumber=0;
  year=0;
  vauthor.clear();
  vcorrespondingauthor.clear();
  title="";
  journal="";
  arxiv="";
  supplementary="";
  supplementary_url="";
  place="";
  date="";
  link="";
  type="";
  _isinvited=FALSE;
  _isonline=FALSE;
  host="";
  abstract="";
  pdf="";
  doi="";
  bibtex="";
  bibtex_journal="";
  bibtex_volume="";
  bibtex_issue="";
  bibtex_pages="";
  bibtex_year="";
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

void _outreach::free() {
}

void _outreach::copy(const _outreach& b) {
  // const _outreach& _outreach::operator=(const _outreach& b)       // operator=
  free();
  newflag=b.newflag;
  wnumber=b.wnumber;
  year=b.year;
  vauthor.clear();for(uint i=0;i<b.vauthor.size();i++) vauthor.push_back(b.vauthor.at(i));
  vcorrespondingauthor.clear();for(uint i=0;i<b.vcorrespondingauthor.size();i++) vcorrespondingauthor.push_back(b.vcorrespondingauthor.at(i));
  title=b.title;
  journal=b.journal;
  arxiv=b.arxiv;
  supplementary=b.supplementary;
  supplementary_url=b.supplementary_url;
  place=b.place;
  date=b.date;
  link=b.link;
  type=b.type;
  _isinvited=b._isinvited;
  _isonline=b._isonline;
  host=b.host;
  abstract=b.abstract;
  pdf=b.pdf;
  doi=b.doi;
  bibtex=b.bibtex;
  bibtex_journal=b.bibtex_journal;
  bibtex_volume=b.bibtex_volume;
  bibtex_issue=b.bibtex_issue;
  bibtex_pages=b.bibtex_pages;
  bibtex_year=b.bibtex_year;
  vextra_html.clear();for(uint i=0;i<b.vextra_html.size();i++) vextra_html.push_back(b.vextra_html.at(i));
  vextra_latex.clear();for(uint i=0;i<b.vextra_latex.size();i++) vextra_latex.push_back(b.vextra_latex.at(i));
  vkeyword.clear();for(uint i=0;i<b.vkeyword.size();i++) vkeyword.push_back(b.vkeyword.at(i));
  vsponsor.clear();for(uint i=0;i<b.vsponsor.size();i++) vsponsor.push_back(b.vsponsor.at(i));
  valloy.clear();for(uint i=0;i<b.valloy.size();i++) valloy.push_back(b.valloy.at(i));
}

// copy
_outreach::_outreach(const _outreach& b) {
  // [OBSOLETE]   free();
  // [OBSOLETE]  *this=b;
  copy(b);
}

// copies xtructures: b=a
const _outreach& _outreach::operator=(const _outreach& b) {  // operator=
  if(this!=&b) {
    free();
    copy(b);
  }
  return *this;
}

void _outreach::clear() {
  _outreach outreach_temp;
  copy(outreach_temp);
}

string _outreach::print_string() {
  stringstream z; z << *this;
  return z.str();
}

// ******************************************************************************************************************************************************
// ******************************************************************************************************************************************************
// ARTICLES/PRESENTATIONS
// ******************************************************************************************************************************************************
// ******************************************************************************************************************************************************

uint fixlabel(const vector<string>& vlabel,string& labelfixed) {
  for(uint j=0;j<vlabel.size();j+=2) 
    //    if(labelfixed==vlabel.at(j))
    aurostd::StringSubst(labelfixed,vlabel.at(j),vlabel.at(j+1));
  aurostd::StringSubst(labelfixed,"_WEB_PDF_",WEB_PDF);
  aurostd::StringSubst(labelfixed,"_WEB_DOI_",WEB_DOI);
  return TRUE;
}

uint fixlabel(const vector<string>& vlabel,vector<string>& vlabelfixed) {
  for(uint i=0;i<vlabelfixed.size();i++) {
    for(uint j=0;j<vlabel.size();j+=2) 
      //      if(vlabelfixed.at(i)==vlabel.at(j))
      aurostd::StringSubst(vlabelfixed.at(i),vlabel.at(j),vlabel.at(j+1));
    aurostd::StringSubst(vlabelfixed.at(i),"_WEB_PDF_",WEB_PDF);
    aurostd::StringSubst(vlabelfixed.at(i),"_WEB_DOI_",WEB_DOI);
  }
  return vlabelfixed.size();
}


uint bibtex2file(string bibtex,string _authors,string _title,string journal,string volume,string issue,string pages,string year,string abstract,string doi,string bibfile) { //SC20201229
  stringstream bibcontent;
  // FIX THE AUTHORS
  string authors=_authors;
  vector<string> vfix;
  aurostd::string2tokens("Buongiorno Nardelli,van Roekeghem,Aspuru-Guzik,Hattrick-Simpers,DeCost,de Coss,De Santo,De Gennaro,Al Rahal Al Orabi,de Jong,D'Amico,van der Zwaag,van de Walle,Di Stefano,Ojeda Mota,Simmons Jr.,Mattos Jr.",vfix,",");
  for(size_t i=0;i<vfix.size();i++) aurostd::StringSubst(authors,vfix[i],string("{"+vfix[i]+"}")); // FIX
  aurostd::StringSubst(authors,".",".~"); 
  aurostd::StringSubst(authors,"~ "," "); 
  aurostd::StringSubst(authors,"~ "," "); 
  aurostd::StringSubst(authors,"~-","-"); 
  aurostd::StringSubst(authors,","," and");
  aurostd::StringSubst(authors,"and and","and");
  authors=aurostd::html2latex(authors);
  // FIX THE TITLES
  string title="{"+_title+"}";
  //  aurostd::string2tokens("SnSe,SnTe,GeTe,CoSb,FePt,AgPt,LSO:Ce,Eu:SrI,Mo-Ru-Ta-W,BaSnO,Heusler,CsI(Tl),NaI(Tl),LaBr,L1,IrV,RhV,Pt$_{8}$Ti,Mg-B,Fe:Mo:C,Fe-C,LiB,MgB,AFLUX,LUX,ACBN0,PAOFLOW,Ag,Au,Cu,Mg,Bi,In,Sb,Fe,Mo,Kr,Ar,Ne,Cs,Xe",vfix,",");
  for(uint i=0;i<vfix.size();i++) aurostd::StringSubst(title,vfix.at(i),string("{"+vfix.at(i)+"}")); // FIX 
  // NOW build
  bibcontent << "@article{" << bibtex << "," << endl;
  if(authors.size()) bibcontent << " author={" << authors << "}";  // AUTHORS
  if(title.size()) bibcontent << "," << endl << " title={" << aurostd::html2latex(title) << "}"; // TITLE
  if(journal.size()) bibcontent << "," << endl << " journal={" << journal << "}";   // JOURNAL
  if(volume.size()) bibcontent << "," << "volume={" << volume << "}";   // VOLUME
  if(issue.size()) bibcontent << "," << "issue={" << issue << "}";   // ISSUE
  if(pages.size()) bibcontent << "," << "pages={" << pages << "}";   // PAGES
  if(year.size()) bibcontent << "," << "year={" << year << "}";   // YEAR
  if(abstract.size()) bibcontent << "," << endl << " abstract={" << abstract << "}";   // ABSTRACT
  if(doi.size()) bibcontent << "," << "doi={" << doi << "}";   // DOI
  bibcontent << endl << "}" << endl;
  bibcontent << "% Automatically generated - AFLOW " << AFLOW_VERSION << endl;

  // save
  aurostd::stringstream2file(bibcontent,bibfile);
  return bibcontent.str().size();
}


ostream& operator<<(ostream& oss,const _outreach& outreach) {
  // ******************************************************************************************************************************************************
  bool compact=TRUE;
  stringstream newline;newline << endl;

  // ***************************************************************************
  // ARTICLE
  if(aurostd::substring2bool(outreach.type,"ARTICLE")) {

    // generate authors_json
    string authors_json;
    authors_json="[";
    for(size_t iauth=0;iauth<outreach.vauthor.size();iauth++) {
      authors_json+="\""+outreach.vauthor[iauth]+"\"";
      if(iauth!=outreach.vauthor.size()-1) authors_json+=",";
    }
    authors_json+="]";

    // generate authors_txt
    string authors_txt;
    authors_txt="";
    for(size_t iauth=0;iauth<outreach.vauthor.size();iauth++) {
        authors_txt+=outreach.vauthor[iauth];
        if(outreach.vauthor.size()==2 && iauth==outreach.vauthor.size()-2) authors_txt+=" and ";
        if(outreach.vauthor.size()!=2 && iauth==outreach.vauthor.size()-2) authors_txt+=", and ";
        if(iauth!=outreach.vauthor.size()-2 && iauth!=outreach.vauthor.size()-1) authors_txt+=", ";
      }

    string wnumber=aurostd::utype2string(outreach.wnumber);
    if(wnumber.length()<2) wnumber="0"+wnumber;

    //  generate journal HTML style
    string journal=""; 
    if(!outreach.bibtex_journal.size()) { // with journal
      journal=outreach.journal;           // with journal
    } else {
      journal=outreach.bibtex_journal;
      if(outreach.bibtex_volume.size()) journal+=" <b>"+outreach.bibtex_volume+"</b>";
      if(outreach.bibtex_issue.size()) journal+="("+outreach.bibtex_issue+")";
      if(outreach.bibtex_volume.size()) journal+=",";
      if(outreach.bibtex_pages.size()) journal+=" "+outreach.bibtex_pages;
      if(outreach.bibtex_year.size()) journal+=" ("+outreach.bibtex_year+")";
    }

    // ***************************************************************************
    // TXT
    if(XHOST.vflag_control.flag("PRINT_MODE::TXT")) {
      // 3rd line AUTHORS
      aurostd::StringSubst(authors_txt,"Csányi","Csanyi");
      oss << aurostd::html2txt(authors_txt) << ", ";
      // 4th line TITLE
      oss << "\"" << aurostd::html2txt(outreach.title) << "\", ";
      // 5th line JOURNAL with year
      oss << "" << aurostd::html2txt(journal) << ". ";

      //   oss  << endl;
      // EXTRA STUFF FOR PHP WEB
      // 7th line DOI
      if(outreach.doi.size())
        oss << "doi: " << outreach.doi;
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
      //     for(uint ix=0;ix<outreach.vextra_html.size();ix++)
      //       oss << outreach.vextra_html.at(ix) << endl;
      //     oss << "<br>";
      //     oss << endl; not needed in new SC20210616
    }
    // ***************************************************************************
    // JSON
    if(XHOST.vflag_control.flag("PRINT_MODE::JSON")) {
      oss << "{";
      if(outreach.newflag) {
        if(XHOST.vflag_control.flag("PRINT_MODE::NEW") && outreach.newflag) oss << "\"new\":true,";
        if(XHOST.vflag_control.flag("PRINT_MODE::DATE") && outreach.newflag && outreach.date.size()) oss << "\"date\":\"" << outreach.date << "\",";
      }
      if(XHOST.vflag_control.flag("PRINT_MODE::NUMBER")) {
        if(outreach.wnumber==0) {  // SHORTCUT for CHAPTER
          oss << "\"number\":" << outreach.wnumber << ",";
          oss << "\"chapter\":true,";
        } else {
          oss << "\"number\":" << outreach.wnumber << ",";
          oss << "\"article\":true,";
        }
      }
      // 3rd line AUTHORS
      aurostd::StringSubst(authors_json,"Csányi","Csanyi");
      oss << "\"authors\":" <<  authors_json << ",";
      // 4th line TITLE
      oss << "\"title\":\"" << outreach.title << "\",";
      // 5th line JOURNAL with year
      oss << "\"journal\":\"" << journal << "\",";
      // EXTRA STUFF FOR PHP WEB
      // 7th line DOI
      if(XHOST.vflag_control.flag("PRINT_MODE::DOI") && outreach.doi.size())
        oss << "\"doi\":\"" << outreach.doi << "\",";
      // 6th line PDF    
      if(XHOST.vflag_control.flag("PRINT_MODE::PDF") && outreach.pdf.size())
        oss << "\"pdf\":\"" << WEB_PDF << outreach.pdf << "\",";
      // 6th line SUPPLEMENTARY
      if(XHOST.vflag_control.flag("PRINT_MODE::PDF") && outreach.supplementary.size())
        oss << "\"supplementary\":\"" << WEB_PDF << outreach.supplementary << "\",";
     // 6th line SUPPLEMENTARY_URL
      if(XHOST.vflag_control.flag("PRINT_MODE::DOI") && outreach.supplementary_url.size())
        oss << "\"supplementary_url\":\"" << outreach.supplementary_url << "\",";
      // 6th line ARXIV
      if(outreach.arxiv.size())
        oss << "\"arxiv\":\"" << outreach.arxiv << "\",";
      // 6th line LINK
      if(outreach.link.size())
        oss << "\"link\":\"" << outreach.link << "\",";
      // 7th line BIBTEX
      if(XHOST.vflag_control.flag("PRINT_MODE::BIBTEX") && outreach.bibtex.size() && outreach.bibtex_journal.size() && outreach.doi.size() && outreach.title.size() && outreach.vauthor.size() ) { // needs to have bibtex pdf doi title and authors
        bibtex2file(outreach.bibtex,authors_txt,outreach.title,outreach.bibtex_journal,outreach.bibtex_volume,outreach.bibtex_issue,outreach.bibtex_pages,outreach.bibtex_year,outreach.abstract,outreach.doi,string(WEB_BIBTEX_FILE+outreach.bibtex+".txt"));
        oss << "\"bibtex\":\"" << WEB_BIBTEX << outreach.bibtex << ".txt" << "\",";
      }
      // Xth line EXTRA
      if(outreach.vextra_html.size()) {
        oss << "\"vextra_html\":\"" ;
        for(uint ix=0;ix<outreach.vextra_html.size();ix++) {
          string extrahtml = outreach.vextra_html.at(ix);
          aurostd::StringSubst(extrahtml, "\"", "\\\"");
          if(!aurostd::substring2bool(outreach.vextra_html.at(ix),WEB_PDF) && !aurostd::substring2bool(outreach.vextra_html.at(ix),WEB_DOI)) {
            oss << extrahtml << (compact?" ":"<br>");
          } else {
            if(aurostd::substring2bool(outreach.vextra_html.at(ix),WEB_DOI) && XHOST.vflag_control.flag("PRINT_MODE::DOI")) {
              oss << extrahtml << (compact?" ":"<br>");
            } else {
              if(aurostd::substring2bool(outreach.vextra_html.at(ix),WEB_PDF) && XHOST.vflag_control.flag("PRINT_MODE::PDF")) {
                oss << extrahtml << (compact?" ":"<br>");
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
    if(XHOST.vflag_control.flag("PRINT_MODE::HTML")) {
      oss << "<li>";
      if(outreach.newflag) {
        if(XHOST.vflag_control.flag("PRINT_MODE::NEW") && outreach.newflag) oss << "<blink><b>NEW </b></blink>";
        if(XHOST.vflag_control.flag("PRINT_MODE::DATE") && outreach.newflag && outreach.date.size()) oss << "<b>" << outreach.date << "</b>";
        oss << " "; //endl;
      }
      if(XHOST.vflag_control.flag("PRINT_MODE::NUMBER")) {
        if(outreach.wnumber==0) {
          oss << "<span class=\"pubNumber\">" << "Chapter" << ". </span>";
        } else {
          oss << "<span class=\"pubNumber\">" << outreach.wnumber << ". </span>";
        }
      }
      // 3rd line AUTHORS
      aurostd::StringSubst(authors_txt,"Csányi","Csanyi");
      oss << "<span class=\"pubAuthors\">" <<  authors_txt << "</span>, " << (compact?" ":newline.str());
      // 4th line TITLE
      oss << "<span class=\"pubTitle\">" << outreach.title << "</span>, " << (compact?" ":newline.str());
      // 5th line JOURNAL with year
      oss << "<span class=\"pubJournal\">" << journal << "</span>. " << (compact?" ":newline.str());
      // EXTRA STUFF FOR PHP WEB
      // 7th line DOI
      if(XHOST.vflag_control.flag("PRINT_MODE::DOI") && outreach.doi.size())
        oss << "[<a href="+WEB_DOI+outreach.doi << "><b>doi" << "=" << outreach.doi << "</b></a>] " << (compact?" ":newline.str());
      // 6th line PDF    
      if(XHOST.vflag_control.flag("PRINT_MODE::PDF") && outreach.pdf.size())
        oss << "[<a href="+WEB_PDF+outreach.pdf << "><b>pdf</b></a>] " << (compact?" ":newline.str());
      // 6th line SUPPLEMENTARY
      if(XHOST.vflag_control.flag("PRINT_MODE::PDF") && outreach.supplementary.size())
        oss << "[<a href="+WEB_PDF+outreach.supplementary << "><b>suppl</b></a>] " << (compact?" ":newline.str());
      // 6th line SUPPLEMENTARY_URL
      if(XHOST.vflag_control.flag("PRINT_MODE::DOI") && outreach.supplementary_url.size())
        oss << "[<a href="+outreach.supplementary_url << "><b>suppl</b></a>] " << (compact?" ":newline.str());
      // 6th line ARXIV
      if(outreach.arxiv.size())
        oss << "[<a href="+outreach.arxiv << "><b>arxiv</b></a>] " << (compact?" ":newline.str());
      // 6th line LINK
      if(outreach.link.size())
        oss << "[<a href="+outreach.link << "><b>link</b></a>] " << (compact?" ":newline.str());
      // 7th line BIBTEX
      if(XHOST.vflag_control.flag("PRINT_MODE::BIBTEX") && outreach.bibtex.size() && outreach.bibtex_journal.size() && outreach.doi.size() && outreach.title.size() && outreach.vauthor.size() ) { // needs to have bibtex pdf doi title and authors
        bibtex2file(outreach.bibtex,authors_txt,outreach.title,outreach.bibtex_journal,outreach.bibtex_volume,outreach.bibtex_issue,outreach.bibtex_pages,outreach.bibtex_year,outreach.abstract,outreach.doi,string(WEB_BIBTEX_FILE+outreach.bibtex+".txt"));
        oss << "[<a href=" << WEB_BIBTEX << outreach.bibtex << ".txt" << "><b>bibtex</b></a>] " << (compact?" ":newline.str());
      }
      // Xth line EXTRA
      for(uint ix=0;ix<outreach.vextra_html.size();ix++) {
        if(!aurostd::substring2bool(outreach.vextra_html.at(ix),WEB_PDF) && !aurostd::substring2bool(outreach.vextra_html.at(ix),WEB_DOI)) {
          oss << outreach.vextra_html.at(ix) << (compact?" ":newline.str());
        } else {
          if(aurostd::substring2bool(outreach.vextra_html.at(ix),WEB_DOI) && XHOST.vflag_control.flag("PRINT_MODE::DOI")) {
            oss << outreach.vextra_html.at(ix) << (compact?" ":newline.str());
          } else {
            if(aurostd::substring2bool(outreach.vextra_html.at(ix),WEB_PDF) && XHOST.vflag_control.flag("PRINT_MODE::PDF")) {
              oss << outreach.vextra_html.at(ix) << (compact?" ":newline.str());
            }
          }
        }
      }
      oss << "</li>";
    }
    // ***************************************************************************
    // LATEX
    if(XHOST.vflag_control.flag("PRINT_MODE::LATEX")) {
      // 3rd line AUTHORS
      authors_txt=aurostd::html2latex(authors_txt)+", ";
      oss << " " << authors_txt;
      // 4th line TITLE
      string title="\\textit{"+aurostd::html2latex(outreach.title)+"}, ";
      oss << " " << title;
      // 5th line JOURNAL with year
      oss << " " << aurostd::html2latex(journal) << ". ";
      // EXTRA STUFF FOR LATEX
      bool link=FALSE;
      // 7th line DOI
      if(XHOST.vflag_control.flag("PRINT_MODE::DOI") && outreach.doi.size()) {
        string doi="";
        doi="\\ifthenelse {\\equal{\\hyperlinks}{true}}{";
        doi+="{\\newline \\textsf{\\href{http://dx.doi.org/"+outreach.doi+"}{DOI: "+aurostd::html2latex(outreach.doi)+"}}}";
        doi+="}{{\\newline \\textsf{DOI: "+aurostd::html2latex(outreach.doi)+"}}}";
        oss << " " << doi;
        link=TRUE;
      } else {
        if(outreach.wnumber!=0 && outreach.wnumber!=27 && outreach.wnumber!=14 &&
            outreach.wnumber!=11 && outreach.wnumber!=8 && outreach.wnumber!=2 &&
            outreach.wnumber!=1) {
          string doi=""; // ="{\\newline \\textsf{\\href{http://dx.doi.org/}{DOI: N/A}}}";
          oss << " " << doi;
          link=TRUE;
        }
      }
      // 7th line PDF
      if((XHOST.vflag_control.flag("PRINT_MODE::PDF") || (!XHOST.vflag_control.flag("PRINT_MODE::PDF") && XHOST.vflag_control.flag("PRINT_MODE::DOI"))) && outreach.doi.length()==0 && outreach.pdf.size()) {
        string pdf="";
        pdf="\\ifthenelse {\\equal{\\hyperlinks}{true}}{\\textsf{[\\href{"+WEB_PDF+outreach.pdf+"}{pdf}]}}{}";
        oss << " " << pdf;
        link=TRUE;
      }
      if(link==FALSE) {;}; // if(!XHOST.QUIET_CERR) cerr  << wnumber << endl;
      // LATEX EXTRA
      for(uint ix=0;ix<outreach.vextra_latex.size();ix++) {
        if(!aurostd::substring2bool(outreach.vextra_latex.at(ix),WEB_PDF) && !aurostd::substring2bool(outreach.vextra_latex.at(ix),WEB_DOI)) {
          oss << " " << outreach.vextra_latex.at(ix);
        } else {
          if(aurostd::substring2bool(outreach.vextra_latex.at(ix),WEB_DOI) && XHOST.vflag_control.flag("PRINT_MODE::DOI")) {
            oss << " " << outreach.vextra_latex.at(ix);
          } else {
            if(aurostd::substring2bool(outreach.vextra_latex.at(ix),WEB_PDF) && XHOST.vflag_control.flag("PRINT_MODE::PDF")) {
              oss << " " <<  outreach.vextra_latex.at(ix);
            }
          }
        }
      }
      // oss << " % ART" << wnumber << " - aflow " << string(AFLOW_VERSION) << endl; OLD NEEDS endls SC20210616
      oss << " % ART" << wnumber << " - aflow " << string(AFLOW_VERSION);// << endl; OLD NEEDS endls SC20210616
    }
  }
  // ***************************************************************************
  // PRESENTATION
  if(aurostd::substring2bool(outreach.type,"PRESENTATION")) {
    if(XHOST.vflag_control.flag("PRINT_MODE::TXT")) {
      if(outreach._isinvited && outreach.type=="PRESENTATION_TALK") oss << "Invited talk";
      if(outreach._isinvited && outreach.type=="PRESENTATION_SEMINAR") oss << "Invited seminar";
      if(outreach._isinvited && outreach.type=="PRESENTATION_COLLOQUIUM") oss << "Invited colloquium";
      if(outreach._isinvited && outreach.type=="PRESENTATION_PLENARY") oss << "Plenary Speaker";
      if(outreach._isinvited && outreach.type=="PRESENTATION_KEYNOTE") oss << "Keynote Speaker";
      if(outreach._isinvited && outreach.type=="PRESENTATION_TUTORIAL") oss << "Tutorial";
      if(outreach._isinvited && outreach.type=="PRESENTATION_PANELIST") oss << "Invited panelist";
      if(outreach._isonline) { oss << " (online):"; } else {  oss << ":"; }
      oss << " " << aurostd::latex2txt(outreach.title) << "; ";
      oss << "" << aurostd::latex2txt(outreach.place) << ", ";
      oss << "" << aurostd::latex2txt(outreach.date) << ". ";
      if(outreach.link.size()) {
        oss << "link: " << outreach.link;
      }
    }
    if(XHOST.vflag_control.flag("PRINT_MODE::HTML")) {
      if(outreach._isinvited && outreach.type=="PRESENTATION_TALK") oss << "<b> Invited talk";
      if(outreach._isinvited && outreach.type=="PRESENTATION_SEMINAR") oss << "<b> Invited seminar";
      if(outreach._isinvited && outreach.type=="PRESENTATION_COLLOQUIUM") oss << "<b> Invited colloquium";
      if(outreach._isinvited && outreach.type=="PRESENTATION_PLENARY") oss << "<b> Plenary Speaker";
      if(outreach._isinvited && outreach.type=="PRESENTATION_KEYNOTE") oss << "<b> Keynote Speaker";
      if(outreach._isinvited && outreach.type=="PRESENTATION_TUTORIAL") oss << "<b> Tutorial";
      if(outreach._isinvited && outreach.type=="PRESENTATION_PANELIST") oss << "<b> Invited panelist";
      if(outreach._isonline) { oss << " (online):</b>"; } else {  oss << ":</b>"; }
      //    oss << "<br>" << endl;
      oss << "<i> " << aurostd::latex2html(outreach.title) << "</i>; " << endl;
      oss << "" << aurostd::latex2html(outreach.place) << ", " << endl;
      oss << "" << aurostd::latex2html(outreach.date) << ". ";// << endl;
      oss << "<br>" << endl;
      if(outreach.link.size()) {
        oss << "[<a href="+outreach.link << "><b>link" << "=" << outreach.link << "</b></a>] " << (compact?" ":newline.str());
        oss << "<br>" << endl;
      }
    }
    if(XHOST.vflag_control.flag("PRINT_MODE::LATEX")) {
      if(outreach._isinvited && outreach.type=="PRESENTATION_TALK") oss << "\\textbf{Invited talk";
      if(outreach._isinvited && outreach.type=="PRESENTATION_SEMINAR") oss << "\\textbf{Invited seminar";
      if(outreach._isinvited && outreach.type=="PRESENTATION_COLLOQUIUM") oss << "\\textbf{Invited colloquium";
      if(outreach._isinvited && outreach.type=="PRESENTATION_PLENARY") oss << "\\textbf{Plenary Speaker";
      if(outreach._isinvited && outreach.type=="PRESENTATION_KEYNOTE") oss << "\\textbf{Keynote Speaker";
      if(outreach._isinvited && outreach.type=="PRESENTATION_TUTORIAL") oss << "\\textbf{Tutorial";
      if(outreach._isinvited && outreach.type=="PRESENTATION_PANELIST") oss << "\\textbf{Invited panelist";
      if(outreach._isonline) { oss << " (online):}"; } else {  oss << ":}"; }
      oss << endl;
      // 3rd line AUTHORS
      if(0) {
        for(uint iauth=0;iauth<outreach.vauthor.size();iauth++)
          oss << outreach.vauthor.at(iauth) << ", ";
        oss << endl;
      }
      oss << "\\textit{" << outreach.title << "}; " << endl;
      oss << "" << outreach.place << ", " << endl;
      oss << "" << outreach.date << ". ";// << endl;
      if(outreach.link.size()) {
        oss << endl;
        oss << "\\ifthenelse {\\equal{\\hyperlinks}{true}}{";
        oss << "{\\newline \\textsf{\\href{" << outreach.link << "}{LINK: " << aurostd::html2latex(outreach.link) << "}}}";
        oss << "}{{\\newline \\textsf{LINK: " << aurostd::html2latex(outreach.link) << "}}}";
      }
    }
  }
  return oss;
}

// ******************************************************************************************************************************************************
void voutreach_print_publications_EXTRA(ostream& oss);

vector<_outreach> voutreach_global_list;
uint voutreach_global_max_year=0,voutreach_global_min_year=9999;


// ******************************************************************************************************************************************************
// ******************************************************************************************************************************************************
// ARTICLES
// ******************************************************************************************************************************************************
// ******************************************************************************************************************************************************

uint voutreach_remove_duplicate(vector<_outreach>& voutreach) {
  bool LOCAL_VERBOSE=1;//FALSE;
  if(LOCAL_VERBOSE && !XHOST.QUIET_CERR) cerr << "voutreach_remove_duplicate: voutreach.size()=" <<  voutreach.size() << endl;
  for(uint i=0;i<voutreach.size();i++)
    for(uint j=i+1;j<voutreach.size();j++)
      if(i!=j && voutreach.at(j).year==voutreach.at(i).year) { // same year
        // if(LOCAL_VERBOSE && !XHOST.QUIET_CERR) cerr << "voutreach_remove_duplicate: found same year: (i,j)=(" << i << "," << j << ")  " << voutreach.size() << endl;
        if(voutreach.at(j).vauthor.size()==voutreach.at(i).vauthor.size()) { // same number of authors
          // if(LOCAL_VERBOSE && !XHOST.QUIET_CERR) cerr << "voutreach_remove_duplicate: found same vauthor.size: (i,j)=(" << i << "," << j << ")  " << voutreach.size() << endl;
          if(voutreach.at(j).vauthor.at(0)==voutreach.at(i).vauthor.at(0)) { // same first author
            // if(LOCAL_VERBOSE && !XHOST.QUIET_CERR) cerr << "voutreach_remove_duplicate: found same vauthor.at(0): (i,j)=(" << i << "," << j << ")  " << voutreach.size() << endl;
            if(voutreach.at(j).vauthor.at(voutreach.at(j).vauthor.size()-1)==voutreach.at(i).vauthor.at(voutreach.at(i).vauthor.size()-1)) { // same last author
              // if(LOCAL_VERBOSE && !XHOST.QUIET_CERR) cerr << "voutreach_remove_duplicate: found same vauthor.at(N): (i,j)=(" << i << "," << j << ")  " << voutreach.size() << endl;
              vector<string> vwordi;aurostd::string2tokens(voutreach.at(i).title,vwordi," ");
              vector<string> vwordj;aurostd::string2tokens(voutreach.at(j).title,vwordj," ");
              if(vwordj.size()==vwordi.size()) { // same spaces
                // if(LOCAL_VERBOSE && !XHOST.QUIET_CERR) cerr << "voutreach_remove_duplicate: found same number of title words: (i,j)=(" << i << "," << j << ")  " << vwordj.size() << endl;
                if(aurostd::abs(voutreach.at(i).title.length()-voutreach.at(j).title.length())<3) {
                  if(LOCAL_VERBOSE && !XHOST.QUIET_CERR) cerr << "voutreach_remove_duplicate: found similar title length: (i,j)=(" << i << "," << j << ")  " << voutreach.at(i).title.length() << "," << voutreach.at(j).title.length() << endl;
                  { // ACT
                    if(LOCAL_VERBOSE && !XHOST.QUIET_CERR) cerr << "voutreach_remove_duplicate: found same year: (i,j)=(" << i << "," << j << ")  " << voutreach.size() << endl;
                    if(LOCAL_VERBOSE && !XHOST.QUIET_CERR) cerr << "voutreach_remove_duplicate: found same vauthor.size: (i,j)=(" << i << "," << j << ")  " << voutreach.size() << endl;
                    if(LOCAL_VERBOSE && !XHOST.QUIET_CERR) cerr << "voutreach_remove_duplicate: found same vauthor.at(0): (i,j)=(" << i << "," << j << ")  " << voutreach.size() << endl;
                    if(LOCAL_VERBOSE && !XHOST.QUIET_CERR) cerr << "voutreach_remove_duplicate: found same vauthor.at(N): (i,j)=(" << i << "," << j << ")  " << voutreach.size() << endl;
                    if(LOCAL_VERBOSE && !XHOST.QUIET_CERR) cerr << "voutreach_remove_duplicate: found same number of title words: (i,j)=(" << i << "," << j << ")  " << vwordj.size() << endl;
                    if(LOCAL_VERBOSE && !XHOST.QUIET_CERR) cerr << "voutreach_remove_duplicate: found similar title length: (i,j)=(" << i << "," << j << ")  " << voutreach.at(i).title.length() << "," << voutreach.at(j).title.length() << endl;
                    if(LOCAL_VERBOSE && !XHOST.QUIET_CERR) cerr << "voutreach_remove_duplicate: removing (i,j)=(" << i << "," << j << ")  title1=" << voutreach.at(j).title << " | " << voutreach.at(i).title << endl;
                    if(i!=j) {voutreach.erase(voutreach.begin()+j);i=0;j=0;}
                    if(LOCAL_VERBOSE && !XHOST.QUIET_CERR) cerr << "voutreach_remove_duplicate: Article: found duplicate (i,j)=(" << i << "," << j << ")  " << voutreach.size() << endl;
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
  if(1) { // PATENTS
    oss << "<li><span class=\"pubYears\">Patents</span></li>" << endl;
    oss << "<li>" << endl;
    oss << "<span class=\"pubAuthors\">G. Ceder, C. Fischer, K. Tibbetts, D. Morgan, and S. Curtarolo</span>, " << endl;
    oss << "<span class=\"pubTitle\">Systems and Methods for predicting materials properties</span>,"  << endl;
    oss << "<span class=\"pubJournal\">US Patent #7292958 (2007)</span>. [<a HREF=" << WEB_PDF << "PAT7292958.pdf>pdf</a>]" << endl;
    oss << "</li>" << endl;
  }
  // **********************************************************************************************
  if(1) { // THESES
    oss << "<li><span class=\"pubYears\">Theses</span></li>" << endl;
    oss << "<li>" << endl;
    oss << "<span class=\"pubNumber\">4. </span><span class=\"pubAuthors\">S. Curtarolo</span>," << endl;
    oss << "<span class=\"pubTitle\">Coarse-Graining and Data Mining Approaches to the Prediction of Structures and their Dynamics</span>, " << endl;
    oss << "<span class=\"pubJournal\">Ph.D., Massachusetts Institute of Technology (2003)</span>." << endl;
    oss << "[<a href=" << WEB_PDF << "scurtarolo4.pdf>pdf</a>]<br>" << endl;
    oss << "</li>" << endl;
    oss << "<li>" << endl;
    oss << "<span class=\"pubNumber\">3. </span><span class=\"pubAuthors\">S. Curtarolo</span>," << endl;
    oss << "<span class=\"pubTitle\">Adsorption Problems investigated with Computer Simulation</span>, " << endl;
    oss << "<span class=\"pubJournal\">M.S., Pennsylvania State University (1999)</span>." << endl;
    oss << "[<a href=" << WEB_PDF << "scurtarolo3.pdf>pdf</a>]<br>" << endl;
    oss << "</li>" << endl;
    oss << "<li>" << endl;
    oss << "<span class=\"pubNumber\">2. </span><span class=\"pubAuthors\">S. Curtarolo</span>," << endl;
    oss << "<span class=\"pubTitle\">Influenza della rugosita' sul prewetting di Neon su Magnesio</span>, " << endl;
    oss << "<span class=\"pubJournal\">Laurea in Fisica, Universita` di Padova (1998)</span>." << endl;
    oss << "[<a href=" << WEB_PDF << "scurtarolo2.pdf>pdf</a>] (in Italian). <br>" << endl;
    oss << "</li>" << endl;
    oss << "<li>" << endl;
    oss << "<span class=\"pubNumber\">1. </span><span class=\"pubAuthors\">S. Curtarolo</span>," << endl;
    oss << "<span class=\"pubTitle\">Approccio analitico e numerico allo studio degli adattatori dielettrici</span>, " << endl;
    oss << "<span class=\"pubJournal\">Laurea in Ingegneria Elettronica, Universita` di Padova (1995)</span>." << endl;
    oss << "[<!A href=" << WEB_PDF << "scurtarolo1.pdf>pdf, not available<!/a>] (in Italian). <br>" << endl;
    oss << "</li>" << endl;
  }
}

// ***************************************************************************
void SystemReferences(const string& system_in,vector<uint>& voutreach_wnumber) {  // ADD REFERENCES
  voutreach_wnumber.clear();

  string system=KBIN::VASP_PseudoPotential_CleanName(system_in);
  vector<_outreach> voutreach;
  voutreach_load(voutreach,"PUBLICATIONS");
  // if(!XHOST.QUIET_CERR) cerr << "LOADED " << voutreach.size() << " articles" << endl;
  for(uint iart=0;iart<voutreach.size();iart++)
    if(voutreach.at(iart).valloy.size()>0) 
      for(uint itoken=0;itoken<voutreach.at(iart).valloy.size();itoken++) 
        if(system==voutreach.at(iart).valloy.at(itoken)) 
          voutreach_wnumber.push_back(voutreach.at(iart).wnumber);
}

// ***************************************************************************

void SystemReferences(const string& system_in,vector<string>& vrefs,bool AUTHOR_ETAL) {  // ADD REFERENCES
  vrefs.clear();
  vector<uint> voutreach_wnumber;
  SystemReferences(system_in,voutreach_wnumber);
  vector<_outreach> voutreach;
  voutreach_load(voutreach,"PUBLICATIONS");

  stringstream oss;

  for(uint iarticle=0;iarticle<voutreach_wnumber.size();iarticle++) {
    for(uint i=0;i<voutreach.size();i++) {
      if(voutreach_wnumber.at(iarticle)==voutreach.at(i).wnumber) {
        oss.clear();oss.str(std::string());
        // GOT THE ARTICLE
        //	if(iarticle+1<10)  oss << "<sup>&nbsp;&nbsp;" << iarticle+1 << "</sup> "; // LABEL
        //	if(iarticle+1>=10) oss << "<sup>" << iarticle+1 << "</sup> "; // LABEL
        // AUTHORS
        string authors;
        if(!AUTHOR_ETAL || voutreach.at(i).vauthor.size()<=4) { // all authors OR <=4
          for(uint iauth=0;iauth<voutreach.at(i).vauthor.size();iauth++) {
            authors+=voutreach.at(i).vauthor.at(iauth);
            if(voutreach.at(i).vauthor.size()==2 && iauth==voutreach.at(i).vauthor.size()-2) authors+=" and ";
            if(voutreach.at(i).vauthor.size()!=2 && iauth==voutreach.at(i).vauthor.size()-2) authors+=", and ";
            if(iauth!=voutreach.at(i).vauthor.size()-2) authors+=", ";
          }
          //	aurostd::StringSubst(authors,"Csányi","Csanyi");
        } else { // first et al..
          authors+=voutreach.at(i).vauthor.at(0)+" et al., ";
        }
        oss << authors;// << endl;
        // TITLE
        string title=voutreach.at(i).title;
        title=""+title+", ";
        oss << title;
        // JOURNAL with year
        string journal=voutreach.at(i).journal;
        oss << "" << journal << ". ";
        vrefs.push_back(aurostd::html2txt(oss.str()));
        oss.clear();oss.str(std::string());
      }
    }
  }
}

// ******************************************************************************************************************************************************

bool SystemInSystems(const string& system,const string& systems) {
  vector<string> tokens;
  aurostd::string2tokens(systems,tokens,",");
  for(uint i=0;i<tokens.size();i++)
    if(system==tokens.at(i))
      return TRUE;
  return FALSE;
}

bool SystemInSystems(const string& system,const vector<string>& vsystems) {
  for(uint i=0;i<vsystems.size();i++)
    if(system==vsystems.at(i))
      return TRUE;
  return FALSE;
}

// ******************************************************************************************************************************************************

bool AlloyAlphabeticLIBRARY(const string& s) {
  // never anymore this problem
  if(s=="MoMg"||s=="NaMg"||s=="NbMg"||s=="OsMg"||s=="PbMg"||s=="RbMg"||s=="RbPd"||s=="RbPt"||s=="ReMg"||s=="RhMg"||s=="RuMg"||s=="ScMg"||s=="SiPd") return FALSE;
  if(s=="SiPt"||s=="SnMg"||s=="SnPd"||s=="SnPt"||s=="SrMg"||s=="SrPd"||s=="SrPt"||s=="TaMg"||s=="TiMg"||s=="VMg"||s=="WMg"||s=="YMg"||s=="ZnMg"||s=="ZrMg") return FALSE;
  return TRUE;
}

// ******************************************************************************************************************************************************
// ******************************************************************************************************************************************************
// ARTICLES/PRESENTATIONS
// ******************************************************************************************************************************************************
// ******************************************************************************************************************************************************

void voutreach_print(uint _mode,ostream& oss,string what2print) {
  string message = "";
  aurostd::xoption vflag=XHOST.vflag_outreach;
  bool LDEBUG=(FALSE || XHOST.DEBUG);
  uint mode=_mode;

  // ******************************************************************************************************************************************************
  // ARTICLES
  if(aurostd::substring2bool(what2print,"ARTICLE") || aurostd::substring2bool(what2print,"PUBLICATION")) {
    vector<_outreach> voutreach;
    voutreach_load(voutreach,"PUBLICATIONS");
    bool LOCAL_VERBOSE=1;//FALSE;
    if(!XHOST.QUIET_CERR) cerr << "LOADED " << voutreach.size() << " articles" << endl;
    bool HTRESOURCE_MODE_PHP_PREAMBLE=FALSE;

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
    vector<string> vauthor,vkeyword,valloy;
    //  vector<string> vargument;
    bool flag_simple=FALSE;

    if(LOCAL_VERBOSE && !XHOST.QUIET_CERR) cerr << "voutreach_print [begin]" << endl;
    //   LDEBUG=TRUE;

    if(LDEBUG) oss << "XHOST.vflag_control.flags" << endl;
    if(LDEBUG) oss << "XHOST.vflag_control.flag(\"PRINT_MODE::JSON\")=" << XHOST.vflag_control.flag("PRINT_MODE::JSON") << endl;
    if(LDEBUG) oss << "XHOST.vflag_control.flag(\"PRINT_MODE::HTML\")=" << XHOST.vflag_control.flag("PRINT_MODE::HTML") << endl;
    if(LDEBUG) oss << "XHOST.vflag_control.flag(\"PRINT_MODE::TXT\")=" << XHOST.vflag_control.flag("PRINT_MODE::TXT") << endl;
    if(LDEBUG) oss << "XHOST.vflag_control.flag(\"PRINT_MODE::LATEX\")=" << XHOST.vflag_control.flag("PRINT_MODE::LATEX") << endl;
    if(LDEBUG) oss << "XHOST.vflag_control.flag(\"PRINT_MODE::YEAR\")=" << XHOST.vflag_control.flag("PRINT_MODE::YEAR") << endl;
    if(LDEBUG) oss << "XHOST.vflag_control.flag(\"PRINT_MODE::DOI\")=" << XHOST.vflag_control.flag("PRINT_MODE::DOI") << endl;
    if(LDEBUG) oss << "XHOST.vflag_control.flag(\"PRINT_MODE::BIBTEX\")=" << XHOST.vflag_control.flag("PRINT_MODE::BIBTEX") << endl;
    if(LDEBUG) oss << "XHOST.vflag_control.flag(\"PRINT_MODE::EXTRA\")=" << XHOST.vflag_control.flag("PRINT_MODE::EXTRA") << endl;
    if(LDEBUG) oss << "XHOST.vflag_control.flag(\"PRINT_MODE::NUMBER\")=" << XHOST.vflag_control.flag("PRINT_MODE::NUMBER") << endl;
    if(LDEBUG) oss << "XHOST.vflag_control.flag(\"PRINT_MODE::HYPERLINKS\")=" << XHOST.vflag_control.flag("PRINT_MODE::HYPERLINKS") << endl;
    if(LDEBUG) oss << "XHOST.vflag_control.flag(\"PRINT_MODE::NOTE\")=" << XHOST.vflag_control.flag("PRINT_MODE::NOTE") << endl;
    if(LDEBUG) oss << "XHOST.vflag_control.flag(\"PRINT_MODE::NEW\")=" << XHOST.vflag_control.flag("PRINT_MODE::NEW") << endl;
    if(LDEBUG) oss << "XHOST.vflag_control.flag(\"PRINT_MODE::DATE\")=" << XHOST.vflag_control.flag("PRINT_MODE::DATE") << endl;

    // big hack in chinatown ! if Frisco gives me nothing...
    if(XHOST.vflag_control.flag("PHP::PUBS_ALLOY")) {
      if(XHOST.vflag_control.getattachedscheme("PHP::PUBS_ALLOY")=="--print=html") {
        mode=HTRESOURCE_MODE_PHP_ALLOY;
        XHOST.vflag_control.flag("PRINT_MODE::LATEX",FALSE);
        XHOST.vflag_control.flag("PRINT_MODE::TXT",FALSE);
        XHOST.vflag_control.flag("PRINT_MODE::JSON",FALSE);
        XHOST.vflag_control.flag("PRINT_MODE::HTML",TRUE); // DEFAULT
        valloy.push_back("AFLOW");
        valloy.push_back("AFLOWLIB");
        valloy.push_back("nmatHT");
      }
    } 

    if(mode==HTRESOURCE_MODE_NONE) {mode=HTRESOURCE_MODE_PHP_AUTHOR;} // by default

    if(XHOST.vflag_control.flag("CV::AUTHOR") || mode==HTRESOURCE_MODE_PHP_AUTHOR)  {
      mode=HTRESOURCE_MODE_PHP_AUTHOR;
      if(LOCAL_VERBOSE && !XHOST.QUIET_CERR) cerr << "voutreach_print: [" << XHOST.vflag_control.flag("CV::AUTHOR") << "]" << endl;
      if(XHOST.vflag_control.flag("CV::AUTHOR")) {
        if(LOCAL_VERBOSE && !XHOST.QUIET_CERR) cerr << "voutreach_print: [" << XHOST.vflag_control.getattachedscheme("CV::AUTHOR") << "]" << endl;
        aurostd::string2tokens(XHOST.vflag_control.getattachedscheme("CV::AUTHOR"),vauthor,",");
        if(vauthor.size()==0) {
          message = "No authors specified.";
          throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message);
        }
      }
    }

    if(LOCAL_VERBOSE && !XHOST.QUIET_CERR) cerr << "voutreach_print: vauthor.size()=" << vauthor.size() << endl;

    if(XHOST.vflag_control.flag("PHP::PUBS_ALLOY") || mode==HTRESOURCE_MODE_PHP_ALLOY)  {
      if(XHOST.vflag_control.getattachedscheme("PHP::PUBS_ALLOY")!="--print=html") {
        mode=HTRESOURCE_MODE_PHP_ALLOY;
        aurostd::string2tokens(XHOST.vflag_control.getattachedscheme("PHP::PUBS_ALLOY"),valloy,",");
        if(valloy.size()==0) {
          message = "valloy.size() == 0";
          throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message);
        }
      }
    }
    if(XHOST.vflag_control.flag("PHP::PUBS_KEYWORD") || mode==HTRESOURCE_MODE_PHP_THRUST)  {
      mode=HTRESOURCE_MODE_PHP_THRUST;
      aurostd::string2tokens(XHOST.vflag_control.getattachedscheme("PHP::PUBS_KEYWORD"),vkeyword,",");
      if(vkeyword.size()==0) {
        message = "No keywords specified.";
        throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message);
      }
      for(uint i=0;i<vkeyword.size();i++) if(!XHOST.QUIET_CERR) cerr << "vkeyword.at(" << i << ")=" <<  vkeyword.at(i) << endl;
    }

    if(!XHOST.QUIET_CERR) cerr << "voutreach_print: vauthor.size()=" <<  vauthor.size() << endl;
    if(!XHOST.QUIET_CERR) cerr << "voutreach_print: valloy.size()=" << valloy.size() << endl;
    if(!XHOST.QUIET_CERR) cerr << "voutreach_print: vkeyword.size()=" << vkeyword.size() << endl;
    if(!XHOST.QUIET_CERR) cerr << "voutreach_print: mode=" << mode << endl;
    // ************************************************************************************************************************
    // ************************************************************************************************************************
    // OPERATION ON MULTIPLE VAUTHOR
    if(vauthor.size()>0) {
      // if(!XHOST.QUIET_CERR) cerr << "OPERATION ON MULTIPLE VAUTHOR" << endl;
      uint recent_articles=0;
      for(uint iauthor=0;iauthor<vauthor.size();iauthor++) {
        voutreach_local.clear();
        if(!XHOST.vflag_control.flag("PRINT_MODE::JSON")) oss << endl; // just some  leeway but not in json SC20210616
        if(HTRESOURCE_MODE_PHP_PREAMBLE==TRUE && mode==HTRESOURCE_MODE_PHP_AUTHOR) {
          oss << "<?php" << endl;
          if(flag_simple==FALSE)
            oss << "if($author==\"" << vauthor.at(iauthor) << "\")" << endl;
          oss << "{" << endl;
          oss << "?>" << endl;
        }
        if(XHOST.vflag_control.flag("PRINT_MODE::LATEX")) {
          oss << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl;
          oss << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl;
          // oss << "% PUBLICATIONS" << endl;
          oss << "" << endl;
          oss << "\\section{Publications} \\label{publications}" << endl;
          oss << "[Articles can be accessed through the embedded link. Only published, submitted and ``in press'' articles are listed. The list might be slightly out of order.]" << endl;

          // oss << "\\begin{list}{\\labelitemi}{\\leftmargin=1em}" << endl;
        }

        if(LOCAL_VERBOSE && !XHOST.QUIET_CERR) cerr << "voutreach_print: LOADED " << voutreach.size() << " articles" << endl;

        for(uint i=0;i<voutreach.size();i++)
          if(aurostd::substring2bool(voutreach.at(i).vauthor,vauthor.at(iauthor)))
            voutreach_local.push_back(voutreach.at(i));  // push authors

        if(LOCAL_VERBOSE && !XHOST.QUIET_CERR) cerr << "voutreach_print: LOADED_LOCAL " << voutreach_local.size() << " articles" << endl;
        voutreach_remove_duplicate(voutreach_local);
        if(LOCAL_VERBOSE && !XHOST.QUIET_CERR) cerr << "voutreach_print: LOADED_DUPLICATE " << voutreach_local.size() << " articles" << endl;
        voutreach_rsort_wnumber(voutreach_local); // sort TOP numbers first
        if(LOCAL_VERBOSE && !XHOST.QUIET_CERR) cerr << "voutreach_print: LOADED_SORTED " << voutreach_local.size() << " articles" << endl;

        if(mode==HTRESOURCE_MODE_PHP_AUTHOR && !XHOST.vflag_control.flag("PRINT_MODE::LATEX") && !XHOST.vflag_control.flag("PRINT_MODE::TXT") && !XHOST.vflag_control.flag("PRINT_MODE::JSON")) {oss << "<br><br>" << endl;}
        if(XHOST.vflag_control.flag("PRINT_MODE::LATEX")) {oss << endl;}

        // ************************************************************************************************************************
        // HTML/JSON MODE
        if(XHOST.vflag_control.flag("PRINT_MODE::HTML") || XHOST.vflag_control.flag("PRINT_MODE::JSON")) {
          //	if(mode==HTRESOURCE_MODE_PHP_AUTHOR && !XHOST.vflag_control.flag("PRINT_MODE::LATEX")&& !XHOST.vflag_control.flag("PRINT_MODE::TXT"))
          if(XHOST.vflag_control.flag("PRINT_MODE::HTML")) {
            voss.push_back("<li><span class=\"aflow\"> AFLOW V"+string(AFLOW_VERSION)+" </span></li>");
          }
          if(XHOST.vflag_control.flag("PRINT_MODE::JSON")) {
            voss.push_back("[");
          }
          for(uint year=voutreach_global_max_year;year>=voutreach_global_min_year;year--) {
            bool flag_year=FALSE;
            for(uint i=0;i<voutreach_local.size()&&!flag_year;i++)
              if(year==voutreach_local.at(i).year)
                flag_year=TRUE;
            if(flag_year && recent_articles<AUTHOR_RECENT_ARTICLES && year>=voutreach_global_max_year-AUTHOR_RECENT_YEARS && year<=voutreach_global_max_year) {
              if(XHOST.vflag_control.flag("PRINT_MODE::YEAR")) {
                if(XHOST.vflag_control.flag("PRINT_MODE::HTML")) {
                  voss.push_back("<li><span class=\"pubYears\">"+aurostd::utype2string(year)+"</span></li>");
                }
              }
              for(uint i=0;i<voutreach_local.size();i++) {
                string wnumber=aurostd::utype2string(voutreach_local.at(i).wnumber);
                if(wnumber.length()<3) wnumber="0"+wnumber;
                if(wnumber.length()<2) wnumber="0"+wnumber;
                if(voutreach_local.at(i).year==year && recent_articles<AUTHOR_RECENT_ARTICLES && year>=voutreach_global_max_year-AUTHOR_RECENT_YEARS && year<=voutreach_global_max_year) {
                  voss.push_back(voutreach_local.at(i).print_string());
                  recent_articles++;
                } 
              }
            }
          }
          if(XHOST.vflag_control.flag("PRINT_MODE::HTML"))  { // this is only for HTML, json has all preables elsewhere
            if(HTRESOURCE_MODE_PHP_PREAMBLE==TRUE && mode==HTRESOURCE_MODE_PHP_AUTHOR) {
              voss.push_back("\n <?php");
              voss.push_back("}");
              voss.push_back("?>");
            }
            for(uint j=0;j<voss.size();j++) oss << voss.at(j) << endl;
          }
          if(XHOST.vflag_control.flag("PRINT_MODE::JSON")) {
            voss.push_back("]");
            for(uint j=0;j<voss.size();j++) {
              oss << voss.at(j);
              if(j>0 && j<voss.size()-2) oss << ",";  
              //oss << endl;
            }
            oss << endl;
          }
          if(XHOST.vflag_control.flag("PRINT_MODE::EXTRA")) voutreach_print_publications_EXTRA(oss);
        } // HTML/JSON MODE
        // ************************************************************************************************************************

        // ************************************************************************************************************************
        /* OLD MODE
        // LATEX MODE
        if(XHOST.vflag_control.flag("PRINT_MODE::LATEX")) {
        if(!XHOST.vflag_control.flag("PRINT_MODE::YEAR")) oss << "\\begin{itemize}" << endl;

        for(uint year=voutreach_global_max_year;year>=voutreach_global_min_year;year--) {
        bool flag_year=FALSE;
        for(uint i=0;i<voutreach_local.size()&&!flag_year;i++)
        if(year==voutreach_local.at(i).year)
        flag_year=TRUE;
        if(flag_year && year<=voutreach_global_max_year) {
        oss << endl;
        if(XHOST.vflag_control.flag("PRINT_MODE::YEAR")) oss << "\\ \\\\  \\textbf{ \\color{blue}{[" << year << "]}}" << endl;
        if(XHOST.vflag_control.flag("PRINT_MODE::YEAR")) oss << "\\begin{itemize}" << endl;
        for(uint i=0;i<voutreach_local.size();i++) {
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
        if(XHOST.vflag_control.flag("PRINT_MODE::LATEX")) {
          if(!XHOST.vflag_control.flag("PRINT_MODE::YEAR")) voss.push_back("\\begin{itemize}");
          for(uint year=voutreach_global_max_year;year>=voutreach_global_min_year;year--) {
            bool flag_year=FALSE;
            for(uint i=0;i<voutreach_local.size()&&!flag_year;i++)
              if(year==voutreach_local.at(i).year)
                flag_year=TRUE;
            if(flag_year && year<=voutreach_global_max_year) {
              voss.push_back("");
              if(XHOST.vflag_control.flag("PRINT_MODE::YEAR")) voss.push_back("\\ \\\\  \\textbf{ \\color{blue}{["+aurostd::utype2string(year)+"]}}");
              if(XHOST.vflag_control.flag("PRINT_MODE::YEAR")) voss.push_back("\\begin{itemize}");
              for(uint i=0;i<voutreach_local.size();i++) {
                string wnumber=aurostd::utype2string(voutreach_local.at(i).wnumber);
                if(wnumber.length()<3) wnumber="0"+wnumber;
                if(wnumber.length()<2) wnumber="0"+wnumber;
                if(voutreach_local.at(i).year==year && year>=voutreach_global_max_year-AUTHOR_RECENT_YEARS && year<=voutreach_global_max_year) {
                  // print 1 article
                  if(!XHOST.vflag_control.flag("PRINT_MODE::NUMBER")) {
                    voss.push_back("\\item[$\\bullet$]{}"+voutreach_local.at(i).print_string());
                  } else {
                    if(voutreach_local.at(i).wnumber==0) {
                      voss.push_back("\\item[\\textbf{Chapter.\\,}]{}"+voutreach_local.at(i).print_string());
                    } else {
                      voss.push_back("\\item[\\textbf{"+aurostd::utype2string(voutreach_local.at(i).wnumber)+".\\,}]{}"+voutreach_local.at(i).print_string());
                    }
                  }
                } 
              }
              if(XHOST.vflag_control.flag("PRINT_MODE::YEAR")) voss.push_back("\\end{itemize}");
            }
          }
          voss.push_back("");
          if(!XHOST.vflag_control.flag("PRINT_MODE::YEAR")) voss.push_back("\\end{itemize}");
          //    voss.push_back("\\end{list}");
          voss.push_back("");
          voss.push_back("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%");
          voss.push_back("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%");
          for(uint j=0;j<voss.size();j++) oss << voss.at(j) << endl;
        } // LATEX MODE
        // ************************************************************************************************************************

        // ************************************************************************************************************************
        // OLD MODE
        /* TXT MODE
           if(XHOST.vflag_control.flag("PRINT_MODE::TXT")) {
           for(uint year=voutreach_global_max_year;year>=voutreach_global_min_year;year--) {
           bool flag_year=FALSE;
           for(uint i=0;i<voutreach_local.size()&&!flag_year;i++)
           if(year==voutreach_local.at(i).year)
           flag_year=TRUE;
           if(flag_year && year<=voutreach_global_max_year) {
           oss << endl;
           if(XHOST.vflag_control.flag("PRINT_MODE::YEAR")) oss << "[" << year << "]" << endl;
           for(uint i=0;i<voutreach_local.size();i++) {
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
        if(XHOST.vflag_control.flag("PRINT_MODE::TXT")) {
          for(uint year=voutreach_global_max_year;year>=voutreach_global_min_year;year--) {
            bool flag_year=FALSE;
            for(uint i=0;i<voutreach_local.size()&&!flag_year;i++)
              if(year==voutreach_local.at(i).year)
                flag_year=TRUE;
            if(flag_year && year<=voutreach_global_max_year) {
              if(XHOST.vflag_control.flag("PRINT_MODE::YEAR")){  voss.push_back("");voss.push_back("["+aurostd::utype2string(year)+"]");}
              for(uint i=0;i<voutreach_local.size();i++) {
                string wnumber=aurostd::utype2string(voutreach_local.at(i).wnumber);
                if(wnumber.length()<3) wnumber="0"+wnumber;
                if(wnumber.length()<2) wnumber="0"+wnumber;
                if(voutreach_local.at(i).year==year && year>=voutreach_global_max_year-AUTHOR_RECENT_YEARS && year<=voutreach_global_max_year) {
                  // print 1 article
                  if(XHOST.vflag_control.flag("PRINT_MODE::NUMBER")) {
                    voss.push_back("["+aurostd::utype2string(voutreach_local.at(i).wnumber)+".] "+voutreach_local.at(i).print_string());
                  } else {
                    voss.push_back(voutreach_local.at(i).print_string());
                  }
                }
              }
            }
          }
          voss.push_back("");
          voss.push_back("");
          voss.push_back("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%");
          voss.push_back("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%");
          for(uint j=0;j<voss.size();j++) oss << voss.at(j) << endl;
        } // TXT MODE
        // ************************************************************************************************************************
      } // iauthor
    }
    // ************************************************************************************************************************
    // ************************************************************************************************************************
    // OPERATION ON MULTIPLE VTHRUST
    if(vkeyword.size()>0) {
      // if(!XHOST.QUIET_CERR) cerr << "OPERATION ON MULTIPLE VTHRUST" << endl;
      uint recent_articles=0;
      voutreach_local.clear();
      for(uint ithrust=0;ithrust<vkeyword.size();ithrust++) {
        oss << endl;
        if(HTRESOURCE_MODE_PHP_PREAMBLE==TRUE && mode==HTRESOURCE_MODE_PHP_THRUST) {
          oss << "<?php" << endl;
          if(flag_simple==FALSE)
            oss << "if($thrust==\"" << vkeyword.at(ithrust) << "\")" << endl;
          oss << "{" << endl;
          oss << "?>" << endl;
        }
        for(uint i=0;i<voutreach.size();i++) {
          // if(!XHOST.QUIET_CERR) cerr << "voutreach.at(" << i << ").vkeyword.size()=" << voutreach.at(i).vkeyword.size() << endl;
          for(uint j=0;j<voutreach.at(i).vkeyword.size();j++)
            if(voutreach.at(i).vkeyword.at(j)==vkeyword.at(ithrust))
              voutreach_local.push_back(voutreach.at(i));  // push thrust
        }
      } // ithrust
      voutreach_remove_duplicate(voutreach_local);
      voutreach_sort_year(voutreach_local); // sort TOP numbers first

      //    if(mode==HTRESOURCE_MODE_PHP_THRUST) {oss << "<!br><!br>" << endl;}
      // ************************************************************************************************************************
      // HTML MODE
      if(XHOST.vflag_control.flag("PRINT_MODE::HTML")) {
        if(mode==HTRESOURCE_MODE_PHP_THRUST) {
          if(LOCAL_VERBOSE && !XHOST.QUIET_CERR)
            oss << "<li><span class=\"aflow\">" << "AFLOW V" << string(AFLOW_VERSION) << " </span></li>" << endl;
          for(uint year=voutreach_global_max_year;year>=voutreach_global_min_year;year--) {
            bool flag_year=FALSE;
            for(uint i=0;i<voutreach_local.size()&&!flag_year;i++)
              if(year==voutreach_local.at(i).year)
                flag_year=TRUE;
            if(flag_year && recent_articles<THRUST_RECENT_ARTICLES && year>=voutreach_global_max_year-THRUST_RECENT_YEARS && year<=voutreach_global_max_year) {
              for(uint i=0;i<voutreach_local.size();i++) {
                string wnumber=aurostd::utype2string(voutreach_local.at(i).wnumber);
                if(wnumber.length()<3) wnumber="0"+wnumber;
                if(wnumber.length()<2) wnumber="0"+wnumber;
                if(voutreach_local.at(i).year==year && recent_articles<THRUST_RECENT_ARTICLES && year>=voutreach_global_max_year-THRUST_RECENT_YEARS && year<=voutreach_global_max_year) {
                  XHOST.vflag_control.flag("PRINT_MODE::NUMBER",FALSE);
                  oss << voutreach_local.at(i) << endl;
                  recent_articles++;
                }
              }
            }
          }
          if(HTRESOURCE_MODE_PHP_PREAMBLE==TRUE && mode==HTRESOURCE_MODE_PHP_THRUST) {
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
    if(valloy.size()>0) {
      // if(!XHOST.QUIET_CERR) cerr << "OPERATION ON MULTIPLE VALLOY" << endl;
      //   uint recent_articles=0;
      voutreach_local.clear();
      for(uint ialloy=0;ialloy<valloy.size();ialloy++) {
        vector<uint> voutreach_wnumber;
        SystemReferences(KBIN::VASP_PseudoPotential_CleanName(valloy.at(ialloy)),voutreach_wnumber);
        for(uint iwnumber=0;iwnumber<voutreach_wnumber.size();iwnumber++) 
          for(uint iarticle=0;iarticle<voutreach.size();iarticle++)
            if(voutreach.at(iarticle).wnumber==voutreach_wnumber.at(iwnumber))
              //	  if(voutreach.at(iart)!=65 && voutreach.at(iart)!=75)     // remove aflow and aflowlib papers
              voutreach_local.push_back(voutreach.at(iarticle));                     // remove aflow and aflowlib papers
        // if(!XHOST.QUIET_CERR) cerr << voutreach_local.size() << endl;

      } // ialloy
      voutreach_remove_duplicate(voutreach_local);
      voutreach_sort_year(voutreach_local); // sort TOP numbers first

      //    if(mode==HTRESOURCE_MODE_PHP_ALLOY) {oss << "<!br><!br>" << endl;}
      // ************************************************************************************************************************
      // PHP AND WEB MODE
      if(mode==HTRESOURCE_MODE_PHP_ALLOY) {
        if(LOCAL_VERBOSE && !XHOST.QUIET_CERR)
          if(XHOST.vflag_control.flag("PRINT_MODE::HTML") || mode==HTRESOURCE_MODE_PHP_AUTHOR)
            oss << "<li><span class=\"aflow\">" << "AFLOW V" << string(AFLOW_VERSION) << " </span></li>" << endl;
        oss << "<ol class=\"reference\">" << endl;
        for(uint i=0;i<voutreach_local.size();i++) {
          //	string wnumber=voutreach_local.at(i).wnumber;
          //	if(wnumber.length()<2) wnumber="0"+wnumber;
          XHOST.vflag_control.flag("PRINT_MODE::NUMBER",FALSE);
          voutreach_local.at(i).vextra_html.clear(); // delete it
          oss << voutreach_local.at(i) << endl;
          // recent_articles++;
        }
        oss << "</ol>" << endl;
        if(HTRESOURCE_MODE_PHP_PREAMBLE==TRUE && mode==HTRESOURCE_MODE_PHP_ALLOY) {
          oss << "<?php" << endl;
          oss << "}" << endl;
          oss << "?>" << endl;
        }
      } // PHP WEB
      // ************************************************************************************************************************
    }

    //   for(uint i=0;i<voutreach_local.size();i++)
    //   oss << voutreach_local.at(i) << endl;
  }

  // ******************************************************************************************************************************************************
  // PRESENTATIOS
  if(aurostd::substring2bool(what2print,"PRESENTATION")) {
    vector<_outreach> voutreach;
    voutreach_load(voutreach,"PRESENTATIONS");
    if(!XHOST.QUIET_CERR) cerr << "LOADED " << voutreach.size() << " presentations" << endl;

    if(mode==HTRESOURCE_MODE_NONE) {mode=HTRESOURCE_MODE_PHP_AUTHOR;}

    vector<_outreach> voutreach_local;
    vector<string> vauthors;
    vauthors.push_back("Curtarolo");

    uint noutreach=0;

    // ************************************************************************************************************************
    // TXT MODE
    if(XHOST.vflag_control.flag("PRINT_MODE::TXT")) {
      for(uint iauthor=0;iauthor<vauthors.size();iauthor++) {
        voutreach_local.clear();

        oss << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl;
        for(uint i=0;i<voutreach.size();i++)
          if(aurostd::substring2bool(voutreach.at(i).vauthor,vauthors.at(iauthor)))
            voutreach_local.push_back(voutreach.at(i));	
        for(uint year=MAX_YEAR_PRESENTATIONS;year>=1998;year--) {
          bool flag_year=FALSE;
          for(uint i=0;i<voutreach_local.size()&&!flag_year;i++)
            if(year==voutreach_local.at(i).year) flag_year=TRUE;
          if(flag_year) {
            for(uint i=0;i<voutreach_local.size();i++)
              if(voutreach_local.at(i).year==year) {
                // print 1 outreach
                oss << endl;
                oss << voutreach_local.at(i) << endl;
              }
          }
        }
        oss << "" << endl;
      } // iauthor
    } // TXT MODE
    // ************************************************************************************************************************
    // LATEX MODE
    if(XHOST.vflag_control.flag("PRINT_MODE::LATEX")) {
      for(uint iauthor=0;iauthor<vauthors.size();iauthor++) {
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


        for(uint i=0;i<voutreach.size();i++)
          if(aurostd::substring2bool(voutreach.at(i).vauthor,vauthors.at(iauthor)))
            voutreach_local.push_back(voutreach.at(i));

        for(uint year=MAX_YEAR_PRESENTATIONS;year>=1998;year--) {
          bool flag_year=FALSE;
          for(uint i=0;i<voutreach_local.size()&&!flag_year;i++)
            if(year==voutreach_local.at(i).year) flag_year=TRUE;
          if(flag_year) {
            //	oss << endl << "<b><font size=\"3\"><font COLOR=blue><i>" << year << ".</i> </font></b><br><br>" << endl;
            for(uint i=0;i<voutreach_local.size();i++)
              if(voutreach_local.at(i).year==year) {
                // print 1 outreach
                oss << endl;
                if(vauthors.at(iauthor)=="Curtarolo") 
                  oss << "\\item[\\textbf{" << voutreach_local.size()-noutreach++ << ".\\,}]{}" << endl;
                oss << voutreach_local.at(i) << endl;
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
    if(XHOST.vflag_control.flag("PRINT_MODE::HTML")) {
      if(mode==HTRESOURCE_MODE_PHP_AUTHOR) {
        for(uint iauthor=0;iauthor<vauthors.size();iauthor++) {
          voutreach_local.clear();

          for(uint i=0;i<voutreach.size();i++)
            if(aurostd::substring2bool(voutreach.at(i).vauthor,vauthors.at(iauthor)))
              voutreach_local.push_back(voutreach.at(i));

          for(uint year=MAX_YEAR_PRESENTATIONS;year>=1998;year--) {
            bool flag_year=FALSE;
            for(uint i=0;i<voutreach_local.size()&&!flag_year;i++)
              if(year==voutreach_local.at(i).year) flag_year=TRUE;
            if(flag_year) {
              //	oss << endl << "<b><font size=\"3\"><font COLOR=blue><i>" << year << ".</i> </font></b><br><br>" << endl;
              for(uint i=0;i<voutreach_local.size();i++)
                if(voutreach_local.at(i).year==year) {
                  // print 1 outreach
                  oss << endl;
                  if(vauthors.at(iauthor)=="Curtarolo") {
                    oss << "<font COLOR=red><b>" << voutreach_local.size()-noutreach++  << ". </b></font>";
                  }
                  XHOST.vflag_control.flag("PRINT_MODE::LATEX",FALSE);
                  XHOST.vflag_control.flag("PRINT_MODE::TXT",FALSE);
                  XHOST.vflag_control.flag("PRINT_MODE::JSON",FALSE);
                  XHOST.vflag_control.flag("PRINT_MODE::HTML",TRUE);              
                  // [OBSOLETE]		voutreach_local.at(i).print_mode="HTML";
                  oss << voutreach_local.at(i) << endl;
                }
            }
          }
        }
      }
    } // HTML MODE
  }  // PRESENTATIONS
}

// ******************************************************************************************************************************************************
void voutreach_print_everything(ostream& oss,const vector<string>& vitems,string msg1,string msg2,string sectionlabel) {
  bool flag_ITMZ=FALSE;
  bool flag_LIST=TRUE;
  if(!XHOST.QUIET_CERR) cerr << "LOADED " << vitems.size() << " " << msg1 << endl;

  if(XHOST.vflag_control.flag("PRINT_MODE::LATEX")) {
    oss << "% " << msg2 << endl;
    oss << sectionlabel << endl;
    if(flag_ITMZ) oss << "\\begin{itemize}" << endl;
    if(flag_LIST) oss << "\\begin{list}{\\labelitemi}{\\leftmargin=1em}" << endl;
    for(uint i=0;i<vitems.size();i++)
      oss << "\\item{}" << " " << vitems.at(i) << "  % " << msg2 << endl;
    if(flag_ITMZ) oss << "\\end{itemize}" << endl;
    if(flag_LIST) oss << "\\end{list}" << endl;
    //    oss << "" << endl;
  }
}

// ******************************************************************************************************************************************************
vector<_outreach> voutreach_presentations;
vector<_outreach> voutreach_publications;
int voutreach_call=0;
uint voutreach_load(vector<_outreach>& voutreach,string what2print) {
  // if(!XHOST.QUIET_CERR) cerr << voutreach_call++ << endl;
  // check cache
  bool LDEBUG=0;
  voutreach.clear();
  vector<string> valabel;
  vector<string> vjlabel;
  if(aurostd::substring2bool(what2print,"ARTICLE") || aurostd::substring2bool(what2print,"PUBLICATION")) {
    if(voutreach_publications.size()) {
      for(uint i=0;i<voutreach_publications.size();i++)
        voutreach.push_back(voutreach_publications.at(i));
      return voutreach.size();
    }  
  }
  if(aurostd::substring2bool(what2print,"PRESENTATION")) {
    if(voutreach_presentations.size()) {
      for(uint i=0;i<voutreach_presentations.size();i++)
        voutreach.push_back(voutreach_presentations.at(i));
      return voutreach.size();
    }  
  }
  // no cache
  // <b><font size="3"><font COLOR=blue><i>2010-current.</i> </font></b><br><br>
  _outreach ptmp;

  if(LDEBUG) cerr << "voutreach_load [0]" << endl;

  string cv2open="f144468a7ccc2d3a72ba44000715efdb";
  // if(!XHOST.QUIET_CERR) cerr << XHOST.vflag_control.getattachedscheme("CV::AUTHOR") << endl;
  if(aurostd::toupper(XHOST.vflag_control.getattachedscheme("CV::AUTHOR"))=="CURTAROLO" || aurostd::toupper(XHOST.vflag_control.getattachedscheme("CV::AUTHOR"))=="SCURTAROLO") cv2open="f144468a7ccc2d3a72ba44000715efdb";
  // [OBSOLETE] if(aurostd::toupper(XHOST.vflag_control.getattachedscheme("CV::AUTHOR"))=="OSES" || aurostd::toupper(XHOST.vflag_control.getattachedscheme("CV::AUTHOR"))=="COSES") cv2open="d0f1b0e47f178ae627a388d3bf65d2d2";
  // [OBSOLETE] if(aurostd::toupper(XHOST.vflag_control.getattachedscheme("CV::AUTHOR"))=="TOHER" || aurostd::toupper(XHOST.vflag_control.getattachedscheme("CV::AUTHOR"))=="CTOHER") cv2open="decf00ca3ad2fe494eea8e543e929068";

  if(LDEBUG) cerr << "voutreach_load [1]" << endl;
  vector<string> vpres,ktokens,tokens;
  if(aurostd::substring2bool(what2print,"ARTICLE") || aurostd::substring2bool(what2print,"PUBLICATION")) aurostd::string2vectorstring(aurostd::RemoveComments(init::InitGlobalObject(cv2open)),vpres);
  if(aurostd::substring2bool(what2print,"PRESENTATION")) aurostd::string2vectorstring(aurostd::RemoveComments(init::InitGlobalObject(cv2open)),vpres);
  if(aurostd::substring2bool(what2print,"EDUCATION")) aurostd::string2vectorstring(aurostd::RemoveComments(init::InitGlobalObject(cv2open)),vpres);
  if(aurostd::substring2bool(what2print,"RESEARCH")) aurostd::string2vectorstring(aurostd::RemoveComments(init::InitGlobalObject(cv2open)),vpres);
  if(aurostd::substring2bool(what2print,"ACADEMIC")) aurostd::string2vectorstring(aurostd::RemoveComments(init::InitGlobalObject(cv2open)),vpres);
  if(aurostd::substring2bool(what2print,"SERVICEOUTSIDE")) aurostd::string2vectorstring(aurostd::RemoveComments(init::InitGlobalObject(cv2open)),vpres);
  if(aurostd::substring2bool(what2print,"SERVICEINSIDE")) aurostd::string2vectorstring(aurostd::RemoveComments(init::InitGlobalObject(cv2open)),vpres);
  if(aurostd::substring2bool(what2print,"TEACHING")) aurostd::string2vectorstring(aurostd::RemoveComments(init::InitGlobalObject(cv2open)),vpres);
  if(aurostd::substring2bool(what2print,"ADVISING")) aurostd::string2vectorstring(aurostd::RemoveComments(init::InitGlobalObject(cv2open)),vpres);
  if(aurostd::substring2bool(what2print,"PATENTS")) aurostd::string2vectorstring(aurostd::RemoveComments(init::InitGlobalObject(cv2open)),vpres);
  if(aurostd::substring2bool(what2print,"PRESS")) aurostd::string2vectorstring(aurostd::RemoveComments(init::InitGlobalObject(cv2open)),vpres);
  if(aurostd::substring2bool(what2print,"AWARDS")) aurostd::string2vectorstring(aurostd::RemoveComments(init::InitGlobalObject(cv2open)),vpres);

   if(LDEBUG) cerr << "voutreach_load [2]" << endl;

  string iline,jline,kline;
  uint i,j,k,l;

  for(i=0;i<vpres.size();i++) {
    //  cout << vpres.at(i) << endl;

    iline=vpres.at(i);

    if(LDEBUG) cerr << "voutreach_load [2.line]=" << iline << endl;

    iline=aurostd::StringSubst(iline," ","");
    if(iline=="OBJECT={") { // found an object
      for(j=i;i<vpres.size();j++) {
        jline=vpres.at(j);
        jline=aurostd::StringSubst(jline," ","");
        if(jline.at(0)=='}') break;
      }
      // now I have an object between i+1 and j-1
      ptmp.clear();
      for(k=i+1;k<j;k++) {
        kline=vpres.at(k);
        kline=aurostd::StringSubst(kline," ","");
        aurostd::string2tokens(kline,ktokens,"=");
        aurostd::string2tokens(vpres.at(k),tokens,"=");
        // if(!XHOST.QUIET_CERR) cerr << tokens.at(1) << endl;
        string object;
        if(tokens.size()>0 && ktokens.size()>0) {
          object=tokens.at(1);
          for(l=2;l<tokens.size();l++)
            object+="="+tokens.at(l);
          // if(!XHOST.QUIET_CERR) cerr << object << endl;
          for(l=0;l<5&&object.size();l++) { // cleanup
            if(object.at(0)==' ' || object.at(0)=='\"')
              object=object.substr(1,object.length());
            if(object.at(object.length()-1)==' ' || object.at(object.length()-1)==';' || 
                object.at(object.length()-1)=='"') object=object.substr(0,object.length()-1);
          }
          // if(!XHOST.QUIET_CERR) cerr << object << endl;
          if(ktokens.at(0)=="YEAR" || ktokens.at(0)=="year") ptmp.year=aurostd::string2utype<uint>(ktokens.at(1)); // check year
          if(ktokens.at(0)=="TYPE" || ktokens.at(0)=="type") { // check type
            if(object=="TALK" || object=="PRESENTATION_TALK") {ptmp.type="PRESENTATION_TALK";ptmp._isinvited=TRUE;}
            if(object=="SEMINAR" || object=="PRESENTATION_SEMINAR") {ptmp.type="PRESENTATION_SEMINAR";ptmp._isinvited=TRUE;}
            if(object=="COLLOQUIUM" || object=="PRESENTATION_COLLOQUIUM") {ptmp.type="PRESENTATION_COLLOQUIUM";ptmp._isinvited=TRUE;}
            if(object=="KEYNOTE" || object=="PRESENTATION_KEYNOTE") {ptmp.type="PRESENTATION_KEYNOTE";ptmp._isinvited=TRUE;}
            if(object=="PLENARY" || object=="PRESENTATION_PLENARY") {ptmp.type="PRESENTATION_PLENARY";ptmp._isinvited=TRUE;}
            if(object=="TUTORIAL" || object=="PRESENTATION_TUTORIAL") {ptmp.type="PRESENTATION_TUTORIAL";ptmp._isinvited=TRUE;}
            if(object=="PANELIST" || object=="PRESENTATION_PANELIST") {ptmp.type="PRESENTATION_PANELIST";ptmp._isinvited=TRUE;}
            if(object=="CONTRIBUTED" || object=="PRESENTATION_CONTRIBUTED") ptmp.type="PRESENTATION_CONTRIBUTED";
            if(object=="POSTER" || object=="PRESENTATION_POSTER") ptmp.type="PRESENTATION_POSTER";
            if(object=="ARTICLE" || object=="article") ptmp.type="ARTICLE";
            if(object=="ALABEL" || object=="alabel") ptmp.type="ALABEL";
            if(object=="JLABEL" || object=="jlabel") ptmp.type="JLABEL";
            if(object=="PUBLICATION" || object=="publication") ptmp.type="ARTICLE";
            if(object=="EDUCATION" || object=="education") ptmp.type="EDUCATION";
            if(object=="RESEARCH" || object=="research") ptmp.type="RESEARCH";
            if(object=="ACADEMIC" || object=="academic") ptmp.type="ACADEMIC";
            if(object=="SERVICEOUTSIDE" || object=="serviceoutside") ptmp.type="SERVICEOUTSIDE";
            if(object=="SERVICEINSIDE" || object=="serviceinside") ptmp.type="SERVICEINSIDE";
            if(object=="TEACHING" || object=="teaching") ptmp.type="TEACHING";
            if(object=="ADVISING" || object=="advising") ptmp.type="ADVISING";
            if(object=="PATENTS" || object=="patents") ptmp.type="PATENTS";
            if(object=="PRESS" || object=="press") ptmp.type="PRESS";
            if(object=="AWARDS" || object=="awards") ptmp.type="AWARDS";
          }
          if(ktokens.at(0)=="ONLINE" || ktokens.at(0)=="online") { // check ONLINESS
            if(object=="YES" || object=="yes" || object=="TRUE" || object=="true" || object=="1") {ptmp._isonline=TRUE;}
            if(object=="NO" || object=="no" || object=="FALSE" || object=="false" || object=="0") {ptmp._isonline=FALSE;}
          }

          if(ktokens.at(0)=="TITLE" || ktokens.at(0)=="title") ptmp.title=object; // check title
          if(ktokens.at(0)=="JOURNAL" || ktokens.at(0)=="journal") ptmp.journal=object; // check journal
          if(ktokens.at(0)=="LINK" || ktokens.at(0)=="link") ptmp.link=object; // check link
          if(ktokens.at(0)=="ARXIV" || ktokens.at(0)=="arxiv") ptmp.arxiv=object; // check arxiv
          if(ktokens.at(0)=="SUPPLEMENTARY" || ktokens.at(0)=="supplementary") ptmp.supplementary=object; // check supplementary
	  if(ktokens.at(0)=="SUPPLEMENTARY_URL" || ktokens.at(0)=="supplementary_url") ptmp.supplementary_url=object; // check supplementary_url
          if(ktokens.at(0)=="BIBTEX" || ktokens.at(0)=="bibtex") ptmp.bibtex=object; // check bibtex
          if(ktokens.at(0)=="BIBTEX_JOURNAL" || ktokens.at(0)=="bibtex_journal") ptmp.bibtex_journal=object; // check bibtex_journal
          if(ktokens.at(0)=="BIBTEX_VOLUME" || ktokens.at(0)=="bibtex_volume") ptmp.bibtex_volume=object; // check bibtex_volume
          if(ktokens.at(0)=="BIBTEX_ISSUE" || ktokens.at(0)=="bibtex_issue") ptmp.bibtex_issue=object; // check bibtex_issue
          if(ktokens.at(0)=="BIBTEX_PAGES" || ktokens.at(0)=="bibtex_pages") ptmp.bibtex_pages=object; // check bibtex_pages
          if(ktokens.at(0)=="BIBTEX_YEAR" || ktokens.at(0)=="bibtex_year") ptmp.bibtex_year=object; // check bibtex_year
          if(ktokens.at(0)=="PLACE" || ktokens.at(0)=="place") ptmp.place=object; // check place
          if(ktokens.at(0)=="DATE" || ktokens.at(0)=="date") ptmp.date=object; // check date
          if(ktokens.at(0)=="HOST" || ktokens.at(0)=="host") ptmp.host=object; // check host
          if(ktokens.at(0)=="ABSTRACT" || ktokens.at(0)=="abstract") ptmp.abstract=object; // check abstract
          if(ktokens.at(0)=="PDF" || ktokens.at(0)=="pdf") ptmp.pdf=object; // check pdf
          if(ktokens.at(0)=="DOI" || ktokens.at(0)=="doi") ptmp.doi=object; // check doi
          if(ktokens.at(0)=="NEW" || ktokens.at(0)=="new") ptmp.newflag=TRUE; // check doi
          if(ktokens.at(0)=="NEWFLAG" || ktokens.at(0)=="newflag") ptmp.newflag=TRUE; // check doi
          if(ktokens.at(0)=="WNUMBER" || ktokens.at(0)=="wnumber") ptmp.wnumber=aurostd::string2utype<uint>(object); // check wnumber  
          if(ktokens.at(0)=="AUTHOR" || ktokens.at(0)=="author") aurostd::string2tokensAdd(object,ptmp.vauthor,",");
          if(ktokens.at(0)=="EXTRA_HTML" || ktokens.at(0)=="extra_html") ptmp.vextra_html.push_back(object);
          if(ktokens.at(0)=="EXTRA_LATEX" || ktokens.at(0)=="extra_latex") ptmp.vextra_latex.push_back(object);
          if(ktokens.at(0)=="KEYWORD" || ktokens.at(0)=="keyword") aurostd::string2tokensAdd(object,ptmp.vkeyword,",");
          if(ktokens.at(0)=="SPONSOR" || ktokens.at(0)=="sponsor") aurostd::string2tokensAdd(object,ptmp.vsponsor,",");
          if(ktokens.at(0)=="ALLOY" || ktokens.at(0)=="alloy") aurostd::string2tokensAdd(object,ptmp.valloy,",");
          if(ktokens.at(0)=="ALABEL" || ktokens.at(0)=="alabel") aurostd::string2tokensAdd(object,valabel,",");
          if(ktokens.at(0)=="JLABEL" || ktokens.at(0)=="jlabel") aurostd::string2tokensAdd(object,vjlabel,",");
        }
      }

      if(ptmp.bibtex=="" && ptmp.bibtex_volume!="") { stringstream oss;oss << "curtarolo:art" << ptmp.wnumber;ptmp.bibtex=oss.str(); }

      if(ptmp.type=="") {
        for(k=i;k<=j;k++) 
          if(!XHOST.QUIET_CERR) cerr << "entry=[" << vpres.at(k) << "]" << endl;
      } else {
        if((aurostd::substring2bool(ptmp.type,"ARTICLE") || aurostd::substring2bool(ptmp.type,"PUBLICATION")) && (aurostd::substring2bool(what2print,"ARTICLE") || aurostd::substring2bool(what2print,"PUBLICATION")))
          voutreach.push_back(ptmp);
        if(aurostd::substring2bool(ptmp.type,"PRESENTATION") && aurostd::substring2bool(what2print,"PRESENTATION")) {
          if(ptmp.vauthor.size()==0) ptmp.vauthor.push_back("S. Curtarolo");  // for SC CV
          voutreach.push_back(ptmp);
        }
        if(aurostd::substring2bool(ptmp.type,"EDUCATION") && aurostd::substring2bool(what2print,"EDUCATION")) voutreach.push_back(ptmp);
        if(aurostd::substring2bool(ptmp.type,"RESEARCH") && aurostd::substring2bool(what2print,"RESEARCH")) voutreach.push_back(ptmp);
        if(aurostd::substring2bool(ptmp.type,"ACADEMIC") && aurostd::substring2bool(what2print,"ACADEMIC")) voutreach.push_back(ptmp);
        if(aurostd::substring2bool(ptmp.type,"SERVICEOUTSIDE") && aurostd::substring2bool(what2print,"SERVICEOUTSIDE")) voutreach.push_back(ptmp);
        if(aurostd::substring2bool(ptmp.type,"SERVICEINSIDE") && aurostd::substring2bool(what2print,"SERVICEINSIDE")) voutreach.push_back(ptmp);
        if(aurostd::substring2bool(ptmp.type,"TEACHING") && aurostd::substring2bool(what2print,"TEACHING")) voutreach.push_back(ptmp);
        if(aurostd::substring2bool(ptmp.type,"ADVISING") && aurostd::substring2bool(what2print,"ADVISING")) voutreach.push_back(ptmp);
        if(aurostd::substring2bool(ptmp.type,"PATENTS") && aurostd::substring2bool(what2print,"PATENTS")) voutreach.push_back(ptmp);
        if(aurostd::substring2bool(ptmp.type,"PRESS") && aurostd::substring2bool(what2print,"PRESS")) voutreach.push_back(ptmp);
        if(aurostd::substring2bool(ptmp.type,"AWARDS") && aurostd::substring2bool(what2print,"AWARDS")) voutreach.push_back(ptmp);
      }
    }
  }

  if(LDEBUG) cerr << "voutreach_load [3]" << endl;

  for(uint i=0;i<voutreach.size();i++) {
    fixlabel(valabel,voutreach.at(i).vauthor);
    //   cout << voutreach.at(i).vauthor.size() << endl;
    for(uint j=0;j<voutreach.at(i).vauthor.size();j++) { // fix corresponding  authors
      aurostd::StringSubst(voutreach.at(i).vauthor.at(j),"*","");  //fix corresponding authors
      if(0)   {
        vector<string> tokens;
        aurostd::string2tokens(voutreach.at(i).vauthor.at(j),tokens," ");
        if(tokens.size()>2) {
          cout << tokens.size() << " [" << voutreach.at(i).vauthor.at(j) << "] ";
          for (k=0;k< tokens.size();k++)
            cout << tokens.at(k) << " ";
          cout << endl;
        }
      }
    }
  }

  if(LDEBUG) cerr << "voutreach_load [4]" << endl;

  for(uint i=0;i<voutreach.size();i++) fixlabel(valabel,voutreach.at(i).vextra_html);
  for(uint i=0;i<voutreach.size();i++) fixlabel(valabel,voutreach.at(i).vextra_latex);
  for(uint i=0;i<voutreach.size();i++) fixlabel(vjlabel,voutreach.at(i).journal);
  for(uint i=0;i<voutreach.size();i++) fixlabel(vjlabel,voutreach.at(i).bibtex_journal);
  for(uint i=0;i<voutreach.size();i++) fixlabel(vjlabel,voutreach.at(i).vextra_html);
  for(uint i=0;i<voutreach.size();i++) fixlabel(vjlabel,voutreach.at(i).vextra_latex);
  if(aurostd::substring2bool(ptmp.type,"ARTICLE") || aurostd::substring2bool(ptmp.type,"PUBLICATION"))  
    voutreach_remove_duplicate(voutreach);  // talks can be duplicated
  // SAVE the global one
  if(voutreach_global_list.size()==0) voutreach_global_list=voutreach;
  // CHECK MAX YEAR
  for(uint i=0;i<voutreach.size();i++) if(voutreach_global_max_year<voutreach.at(i).year) voutreach_global_max_year=voutreach.at(i).year;
  // CHECK MIN YEAR
  for(uint i=0;i<voutreach.size();i++) if(voutreach_global_min_year>voutreach.at(i).year) voutreach_global_min_year=voutreach.at(i).year;
  // DONE

  if(LDEBUG) cerr << "voutreach_load [5]" << endl;

  // save cache
  if(aurostd::substring2bool(ptmp.type,"ARTICLE") || aurostd::substring2bool(ptmp.type,"PUBLICATION")) {
    voutreach_publications.clear();
    for(uint i=0;i<voutreach.size();i++) {
      voutreach_publications.push_back(voutreach.at(i));
      // if(!XHOST.QUIET_CERR) cerr << voutreach.at(i) << endl;
    }
  } 
  if(LDEBUG) cerr << "voutreach_load [6]" << endl;

  if(aurostd::substring2bool(ptmp.type,"PRESENTATION")) {
    voutreach_presentations.clear();
    for(uint i=0;i<voutreach.size();i++)
      voutreach_presentations.push_back(voutreach.at(i));
  }

  if(LDEBUG) cerr << "voutreach_load [9]" << endl;

  return voutreach.size();
}


// ******************************************************************************************************************************************************
void HT_CHECK_GRANTS(ostream& oss) {//,const vector<string>& vitems,string msg1,string msg2,string sectionlabel)
  aurostd::xoption vflag=XHOST.vflag_outreach;
  vector<_outreach> voutreach;
  voutreach_load(voutreach,"PUBLICATIONS");
  if(!vflag.flag("GRANTS")) {
    string message = "No grants.";
    throw aurostd::xerror(__AFLOW_FILE__, __AFLOW_FUNC__, message);
  }
  if(!XHOST.QUIET_CERR) cerr << "LOADED " << voutreach.size() << " " << endl;
  string grant=vflag.getattachedscheme("GRANTS");
  //  if(mode==HTRESOURCE_MODE_NONE) {mode=HTRESOURCE_MODE_PHP_AUTHOR;} // by default

  if(XHOST.vflag_control.flag("PRINT_MODE::LATEX")) {
    for(uint i=0;i<voutreach.size();i++) {
      bool found=FALSE;
      for(uint j=0;j<voutreach.at(i).vsponsor.size()&&found==FALSE;j++)
        if(aurostd::substring2bool(voutreach.at(i).vsponsor.at(j),grant)) found=TRUE;
      if(found) {
        XHOST.vflag_control.flag("PRINT_MODE::TXT",TRUE);
        // [OBSOLETE]	voutreach.at(i).print_mode="TXT";
        oss << voutreach.at(i) << endl;// << endl;
      }
    } 
  }
}

// ******************************************************************************************************************************************************
bool ProcessPhpLatexCv(void) {
  aurostd::xoption vflag=XHOST.vflag_outreach;

  bool Arun=FALSE;
  vector<_outreach> voutreach;
  if(!XHOST.vflag_control.flag("PRINT_MODE::JSON") &&
      !XHOST.vflag_control.flag("PRINT_MODE::HTML") && 
      !XHOST.vflag_control.flag("PRINT_MODE::LATEX") && 
      !XHOST.vflag_control.flag("PRINT_MODE::TXT")) {
    XHOST.vflag_control.flag("PRINT_MODE::HTML",TRUE);
  }
  if(XHOST.vflag_control.flag("PRINT_MODE::JSON")) {XHOST.vflag_control.flag("PRINT_MODE::LATEX",FALSE);XHOST.vflag_control.flag("PRINT_MODE::TXT",FALSE);XHOST.vflag_control.flag("PRINT_MODE::HTML",FALSE);}
  if(XHOST.vflag_control.flag("PRINT_MODE::HTML")) {XHOST.vflag_control.flag("PRINT_MODE::LATEX",FALSE);XHOST.vflag_control.flag("PRINT_MODE::TXT",FALSE);XHOST.vflag_control.flag("PRINT_MODE::JSON",FALSE);}
  if(XHOST.vflag_control.flag("PRINT_MODE::LATEX")) {XHOST.vflag_control.flag("PRINT_MODE::TXT",FALSE);XHOST.vflag_control.flag("PRINT_MODE::HTML",FALSE);XHOST.vflag_control.flag("PRINT_MODE::JSON",FALSE);}
  if(XHOST.vflag_control.flag("PRINT_MODE::TXT")) {XHOST.vflag_control.flag("PRINT_MODE::HTML",FALSE);XHOST.vflag_control.flag("PRINT_MODE::LATEX",FALSE);XHOST.vflag_control.flag("PRINT_MODE::JSON",FALSE);}

  if(!Arun && XHOST.vflag_control.flag("PHP::CENTER_MISSION")) {
    center_print(HTRESOURCE_MODE_PHP_AUTHOR,cout);
    Arun=TRUE;
  }
  if(!Arun && XHOST.vflag_control.flag("PHP::PUBS_ALLOY")) {
    voutreach_print(HTRESOURCE_MODE_PHP_ALLOY,cout,"ARTICLES");
    Arun=TRUE;
  }
  if(!Arun && XHOST.vflag_control.flag("PHP::PUBS_KEYWORD")) {
    voutreach_print(HTRESOURCE_MODE_PHP_THRUST,cout,"ARTICLES");
    Arun=TRUE;
  }
  if(!Arun && XHOST.vflag_control.flag("PHP::ITALKS")) {
    voutreach_print(HTRESOURCE_MODE_PHP_AUTHOR,cout,"PRESENTATIONS");
    Arun=TRUE;
  }
  if(!Arun && XHOST.vflag_control.flag("CV::PUBS")) {
    voutreach_print(HTRESOURCE_MODE_PHP_AUTHOR,cout,"ARTICLES");
    Arun=TRUE;
  }
  if(!Arun && vflag.flag("GRANTS")) {
    HT_CHECK_GRANTS(cout);
    Arun=TRUE;
  } 
  if(!Arun && XHOST.vflag_control.flag("CV::ITALKS")) {
    voutreach_print(HTRESOURCE_MODE_PHP_AUTHOR,cout,"PRESENTATIONS");
    Arun=TRUE;
  }
  if(!Arun && XHOST.vflag_control.flag("CV::ACADEMIC")) {
    voutreach_load(voutreach,"ACADEMIC");
    voutreach_print_everything(cout,voutreach.at(0).vextra_latex,"Academic","ACADEMIC","\\section{Academic Positions}\\label{academic_positions}");
    Arun=TRUE;
  }
  if(!Arun && XHOST.vflag_control.flag("CV::RESEARCH")) {
    voutreach_load(voutreach,"RESEARCH");
    voutreach_print_everything(cout,voutreach.at(0).vextra_latex,"Research","RESEARCH","\\section{Research Experience}\\label{research_experience}");
    Arun=TRUE;
  }
  if(!Arun && XHOST.vflag_control.flag("CV::EDUCATION")) {
    voutreach_load(voutreach,"EDUCATION");
    voutreach_print_everything(cout,voutreach.at(0).vextra_latex,"Education","EDUCATION","\\section{Education}\\label{education}");
    Arun=TRUE;
  }
  if(!Arun && XHOST.vflag_control.flag("CV::TEACHING")) {
    voutreach_load(voutreach,"TEACHING");
    voutreach_print_everything(cout,voutreach.at(0).vextra_latex,"Teaching","TEACHING","\\section{Teaching Experience}\\label{teaching_experience}");
    Arun=TRUE;
  }
  if(!Arun && XHOST.vflag_control.flag("CV::ADVISING")) {
    voutreach_load(voutreach,"ADVISING");
    voutreach_print_everything(cout,voutreach.at(0).vextra_latex,"Advising","ADVISING","\\section{Advising Experience}\\label{advising_experience}");
    Arun=TRUE;
  }
  if(!Arun && XHOST.vflag_control.flag("CV::AWARDS")) {
    voutreach_load(voutreach,"AWARDS");
    voutreach_print_everything(cout,voutreach.at(0).vextra_latex,"Awards","AWARDS","\\section{Awards and Honors}\\label{awards}");
    Arun=TRUE;
  }
  if(!Arun && XHOST.vflag_control.flag("CV::PRESS")) {
    voutreach_load(voutreach,"PRESS");
    voutreach_print_everything(cout,voutreach.at(0).vextra_latex,"Press","PRESS","\\section{Press and news releases}\\label{pressreleases}");
    Arun=TRUE;
  }
  if(!Arun && XHOST.vflag_control.flag("CV::PATENTS")) {
    voutreach_load(voutreach,"PATENTS");
    voutreach_print_everything(cout,voutreach.at(0).vextra_latex,"Patents","PATENTS","\\section{Patents}\\label{patents}");
    Arun=TRUE;
  }
  if(!Arun && XHOST.vflag_control.flag("CV::SERVICE_OUTSIDE")) {
    voutreach_load(voutreach,"SERVICEOUTSIDE");
    voutreach_print_everything(cout,voutreach.at(0).vextra_latex,"ServiceOutside","SERVICE OUTSIDE","\\section{Outreach and Professional Activities}\\label{outreach}");
    Arun=TRUE;
  }
  if(!Arun && XHOST.vflag_control.flag("CV::SERVICE_INSIDE")) {
    voutreach_load(voutreach,"SERVICEINSIDE");
    voutreach_print_everything(cout,voutreach.at(0).vextra_latex,"ServiceInside","SERVICE INSIDE","\\section{Duke University - Academic Service Activities}\\label{duke_service}");
    Arun=TRUE;
  }
  if(!Arun && XHOST.vflag_control.flag("CV::AUTHOR")) {
    // something need to be specified
    // [OBSOLETE]   XHOST.vflag_control.flag("PRINT_MODE::HTML",TRUE);  // override
    // [OBSOLETE]   XHOST.vflag_control.flag("PRINT_MODE::JSON",FALSE);  // override
    // [OBSOLETE]    XHOST.vflag_control.flag("PRINT_MODE::HTML",FALSE); // override
    // [OBSOLETE]    XHOST.vflag_control.flag("PRINT_MODE::LATEX",FALSE); // override
    // [OBSOLETE]  print_mode="HTML"; // override
    voutreach_print(HTRESOURCE_MODE_PHP_AUTHOR,cout,"ARTICLES");
    Arun=TRUE;
  }

  return Arun;
}


// ******************************************************************************************************************************************************

uint voutreach_sort_year(vector<_outreach>& voutreach) {
  if(voutreach.size()>0)
    if(voutreach.at(0).year!=0) 
      sort(voutreach.begin(),voutreach.end(),_sort_outreach_outreach_year_());
  return voutreach.size();
}
uint voutreach_rsort_year(vector<_outreach>& voutreach) {
  if(voutreach.size()>0)
    if(voutreach.at(0).year!=0) 
      sort(voutreach.begin(),voutreach.end(),_rsort_outreach_outreach_year_());
  return voutreach.size();
}
uint voutreach_sort_wnumber(vector<_outreach>& voutreach) {
  if(voutreach.size()>0)
    if(voutreach.at(0).wnumber!=0) 
      sort(voutreach.begin(),voutreach.end(),_sort_outreach_outreach_wnumber_());
  return voutreach.size();
}
uint voutreach_rsort_wnumber(vector<_outreach>& voutreach) {
  if(voutreach.size()>0) 
    if(voutreach.at(0).wnumber!=0) 
      sort(voutreach.begin(),voutreach.end(),_rsort_outreach_outreach_wnumber_());
  return voutreach.size();
}

// ******************************************************************************************************************************************************
void center_print(uint mode, ostream& oss) {
  aurostd::xoption vflag=XHOST.vflag_outreach;

  if(mode) {;} // dummy load

  if(XHOST.vflag_control.flag("PHP::CENTER_MISSION")) {
    //  if(mode=print_mode && mode=HTRESOURCE_MODE_LATEX)
    oss << "  <br>" << endl;
  }
}

#endif //   _AFLOW_WEB_INTERFACE_CPP_

// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *                                                                         *
// ***************************************************************************
