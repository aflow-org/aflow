// **************************************************************************
// *                                                                        *
// *   VASPIN STEFANO CURTAROLO MASSACHUSETTS INSTITUTE OF TECHNOLOGY 2002  *
// *                                                                        *
// **************************************************************************


#include<fstream>
#include<sstream>
#include<string>
//#include<stringstream>

using std::string;
using std::ifstream;
using std::ofstream;
using std::ostringstream;
using std::ios_base;

#define BUFFER_MAXLEN 1024
#define AAmax 20
#define BBmax 20

void LatexSideMarginThesisON(void){
  cout << "\\oddsidemargin -0.25 in " << endl;
  cout << "\\evensidemargin -0.25 in " << endl;
}

void LatexSideMarginThesisOFF(void){
  cout << "\\oddsidemargin 0.2 in " << endl;
  cout << "\\evensidemargin 0.2 in " << endl;
}

void ZungerConversionChart(void);

void Rational(ostream &os,double x) {
  if(x<0) os << "-";
  x=abs(x);
  if(isequal(x,0.50000000000000)) {os << "1/2"; return;};
  if(isequal(x,0.08333333333333)) {os << "1/12"; return;};
  if(isequal(x,0.41666666666667)) {os << "5/12"; return;};
  if(isequal(x,0.58333333333333)) {os << "7/12"; return;};
  if(isequal(x,0.91666666666667)) {os << "11/12"; return;};
  if(isequal(x,1.50000000000000)) {os << "3/2"; return;};
  if(isequal(x,2.50000000000000)) {os << "5/2"; return;};
  if(isequal(x,0.33333333333333)) {os << "1/3"; return;};
  if(isequal(x,0.66666666666666)) {os << "2/3"; return;};
  if(isequal(x,1.33333333333333)) {os << "4/3"; return;};
  if(isequal(x,1.66666666666666)) {os << "5/3"; return;};
  if(isequal(x,2.33333333333333)) {os << "7/3"; return;};
  if(isequal(x,2.66666666666666)) {os << "8/3"; return;};
  if(isequal(x,0.12500000000000)) {os << "1/8"; return;};
  if(isequal(x,0.25000000000000)) {os << "1/4"; return;};
  if(isequal(x,0.37500000000000)) {os << "3/8"; return;};
  if(isequal(x,0.62500000000000)) {os << "5/8"; return;};
  if(isequal(x,0.75000000000000)) {os << "3/4"; return;};
  if(isequal(x,0.87500000000000)) {os << "7/8"; return;};
  if(isequal(x,0.11111111111111)) {os << "1/9"; return;};
  if(isequal(x,0.22222222222222)) {os << "2/9"; return;};
  if(isequal(x,0.44444444444444)) {os << "4/9"; return;};
  if(isequal(x,0.55555555555555)) {os << "5/9"; return;};
  if(isequal(x,0.77777777777777)) {os << "7/9"; return;};
  if(isequal(x,0.88888888888888)) {os << "8/9"; return;};
  if(isequal(x,1.73205080756888)) {os << "\\sqrt{3}"; return;};
  if(isequal(x,2.59807621135332)) {os << "3\\sqrt{3}/2"; return;};
  if(isequal(x,1.63299316185545)) {os << "\\sqrt{8/3}"; return;};
  if(isequal(x,3.26598632371090)) {os << "2\\sqrt{8/3}"; return;};
  if(isequal(x,0.86602540378444)) {os << "\\sqrt{3}/2"; return;};
  if(isequal(x,0.16666666666667)) {os << "1/6"; return;};
  if(isequal(x,0.83333333333333)) {os << "5/6"; return;};
  if(isequal(x,0.10000000000000)) {os << "1/10"; return;};
  if(isequal(x,0.20000000000000)) {os << "1/5"; return;};
  if(isequal(x,0.30000000000000)) {os << "3/10"; return;};
  if(isequal(x,0.40000000000000)) {os << "4/10"; return;};
  if(isequal(x,0.60000000000000)) {os << "3/5"; return;};
  if(isequal(x,0.70000000000000)) {os << "7/10"; return;};
  if(isequal(x,0.80000000000000)) {os << "4/5"; return;};
  if(isequal(x,0.90000000000000)) {os << "9/10"; return;};
  if(isequal(x,0.05555555555556)) {os << "1/18"; return;};
  if(isequal(x,0.27777777777778)) {os << "5/18"; return;};
  if(isequal(x,0.38888888888889)) {os << "7/18"; return;};
  if(isequal(x,0.61111111111111)) {os << "11/18"; return;};
  if(isequal(x,0.72222222222222)) {os << "13/18"; return;};
  if(isequal(x,0.94444444444444)) {os << "17/18"; return;};
  if(isequal(x,1.41421356237310)) {os << "\\sqrt{2}"; return;};
  if(isequal(x,2.23606797749979)) {os << "\\sqrt{5}"; return;};
  if(isequal(x,2.64575131106459)) {os << "\\sqrt{7}"; return;};
  if(isequal(x,3.16227766016838)) {os << "\\sqrt{10}"; return;};
  if(isequal(x,3.31662479035540)) {os << "\\sqrt{11}"; return;};
  if(isequal(x,3.46410161513775)) {os << "\\sqrt{12}"; return;};
  if(isequal(x,3.60555127546399)) {os << "\\sqrt{13}"; return;};
  if(isequal(x,3.74165738677394)) {os << "\\sqrt{14}"; return;};
  if(isequal(x,3.87298334620742)) {os << "\\sqrt{15}"; return;};
  os << x ; return;
}

void read_VASPIN(ostringstream* l,
    int number,
    const char* _labels,
    const char* _strukturbericht,
    const char* _superstr,
    const char* _superlatt,
    const char* _proto,
    bool VERB) {
  ifstream File,FileResult;
  ofstream FileOut;
  ostringstream aus0,aus1;
  char c=0,buffer[BUFFER_MAXLEN];//c=` `;
  string spacegroup,latticetype,pearson;
  string FileName("VASPIN/vasp.in.");
  string FileNameResult("/tmp/xresult");
  int AA,BB,pearsonN;
  double volume,volumeP,Cb,a_lat,b_lat,c_lat,alpha_lat,beta_lat,gamma_lat;
  xmatrix<double> A(3,3);
  xmatrix<double> Apos(AAmax,3),Bpos(BBmax,3);
  int j=1;
  l[j++] << "& " << _strukturbericht << "\t";;  // STRUKTURBERICHT
  l[j++] << "& " << _proto << "\t";;  // PROTOTYPE
  l[j++] << "& " << _labels << "\t";;  // LABEL
  { // OPEN VASPIN FILE
    for(int i=0;i<BUFFER_MAXLEN;i++)aus0<<c;aus0.seekp(0,ios_base::beg);
    aus0 << number;
    if(number==17) aus0 << ".unrelaxed";
    if(number==120) aus0 << ".unrelaxed";
    FileName += aus0.str();
    if(VERB) cerr << "FILE=[" << FileName << "]" << "  ";
    File.open(FileName.c_str(),std::ios::in);
    if(!File) xerror("File not found",FileName.c_str());
  }
  { // GET AA and BB NUMBERS
    for(int i=0;i<BUFFER_MAXLEN;i++)aus1<<c;aus1.seekp(0,ios_base::beg);  // RESET
    aus1 << "grep -c BB " << FileName.c_str() << " > " << FileNameResult.c_str() << endl;
    system(aus1.str().c_str());
    FileResult.open(FileNameResult.c_str(),std::ios::in);
    FileResult.getline(buffer,BUFFER_MAXLEN,'\n');
    FileResult.close();
    BB=atoi(buffer);
    aus1.seekp(0,ios_base::beg);  // RESET
    aus1 << "grep -c AA " << FileName.c_str() << " > " << FileNameResult.c_str() << endl;
    system(aus1.str().c_str());
    FileResult.open(FileNameResult.c_str(),std::ios::in);
    FileResult.getline(buffer,BUFFER_MAXLEN,'\n');
    FileResult.close();
    AA=atoi(buffer);
    if(VERB) cerr << "AA=" << AA << "\t";
    if(VERB) cerr << "BB=" << BB << "\t";
    volume = 18.1992*double(AA)+22.9326*double(BB);
    if(VERB) cerr << "volume=" << volume << "\t";
  }
  if(AA>AAmax) xerror("AA>Amax: ",AA,">",AAmax);
  if(BB>BBmax) xerror("BB>Bmax: ",BB,">",BBmax);
  if(BB>0 && BB<AA) xerror("BB<AA");
  Cb=(double) BB/(AA+BB);
  { // FORMULA
    l[j] << "& ";
    if(AA==1) l[j] << "A";
    if(AA>1)  l[j] << "A$_{"<< AA << "}$";
    if(BB==1) l[j] << "B";
    if(BB>1)  l[j] << "B$_{"<< BB << "}$";
    l[j]<<"\t"; j++;
  }
  {
    l[j++] << "& " << _superstr << "\t";;  // superstructure
  }
  { // C_b
    //  l[j] << "& ";Rational(l[j],Cb);l[j]<<"\t";j++;
    l[j].precision(4);
    l[j] << "& " << Cb <<"\t";j++;
  }
  { // CREATE VASPIN WITH PARAMETERS
    aus0.seekp(0,ios_base::beg);for(int i=0;i<BUFFER_MAXLEN;i++)aus0<<c;aus0.seekp(0,ios_base::beg);  // RESET
    aus0 << "cat " << FileName.c_str() << " | sed \"s/AA/Au/g\" | sed \"s/BB/Cd/g\" | sed \"s/SCALExxx/-" << volume << "/g\"> /tmp/vasp.in" << endl;
  }
  { // GET SPACEGROUP
    aus0.seekp(0,ios_base::beg);for(int i=0;i<BUFFER_MAXLEN;i++) aus0<<c;aus0.seekp(0,ios_base::beg);  // RESET
    aus0 << "cd /tmp && /home/auro/work/gndstate/GND_SCRIPTS_AB/BIN/xplato_spacegroup | sed \"s/#/\\\\\\\\#/g\" "<<" > "<<FileNameResult.c_str()<<endl;
    system(aus0.str().c_str());
    // if(VERB) cerr << aus0.str() << "\t";
    FileResult.open(FileNameResult.c_str(),std::ios::in);
    FileResult.getline(buffer,BUFFER_MAXLEN,'\n');
    FileResult.close();
    spacegroup += buffer;
    if(VERB) cerr << "SG=" << spacegroup << "  ";
    l[j] << "&{\\small " << spacegroup << "}\t";j++;
    //   l[j] << "&{\\footnotesize " << spacegroup << "}\t";j++;
  }
  { // GET PEARSON_NUMBER
    aus0.seekp(0,ios_base::beg);for(int i=0;i<BUFFER_MAXLEN;i++) aus0<<c;aus0.seekp(0,ios_base::beg);  // RESET
    aus0 << "cd /tmp && /home/auro/work/gndstate/GND_SCRIPTS_AB/BIN/xplato_pearsonN" << " > " << FileNameResult.c_str() << endl;
    system(aus0.str().c_str());
    FileResult.open(FileNameResult.c_str(),std::ios::in);
    FileResult.getline(buffer,BUFFER_MAXLEN,'\n');
    FileResult.close();
    volumeP=atof(buffer+58);
    if(VERB) cerr << "VOLUMEP=" << volumeP << " ";
    pearsonN=(int) (((double) AA+BB)*((double)volumeP/volume+0.1));
  }
  { // GET PEARSON_SYMBOL    
    aus0.seekp(0,ios_base::beg);for(int i=0;i<BUFFER_MAXLEN;i++) aus0<<c;aus0.seekp(0,ios_base::beg);  // RESET
    aus0 << "cd /tmp && /home/auro/work/gndstate/GND_SCRIPTS_AB/BIN/xplato_pearson" << " > " << FileNameResult.c_str() << endl;
    system(aus0.str().c_str());
    // if(VERB) cerr << aus0.str() << "\t";
    FileResult.open(FileNameResult.c_str(),std::ios::in);
    FileResult.getline(buffer,BUFFER_MAXLEN,'\n');
    FileResult.close();
    pearson += buffer;
    if(VERB) cerr << "PEARSON=" << pearson << pearsonN << "\t";
    l[j] << "& " << pearson << pearsonN  << "\t";j++;
  }
  { // GET LATTICETYPE
    aus0.seekp(0,ios_base::beg);for(int i=0;i<BUFFER_MAXLEN;i++) aus0<<c;aus0.seekp(0,ios_base::beg);  // RESET
    aus0 << "cd /tmp && /home/auro/work/gndstate/GND_SCRIPTS_AB/BIN/xplato_latticetype" << " > " << FileNameResult.c_str() << endl;
    system(aus0.str().c_str());
    // if(VERB) cerr << aus0.str() << "\t";
    FileResult.open(FileNameResult.c_str(),std::ios::in);
    FileResult.getline(buffer,BUFFER_MAXLEN,'\n');
    FileResult.close();
    latticetype += buffer;
    if(VERB) cerr << "LATTICE=" << latticetype << " ";
    //   l[j] << "& " << latticetype << "\t";j++;
    l[j] << "&{\\small  " << latticetype << "}\t";j++;
    aus0.seekp(0,ios_base::beg);for(int i=0;i<BUFFER_MAXLEN;i++) aus0<<c;aus0.seekp(0,ios_base::beg);  // RESET
    aus0 << "cd /tmp && /home/auro/work/gndstate/GND_SCRIPTS_AB/BIN/xplato_latticetype_abc" << " > " << FileNameResult.c_str() << endl;
    system(aus0.str().c_str());
    // if(VERB) cerr << aus0.str() << "\t";
    FileResult.open(FileNameResult.c_str(),std::ios::in);
    FileResult.getline(buffer,BUFFER_MAXLEN,'\n');a_lat=atof(buffer);//if(VERB) cerr << "a_lat=" << a_lat << endl;
    FileResult.getline(buffer,BUFFER_MAXLEN,'\n');b_lat=atof(buffer);//if(VERB) cerr << "b_lat=" << b_lat << endl;
    FileResult.getline(buffer,BUFFER_MAXLEN,'\n');c_lat=atof(buffer);//if(VERB) cerr << "c_lat=" << c_lat << endl;
    FileResult.getline(buffer,BUFFER_MAXLEN,'\n');alpha_lat=atof(buffer);//if(VERB) cerr << "alpha_lat=" << alpha_lat << endl;
    FileResult.getline(buffer,BUFFER_MAXLEN,'\n');beta_lat=atof(buffer); //if(VERB) cerr << "beta_lat=" << beta_lat << endl;
    FileResult.getline(buffer,BUFFER_MAXLEN,'\n');gamma_lat=atof(buffer);//if(VERB) cerr << "gamma_lat=" << gamma_lat << endl;
    FileResult.close();
    l[j].precision(4);l[j] << "&$(" << a_lat << "," << b_lat << "," << c_lat << ")$\t";j++;
    l[j].precision(4);l[j] << "&$(" << alpha_lat << "," << beta_lat << "," << gamma_lat << ")$\t";j++;
  }
  { // EQUIVALENT SUPERLATTIVE
    l[j++] << "& " << _superlatt << "\t";;  // LABEL
  }
  { // GET A1 A2 A3
    aus0.seekp(0,ios_base::beg);for(int i=0;i<BUFFER_MAXLEN;i++) aus0<<c;aus0.seekp(0,ios_base::beg);  // RESET
    aus0 << "cd /tmp && /home/auro/work/gndstate/GND_SCRIPTS_AB/BIN/xplato_primitive " << " > " << FileNameResult.c_str() << endl;
    system(aus0.str().c_str());
    // if(VERB) cerr << aus0.str() << "\t";
    FileResult.open(FileNameResult.c_str(),std::ios::in);
    FileResult.getline(buffer,BUFFER_MAXLEN,'\n');A(1,1)=atof(buffer);//if(VERB) cerr << "A(1,1)=" << A(1,1) << endl;
    FileResult.getline(buffer,BUFFER_MAXLEN,'\n');A(1,2)=atof(buffer);//if(VERB) cerr << "A(1,2)=" << A(1,2) << endl;
    FileResult.getline(buffer,BUFFER_MAXLEN,'\n');A(1,3)=atof(buffer);//if(VERB) cerr << "A(1,3)=" << A(1,3) << endl;
    FileResult.getline(buffer,BUFFER_MAXLEN,'\n');A(2,1)=atof(buffer);//if(VERB) cerr << "A(2,1)=" << A(2,1) << endl;
    FileResult.getline(buffer,BUFFER_MAXLEN,'\n');A(2,2)=atof(buffer);//if(VERB) cerr << "A(2,2)=" << A(2,2) << endl;
    FileResult.getline(buffer,BUFFER_MAXLEN,'\n');A(2,3)=atof(buffer);//if(VERB) cerr << "A(2,3)=" << A(2,3) << endl;
    FileResult.getline(buffer,BUFFER_MAXLEN,'\n');A(3,1)=atof(buffer);//if(VERB) cerr << "A(3,1)=" << A(3,1) << endl;
    FileResult.getline(buffer,BUFFER_MAXLEN,'\n');A(3,2)=atof(buffer);//if(VERB) cerr << "A(3,2)=" << A(3,2) << endl;
    FileResult.getline(buffer,BUFFER_MAXLEN,'\n');A(3,3)=atof(buffer);//if(VERB) cerr << "A(3,3)=" << A(3,3) << endl;
    FileResult.close();
    l[j] << "&\t\t";j++;  // primitive vectors
    l[j]<<"&$(";Rational(l[j],A(1,1));l[j]<<",";Rational(l[j],A(1,2));l[j]<<",";Rational(l[j],A(1,3));l[j]<<")$";j++;
    l[j]<<"&$(";Rational(l[j],A(2,1));l[j]<<",";Rational(l[j],A(2,2));l[j]<<",";Rational(l[j],A(2,3));l[j]<<")$";j++;
    l[j]<<"&$(";Rational(l[j],A(3,1));l[j]<<",";Rational(l[j],A(3,2));l[j]<<",";Rational(l[j],A(3,3));l[j]<<")$";j++;
    l[j] << "&\t\t";j++;  // cartesian
  }  
  { // GET A1 A2 A3
    aus0.seekp(0,ios_base::beg);for(int i=0;i<BUFFER_MAXLEN;i++) aus0<<c;aus0.seekp(0,ios_base::beg);  // RESET
    aus0 << "cd /tmp && /home/auro/work/gndstate/GND_SCRIPTS_AB/BIN/xplato_positions " << " > " << FileNameResult.c_str() << endl;
    system(aus0.str().c_str());
    // if(VERB) cerr << aus0.str() << "\t";
    FileResult.open(FileNameResult.c_str(),std::ios::in);
    l[j] << "&\t\t";j++;  // atomic positions
    for(int i=1;i<=AAmax;i++) {
      if(i<=AA) {
        FileResult.getline(buffer,BUFFER_MAXLEN,'\n');Apos(i,1)=atof(buffer);if(VERB) cerr << "Apos("<<i<<",1)=" << Apos(i,1) << endl;
        FileResult.getline(buffer,BUFFER_MAXLEN,'\n');Apos(i,2)=atof(buffer);if(VERB) cerr << "Apos("<<i<<",2)=" << Apos(i,2) << endl;
        FileResult.getline(buffer,BUFFER_MAXLEN,'\n');Apos(i,3)=atof(buffer);if(VERB) cerr << "Apos("<<i<<",3)=" << Apos(i,3) << endl;
        //l[j]<<"& [A]$(";Rational(l[j],Apos(i,1));l[j]<<",";Rational(l[j],Apos(i,2));l[j]<<",";Rational(l[j],Apos(i,3));l[j]<<")$";j++;
        l[j]<<"& $(";Rational(l[j],Apos(i,1));l[j]<<",";Rational(l[j],Apos(i,2));l[j]<<",";Rational(l[j],Apos(i,3));l[j]<<")$";j++;
      } else {
        //l[j] << "& [A]    $-$     ";j++;
        l[j] << "&     $-$     ";j++;
      }
    }
    for(int i=1;i<=BBmax;i++) {
      if(i<=BB) {
        FileResult.getline(buffer,BUFFER_MAXLEN,'\n');Bpos(i,1)=atof(buffer);if(VERB) cerr << "Bpos("<<i<<",1)=" << Bpos(i,1) << endl;
        FileResult.getline(buffer,BUFFER_MAXLEN,'\n');Bpos(i,2)=atof(buffer);if(VERB) cerr << "Bpos("<<i<<",2)=" << Bpos(i,2) << endl;
        FileResult.getline(buffer,BUFFER_MAXLEN,'\n');Bpos(i,3)=atof(buffer);if(VERB) cerr << "Bpos("<<i<<",3)=" << Bpos(i,3) << endl;
        //l[j]<<"& [B]$(";Rational(l[j],Bpos(i,1));l[j]<<",";Rational(l[j],Bpos(i,2));l[j]<<",";Rational(l[j],Bpos(i,3));l[j]<<")$";j++;
        l[j]<<"& $(";Rational(l[j],Bpos(i,1));l[j]<<",";Rational(l[j],Bpos(i,2));l[j]<<",";Rational(l[j],Bpos(i,3));l[j]<<")$";j++;
      } else {
        //l[j] << "& [B]    $-$     ";j++;
        l[j] << "&     $-$     ";j++;
      }
    }
    l[j] << "&\t\t";j++;  // fractional
  }

  //system(aus0.str().c_str());
  //if(VERB) cerr << aus0.str() << "\t";
  //  system(aus0.str().c_str());
  if(VERB) cerr << endl;
}


#define MAXLINESOSTREAM 64

void MakeVaspIn(bool VERB) {
  cerr << "% VASP IN PROCEDURE " << endl;
  int i,j;
  int i1=17,i2=8,i3=11,i4=12,i5=12,i6=13,i7=17; //some parameters
  bool _FCC=0,_BCC=0,_HCP=0,_EXT=0;

  _FCC=1,_BCC=1,_HCP=1,_EXT=1;

  cerr << "% VASP IN PROCEDURE " << endl;
  //cout << "\\documentclass" << APENNSY_LATEX_DOC_TYPE << endl;
  //cout << "\\usepackage{epsfig}" << endl;
  //cout << "\\begin{document}" << endl;
  //cout << "\\voffset 15mm" << endl;
  //cout << "\\hoffset -5mm" << endl;
  //cout << "\\textheight 240mm" << endl;
  //cout << "\\textwidth 180mm" << endl;
  //ZungerConversionChart();
  cout << "% VASP IN PROCEDURE " << endl;
  cout << "\\def\\baselinestretch{1}" << endl;
  cout << "\\section{{Tables of the structure prototypes}} " << endl;
  cout << "\\label{appendix.LIST_superstructures} " << endl;
  if(VERB) cerr << " *********************************************************************** FCC " << endl;
  if(VERB) cerr << " *********************************************************************** FCC " << endl;
  if(VERB) cerr << " *********************************************************************** FCC " << endl;
  if(VERB) cerr << " *********************************************************************** FCC " << endl;
  if(_FCC) {
    { // FCC1a TABLE
      cerr << "FCC1a TABLE" << endl;// FCC1 TABLE
      cout << "\\subsection{{FCC superstructures}} " << endl;
      cout << "\\label{appendix.FCC_superstructures} " << endl;
      LatexSideMarginThesisON();
      ostringstream lines[MAXLINESOSTREAM];
      j=1;
      lines[j++] << "Strukturbericht";lines[j++] << "Prototype      ";
      lines[j++] << "Label          ";lines[j++] << "Formula        ";
      lines[j++] << "Parent lattice ";lines[j++] << "C$_b$          ";
      lines[j++] << "Space Group    ";lines[j++] << "Pearson Symbol ";
      lines[j++] << "Lattice type   ";lines[j++] << "$(a,b,c)$ (\\AA)";lines[j++] << "$(\\alpha,\\beta,\\gamma)$ degrees";
      lines[j++] << "Superlattice   ";
      lines[j++]<<"Primitive vectors";lines[j++]<<"{\\bf a$_1$}/a";lines[j++]<<"{\\bf a$_2$}/a";lines[j++]<<"{\\bf a$_3$}/a";lines[j++]<<"(cartesian)";
      lines[j++]<<"Atomic positions";
      for(i=1;i<=AAmax;i++) lines[j++] <<"{\\bf A$"<<i<<"$}     ";
      for(i=1;i<=BBmax;i++) lines[j++] <<"{\\bf B$"<<i<<"$}     ";
      lines[j++] << "(fractional)   ";
      read_VASPIN(lines,2,"1,2 ","A1"                 ,"FCC","                ","Cu      ",VERB);
      read_VASPIN(lines,3,"3   ","L1$_0$ "            ,"FCC","AB along [001]  ","AuCu    ",VERB);
      read_VASPIN(lines,4,"4   ","L1$_1$"             ,"FCC","AB along [111]  ","CuPt     ",VERB);
      read_VASPIN(lines,5,"5,6 ","$\\beta1/\\beta2$-FCC$_{AB2}^{[100]}$"  ,"FCC","AB2 along [100] ","        ",VERB);
      for(i=1;i<j;i++) lines[i] << "\t \\\\ ";
      cout << "{\\begin{center}" << endl;
      cout << "{\\footnotesize" << endl;
      cout << "\\begin{tabular}{||c||c|c|c|c|c||}\\hline\\hline   " << endl;
      for(i=1;i<j;i++) {
        string aus=lines[i].str();
        if(i<=i1 || (i>i1 && (((int) aus.find(")"))>0)) || ((((int) aus.find("posi"))>0)) ) {
          if(i!=i5 || (int) aus.find("A")>0) {
            cout << lines[i].str();
            if(i<=i2 || i==i1 || i==i3 || i==i4 ) cout << "\\hline";
            cout << endl;}}}
      //     cout << "\\hline \\hline" << endl << "\\end{tabular}}" << endl << "\\end{center}}" << endl;
      cout << "\\hline \\hline \\hline \\hline " << endl;
    }
    { // FCC1b TABLE
      cerr << "FCC1b TABLE" << endl; // FCC1b TABLE
      ostringstream lines[MAXLINESOSTREAM];
      j=1;
      lines[j++] << "Strukturbericht";lines[j++] << "Prototype      ";
      lines[j++] << "Label          ";lines[j++] << "Formula        ";
      lines[j++] << "Parent lattice ";lines[j++] << "C$_b$          ";
      lines[j++] << "Space Group    ";lines[j++] << "Pearson Symbol ";
      lines[j++] << "Lattice type   ";lines[j++] << "$(a,b,c)$ (\\AA)";lines[j++] << "$(\\alpha,\\beta,\\gamma)$ degrees";
      lines[j++] << "Superlattice   ";
      lines[j++]<<"Primitive vectors";lines[j++]<<"{\\bf a$_1$}/a";lines[j++]<<"{\\bf a$_2$}/a";lines[j++]<<"{\\bf a$_3$}/a";lines[j++]<<"(cartesian)";
      lines[j++]<<"Atomic positions";
      for(i=1;i<=AAmax;i++) lines[j++] <<"{\\bf A$"<<i<<"$}     ";
      for(i=1;i<=BBmax;i++) lines[j++] <<"{\\bf B$"<<i<<"$}     ";
      lines[j++] << "(fractional)   ";
      read_VASPIN(lines,7,"7,8 "," - "                ,"FCC","AB2 along [011] ","MoPt$_2$",VERB);
      read_VASPIN(lines,9,"9,10","$\\alpha1/\\alpha2$-FCC$_{AB2}^{[111]}$","FCC","AB2 along [111] ","        ",VERB);
      read_VASPIN(lines,11,"11,12"," -   ","FCC","          ","       ",VERB);
      read_VASPIN(lines,13,"13,15","Z1/Z3-FCC$_{AB3}^{[001]}$","FCC","AB3 along [001] "," ",VERB);
      for(i=1;i<j;i++) lines[i] << "\t \\\\ ";
      // cout << "{\\begin{center}" << endl;
      // cout << "{\\footnotesize" << endl;
      // cout << "\\begin{tabular}{||c||c|c|c|c|c|c||}\\hline\\hline   " << endl;
      for(i=1;i<j;i++) {
        string aus=lines[i].str();
        if(i<=i1 || (i>i1 && (((int) aus.find(")"))>0)) || ((((int) aus.find("posi"))>0)) ) {
          if(i!=i5 || (int) aus.find("A")>0) {
            cout << lines[i].str();
            if(i<=i2 || i==i1 || i==i3 || i==i4 ) cout << "\\hline";
            cout << endl;}}}
      cout << "\\hline\\hline" << endl << "\\end{tabular}}" << endl << "\\end{center}}" << endl;
    }
  }
  cout << "\\newpage" << endl;
  if(_FCC) {
    { // FCC2a TABLE
      cerr << "FCC2a TABLE" << endl; // FCC2 TABLE
      //    cout << "\\subsection{FCC2 superstructures} " << endl;
      ostringstream lines[MAXLINESOSTREAM];
      j=1;
      lines[j++] << "Strukturbericht";lines[j++] << "Prototype      ";
      lines[j++] << "Label          ";lines[j++] << "Formula        ";
      lines[j++] << "Parent lattice ";lines[j++] << "C$_b$          ";
      lines[j++] << "Space Group    ";lines[j++] << "Pearson Symbol ";
      lines[j++] << "Lattice type   ";lines[j++] << "$(a,b,c)$ (\\AA)";lines[j++] << "$(\\alpha,\\beta,\\gamma)$ degrees";
      lines[j++] << "Superlattice   ";
      lines[j++]<<"Primitive vectors";lines[j++]<<"{\\bf a$_1$}/a";lines[j++]<<"{\\bf a$_2$}/a";lines[j++]<<"{\\bf a$_3$}/a";lines[j++]<<"(cartesian)";
      lines[j++]<<"Atomic positions";
      for(i=1;i<=AAmax;i++) lines[j++] <<"{\\bf A$"<<i<<"$}     ";
      for(i=1;i<=BBmax;i++) lines[j++] <<"{\\bf B$"<<i<<"$}     ";
      lines[j++] << "(fractional)   ";
      read_VASPIN(lines,14,"14"   ,"Z2-FCC$_{A2B2}^{[001]}$","FCC","A2B2 along [001]"," ",VERB);
      read_VASPIN(lines,16,"16,18","     ","FCC","        ","         ",VERB);
      read_VASPIN(lines,17,"17"   ,"W2-FCC$_{A2B2}^{[311]}$","FCC","A2B2 along [311]"," ",VERB);
      read_VASPIN(lines,19,"19,21","Y1/Y3-FCC$_{AB3}^{[011]}$","FCC","AB3  along [011]"," ",VERB);// 19,21   oP?  Pmmm   Y1/3  AB3   #47  AB3  along [011]z
      for(i=1;i<j;i++) lines[i] << "\t \\\\ ";
      cout << "{\\begin{center}" << endl;
      cout << "{\\footnotesize" << endl;
      cout << "\\begin{tabular}{||c||c|c|c|c|c||}\\hline\\hline   " << endl;
      for(i=1;i<j;i++) {
        string aus=lines[i].str();
        if(i<=i1 || (i>i1 && (((int) aus.find(")"))>0)) || ((((int) aus.find("posi"))>0)) ) {
          if(i!=i5 || (int) aus.find("A")>0) {
            cout << lines[i].str();
            if(i<=i2 || i==i1 || i==i3 || i==i4 ) cout << "\\hline";
            cout << endl;}}}
      //     cout << "\\hline \\hline" << endl << "\\end{tabular}}" << endl << "\\end{center}}" << endl;
      cout << "\\hline \\hline \\hline \\hline " << endl;
    }
    { // FCC2b TABLE
      cerr << "FCC2b TABLE" << endl; // FCC4 TABLE
      ostringstream lines[MAXLINESOSTREAM];
      j=1;
      lines[j++] << "Strukturbericht";lines[j++] << "Prototype      ";
      lines[j++] << "Label          ";lines[j++] << "Formula        ";
      lines[j++] << "Parent lattice ";lines[j++] << "C$_b$          ";
      lines[j++] << "Space Group    ";lines[j++] << "Pearson Symbol ";
      lines[j++] << "Lattice type   ";lines[j++] << "$(a,b,c)$ (\\AA)";lines[j++] << "$(\\alpha,\\beta,\\gamma)$ degrees";
      lines[j++] << "Superlattice   ";
      lines[j++]<<"Primitive vectors";lines[j++]<<"{\\bf a$_1$}/a";lines[j++]<<"{\\bf a$_2$}/a";lines[j++]<<"{\\bf a$_3$}/a";lines[j++]<<"(cartesian)";
      lines[j++]<<"Atomic positions";
      for(i=1;i<=AAmax;i++) lines[j++] <<"{\\bf A$"<<i<<"$}     ";
      for(i=1;i<=BBmax;i++) lines[j++] <<"{\\bf B$"<<i<<"$}     ";
      lines[j++] << "(fractional)   ";
      read_VASPIN(lines,20,"20   ","Y2-FCC$_{A2B2}^{[011]}$","FCC","A2B2 along [011]","        ",VERB);// 20    oP? Pmmm     Y2      A2B2   #47  A2B2 along [011]
      read_VASPIN(lines,22,"22,24","D0$_{22}$   ","FCC","AB4  along [201]","Al$_3$Ti",VERB);// 22,24 tI8 I4/mmm   D0_22   Al3Ti  #139 AB4  along [201]
      read_VASPIN(lines,23,"23   ","CH-FCC$_{A2B2}^{[201]}$","FCC","A2B2 along [201]","NbP     ",VERB);// 23    tI8 I41/amd  CH "40" NbP    #141 A2B2 along [201]
      read_VASPIN(lines,25,"25,26","L1$_2$ ","FCC","                ","AuCu$_3$",VERB);// 25,26 cP4 Pm3(bar)m  L1_2    AuCu3  #221  
      for(i=1;i<j;i++) lines[i] << "\t \\\\ ";
      // cout << "{\\begin{center}" << endl;
      // cout << "{\\footnotesize" << endl;
      // cout << "\\begin{tabular}{||c||c|c|c|c|c|c||}\\hline\\hline   " << endl;
      for(i=1;i<j;i++) {
        string aus=lines[i].str();
        if(i<=i1 || (i>i1 && (((int) aus.find(")"))>0)) || ((((int) aus.find("posi"))>0)) ) {
          if(i!=i5 || (int) aus.find("A")>0) {
            cout << lines[i].str();
            if(i<=i2 || i==i1 || i==i3 || i==i4 ) cout << "\\hline";
            cout << endl;}}}
      cout << "\\hline\\hline" << endl << "\\end{tabular}}" << endl << "\\end{center}}" << endl;
    }
  }
  cout << "\\newpage" << endl;
  if(_FCC) {
    { // FCC2c TABLE
      cerr << "FCC2c TABLE" << endl; // FCC2 TABLE
      //    cout << "\\subsection{FCC2 superstructures} " << endl;
      ostringstream lines[MAXLINESOSTREAM];
      j=1;
      lines[j++] << "Strukturbericht";lines[j++] << "Prototype      ";
      lines[j++] << "Label          ";lines[j++] << "Formula        ";
      lines[j++] << "Parent lattice ";lines[j++] << "C$_b$          ";
      lines[j++] << "Space Group    ";lines[j++] << "Pearson Symbol ";
      lines[j++] << "Lattice type   ";lines[j++] << "$(a,b,c)$ (\\AA)";lines[j++] << "$(\\alpha,\\beta,\\gamma)$ degrees";
      lines[j++] << "Superlattice   ";
      lines[j++]<<"Primitive vectors";lines[j++]<<"{\\bf a$_1$}/a";lines[j++]<<"{\\bf a$_2$}/a";lines[j++]<<"{\\bf a$_3$}/a";lines[j++]<<"(cartesian)";
      lines[j++]<<"Atomic positions";
      for(i=1;i<=AAmax;i++) lines[j++] <<"{\\bf A$"<<i<<"$}     ";
      for(i=1;i<=BBmax;i++) lines[j++] <<"{\\bf B$"<<i<<"$}     ";
      lines[j++] << "(fractional)   ";
      read_VASPIN(lines,27,"27,29","V1/V3-FCC$_{AB3}^{[111]}$","FCC","AB3  along [111]"," ",VERB);// 27,29 hR? R-3m  V1/3 AB3    #166 AB3  along [111]
      read_VASPIN(lines,28,"28   ","V2-FCC$_{A2B2}^{[111]}$","FCC","A2B2 along [111]","  ",VERB);// 28  hR? R-3m    V2    A2B2   #166 A2B2 along [111]
      read_VASPIN(lines,30,"30/35","       ","FCC","                ","        ",VERB);//                                                          
      read_VASPIN(lines,36,"36/39","       ","FCC","                ","        ",VERB);//
      for(i=1;i<j;i++) lines[i] << "\t \\\\ ";
      cout << "{\\begin{center}" << endl;
      cout << "{\\footnotesize" << endl;
      cout << "\\begin{tabular}{||c||c|c|c|c|c||}\\hline\\hline   " << endl;
      for(i=1;i<j;i++) {
        string aus=lines[i].str();
        if(i<=i1 || (i>i1 && (((int) aus.find(")"))>0)) || ((((int) aus.find("posi"))>0)) ) {
          if(i!=i5 || (int) aus.find("A")>0) {
            cout << lines[i].str();
            if(i<=i2 || i==i1 || i==i3 || i==i4 ) cout << "\\hline";
            cout << endl;}}}
      cout << "\\hline\\hline" << endl << "\\end{tabular}}" << endl << "\\end{center}}" << endl;
    }
  }
  LatexSideMarginThesisOFF();
  cout << "\\newpage" << endl;
  if(VERB) cerr << " *********************************************************************** BCC " << endl;
  if(VERB) cerr << " *********************************************************************** BCC " << endl;
  if(VERB) cerr << " *********************************************************************** BCC " << endl;
  if(VERB) cerr << " *********************************************************************** BCC " << endl;
  if(_BCC) {
    {// BCC1a TABLE
      cerr << "BCC1a TABLE" << endl; // BCC1 TABLE
      cout << "\\subsection{{BCC superstructures}} " << endl;
      cout << "\\label{appendix.BCC_superstructures} " << endl;
      LatexSideMarginThesisON();
      ostringstream lines[MAXLINESOSTREAM];
      j=1;
      lines[j++] << "Strukturbericht";lines[j++] << "Prototype      ";
      lines[j++] << "Label          ";lines[j++] << "Formula        ";
      lines[j++] << "Parent lattice ";lines[j++] << "C$_b$          ";
      lines[j++] << "Space Group    ";lines[j++] << "Pearson Symbol ";
      lines[j++] << "Lattice type   ";lines[j++] << "$(a,b,c)$ (\\AA)";lines[j++] << "$(\\alpha,\\beta,\\gamma)$ degrees";
      lines[j++] << "Superlattice   ";
      lines[j++]<<"Primitive vectors";lines[j++]<<"{\\bf a$_1$}/a";lines[j++]<<"{\\bf a$_2$}/a";lines[j++]<<"{\\bf a$_3$}/a";lines[j++]<<"(cartesian)";
      lines[j++]<<"Atomic positions";
      for(i=1;i<=AAmax;i++) lines[j++] <<"{\\bf A$"<<i<<"$}     ";
      for(i=1;i<=BBmax;i++) lines[j++] <<"{\\bf B$"<<i<<"$}     ";
      lines[j++] << "(fractional)   ";
      read_VASPIN(lines,58,"58,59","A2       ","BCC","                ","W      ",VERB);   // 58,59  bcc/cI2 Im-3m  A2    W       #229  
      read_VASPIN(lines,60,"60","","BCC","AB along [101]","$\\gamma$-IrV",VERB);// 60   oC8   Cmmm   gamma-IrV #65 AB along [101]
      read_VASPIN(lines,61,"61","B2","BCC"," AB along [001] ","CsCl    ",VERB);            // 61     cP2   Pm-3m    B2    CsCl    #221  AB along [001]
      read_VASPIN(lines,62,"62,63","   ","BCC"," AB2 along [211]","   ",VERB);             // 62,63  hP?   P-3m1    ?     AB2     #164  AB2 along [211]
      for(i=1;i<j;i++) lines[i] << "\t \\\\ ";
      cout << "{\\begin{center}" << endl;
      cout << "{\\footnotesize" << endl;
      cout << "\\begin{tabular}{||c||c|c|c|c|c||}\\hline\\hline   " << endl;
      for(i=1;i<j;i++) {
        string aus=lines[i].str();
        if(i<=i1 || (i>i1 && (((int) aus.find(")"))>0)) || ((((int) aus.find("posi"))>0)) ) {
          if(i!=i5 || (int) aus.find("A")>0) {
            cout << lines[i].str();
            if(i<=i2 || i==i1 || i==i3 || i==i4 ) cout << "\\hline";
            cout << endl;}}}
      //     cout << "\\hline \\hline" << endl << "\\end{tabular}}" << endl << "\\end{center}}" << endl;
      cout << "\\hline \\hline \\hline \\hline " << endl;
    }
    { // BCC1b TABLE
      cerr << "BCC1b TABLE" << endl; // BCC2 TABLE
      ostringstream lines[MAXLINESOSTREAM];
      j=1;
      lines[j++] << "Strukturbericht";lines[j++] << "Prototype      ";
      lines[j++] << "Label          ";lines[j++] << "Formula        ";
      lines[j++] << "Parent lattice ";lines[j++] << "C$_b$          ";
      lines[j++] << "Space Group    ";lines[j++] << "Pearson Symbol ";
      lines[j++] << "Lattice type   ";lines[j++] << "$(a,b,c)$ (\\AA)";lines[j++] << "$(\\alpha,\\beta,\\gamma)$ degrees";
      lines[j++] << "Superlattice   ";
      lines[j++]<<"Primitive vectors";lines[j++]<<"{\\bf a$_1$}/a";lines[j++]<<"{\\bf a$_2$}/a";lines[j++]<<"{\\bf a$_3$}/a";lines[j++]<<"(cartesian)";
      lines[j++]<<"Atomic positions";
      for(i=1;i<=AAmax;i++) lines[j++] <<"{\\bf A$"<<i<<"$}     ";
      for(i=1;i<=BBmax;i++) lines[j++] <<"{\\bf B$"<<i<<"$}     ";
      lines[j++] << "(fractional)   ";
      read_VASPIN(lines,64,"64,65","   ","BCC"," AB2 along [011]","  ",VERB);              // 64,65  oF?   Fmmm     ?     AB2     #69   AB2 along [011]
      read_VASPIN(lines,66,"66,67","C11$_b$ ","BCC"," ","MoSi$_2$",VERB);// 66,67  tI6   I4/mmm   C11_b MoSi2   #139  
      read_VASPIN(lines,68,"68,69","A2      ","BCC"," "," ",VERB);// 68,69   hR?       R-3m     ?   AB3       #166   UNDEF
      read_VASPIN(lines,70,"70,72","        ","BCC"," "," ",VERB);// 70,72   oC?       Cmmm     ?   AB3       #65    UNDEF  
      for(i=1;i<j;i++) lines[i] << "\t \\\\ ";
      // cout << "{\\begin{center}" << endl;
      // cout << "{\\footnotesize" << endl;
      // cout << "\\begin{tabular}{||c||c|c|c|c|c|c||}\\hline\\hline   " << endl;
      for(i=1;i<j;i++) {
        string aus=lines[i].str();
        if(i<=i1 || (i>i1 && (((int) aus.find(")"))>0)) || ((((int) aus.find("posi"))>0)) ) {
          if(i!=i5 || (int) aus.find("A")>0) {
            cout << lines[i].str();
            if(i<=i2 || i==i1 || i==i3 || i==i4 ) cout << "\\hline";
            cout << endl;}}}
      cout << "\\hline\\hline" << endl << "\\end{tabular}}" << endl << "\\end{center}}" << endl;
    }
  }
  cout << "\\newpage" << endl;
  if(_BCC) {
    { // BCC2 TABLE
      cerr << "BCC2a TABLE" << endl; // BCC2 TABLE
      // cout << "\\subsection{BCC2a superstructures} " << endl;
      ostringstream lines[MAXLINESOSTREAM];
      j=1;
      lines[j++] << "Strukturbericht";lines[j++] << "Prototype      ";
      lines[j++] << "Label          ";lines[j++] << "Formula        ";
      lines[j++] << "Parent lattice ";lines[j++] << "C$_b$          ";
      lines[j++] << "Space Group    ";lines[j++] << "Pearson Symbol ";
      lines[j++] << "Lattice type   ";lines[j++] << "$(a,b,c)$ (\\AA)";lines[j++] << "$(\\alpha,\\beta,\\gamma)$ degrees";
      lines[j++] << "Superlattice   ";
      lines[j++]<<"Primitive vectors";lines[j++]<<"{\\bf a$_1$}/a";lines[j++]<<"{\\bf a$_2$}/a";lines[j++]<<"{\\bf a$_3$}/a";lines[j++]<<"(cartesian)";
      lines[j++]<<"Atomic positions";
      for(i=1;i<=AAmax;i++) lines[j++] <<"{\\bf A$"<<i<<"$}     ";
      for(i=1;i<=BBmax;i++) lines[j++] <<"{\\bf B$"<<i<<"$}     ";
      lines[j++] << "(fractional)   ";
      read_VASPIN(lines,71,"71   ","B2      ","BCC"," "," ",VERB);// 71      oC?       Cmma     ?   A2B2      #67    UNDEF          
      read_VASPIN(lines,73,"73,75","        ","BCC"," "," ",VERB);// 73,75   mP?       P2/m     ?   AB3       #10    UNDEF  
      read_VASPIN(lines,74,"74   ","          ","BCC"," "," ",VERB);// 74      mP?       P21/m    ?   A2B2      #11    UNDEF  
      read_VASPIN(lines,76,"76,78","C11$_b$   ","BCC"," "," ",VERB);// 76,78   tP?       P4/mmm   ?   AB3       #123   UNDEF  
      for(i=1;i<j;i++) lines[i] << "\t \\\\ ";
      cout << "{\\begin{center}" << endl;
      cout << "{\\footnotesize" << endl;
      cout << "\\begin{tabular}{||c||c|c|c|c|c|c||}\\hline\\hline   " << endl;
      for(i=1;i<j;i++) {
        string aus=lines[i].str();
        if(i<=i1 || (i>i1 && (((int) aus.find(")"))>0)) || ((((int) aus.find("posi"))>0)) ) {
          if(i!=i5 || (int) aus.find("A")>0) {
            cout << lines[i].str();
            if(i<=i2 || i==i1 || i==i3 || i==i4 ) cout << "\\hline";
            cout << endl;}}}
      //     cout << "\\hline \\hline" << endl << "\\end{tabular}}" << endl << "\\end{center}}" << endl;
      cout << "\\hline \\hline \\hline \\hline " << endl;
    }
    { // BCC2b TABLE
      cerr << "BCC2b TABLE" << endl; // BCC4 TABLE
      ostringstream lines[MAXLINESOSTREAM];
      j=1;
      lines[j++] << "Strukturbericht";lines[j++] << "Prototype      ";
      lines[j++] << "Label          ";lines[j++] << "Formula        ";
      lines[j++] << "Parent lattice ";lines[j++] << "C$_b$          ";
      lines[j++] << "Space Group    ";lines[j++] << "Pearson Symbol ";
      lines[j++] << "Lattice type   ";lines[j++] << "$(a,b,c)$ (\\AA)";lines[j++] << "$(\\alpha,\\beta,\\gamma)$ degrees";
      lines[j++] << "Superlattice   ";
      lines[j++]<<"Primitive vectors";lines[j++]<<"{\\bf a$_1$}/a";lines[j++]<<"{\\bf a$_2$}/a";lines[j++]<<"{\\bf a$_3$}/a";lines[j++]<<"(cartesian)";
      lines[j++]<<"Atomic positions";
      for(i=1;i<=AAmax;i++) lines[j++] <<"{\\bf A$"<<i<<"$}     ";
      for(i=1;i<=BBmax;i++) lines[j++] <<"{\\bf B$"<<i<<"$}     ";
      lines[j++] << "(fractional)   ";
      read_VASPIN(lines,77,"77   ","B11 ","BCC"," ","$\\gamma$-CuTi",VERB); // 77  tP4   P4/nmm  B11 gamma-CuTi #129   UNDEF  
      read_VASPIN(lines,79,"79,81","          ","BCC"," "," ",VERB);// 79,81   oI?       Immm    ?    AB3       #71    UNDEF
      read_VASPIN(lines,80,"80   ","          ","BCC"," "," ",VERB);// 80      oI?       Imma    ?    A2B2      #74    UNDEF  
      read_VASPIN(lines,82,"82,83","L6$_0$    ","BCC"," ","CuTi$_3$",VERB); // 82,83   tP4    P4/mmm     L6_0    CuTi3     #123    
      for(i=1;i<j;i++) lines[i] << "\t \\\\ ";
      // cout << "{\\begin{center}" << endl;
      // cout << "{\\footnotesize" << endl;
      // cout << "\\begin{tabular}{||c||c|c|c|c|c|c||}\\hline\\hline   " << endl;
      for(i=1;i<j;i++) {
        string aus=lines[i].str();
        if(i<=i1 || (i>i1 && (((int) aus.find(")"))>0)) || ((((int) aus.find("posi"))>0)) ) {
          if(i!=i5 || (int) aus.find("A")>0) {
            cout << lines[i].str();
            if(i<=i2 || i==i1 || i==i3 || i==i4 ) cout << "\\hline";
            cout << endl;}}}
      cout << "\\hline\\hline" << endl << "\\end{tabular}}" << endl << "\\end{center}}" << endl;
    }
  }
  cout << "\\newpage" << endl;
  if(_BCC) {
    { // BCC2c TABLE
      cerr << "BCC2c TABLE" << endl; // BCC2 TABLE
      //    cout << "\\subsection{BCC2 superstructures} " << endl;
      ostringstream lines[MAXLINESOSTREAM];
      j=1;
      lines[j++] << "Strukturbericht";lines[j++] << "Prototype      ";
      lines[j++] << "Label          ";lines[j++] << "Formula        ";
      lines[j++] << "Parent lattice ";lines[j++] << "C$_b$          ";
      lines[j++] << "Space Group    ";lines[j++] << "Pearson Symbol ";
      lines[j++] << "Lattice type   ";lines[j++] << "$(a,b,c)$ (\\AA)";lines[j++] << "$(\\alpha,\\beta,\\gamma)$ degrees";
      lines[j++] << "Superlattice   ";
      lines[j++]<<"Primitive vectors";lines[j++]<<"{\\bf a$_1$}/a";lines[j++]<<"{\\bf a$_2$}/a";lines[j++]<<"{\\bf a$_3$}/a";lines[j++]<<"(cartesian)";
      lines[j++]<<"Atomic positions";
      for(i=1;i<=AAmax;i++) lines[j++] <<"{\\bf A$"<<i<<"$}     ";
      for(i=1;i<=BBmax;i++) lines[j++] <<"{\\bf B$"<<i<<"$}     ";
      lines[j++] << "(fractional)   ";
      read_VASPIN(lines,84,"84,86","D0$_3$    ","BCC"," ","AlFe$_3$",VERB); // 84,86   cF16   Fm-3m      D0_3    AlFe3     #225
      read_VASPIN(lines,85,"85   ","B32       ","BCC"," ","NaTl    ",VERB); // 85      cF16   Fd-3m      B32     NaTl      #227  
      read_VASPIN(lines,87,"87/92","          ","BCC"," "," ",VERB);
      read_VASPIN(lines,93,"93/98","          ","BCC"," "," ",VERB);
      for(i=1;i<j;i++) lines[i] << "\t \\\\ ";
      cout << "{\\begin{center}" << endl;
      cout << "{\\footnotesize" << endl;
      cout << "\\begin{tabular}{||c||c|c|c|c|c||}\\hline\\hline   " << endl;
      for(i=1;i<j;i++) {
        string aus=lines[i].str();
        if(i<=i1 || (i>i1 && (((int) aus.find(")"))>0)) || ((((int) aus.find("posi"))>0)) ) {
          if(i!=i5 || (int) aus.find("A")>0) {
            cout << lines[i].str();
            if(i<=i2 || i==i1 || i==i3 || i==i4 ) cout << "\\hline";
            cout << endl;}}}
      cout << "\\hline\\hline" << endl << "\\end{tabular}}" << endl << "\\end{center}}" << endl;
    }
  }
  LatexSideMarginThesisOFF();
  cout << "\\newpage" << endl;
  if(VERB) cerr << " *********************************************************************** HCP " << endl;
  if(VERB) cerr << " *********************************************************************** HCP " << endl;
  if(VERB) cerr << " *********************************************************************** HCP " << endl;
  if(VERB) cerr << " *********************************************************************** HCP " << endl;
  if(_HCP) {
    { // HCP1a TABLE
      cout << "\\subsection{{HCP superstructures}} " << endl;
      cout << "\\label{appendix.HCP_superstructures} " << endl;
      LatexSideMarginThesisON();
      cerr << "HCP1a TABLE" << endl; // HCP1 TABLE
      ostringstream lines[MAXLINESOSTREAM];
      j=1;
      lines[j++] << "Strukturbericht";lines[j++] << "Prototype      ";
      lines[j++] << "Label          ";lines[j++] << "Formula        ";
      lines[j++] << "Parent lattice ";lines[j++] << "C$_b$          ";
      lines[j++] << "Space Group    ";lines[j++] << "Pearson Symbol ";
      lines[j++] << "Lattice type   ";lines[j++] << "$(a,b,c)$ (\\AA)";lines[j++] << "$(\\alpha,\\beta,\\gamma)$ degrees";
      lines[j++] << "Superlattice   ";
      lines[j++]<<"Primitive vectors";lines[j++]<<"{\\bf a$_1$}/a";lines[j++]<<"{\\bf a$_2$}/a";lines[j++]<<"{\\bf a$_3$}/a";lines[j++]<<"(cartesian)";
      lines[j++]<<"Atomic positions";
      for(i=1;i<=AAmax;i++) lines[j++] <<"{\\bf A$"<<i<<"$}     ";
      for(i=1;i<=BBmax;i++) lines[j++] <<"{\\bf B$"<<i<<"$}     ";
      lines[j++] << "(fractional)   ";
      read_VASPIN(lines,115,"115,117","A3        ","HCP"," ","Mg  ",VERB);// 115,117 hcp/hcp2  P6_3/mmc      A3      Mg        #194
      read_VASPIN(lines,116,"116    ","B$_h$     ","HCP"," ","WC  ",VERB);// 116     hP2       P_6m2         B_h     WC        #187          
      read_VASPIN(lines,118,"118,121","          ","HCP"," ","    ",VERB);// 118,121 oP?       Pmm2          ?       AB3       #25    UNDEF  
      read_VASPIN(lines,119,"119    ","          ","HCP"," ","CuTe",VERB);// 119     oP4       Pmmn          -       CuTe      #59    
      for(i=1;i<j;i++) lines[i] << "\t \\\\ ";
      cout << "{\\begin{center}" << endl;
      cout << "{\\footnotesize" << endl;
      cout << "\\begin{tabular}{||c||c|c|c|c||}\\hline\\hline   " << endl;
      for(i=1;i<j;i++) {
        string aus=lines[i].str();
        if(i<=i1 || (i>i1 && (((int) aus.find(")"))>0)) || ((((int) aus.find("posi"))>0)) ) {
          if(i!=i5 || (int) aus.find("A")>0) {
            cout << lines[i].str();
            if(i<=i2 || i==i1 || i==i3 || i==i4 ) cout << "\\hline";
            cout << endl;}}}
      //     cout << "\\hline \\hline" << endl << "\\end{tabular}}" << endl << "\\end{center}}" << endl;
      cout << "\\hline \\hline \\hline \\hline " << endl;
    }
    { // HCP1b TABLE
      cerr << "HCP1b TABLE" << endl; // HCP2 TABLE
      ostringstream lines[MAXLINESOSTREAM];
      j=1;
      lines[j++] << "Strukturbericht";lines[j++] << "Prototype      ";
      lines[j++] << "Label          ";lines[j++] << "Formula        ";
      lines[j++] << "Parent lattice ";lines[j++] << "C$_b$          ";
      lines[j++] << "Space Group    ";lines[j++] << "Pearson Symbol ";
      lines[j++] << "Lattice type   ";lines[j++] << "$(a,b,c)$ (\\AA)";lines[j++] << "$(\\alpha,\\beta,\\gamma)$ degrees";
      lines[j++] << "Superlattice   ";
      lines[j++]<<"Primitive vectors";lines[j++]<<"{\\bf a$_1$}/a";lines[j++]<<"{\\bf a$_2$}/a";lines[j++]<<"{\\bf a$_3$}/a";lines[j++]<<"(cartesian)";
      lines[j++]<<"Atomic positions";
      for(i=1;i<=AAmax;i++) lines[j++] <<"{\\bf A$"<<i<<"$}     ";
      for(i=1;i<=BBmax;i++) lines[j++] <<"{\\bf B$"<<i<<"$}     ";
      lines[j++] << "(fractional)   ";
      read_VASPIN(lines,120,"120    ","B19  ","HCP"," ","AuCd",VERB);// 120     oP4       Pmma          B19     AuCd      #51
      read_VASPIN(lines,122,"122,124","     ","HCP"," ","    ",VERB);// 122,124 oI?       Imm2          ?       AB3       #44    UNDEF  
      read_VASPIN(lines,123,"123    ","     ","HCP"," "," ",VERB); // 123     mC?       C2/m          ?       A2B2      #12    UNDEF  
      read_VASPIN(lines,125,"125,127","     ","HCP"," "," ",VERB); // 125,127 hP?       P-6m2         ?       AB3       #187   UNDEF
      for(i=1;i<j;i++) lines[i] << "\t \\\\ ";
      // cout << "{\\begin{center}" << endl;
      // cout << "{\\footnotesize" << endl;
      // cout << "\\begin{tabular}{||c||c|c|c|c||}\\hline\\hline   " << endl;
      for(i=1;i<j;i++) {
        string aus=lines[i].str();
        if(i<=i1 || (i>i1 && (((int) aus.find(")"))>0)) || ((((int) aus.find("posi"))>0)) ) {
          if(i!=i5 || (int) aus.find("A")>0) {
            cout << lines[i].str();
            if(i<=i2 || i==i1 || i==i3 || i==i4 ) cout << "\\hline";
            cout << endl;}}}
      cout << "\\hline\\hline" << endl << "\\end{tabular}}" << endl << "\\end{center}}" << endl;
    }
  }
  cout << "\\newpage" << endl;
  if(_HCP) {
    { // HCP2a TABLE
      cerr << "HCP2a TABLE" << endl; // HCP3 TABLE
      ostringstream lines[MAXLINESOSTREAM];
      j=1;
      lines[j++] << "Strukturbericht";lines[j++] << "Prototype      ";
      lines[j++] << "Label          ";lines[j++] << "Formula        ";
      lines[j++] << "Parent lattice ";lines[j++] << "C$_b$          ";
      lines[j++] << "Space Group    ";lines[j++] << "Pearson Symbol ";
      lines[j++] << "Lattice type   ";lines[j++] << "$(a,b,c)$ (\\AA)";lines[j++] << "$(\\alpha,\\beta,\\gamma)$ degrees";
      lines[j++] << "Superlattice   ";
      lines[j++]<<"Primitive vectors";lines[j++]<<"{\\bf a$_1$}/a";lines[j++]<<"{\\bf a$_2$}/a";lines[j++]<<"{\\bf a$_3$}/a";lines[j++]<<"(cartesian)";
      lines[j++]<<"Atomic positions";
      for(i=1;i<=AAmax;i++) lines[j++] <<"{\\bf A$"<<i<<"$}     ";
      for(i=1;i<=BBmax;i++) lines[j++] <<"{\\bf B$"<<i<<"$}     ";
      lines[j++] << "(fractional)   ";
      read_VASPIN(lines,126,"126    "," ","HCP"," "," ",VERB); // 126     hP?       P-3m1         ?       A2B2      #164   UNDEF
      read_VASPIN(lines,128,"128,132"," ","HCP"," "," ",VERB); // 128,132 mC?       Cm            ?       AB5       #8     UNDEF  
      read_VASPIN(lines,134,"129,134"," ","HCP"," "," ",VERB); // 129,134 mC?       C2/m          ?       A4B2      #12
      read_VASPIN(lines,137,"130,137"," ","HCP"," "," ",VERB); // 130.137 mC?       C2/m          ?       A4B2      #12
      for(i=1;i<j;i++) lines[i] << "\t \\\\ ";
      cout << "{\\begin{center}" << endl;
      cout << "{\\footnotesize" << endl;
      cout << "\\begin{tabular}{||c||c|c|c|c|c|c||}\\hline\\hline   " << endl;
      for(i=1;i<j;i++) {
        string aus=lines[i].str();
        if(i<=i1 || (i>i1 && (((int) aus.find(")"))>0)) || ((((int) aus.find("posi"))>0)) ) {
          if(i!=i5 || (int) aus.find("A")>0) {
            cout << lines[i].str();
            if(i<=i2 || i==i1 || i==i3 || i==i4 ) cout << "\\hline";
            cout << endl;}}}
      //     cout << "\\hline \\hline" << endl << "\\end{tabular}}" << endl << "\\end{center}}" << endl;
      cout << "\\hline \\hline \\hline \\hline " << endl;
    }
    { // HCP2b TABLE
      cerr << "HCP2b TABLE" << endl; // HCP3 TABLE
      ostringstream lines[MAXLINESOSTREAM];
      j=1;
      lines[j++] << "Strukturbericht";lines[j++] << "Prototype      ";
      lines[j++] << "Label          ";lines[j++] << "Formula        ";
      lines[j++] << "Parent lattice ";lines[j++] << "C$_b$          ";
      lines[j++] << "Space Group    ";lines[j++] << "Pearson Symbol ";
      lines[j++] << "Lattice type   ";lines[j++] << "$(a,b,c)$ (\\AA)";lines[j++] << "$(\\alpha,\\beta,\\gamma)$ degrees";
      lines[j++] << "Superlattice   ";
      lines[j++]<<"Primitive vectors";lines[j++]<<"{\\bf a$_1$}/a";lines[j++]<<"{\\bf a$_2$}/a";lines[j++]<<"{\\bf a$_3$}/a";lines[j++]<<"(cartesian)";
      lines[j++]<<"Atomic positions";
      for(i=1;i<=AAmax;i++) lines[j++] <<"{\\bf A$"<<i<<"$}     ";
      for(i=1;i<=BBmax;i++) lines[j++] <<"{\\bf B$"<<i<<"$}     ";
      lines[j++] << "(fractional)   ";
      read_VASPIN(lines,135,"131,135"," ","HCP"," "," ",VERB); // 131,135 mC?       C2/m          ?       A4B2      #12
      read_VASPIN(lines,133,"133,140"," ","HCP"," "," ",VERB); // 133,140 mC?       Cm            ?       A2B4      #8     UNDEF  
      read_VASPIN(lines,136,"136    "," ","HCP"," "," ",VERB); // 136     mC?       Cm            ?       A3B3      #8     UNDEF  
      read_VASPIN(lines,138,"138,139"," ","HCP"," "," ",VERB); // 138,139 mC?       Cm            ?       A3B3      #8     UNDEF  
      for(i=1;i<j;i++) lines[i] << "\t \\\\ ";
      // cout << "{\\begin{center}" << endl;
      // cout << "{\\footnotesize" << endl;
      // cout << "\\begin{tabular}{||c||c|c|c|c|c|c||}\\hline\\hline   " << endl;
      for(i=1;i<j;i++) {
        string aus=lines[i].str();
        if(i<=i1 || (i>i1 && (((int) aus.find(")"))>0)) || ((((int) aus.find("posi"))>0)) ) {
          if(i!=i5 || (int) aus.find("A")>0) {
            cout << lines[i].str();
            if(i<=i2 || i==i1 || i==i3 || i==i4 ) cout << "\\hline";
            cout << endl;}}}
      cout << "\\hline\\hline" << endl << "\\end{tabular}}" << endl << "\\end{center}}" << endl;
    }
  }
  cout << "\\newpage" << endl;
  if(_HCP) {
    { // HCP3a TABLE
      cerr << "HCP3a TABLE" << endl; // HCP3 TABLE
      ostringstream lines[MAXLINESOSTREAM];
      j=1;
      lines[j++] << "Strukturbericht";lines[j++] << "Prototype      ";
      lines[j++] << "Label          ";lines[j++] << "Formula        ";
      lines[j++] << "Parent lattice ";lines[j++] << "C$_b$          ";
      lines[j++] << "Space Group    ";lines[j++] << "Pearson Symbol ";
      lines[j++] << "Lattice type   ";lines[j++] << "$(a,b,c)$ (\\AA)";lines[j++] << "$(\\alpha,\\beta,\\gamma)$ degrees";
      lines[j++] << "Superlattice   ";
      lines[j++]<<"Primitive vectors";lines[j++]<<"{\\bf a$_1$}/a";lines[j++]<<"{\\bf a$_2$}/a";lines[j++]<<"{\\bf a$_3$}/a";lines[j++]<<"(cartesian)";
      lines[j++]<<"Atomic positions";
      for(i=1;i<=AAmax;i++) lines[j++] <<"{\\bf A$"<<i<<"$}     ";
      for(i=1;i<=BBmax;i++) lines[j++] <<"{\\bf B$"<<i<<"$}     ";
      lines[j++] << "(fractional)   ";
      read_VASPIN(lines,141,"141,145"," ","HCP"," "," ",VERB); // 141,145 oA?  Amm2  ?   AB5       #38
      read_VASPIN(lines,147,"142,147"," ","HCP"," "," ",VERB); // 142,147 oC?  Cmcm  ?   A4B2      #63
      read_VASPIN(lines,148,"144,148"," ","HCP"," "," ",VERB); // 144,148  Cmcm      ?   A4B2      #63           NOT RUN
      read_VASPIN(lines,146,"146,153"," ","HCP"," "," ",VERB); // 146,153  Amm2      ?   A2B4      #38           NOT RUN
      for(i=1;i<j;i++) lines[i] << "\t \\\\ ";
      cout << "{\\begin{center}" << endl;
      cout << "{\\footnotesize" << endl;
      cout << "\\begin{tabular}{||c||c|c|c|c|c|c||}\\hline\\hline   " << endl;
      for(i=1;i<j;i++) {
        string aus=lines[i].str();
        if(i<=i1 || (i>i1 && (((int) aus.find(")"))>0)) || ((((int) aus.find("posi"))>0)) ) {
          if(i!=i5 || (int) aus.find("A")>0) {
            cout << lines[i].str();
            if(i<=i2 || i==i1 || i==i3 || i==i4 ) cout << "\\hline";
            cout << endl;}}}
      cout << "\\hline \\hline" << endl << "\\end{tabular}}" << endl << "\\end{center}}" << endl;
    }
  }
  LatexSideMarginThesisOFF();
  cout << "\\newpage" << endl;  
  if(VERB) cerr << " *********************************************************************** EXT " << endl;
  if(VERB) cerr << " *********************************************************************** EXT " << endl;
  if(VERB) cerr << " *********************************************************************** EXT " << endl;
  if(VERB) cerr << " *********************************************************************** EXT " << endl;
  if(_EXT) {
    cout << "\\subsection{{EXTRA structures}} " << endl;
    cout << "\\label{appendix.EXTRA_superstructures} " << endl;
    LatexSideMarginThesisON();
    { // EXT1a TABLE
      cerr << "EXT1a TABLE" << endl;// EXT1 TABLE
      ostringstream lines[MAXLINESOSTREAM];
      j=1;
      lines[j++] << "Strukturbericht";lines[j++] << "Prototype      ";
      lines[j++] << "Label          ";lines[j++] << "Formula        ";
      lines[j++] << "Parent lattice ";lines[j++] << "C$_b$          ";
      lines[j++] << "Space Group    ";lines[j++] << "Pearson Symbol ";
      lines[j++] << "Lattice type   ";lines[j++] << "$(a,b,c)$ (\\AA)";lines[j++] << "$(\\alpha,\\beta,\\gamma)$ degrees";
      lines[j++] << "Superlattice   ";
      lines[j++]<<"Primitive vectors";lines[j++]<<"{\\bf a$_1$}/a";lines[j++]<<"{\\bf a$_2$}/a";lines[j++]<<"{\\bf a$_3$}/a";lines[j++]<<"(cartesian)";
      lines[j++]<<"Atomic positions";
      for(i=1;i<=AAmax;i++) lines[j++] <<"{\\bf A$"<<i<<"$}     ";
      for(i=1;i<=BBmax;i++) lines[j++] <<"{\\bf B$"<<i<<"$}     ";
      lines[j++] << "(fractional)   ";
      read_VASPIN(lines,178,"178,179","C14       "," "," ","MgZn$_2$   ",VERB);// 178,179  hP12   P6_3/mm   C14    MgZn2 #194 A4B8Hex Laves
      read_VASPIN(lines,180,"180,181","D0$_{11}$ "," "," ","Fe$_3$C    ",VERB);// 180,181  oP16   Pnma      D0_11  Fe3C
      read_VASPIN(lines,183,"182,183","C15       "," "," ","Cu$_2$Mg   ",VERB);// 182,183  cF24   Fd-3m     C15    Cu2Mg Cub Laves Phase
      read_VASPIN(lines,184,"184,185","A15       "," "," ","Cr$_3$Si   ",VERB);// 184,185  cP8    Pm3(bar)n A15    Cr3Si
      for(i=1;i<j;i++) lines[i] << "\t \\\\ ";
      cout << "{\\begin{center}" << endl;
      cout << "{\\footnotesize" << endl;
      cout << "\\begin{tabular}{||c||c|c|c|c||}\\hline\\hline   " << endl;
      for(i=1;i<j;i++) {
        string aus=lines[i].str();
        if(i<=i1 || (i>i1 && (((int) aus.find(")"))>0)) || ((((int) aus.find("posi"))>0)) ) {
          if(i!=i5 || (int) aus.find("A")>0) {
            if((i<i6 || i>i7) && (int) aus.find("aren")<0) {
              cout << lines[i].str();
              if(i<=i2 || i==i1 || i==i3 || i==i4 ) cout << "\\hline";
              cout << endl;}}}}
      //     cout << "\\hline \\hline" << endl << "\\end{tabular}}" << endl << "\\end{center}}" << endl;
      cout << "\\hline \\hline \\hline \\hline " << endl;
    }
    { // EXT1b TABLE
      cerr << "EXT1b TABLE" << endl; // EXT2 TABLE
      ostringstream lines[MAXLINESOSTREAM];
      j=1;
      lines[j++] << "Strukturbericht";lines[j++] << "Prototype      ";
      lines[j++] << "Label          ";lines[j++] << "Formula        ";
      lines[j++] << "Parent lattice ";lines[j++] << "C$_b$          ";
      lines[j++] << "Space Group    ";lines[j++] << "Pearson Symbol ";
      lines[j++] << "Lattice type   ";lines[j++] << "$(a,b,c)$ (\\AA)";lines[j++] << "$(\\alpha,\\beta,\\gamma)$ degrees";
      lines[j++] << "Superlattice   ";
      lines[j++]<<"Primitive vectors";lines[j++]<<"{\\bf a$_1$}/a";lines[j++]<<"{\\bf a$_2$}/a";lines[j++]<<"{\\bf a$_3$}/a";lines[j++]<<"(cartesian)";
      lines[j++]<<"Atomic positions";
      for(i=1;i<=AAmax;i++) lines[j++] <<"{\\bf A$"<<i<<"$}     ";
      for(i=1;i<=BBmax;i++) lines[j++] <<"{\\bf B$"<<i<<"$}     ";
      lines[j++] << "(fractional)   ";
      read_VASPIN(lines,186,"186,187","D0$_{19}$ "," "," ","Ni$_3$Sn        ",VERB);// 186,187  hP8    P6_3/mmc  D0_19  Ni3Sn #194
      read_VASPIN(lines,188,"188,189","C49       "," "," ","ZrSi$_2$        ",VERB);// 188,189  hP6    Cmcm      C49    ZrSi2
      read_VASPIN(lines,190,"190,191"," $\\omega$ Z=1/4"," "," ","$\\omega$ Z=1/4 ",VERB);// 190,191  hP3  P3m1  some omega phase with Z=1/4 !
      read_VASPIN(lines,192,"192,193","B33       "," "," ","CrB             ",VERB);// 192,193  oC8    Cmcm      B33    CrB     #63
      for(i=1;i<j;i++) lines[i] << "\t \\\\ ";
      // cout << "{\\begin{center}" << endl;
      // cout << "{\\footnotesize" << endl;
      // cout << "\\begin{tabular}{||c||c|c|c|c||}\\hline\\hline   " << endl;
      for(i=1;i<j;i++) {
        string aus=lines[i].str();
        if(i<=i1 || (i>i1 && (((int) aus.find(")"))>0)) || ((((int) aus.find("posi"))>0)) ) {
          if(i!=i5 || (int) aus.find("A")>0) {
            if((i<i6 || i>i7) && (int) aus.find("aren")<0) {
              cout << lines[i].str();
              if(i<=i2 || i==i1 || i==i3 || i==i4 ) cout << "\\hline";
              cout << endl;}}}}
      cout << "\\hline\\hline" << endl << "\\end{tabular}}" << endl << "\\end{center}}" << endl;
    }
  }
  cout << "\\newpage" << endl;  
  if(_EXT) {
    { // EXT2a TABLE
      cerr << "EXT2a TABLE" << endl;// EXT2a TABLE
      // cout << "\\subsection{EXT2a superstructures} " << endl;
      ostringstream lines[MAXLINESOSTREAM];
      j=1;
      lines[j++] << "Strukturbericht";lines[j++] << "Prototype      ";
      lines[j++] << "Label          ";lines[j++] << "Formula        ";
      lines[j++] << "Parent lattice ";lines[j++] << "C$_b$          ";
      lines[j++] << "Space Group    ";lines[j++] << "Pearson Symbol ";
      lines[j++] << "Lattice type   ";lines[j++] << "$(a,b,c)$ (\\AA)";lines[j++] << "$(\\alpha,\\beta,\\gamma)$ degrees";
      lines[j++] << "Superlattice   ";
      lines[j++]<<"Primitive vectors";lines[j++]<<"{\\bf a$_1$}/a";lines[j++]<<"{\\bf a$_2$}/a";lines[j++]<<"{\\bf a$_3$}/a";lines[j++]<<"(cartesian)";
      lines[j++]<<"Atomic positions";
      for(i=1;i<=AAmax;i++) lines[j++] <<"{\\bf A$"<<i<<"$}     ";
      for(i=1;i<=BBmax;i++) lines[j++] <<"{\\bf B$"<<i<<"$}     ";
      lines[j++] << "(fractional)   ";
      read_VASPIN(lines,201,"201    ","B1     "," "," ","NaCl         ",VERB);// 201      cF8    Fm-3m     B1     NaCl    #225
      read_VASPIN(lines,202,"202,203","D1$_3$ "," "," ","BaAl$_4$     ",VERB);// 202,203  tI10   I4/mmm    D1_3   BaAl4   #139  NEW
      read_VASPIN(lines,204,"204,205","D2$_d$ "," "," ","CaCu$_5$     ",VERB);// 204,205  hP6    P6/mmm    D2_d   CaCu5   #191
      // read_VASPIN(lines,207,"206,207","D2$_b$ "," "," ","ThMn$_{12}$  ",VERB);// 206,207  tI26   I4/mmm    D2_b   ThMn12  #139      
      read_VASPIN(lines,208,"208,209","C22      "," "," ","Fe$_2$P       ",VERB);// 208,209  hP9    P-62m     C22    Fe2P    #189      
      for(i=1;i<j;i++) lines[i] << "\t \\\\ ";
      cout << "{\\begin{center}" << endl;
      cout << "{\\footnotesize" << endl;
      cout << "\\begin{tabular}{||c||c|c|c|c||}\\hline\\hline   " << endl;
      for(i=1;i<j;i++) {
        string aus=lines[i].str();
        if(i<=i1 || (i>i1 && (((int) aus.find(")"))>0)) || ((((int) aus.find("posi"))>0)) ) {
          if(i!=i5 || (int) aus.find("A")>0) {
            if((i<i6 || i>i7) && (int) aus.find("aren")<0) {
              cout << lines[i].str();
              if(i<=i2 || i==i1 || i==i3 || i==i4 ) cout << "\\hline";
              cout << endl;}}}}
      //     cout << "\\hline \\hline" << endl << "\\end{tabular}}" << endl << "\\end{center}}" << endl;
      cout << "\\hline \\hline \\hline \\hline " << endl;
    }
    { // EXT2b TABLE
      cerr << "EXT2b TABLE" << endl; // EXT4 TABLE
      ostringstream lines[MAXLINESOSTREAM];
      j=1;
      lines[j++] << "Strukturbericht";lines[j++] << "Prototype      ";
      lines[j++] << "Label          ";lines[j++] << "Formula        ";
      lines[j++] << "Parent lattice ";lines[j++] << "C$_b$          ";
      lines[j++] << "Space Group    ";lines[j++] << "Pearson Symbol ";
      lines[j++] << "Lattice type   ";lines[j++] << "$(a,b,c)$ (\\AA)";lines[j++] << "$(\\alpha,\\beta,\\gamma)$ degrees";
      lines[j++] << "Superlattice   ";
      lines[j++]<<"Primitive vectors";lines[j++]<<"{\\bf a$_1$}/a";lines[j++]<<"{\\bf a$_2$}/a";lines[j++]<<"{\\bf a$_3$}/a";lines[j++]<<"(cartesian)";
      lines[j++]<<"Atomic positions";
      for(i=1;i<=AAmax;i++) lines[j++] <<"{\\bf A$"<<i<<"$}     ";
      for(i=1;i<=BBmax;i++) lines[j++] <<"{\\bf B$"<<i<<"$}     ";
      lines[j++] << "(fractional)   ";
      //  read_VASPIN(lines,208,"208,209","C22      "," "," ","Fe$_2$P       ",VERB);// 208,209  hP9    P-62m     C22    Fe2P    #189      
      read_VASPIN(lines,211,"210,211","C37      "," "," ","Co$_2$Si      ",VERB);// 210,211  oP12   Pnma      C37    Co2Si   #62  AL1851.cif    
      //   read_VASPIN(lines,212,"212,213","         "," "," ","Th$_2$Zn$_{17}$",VERB);// 212,213  hR19                    Th2Zn17 #166      
      // read_VASPIN(lines,215,"214,215","D7$_3$   "," "," ","Th$_3$P$_4$   ",VERB);// 214,215  cI28   I-43d     D7_3   Th3P4   #220      
      read_VASPIN(lines,216,"216,217","C32         "," "," ","AlB$_2$          ",VERB);// 216,217  hP3    P6/mmm    C32    AlB2    #191
      read_VASPIN(lines,218,"218    ","B3          "," "," ","ZnS              ",VERB);// 218      cF8    F-43m     B3     ZnS     #216
      read_VASPIN(lines,219,"219    ","B4          "," "," ","ZnS              ",VERB);// 219      hP4    P63/mc    B4     ZnS     #186
      for(i=1;i<j;i++) lines[i] << "\t \\\\ ";
      // cout << "{\\begin{center}" << endl;
      // cout << "{\\footnotesize" << endl;
      // cout << "\\begin{tabular}{||c||c|c|c|c||}\\hline\\hline   " << endl;
      for(i=1;i<j;i++) {
        string aus=lines[i].str();
        if(i<=i1 || (i>i1 && (((int) aus.find(")"))>0)) || ((((int) aus.find("posi"))>0)) ) {
          if(i!=i5 || (int) aus.find("A")>0) {
            if((i<i6 || i>i7) && (int) aus.find("aren")<0) {
              cout << lines[i].str();
              if(i<=i2 || i==i1 || i==i3 || i==i4 ) cout << "\\hline";
              cout << endl;}}}}
      cout << "\\hline\\hline" << endl << "\\end{tabular}}" << endl << "\\end{center}}" << endl;
    }
  }
  cout << "\\newpage" << endl;  
  if(_EXT) {
    {  // EXT3a TABLE
      cerr << "EXT3a TABLE" << endl;// EXT3a TABLE
      ostringstream lines[MAXLINESOSTREAM];
      j=1;
      lines[j++] << "Strukturbericht";lines[j++] << "Prototype      ";
      lines[j++] << "Label          ";lines[j++] << "Formula        ";
      lines[j++] << "Parent lattice ";lines[j++] << "C$_b$          ";
      lines[j++] << "Space Group    ";lines[j++] << "Pearson Symbol ";
      lines[j++] << "Lattice type   ";lines[j++] << "$(a,b,c)$ (\\AA)";lines[j++] << "$(\\alpha,\\beta,\\gamma)$ degrees";
      lines[j++] << "Superlattice   ";
      lines[j++]<<"Primitive vectors";lines[j++]<<"{\\bf a$_1$}/a";lines[j++]<<"{\\bf a$_2$}/a";lines[j++]<<"{\\bf a$_3$}/a";lines[j++]<<"(cartesian)";
      lines[j++]<<"Atomic positions";
      for(i=1;i<=AAmax;i++) lines[j++] <<"{\\bf A$"<<i<<"$}     ";
      for(i=1;i<=BBmax;i++) lines[j++] <<"{\\bf B$"<<i<<"$}     ";
      lines[j++] << "(fractional)   ";
      read_VASPIN(lines,220,"220,221","B8$_1$   "," "," ","NiAs             ",VERB);// 220,221  hP4    P63/mmc   B8_1   NiAs    #194
      read_VASPIN(lines,223,"222,223","D8$_8$   "," "," ","Mn$_5$Si$_3$     ",VERB);// 222,223  hP16   P63/mcm   D8_8   Mn5Si3  #193      
      //read_VASPIN(lines,224,"224,225","         "," "," ","Th$_2$Ni$_{17}$  ",VERB);// 224,225  hP38                    Th2Ni17 #194    
      read_VASPIN(lines,226,"226,227","A15      "," "," ","Cr$_3$Si         ",VERB);// 226,227  cP8    Pm-3n     A15    Cr3Si   #223  
      read_VASPIN(lines,229,"228,229","C38      "," "," ","Cu$_2$Sb         ",VERB);// 228,229  tP6    P4/nmm    C38    Cu2Sb   #129      
      for(i=1;i<j;i++) lines[i] << "\t \\\\ ";
      cout << "{\\begin{center}" << endl;
      cout << "{\\footnotesize" << endl;
      cout << "\\begin{tabular}{||c||c|c|c|c||}\\hline\\hline   " << endl;
      for(i=1;i<j;i++) {
        string aus=lines[i].str();
        if(i<=i1 || (i>i1 && (((int) aus.find(")"))>0)) || ((((int) aus.find("posi"))>0)) ) {
          if(i!=i5 || (int) aus.find("A")>0) {
            if((i<i6 || i>i7) && (int) aus.find("aren")<0) {
              cout << lines[i].str();
              if(i<=i2 || i==i1 || i==i3 || i==i4 ) cout << "\\hline";
              cout << endl;}}}}
      //     cout << "\\hline \\hline" << endl << "\\end{tabular}}" << endl << "\\end{center}}" << endl;
      cout << "\\hline \\hline \\hline \\hline " << endl;
    }
    { // EXT3b TABLE
      cerr << "EXT3b TABLE" << endl; // EXT3b TABLE
      ostringstream lines[MAXLINESOSTREAM];
      j=1;
      lines[j++] << "Strukturbericht";lines[j++] << "Prototype      ";
      lines[j++] << "Label          ";lines[j++] << "Formula        ";
      lines[j++] << "Parent lattice ";lines[j++] << "C$_b$          ";
      lines[j++] << "Space Group    ";lines[j++] << "Pearson Symbol ";
      lines[j++] << "Lattice type   ";lines[j++] << "$(a,b,c)$ (\\AA)";lines[j++] << "$(\\alpha,\\beta,\\gamma)$ degrees";
      lines[j++] << "Superlattice   ";
      lines[j++]<<"Primitive vectors";lines[j++]<<"{\\bf a$_1$}/a";lines[j++]<<"{\\bf a$_2$}/a";lines[j++]<<"{\\bf a$_3$}/a";lines[j++]<<"(cartesian)";
      lines[j++]<<"Atomic positions";
      for(i=1;i<=AAmax;i++) lines[j++] <<"{\\bf A$"<<i<<"$}     ";
      for(i=1;i<=BBmax;i++) lines[j++] <<"{\\bf B$"<<i<<"$}     ";
      lines[j++] << "(fractional)   ";
      read_VASPIN(lines,231,"230,231","        "," "," ","CuAl$_2$ ",VERB);// 230,231  tI12                  CuAl2   #140      
      read_VASPIN(lines,232,"232,233","C11$_b$ "," "," ","MoSi$_2$ ",VERB);// 232,233  tI6    I4/mmm  C11_b  MoSi2   #139  
      read_VASPIN(lines,234,"234,235","C16     "," "," ","Al$_2$Cu ",VERB);// 234,235  tl12   I4/mcm  C16    Al2Cu   #140 130617.cif
      read_VASPIN(lines,246,"246,247","        "," "," ","CaIn$_2$         ",VERB);// 246,247  hP6    P63/mmc          CaIn2   #194 AL5296.cif  
      for(i=1;i<j;i++) lines[i] << "\t \\\\ ";
      // cout << "{\\begin{center}" << endl;
      // cout << "{\\footnotesize" << endl;
      // cout << "\\begin{tabular}{||c||c|c|c|c||}\\hline\\hline   " << endl;
      for(i=1;i<j;i++) {
        string aus=lines[i].str();
        if(i<=i1 || (i>i1 && (((int) aus.find(")"))>0)) || ((((int) aus.find("posi"))>0)) ) {
          if(i!=i5 || (int) aus.find("A")>0) {
            if((i<i6 || i>i7) && (int) aus.find("aren")<0) {
              cout << lines[i].str();
              if(i<=i2 || i==i1 || i==i3 || i==i4 ) cout << "\\hline";
              cout << endl;}}}}
      cout << "\\hline\\hline" << endl << "\\end{tabular}}" << endl << "\\end{center}}" << endl;
    }
  }
  cout << "\\newpage" << endl;  
  if(_EXT) {
    {  // EXT4a TABLE
      cerr << "EXT4a TABLE" << endl;// EXT4a TABLE
      // cout << "\\subsection{EXT7 superstructures} " << endl;
      ostringstream lines[MAXLINESOSTREAM];
      j=1;
      lines[j++] << "Strukturbericht";lines[j++] << "Prototype      ";
      lines[j++] << "Label          ";lines[j++] << "Formula        ";
      lines[j++] << "Parent lattice ";lines[j++] << "C$_b$          ";
      lines[j++] << "Space Group    ";lines[j++] << "Pearson Symbol ";
      lines[j++] << "Lattice type   ";lines[j++] << "$(a,b,c)$ (\\AA)";lines[j++] << "$(\\alpha,\\beta,\\gamma)$ degrees";
      lines[j++] << "Superlattice   ";
      lines[j++]<<"Primitive vectors";lines[j++]<<"{\\bf a$_1$}/a";lines[j++]<<"{\\bf a$_2$}/a";lines[j++]<<"{\\bf a$_3$}/a";lines[j++]<<"(cartesian)";
      lines[j++]<<"Atomic positions";
      for(i=1;i<=AAmax;i++) lines[j++] <<"{\\bf A$"<<i<<"$}     ";
      for(i=1;i<=BBmax;i++) lines[j++] <<"{\\bf B$"<<i<<"$}     ";
      lines[j++] << "(fractional)   ";
      read_VASPIN(lines,238,"238,239","        "," "," ","Co$_3$V  ",VERB);// 238,239  hp24   P-6m2          Co3V    #187 130340.cif
      read_VASPIN(lines,242,"242,243","            "," "," ","NbPd$_3$         ",VERB);// 242,243  oP24   PmnmO2           NbPd3   #59  130379.cif
      read_VASPIN(lines,245,"244,245","D0$_{24}$   "," "," ","Ni$_3$Ti         ",VERB);// 244,245         P63/mmc   D0_24  Ni3Ti   #194 AL5957.cif
      //    read_VASPIN(lines,254,"254,255","            "," "," ","W$_5$Si$_3$      ",VERB);// 254,255  tI32   I4/mcm           W5Si3   #140 AL2650.cif  
      read_VASPIN(lines,256,"256    ","B27         "," "," ","BFe              ",VERB);//  oP8	  Pnma		B27	BFe	  #62	AL1744.cif
      for(i=1;i<j;i++) lines[i] << "\t \\\\ ";
      cout << "{\\begin{center}" << endl;
      cout << "{\\footnotesize" << endl;
      cout << "\\begin{tabular}{||c||c|c|c|c||}\\hline\\hline   " << endl;
      for(i=1;i<j;i++) {
        string aus=lines[i].str();
        if(i<=i1 || (i>i1 && (((int) aus.find(")"))>0)) || ((((int) aus.find("posi"))>0)) ) {
          if(i!=i5 || (int) aus.find("A")>0) {
            if((i<i6 || i>i7) && (int) aus.find("aren")<0 && (int) aus.find("ositions")<0 && (int) aus.find("raction")<0) {
              cout << lines[i].str();
              if(i<=i2 || i==i1 || i==i3 || i==i4 ) cout << "\\hline";
              cout << endl;}}}}
      //     cout << "\\hline \\hline" << endl << "\\end{tabular}}" << endl << "\\end{center}}" << endl;
      cout << "\\hline \\hline \\hline \\hline " << endl;
    }
    { // EXT4b TABLE
      cerr << "EXT4b TABLE" << endl; // EXT4b TABLE
      ostringstream lines[MAXLINESOSTREAM];
      j=1;
      lines[j++] << "Strukturbericht";lines[j++] << "Prototype      ";
      lines[j++] << "Label          ";lines[j++] << "Formula        ";
      lines[j++] << "Parent lattice ";lines[j++] << "C$_b$          ";
      lines[j++] << "Space Group    ";lines[j++] << "Pearson Symbol ";
      lines[j++] << "Lattice type   ";lines[j++] << "$(a,b,c)$ (\\AA)";lines[j++] << "$(\\alpha,\\beta,\\gamma)$ degrees";
      lines[j++] << "Superlattice   ";
      lines[j++]<<"Primitive vectors";lines[j++]<<"{\\bf a$_1$}/a";lines[j++]<<"{\\bf a$_2$}/a";lines[j++]<<"{\\bf a$_3$}/a";lines[j++]<<"(cartesian)";
      lines[j++]<<"Atomic positions";
      for(i=1;i<=AAmax;i++) lines[j++] <<"{\\bf A$"<<i<<"$}     ";
      for(i=1;i<=BBmax;i++) lines[j++] <<"{\\bf B$"<<i<<"$}     ";
      lines[j++] << "(fractional)   ";
      read_VASPIN(lines,249,"248,249","C$_c$       "," "," ","ThSi$_2$         ",VERB);// 248,249  tI12   I41/amd   C_c    ThSi2   #141 AL2272.cif
      read_VASPIN(lines,260,"259,260","C33         "," "," ","Bi$_2$Te$_3$     ",VERB);// 259,260  hR5    R-3m      C33    Bi2Te3  #166 120146.cif  
      read_VASPIN(lines,270,"269,270","C6          "," "," ","CdI$_2$          ",VERB);// 269,270  hP3    P-3m1     C6     CdI2    #164 AL4899.cif
      read_VASPIN(lines,257,"257,258","C18         "," "," ","FeS$_2$          ",VERB);// 257,258  oP6    Pnnm      C18    FeS2    #58  AL1657.cif
      for(i=1;i<j;i++) lines[i] << "\t \\\\ ";
      // cout << "{\\begin{center}" << endl;
      // cout << "{\\footnotesize" << endl;
      // cout << "\\begin{tabular}{||c||c|c|c|c||}\\hline\\hline   " << endl;
      for(i=1;i<j;i++) {
        string aus=lines[i].str();
        if(i<=i1 || (i>i1 && (((int) aus.find(")"))>0)) || ((((int) aus.find("posi"))>0)) ) {
          if(i!=i5 || (int) aus.find("A")>0) {
            if((i<i6 || i>i7) && (int) aus.find("aren")<0 && (int) aus.find("ositions")<0 && (int) aus.find("raction")<0) {
              cout << lines[i].str();
              if(i<=i2 || i==i1 || i==i3 || i==i4 ) cout << "\\hline";
              cout << endl;}}}}
      cout << "\\hline\\hline" << endl << "\\end{tabular}}" << endl << "\\end{center}}" << endl;
    }
  }
  cout << "\\newpage" << endl;  
  if(_EXT) {
    { // EXT5a TABLE
      cerr << "EXT5a TABLE" << endl;// EXT5a TABLE
      // cout << "\\subsection{EXT9 superstructures} " << endl;
      ostringstream lines[MAXLINESOSTREAM];
      j=1;
      lines[j++] << "Strukturbericht";lines[j++] << "Prototype      ";
      lines[j++] << "Label          ";lines[j++] << "Formula        ";
      lines[j++] << "Parent lattice ";lines[j++] << "C$_b$          ";
      lines[j++] << "Space Group    ";lines[j++] << "Pearson Symbol ";
      lines[j++] << "Lattice type   ";lines[j++] << "$(a,b,c)$ (\\AA)";lines[j++] << "$(\\alpha,\\beta,\\gamma)$ degrees";
      lines[j++] << "Superlattice   ";
      lines[j++]<<"Primitive vectors";lines[j++]<<"{\\bf a$_1$}/a";lines[j++]<<"{\\bf a$_2$}/a";lines[j++]<<"{\\bf a$_3$}/a";lines[j++]<<"(cartesian)";
      lines[j++]<<"Atomic positions";
      for(i=1;i<=AAmax;i++) lines[j++] <<"{\\bf A$"<<i<<"$}     ";
      for(i=1;i<=BBmax;i++) lines[j++] <<"{\\bf B$"<<i<<"$}     ";
      lines[j++] << "(fractional)   ";
      read_VASPIN(lines,252,"252,253","C15$_b$     "," "," ","AuBe$_5$         ",VERB);// 252,253  cF24   F-43m     C15_b  AuBe5   #216 AL3669.cif  
      // read_VASPIN(lines,261,"261,262","            "," "," ","Ti$_2$Ni         ",VERB);// 261,262  cF96   Fd-3m            Ti2Ni   #227 IMPOSSIBLE
      read_VASPIN(lines,263,"263,264","            "," "," ","Ti$_3$Cu$_4$     ",VERB);// 263,264  tI14   I4/mmm           Ti3Cu4  #139 120654.cif
      // read_VASPIN(lines,265,"265,266","            "," "," ","Ti$_2$Cu$_3$     ",VERB);// 265,266  tP10                    Ti2Cu3  #129 IMPOSSIBLE
      read_VASPIN(lines,267,"267,268","C2          "," "," ","FeS$_2$          ",VERB);// 267,268  oP6    Pnnm      C2     FeS2    #56  AL1657.cif
      read_VASPIN(lines,252,"252,253","C15$_b$     "," "," ","AuBe$_5$         ",VERB);// 252,253  cF24   F-43m     C15_b  AuBe5   #216 AL3669.cif  
      for(i=1;i<j;i++) lines[i] << "\t \\\\ ";
      cout << "{\\begin{center}" << endl;
      cout << "{\\footnotesize" << endl;
      cout << "\\begin{tabular}{||c||c|c|c|c|c||}\\hline\\hline   " << endl;
      for(i=1;i<j;i++) {
        string aus=lines[i].str();
        if(i<=i1 || (i>i1 && (((int) aus.find(")"))>0)) || ((((int) aus.find("posi"))>0)) ) {
          if(i!=i5 || (int) aus.find("A")>0) {
            if((i<i6 || i>i7) && (int) aus.find("aren")<0) {
              cout << lines[i].str();
              if(i<=i2 || i==i1 || i==i3 || i==i4 ) cout << "\\hline";
              cout << endl;}}}}
      //     cout << "\\hline \\hline" << endl << "\\end{tabular}}" << endl << "\\end{center}}" << endl;
      cout << "\\hline \\hline \\hline \\hline " << endl;
    }
    { // EXT5b TABLE
      cerr << "EXT5b TABLE" << endl; // EXT5b TABLE
      ostringstream lines[MAXLINESOSTREAM];
      j=1;
      lines[j++] << "Strukturbericht";lines[j++] << "Prototype      ";
      lines[j++] << "Label          ";lines[j++] << "Formula        ";
      lines[j++] << "Parent lattice ";lines[j++] << "C$_b$          ";
      lines[j++] << "Space Group    ";lines[j++] << "Pearson Symbol ";
      lines[j++] << "Lattice type   ";lines[j++] << "$(a,b,c)$ (\\AA)";lines[j++] << "$(\\alpha,\\beta,\\gamma)$ degrees";
      lines[j++] << "Superlattice   ";
      lines[j++]<<"Primitive vectors";lines[j++]<<"{\\bf a$_1$}/a";lines[j++]<<"{\\bf a$_2$}/a";lines[j++]<<"{\\bf a$_3$}/a";lines[j++]<<"(cartesian)";
      lines[j++]<<"Atomic positions";
      for(i=1;i<=AAmax;i++) lines[j++] <<"{\\bf A$"<<i<<"$}     ";
      for(i=1;i<=BBmax;i++) lines[j++] <<"{\\bf B$"<<i<<"$}     ";
      lines[j++] << "(fractional)   ";
      read_VASPIN(lines,271,"271,272","            "," "," ","YCd$_3$          ",VERB);// 271,272  oC16   Cmcm             YCd3    #63  18260.cif
      read_VASPIN(lines,274,"273,274","B8$_2$      "," "," ","Ni$_2$In         ",VERB);// 273,274  hP6    P63/mmc   B8_2   Ni2In   #194 28249.cif
      read_VASPIN(lines,276,"275,276","            "," "," ","Ni$_2$Si         ",VERB);// 275,276  hP6    P63/m            Ni2Si   #176  AL5234.cif
      read_VASPIN(lines,278,"277,278","D0$_\\alpha$"," "," ","$\\beta$-Cu$_3$Ti",VERB);// 277,278  oP8    PmmnO1  D0_alpha betaCu3Ti #59 AL1729.cif
      for(i=1;i<j;i++) lines[i] << "\t \\\\ ";
      // cout << "{\\begin{center}" << endl;
      // cout << "{\\footnotesize" << endl;
      // cout << "\\begin{tabular}{||c||c|c|c|c||}\\hline\\hline   " << endl;
      for(i=1;i<j;i++) {
        string aus=lines[i].str();
        if(i<=i1 || (i>i1 && (((int) aus.find(")"))>0)) || ((((int) aus.find("posi"))>0)) ) {
          if(i!=i5 || (int) aus.find("A")>0) {
            if((i<i6 || i>i7) && (int) aus.find("aren")<0) {
              cout << lines[i].str();
              if(i<=i2 || i==i1 || i==i3 || i==i4 ) cout << "\\hline";
              cout << endl;}}}}
      cout << "\\hline\\hline" << endl << "\\end{tabular}}" << endl << "\\end{center}}" << endl;
    }
  }

  if(_EXT) {
    { // EXT6a TABLE
      cerr << "EXT6a TABLE" << endl;// EXT5a TABLE
      // cout << "\\subsection{EXT9 superstructures} " << endl;
      ostringstream lines[MAXLINESOSTREAM];
      j=1;
      lines[j++] << "Strukturbericht";lines[j++] << "Prototype      ";
      lines[j++] << "Label          ";lines[j++] << "Formula        ";
      lines[j++] << "Parent lattice ";lines[j++] << "C$_b$          ";
      lines[j++] << "Space Group    ";lines[j++] << "Pearson Symbol ";
      lines[j++] << "Lattice type   ";lines[j++] << "$(a,b,c)$ (\\AA)";lines[j++] << "$(\\alpha,\\beta,\\gamma)$ degrees";
      lines[j++] << "Superlattice   ";
      lines[j++]<<"Primitive vectors";lines[j++]<<"{\\bf a$_1$}/a";lines[j++]<<"{\\bf a$_2$}/a";lines[j++]<<"{\\bf a$_3$}/a";lines[j++]<<"(cartesian)";
      lines[j++]<<"Atomic positions";
      for(i=1;i<=AAmax;i++) lines[j++] <<"{\\bf A$"<<i<<"$}     ";
      for(i=1;i<=BBmax;i++) lines[j++] <<"{\\bf B$"<<i<<"$}     ";
      lines[j++] << "(fractional)   ";
      read_VASPIN(lines,280,"279,280","D0$_{23}$"," "," ","AL$_3$Zr",VERB);// 279,280 tI16  I4/mmm	D0_23	Al3Zr	  #139	NAVY DATABASE
      for(i=1;i<j;i++) lines[i] << "\t \\\\ ";
      cout << "{\\begin{center}" << endl;
      cout << "{\\footnotesize" << endl;
      cout << "\\begin{tabular}{||c||c|c|c|c|c||}\\hline\\hline   " << endl;
      for(i=1;i<j;i++) {
        string aus=lines[i].str();
        if(i<=i1 || (i>i1 && (((int) aus.find(")"))>0)) || ((((int) aus.find("posi"))>0)) ) {
          if(i!=i5 || (int) aus.find("A")>0) {
            if((i<i6 || i>i7) && (int) aus.find("aren")<0) {
              cout << lines[i].str();
              if(i<=i2 || i==i1 || i==i3 || i==i4 ) cout << "\\hline";
              cout << endl;}}}}
      //   cout << "\\hline \\hline \\hline \\hline " << endl;
      //   }
      //   { // EXT6b TABLE
      //   cerr << "EXT6b TABLE" << endl; // EXT5b TABLE
      //   ostringstream lines[MAXLINESOSTREAM];
      //   j=1;
      //   lines[j++] << "Strukturbericht";lines[j++] << "Prototype      ";
      //   lines[j++] << "Label          ";lines[j++] << "Formula        ";
      //   lines[j++] << "Parent lattice ";lines[j++] << "C$_b$          ";
      //   lines[j++] << "Space Group    ";lines[j++] << "Pearson Symbol ";
      //   lines[j++] << "Lattice type   ";lines[j++] << "$(a,b,c)$ (\\AA)";lines[j++] << "$(\\alpha,\\beta,\\gamma)$ degrees";
      //   lines[j++] << "Superlattice   ";
      //   lines[j++]<<"Primitive vectors";lines[j++]<<"{\\bf a$_1$}/a";lines[j++]<<"{\\bf a$_2$}/a";lines[j++]<<"{\\bf a$_3$}/a";lines[j++]<<"(cartesian)";
      //   lines[j++]<<"Atomic positions";
      //   for(i=1;i<=AAmax;i++) lines[j++] <<"{\\bf A$"<<i<<"$}     ";
      //   for(i=1;i<=BBmax;i++) lines[j++] <<"{\\bf B$"<<i<<"$}     ";
      //   lines[j++] << "(fractional)   ";
      //   for(i=1;i<j;i++) lines[i] << "\t \\\\ ";
      //   for(i=1;i<j;i++) {
      //   string aus=lines[i].str();
      //   if(i<=i1 || (i>i1 && (((int) aus.find(")"))>0)) || ((((int) aus.find("posi"))>0)) ) {
      //   if(i!=i5 || (int) aus.find("A")>0) {
      //   if((i<i6 || i>i7) && (int) aus.find("aren")<0) {
      //   cout << lines[i].str();
      //   if(i<=i2 || i==i1 || i==i3 || i==i4 ) cout << "\\hline";
      //   cout << endl;}}}}
      cout << "\\hline\\hline" << endl << "\\end{tabular}}" << endl << "\\end{center}}" << endl;
  }
  }

  cerr << " *********************************************************************** END " << endl;
  cerr << " *********************************************************************** END " << endl;
  cerr << " *********************************************************************** END " << endl;
  cerr << " *********************************************************************** END " << endl;
  LatexSideMarginThesisOFF();
  //cout << "\\end{document}" << endl;  
}

void ZungerConversionChart(void) {
  cout << " \\subsection{Conversion chart for FCC/BCC/HCP superstructures}                   " << endl;
  cout << "{\\large                                                                       " << endl;
  cout << "  \\begin{center}                                                              " << endl;
  cout << "    DEFINITION:                                                                " << endl;
  cout << "    \\begin{itemize}                                                           " << endl;
  cout << "    \\item{}XXX$_{SSS}^{[DIR]}$                                                " << endl;
  cout << "    \\item{}XXX=parent lattice                                                 " << endl;
  cout << "    \\item{}SSS=stacking                                                       " << endl;
  cout << "    \\item{}DIR=stacking direction                                             " << endl;
  cout << "    \\end{itemize}                                                             " << endl;
  cout << "    \\begin{tabular}{|c|c|c|}\\hline                                           " << endl;
  cout << "      {\\it unusual} name   & plain formalism (s)         & label \\\\ \\hline " << endl;
  cout << "      $\\beta1/\\beta2$      & FCC$_{AB2}^{[100]}$         & 5,6 \\\\ \\hline  " << endl;
  cout << "      $\\alpha1/\\alpha2$    & FCC$_{AB2}^{[111]}$         & 9,10 \\\\ \\hline " << endl;
  cout << "      Z1/Z3                & FCC$_{AB3}^{[001]}$         & 13,15 \\\\ \\hline  " << endl;
  cout << "      Z2                   & FCC$_{A2B2}^{[001]}$        & 14 \\\\ \\hline     " << endl;
  cout << "      W2                   & FCC$_{A2B2}^{[311]}$        & 17 \\\\ \\hline     " << endl;
  cout << "      Y1/Y3                & FCC$_{AB3}^{[011]}$         & 19,21 \\\\ \\hline  " << endl;
  cout << "      Y2                   & FCC$_{A2B2}^{[011]}$        & 20 \\\\ \\hline     " << endl;
  cout << "      CH or 40             & FCC$_{A2B2}^{[201]}$        & 23 \\\\ \\hline     " << endl;
  cout << "      V1/V3                & FCC$_{AB3}^{[111]}$         & 27,29 \\\\ \\hline  " << endl;
  cout << "      V2                   & FCC$_{A2B2}^{[111]}$        & 28 \\\\ \\hline     " << endl;
  cout << "    \\end{tabular}                                                             " << endl;
  cout << "  \\end{center}                                                                " << endl;
  cout << "  }                                                                            " << endl;
  cout << "                                                                               " << endl;
  cout << "  \\newpage                                                                    " << endl;
}

// **************************************************************************
// *                                                                        *
// *   VASPIN STEFANO CURTAROLO MASSACHUSETTS INSTITUTE OF TECHNOLOGY 2002  *
// *                                                                        *
// **************************************************************************

