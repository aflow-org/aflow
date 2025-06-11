// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2023           *
// *                                                                         *
// ***************************************************************************
// Written by Cormac Toher
// cormac.toher@duke.edu

#include "aflow.h"
// [OBSOLETE]  #include "aflow_contrib_cormac.h"

bool   RemoveGnuplotScriptct=true;

// put your stuff here

void debint (double& y, double& Deb);

// ***************************************************************************
// pflow::DEBYE
// ***************************************************************************
namespace pflow {
  void DEBYE(const string& options) {
    // Fit Debye temperature to heat capacity calculated using APL
    // Heat capacity is in file "THERMO"
    // Usage: aflow --debye=THERMO 
    //

    bool LDEBUG=(FALSE || XHOST.DEBUG);
    if(LDEBUG) cerr << __AFLOW_FUNC__ << " BEGIN" << endl;
    vector<string> tokens;
    aurostd::string2tokens(options,tokens,",");
    if(tokens.size()!=1) {
      init::ErrorOption(options,__AFLOW_FUNC__,"aflow --debye=THERMO[.bz2|.gz|.xz]");
    }
    string thermofilename="";
    if(tokens.size()>=1) thermofilename=tokens.at(0);

    string cvfilename = "", debyefilename = "", firstline = "", flword = "", line = "", sysname = "";
    ifstream infilethermo;
    ofstream ofilecv, ofiledebye;
    int cvcolumn = 0, tcolumn = 0, ncolumn = 0, itmin = 0, itmax = 0, id = 0, jstart = 0;
    int initcol = 0, npoints = 0;
    stringstream thermofile;
    double rowvar = 0.0, cvt = 0.0, nkb = 0.0, natoms = 0.0, tol = 0.0, tol2 = 0.0, tinit = 0.0, y = 0.0, Deb = 0.0, fy = 0.0, ya = 0.0, fya = 0.0, yb = 0.0, fyb = 0.0;
    double tdmin = 0.0, tdmax = 0.0, rmsmin = 0.0, tdtrial = 0.0, smsq = 0.0, rms = 0.0, tdbest = 0.0, DEB_min = 0.0, DEB_max = 0.0, Tmin = 0.0, Tmax = 0.0;
    bool natomcalcd = false;

    nkb = 1.0;
    tol = 0.00001;
    tinit = 100.0;
    tol2 = 1E-8;

    // Open THERMO file (APL output) to determine number of temperature points
    if(aurostd::FileExist(thermofilename)) {
      if(thermofilename == "THERMO") {
        infilethermo.open("THERMO");
      } else if(thermofilename == "THERMO.bz2") {
        aurostd::execute("bzip2 -dqf " + thermofilename);
        infilethermo.open("THERMO");
      } else if(thermofilename == "THERMO.gz") {
        aurostd::execute("gzip -dqf " + thermofilename);
        infilethermo.open("THERMO");
      } else if(thermofilename == "THERMO.xz") {
        aurostd::execute("xz -dqf " + thermofilename);
        infilethermo.open("THERMO");
      } 

      getline(infilethermo, firstline);
      while(infilethermo.peek() != EOF) {
        getline(infilethermo, line);
        npoints++;
      }   
      infilethermo.close();
    }
    else if(thermofilename == "THERMO" && aurostd::FileExist("THERMO.bz2")) {
      aurostd::execute("bzip2 -dqf THERMO.bz2");
      infilethermo.open("THERMO");
      getline(infilethermo, firstline);
      while(infilethermo.peek() != EOF) {
        getline(infilethermo, line);
        npoints++;
      }    
      infilethermo.close();
    }
    else if(thermofilename == "THERMO" && aurostd::FileExist("THERMO.gz")) {
      aurostd::execute("gzip -dqf THERMO.gz");
      infilethermo.open("THERMO");
      getline(infilethermo, firstline);
      while(infilethermo.peek() != EOF) {
        getline(infilethermo, line);
        npoints++;
      }    
      infilethermo.close();
    }
    else if(thermofilename == "THERMO" && aurostd::FileExist("THERMO.xz")) {
      aurostd::execute("xz -dqf THERMO.xz");
      infilethermo.open("THERMO");
      getline(infilethermo, firstline);
      while(infilethermo.peek() != EOF) {
        getline(infilethermo, line);
        npoints++;
      }    
      infilethermo.close();
    }
    else if((thermofilename == "THERMO.bz2" || thermofilename == "THERMO.gz" || thermofilename == "THERMO.xz") && aurostd::FileExist("THERMO")) {
      infilethermo.open("THERMO");
      getline(infilethermo, firstline);
      while(infilethermo.peek() != EOF) {
        getline(infilethermo, line);
        npoints++;
      }    
      infilethermo.close();
    }
    else {
      cout << "Error: file " << thermofilename << " does not exist in this directory" << endl;
      return;
    }

    infilethermo.open("THERMO");
    // Determine which columns contain temperature and heat capacity data
    for (int i = 0; i < aurostd::CountWordsinString(firstline); i++) {
      infilethermo >> flword;
      if(flword == "#") {
        initcol = -1;
      }
      if(flword == "T(K)") {
        tcolumn = i + initcol;
      }
      if(flword == "Cv(kB/cell)") {
        cvcolumn = i + initcol;
      }
    }

    infilethermo.close();

    thermofile.clear(); thermofile.str(std::string());
    aurostd::file2stringstream("THERMO", thermofile);

    ncolumn = aurostd::CountWordsinString(firstline) + initcol;

    vector<double> Temp(npoints);
    vector<double> Cv(npoints);
    vector<double> ThetaD(npoints);

    getline(thermofile, line);
    // Read in temperature and heat capacity data from THERMO file
    for (int i = 0; i < npoints; i++) {
      for (int j = 0; j < ncolumn; j++) {
        thermofile >> rowvar;
        if(j == tcolumn) {
          Temp[i] = rowvar;
        }
        if(j == cvcolumn) {
          Cv[i] = rowvar;
        }
      }
    }
    // Determines number of atoms in unit cell by reading POSCAR part of _AFLOWIN_ file if it exists
    if(aurostd::FileExist(_AFLOWIN_)) {
      stringstream iafile;
      iafile.clear(); iafile.str(std::string());
      aurostd::file2stringstream(_AFLOWIN_, iafile);    
      //[SD20220520 - OBSOLETE]if(aurostd::substring2bool(iafile.str(),_VASP_POSCAR_MODE_EXPLICIT_START_) && aurostd::substring2bool(iafile.str(),_VASP_POSCAR_MODE_EXPLICIT_STOP_))
      stringstream POSCAR;
      POSCAR.str(aurostd::substring2string(iafile.str(),_VASP_POSCAR_MODE_EXPLICIT_START_,_VASP_POSCAR_MODE_EXPLICIT_STOP_,-1));
      if(!POSCAR.str().empty()) {
        //[SD20220520 - OBSOLETE]POSCAR.clear();   POSCAR.str(std::string());
        //[SD20220520 - OBSOLETE]aurostd::ExtractToStringstreamEXPLICIT(iafile.str(),POSCAR,_VASP_POSCAR_MODE_EXPLICIT_START_,_VASP_POSCAR_MODE_EXPLICIT_STOP_);
        xstructure xstr = xstructure(POSCAR, IOVASP_AUTO);
        int ntypes = xstr.num_each_type.size();
        natoms = 0;
        for(int i=0; i<ntypes; i++) {
          natoms = natoms + xstr.num_each_type.at(i);
        }
        natomcalcd = true;
        sysname = xstr.title;
      }
    }
    // If POSCAR part of _AFLOWIN_ file does not exist, 
    // gives an error and exits
    if(!natomcalcd) {
      cout << "Error: file " << _AFLOWIN_ << " does not exist in this directory, cannot determine number of atoms" << endl;
      return;
    }
    nkb = nkb * natoms;
    // Skip initial temperature point if it is equal to zero
    if(aurostd::abs(Temp[0]) < tol2) {
      jstart = 1;
    }
    else {
      jstart = 0;
    }    
    // For each temperature value, determines the value of the Debye temperature 
    // which gives the corresponding value of the heat capacity as calculated using APL
    // Uses Debye model for heat capacity expression
    // See Aschcroft & Mermin, Solid State Physics, Eqn 23.26
    for (int j = jstart; j < npoints; j++) {
      cvt = Cv[j] / nkb;
      if(j == 1) {	 
        ThetaD[j] = Temp[j] * tinit;
      }
      else {
        ThetaD[j] = ThetaD[j-1];
      }
      y = ThetaD[j] / Temp[j];
      debint(y, Deb);
      fy = 9.0 * Deb - cvt;

      if(aurostd::abs(fy) < tol) {
        //	  Initial trial theta is correct - do nothing
      }
      else if(fy > 0.0) {
        fya = fy;
        ya = y;
        y = y * 2.0;
        debint(y, Deb);
        fy = 9.0 * Deb - cvt;
        if(fy < fya) {
          while (fy > 0.0) {
            y = y * 2.0;
            debint(y, Deb);
            fy = 9.0 * Deb - cvt;
          }
          yb = y;
        }
        else if(fy > fya) {
          while (fy > 0.0) {
            y = y / 2.0;
            debint(y, Deb);
            fy = 9.0 * Deb - cvt;
          }
          yb = y;
        }  
      }
      else if(fy < 0.0) {
        fyb = fy;
        yb = y;
        y = y / 2.0;
        debint(y, Deb);
        fy = 9.0 * Deb - cvt;
        if(fy > fyb) {
          while (fy < 0.0) {
            y = y / 2.0;
            debint(y, Deb);
            fy = 9.0 * Deb - cvt;
          }
          ya = y;
        }
        else if(fy < fyb) {
          while (fy < 0.0) {
            y = y * 2.0;
            debint(y, Deb);
            fy = 9.0 * Deb - cvt;
          }
          yb = y;
        }
      }

      while (aurostd::abs(fy) > tol) {
        y = (ya + yb) / 2.0;
        debint(y, Deb);
        fy = 9.0 * Deb - cvt;
        if(fy > 0.0) {
          ya = y;
        }
        else if(fy < 0.0) {
          yb = y;
        }
      }
      ThetaD[j] = y * Temp[j];

    }

    if(aurostd::abs(Temp[0]) < tol2) {
      tdmin = ThetaD[1];
      tdmax = ThetaD[1];
    }
    else {
      tdmin = ThetaD[1];
      tdmax = ThetaD[1];
    }

    for (int j = 1; j < npoints; j++) {
      if(ThetaD[j] > tdmax) {
        tdmax = ThetaD[j];
      }

      if(ThetaD[j] < tdmin) {
        tdmin = ThetaD[j];
      }
    }
    DEB_min = tdmin - fmod(tdmin,100.0);
    DEB_max = tdmax - fmod(tdmin,100.0) + 100.0;
    Tmin = Temp[0];
    Tmax = Temp[npoints-1];
    tdmin = tdmin - fmod(tdmin, 1.0);
    tdmax = tdmax - fmod(tdmax, 1.0);

    itmin = tdmin;
    itmax = tdmax;
    rmsmin = 1e30;
    // Determines which Debye temperature gives best fit to all heat capacity data in given range
    for (int i = itmin; i <= itmax; i++)
    {
      id = i;
      tdtrial = id;
      smsq = 0.0;
      for (int j = 1; j < npoints; j++) {
        y = tdtrial / Temp[j];
        debint(y, Deb);
        cvt = 9.0 * nkb * Deb;
        smsq = smsq + pow((cvt - Cv[j]), 2);
      }
      rms = sqrt(smsq);
      if(rms < rmsmin) {
        rmsmin = rms;
        tdbest = tdtrial;
      }
    }

    ofiledebye.open("debye_temperature.dat");  
    if(ofiledebye.is_open()) {
      ofiledebye << "# Debye temperature fitted to Cv values obtained from APL" << endl;
      ofiledebye << "# Debye temperature giving best fit to Cv data = " << tdbest << "K " << endl;
      ofiledebye << "# T(K) " << "\t" << "Debye temperature (K) " << endl; 
      for (int i = jstart; i < npoints; i++) {
        ofiledebye << Temp[i] << "\t" << ThetaD[i] << endl;
      }
    }
    else {
      cout << "Error: Unable to open file debye_temperature.dat" <<  endl;
      return;
    } 
    ofiledebye.close();

    string debyedatafile = "debye_data.dat";
    ofstream fdebin;
    fdebin.open(debyedatafile.c_str());  
    if(fdebin.is_open()) {
      if(aurostd::abs(Temp[0]) < tol2) {
        for (int i = 1; i < npoints; i++) {
          fdebin << Temp[i] << "\t" << ThetaD[i] << endl;
        }
      }
      else {
        for (int i = 0; i < npoints; i++) {
          fdebin << Temp[i] << "\t" << ThetaD[i] << endl;
        }		  
      }
    }
    else {
      cout << "Error: Unable to open file debye_data.dat" <<  endl;
      return;
    } 
    fdebin.close();

    //***********************************GENERATING GNUPLOT SCRIPT**************************************************************
    //Writing Gnuplot Script
    string gnuplotscript = "GNUPLOT_debye.gp";
    ofstream fin;
    fin.open(gnuplotscript.c_str());
    fin << "#Generated by AFLOW (Cormac Toher [cormac.toher@duke.edu], 2013, Duke)" << endl;
    fin << "set term postscript eps enhanced color font \"Times-Roman, 40\" size 18, 10.125" << endl;
    fin << "set output " << "\"" << "debye_temperature" <<".eps" << "\"" << endl;
    fin << "set label '" << AFLOWLIB_CONSORTIUM_STRING << "' at screen 0.75, 0.02 font \"Times-Roman, 32\"" << endl;
    fin << endl;

    fin << "#Debye PLOT" << endl;

    fin << "set title 'Debye Temperature (K) "  << sysname << "'"  << endl;
    fin << "set xtics " << endl;
    fin << "set ytics" << endl;
    fin << "set yrange [" << DEB_min << ":" << DEB_max << "]" << endl;
    fin << "set xrange [" << Tmin << ":"  << Tmax << "]" << endl;
    fin << endl;

    fin << "set xlabel 'Temperature (K)' offset graph 0.00" << endl;
    fin << "set ylabel 'Debye Temperature (K)' offset graph 0.00" << endl;
    fin << "set arrow from 0, 0 to graph 1, first 0 nohead lt 3 lw 1.5" << endl;
    fin << "set arrow from 0, 0 to graph 0, first 0 nohead lt 3 lw 1.5" << endl;
    fin << "set key font \"Times-Roman, 40\"" << endl; 	
    fin << endl;
    fin << "plot[][] \\" << endl;
    fin << "\"" << debyedatafile << "\""  << endl;

    fin << endl;

    fin.close();

    //Call gnuplot to plot the Debye temperature	
    aurostd::execute(XHOST.command("gnuplot")+" " + gnuplotscript);
    aurostd::execute(XHOST.command("convert")+"./debye_temperature.eps ./debye_temperature.png");

    //Postprocess
    if(RemoveGnuplotScriptct) { 
      aurostd::execute("rm -f " + debyedatafile);
      aurostd::execute("rm -f " + gnuplotscript); //Delete the gnuplot script
    }

    aurostd::execute("rm -f debye_temperature.eps");

    if(thermofilename == "THERMO.bz2") { aurostd::execute("bzip2 -9fq THERMO"); }
    if(thermofilename == "THERMO.gz") { aurostd::execute("gzip -9fq THERMO"); }
    if(thermofilename == "THERMO.xz") { aurostd::execute("xz -9fq THERMO"); }
    // END
    if(LDEBUG) cerr << __AFLOW_FUNC__ << " END" << endl;
  }
}


void gauleg(double& x1,double& x2,double x[],double w[], int& n) {
  //     Calculate Gauss-Legendre quadrature for integral
  //.....Increase eps if you don't have this floating precision.
  double eps=3.0e-16;
  int i, j, m, id;
  double xm, xl, p1, p2, p3, pp, z, z1;
  bool converged = false;
  //.....The roots are symmetric in the interval, so we only have to find
  //     half of them.
  m=(n+1)/2;
  xm=0.5*(x2+x1);
  xl=0.5*(x2-x1);
  //.....Loop over the desired roots.
  for ( i = 0; i < m; i++) {
    id = i + 1;
    z=cos(PI*(id-0.25)/(n+0.5));
    converged = false;
    while (!converged) {
      p1=1.0;
      p2=0.0;
      //.........Loop up the recurrence relation to get the Legendre polynomial 
      //         evaluated at z.
      for ( j = 1; j <= n; j++) {
        p3=p2;
        p2=p1;
        p1=((2.0*j-1.0)*z*p2-(j-1.0)*p3)/j;
      }
      //.........p1 is now the desired Legendre polynomial. We next compute pp,
      //         derivative , by a standard relation involving p2, the polyn-
      //         omial of one lower order.
      pp=n*(z*p1-p2)/(z*z-1.0);
      z1=z;
      //.........Newton's method.
      z=z1-p1/pp;
      if(aurostd::abs(z-z1) > eps) {
        converged = false;
      }
      else {
        converged = true;
      }
    }
    //.......Scale the root to the desired interval.
    x[i] = xm - xl * z;
    //.......and put in its symmetric counterpart.
    x[n-id] = xm + xl * z;
    //.......compute the weight.
    w[i]=2.0*xl/((1.0-z*z)*pp*pp);
    //.......and its symmetric counterpart.
    w[n-id]=w[i];
  }
  return;
}

double fdebye(double z) {
  return (pow(z, 4) * exp(z)) / (pow((exp(z) - 1), 2));
}

void debint (double& y, double& Deb) {
  //    Evaluate Debye integral
  double eps=1e-12;
  double cero=0.0;
  int maxnl=100;
  double x[maxnl],w[maxnl];
  double debye, debye0, xabs, sum;
  int i, nl;
  //.....error condition controls
  debye=3.0*pi*pi*pi*pi/y/y/y/15.0;
  if(y <= 250) {
    //.....Loop with increasing number of Legendre points.
    debye0=1e30;
    for (nl = 5; nl <= maxnl; nl = nl +5) {
      gauleg (cero,y,x,w,nl);
      sum=0.0;
      for (i = 0; i < nl; i++) {
        sum=sum+w[i]*fdebye(x[i]);
      }
      debye=sum/y/y/y;
      xabs=aurostd::abs(debye-debye0);
      if(xabs < eps) {
        break;
      }
      else {
        debye0=debye;
      }
    }
  }
  Deb = debye;
  return;
}




// ***************************************************************************
