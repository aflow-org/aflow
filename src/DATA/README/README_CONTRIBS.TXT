****************************************************************************************************
*                                                                                                  *
*                          aflow - Automatic-FLOW for materials discovery                          *
*                aflow.org consortium - High-Throughput ab-initio Computing Project                *
*                                                                                                  *
****************************************************************************************************
*                                                                                                  *
*    AFLOW is free software: you can redistribute it and/or modify it under the terms of the       *
*    GNU General Public License as published by the Free Software Foundation, either version 3     *
*    of the License, or (at your option) any later version.                                        *
*                                                                                                  *
*    This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;     *
*    without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.     *
*    See the GNU General Public License for more details.                                          *
*                                                                                                  *
*    You should have received a copy of the GNU General Public License along with this program.    *
*    If not, see <http://www.gnu.org/licenses/>.                                                   *
*                                                                                                  *
****************************************************************************************************
*                                                                                                  *
*     Use of AFLOW software and repositories welcomes references to the following publications:    *
*                                                                                                  *
*  Esters et al.    Comp. Mat. Sci. 216, 111808 (2023) 10.1016/j.commatsci.2022.111808 (AFLOW.org) *
*  Oses et al.      Comp. Mat. Sci. 217, 111889 (2023) 10.1016/j.commatsci.2022.111889 (aflow++)   *
*  Friedrich et al. npj Comput. Mater. 5, 59 (2019)  10.1038/s41524-019-0192-1       (CCE)         *
*  Hicks et al.     Comp. Mat. Sci. 161, S1 (2019)   10.1016/j.commatsci.2018.10.043 (ANRL proto2) *
*  Oses et al.      J. Chem. Inf. Model. (2018)      10.1021/acs.jcim.8b00393        (AFLOW-CHULL) *
*  Gossett et al.   Comp. Mat. Sci. 152, 134 (2018)  10.1016/j.commatsci.2018.03.075 (AFLOW-ML)    *
*  Hicks et al.     Acta Cryst. A74, 184-203 (2018)  10.1107/S2053273318003066       (AFLOW-SYM)   *
*  MBNardelli et al Comp. Mat. Sci. 143, 462 (2018)  10.1016/j.commatsci.2017.11.034 (PAOFLOW)     *
*  Rose et al.      Comp. Mat. Sci. 137, 362 (2017)  10.1016/j.commatsci.2017.04.036 (AFLUX lang)  *
*  Supka et al.     Comp. Mat. Sci. 136, 76 (2017)   10.1016/j.commatsci.2017.03.055 (AFLOWpi)     *
*  Plata et al.     npj Comput. Mater. 3, 45 (2017)  10.1038/s41524-017-0046-7       (AAPL kappa)  *
*  Toher et al.     Phys. Rev.Mater.1, 015401 (2017) 10.1103/PhysRevMaterials.1.015401 (AEL elast) *
*  Mehl et al.      Comp. Mat. Sci. 136, S1 (2017)   10.1016/j.commatsci.2017.01.017 (ANRL proto1) *
*  Calderon et al.  Comp. Mat. Sci. 108A, 233 (2015) 10.1016/j.commatsci.2015.07.019 (standard)    *
*  Toher et al.     Phys. Rev. B 90, 174107 (2014)   10.1103/PhysRevB.90.174107      (AGL Gibbs)   *
*  Taylor et al.    Comp. Mat. Sci. 93, 178 (2014)   10.1016/j.commatsci.2014.05.014 (REST-API)    *
*  Curtarolo et al. Comp. Mat. Sci. 58, 227 (2012)   10.1016/j.commatsci.2012.02.002 (AFLOW.org)   *
*  Curtarolo et al. Comp. Mat. Sci. 58, 218 (2012)   10.1016/j.commatsci.2012.02.005 (aflow++)     *
*                                                                                                  *
****************************************************************************************************
*                                                                                                  *
*                                  aflow/aflow.org - CONTRIBUTORS                                  *
*  2000-2019 Stefano Curtarolo (aflow); 2002-2004 Dane Morgan (convasp); 2007-2011 Wahyu Setyawan  *
*  (--rsm --edos --kband --icsd*); 2008-2011 Roman Chepulskyy (--edos --kband  surfaces);          *
*  2008 Gus Hart (lattice reductions - prototypes); 2009-2011, Ohad Levy (prototypes);             *
*  2009-2010, Michal Jahnatek (APL); 2010-2013 Shidong Wang (cluster expansion); 2010-2013         *
*  Richard Taylor (surfaces, apennsy); 2010-2013 Junkai Xue (prototyper); 2010-2013 Kesong Yang    *
*  (findsym, frozsl, plotband/dos); 2013-2019 Cormac Toher (AGL Debye-Gruneisen, AEL elastic);     *
*  2013-2019 Frisco Rose (API, Aflux); 2013-2018 Pinku Nath (Quasi-harmonic approximation);        *
*  2013-2017 Jose J. Plata (AAPL, thermal cond.); 2014-2023 David Hicks (symmetry, structure       *
*  comparison, prototypes); 2014-2022 Corey Oses (Egap, bader, chull, APL, pocc); 2018-2022 Marco  *
*  Esters (AAPL, thermal cond.); 2016-2019 Denise Ford (GFA); 2018-2020 Rico Friedrich (CCE);      *
*  2021-2023 Simon Divilov (QCA, fitting); 2021-2023 Hagen Eckert (GFA, entry loader, JSON);       *
****************************************************************************************************
*                                 version 3.2.14 - (C) 2003-2023                                   *
****************************************************************************************************

