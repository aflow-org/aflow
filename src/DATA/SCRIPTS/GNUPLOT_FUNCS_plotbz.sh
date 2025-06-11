#!/bin/ksh
# A script to plot brillouin zone and kpath.
# It needs aflow with --bzplotdata option
# input files: POSCAR.bands
#
# written: 2010 wahyu@alumni.duke.edu

cdir=./
echo plotbz.sh::DIR=  $cdir
kfile=${cdir}/KPOINTS.bands
posfile=${cdir}/POSCAR.bands
datfile=${cdir}/plotbz.dat
recipfile=${cdir}/recip.dat #reciprocal vectors
kptsfile=${cdir}/kpts.dat #kpts
kpathfile=${cdir}/kpath.dat #kpath
nonkpathfile=${cdir}/nonkpath.dat #irbz non kpath
fbzfile=${cdir}/fbz.dat #front bz
bbzfile=${cdir}/bbz.dat #back bz

aflow --bzplotdatauseKPOINTS=$kfile < $posfile > $datfile

#-----reading datfile-----
#---reciprocal lattice vectors b1 b2 b3---
exec 0< $datfile
$IFS read glatt
$IFS read gpath
$IFS read xrotview zrotview srest
read srest
$IFS read Lb1 Lb2 Lb3 srest
$IFS read Tb1 Tb2 Tb3 srest
read srest
$IFS read b11 b12 b13
$IFS read b21 b22 b23
$IFS read b31 b32 b33
((hb11=b11/2.0)); ((hb12=b12/2.0)); ((hb13=b13/2.0)); 
((hb21=b21/2.0)); ((hb22=b22/2.0)); ((hb23=b23/2.0)); 
((hb31=b31/2.0)); ((hb32=b32/2.0)); ((hb33=b33/2.0)); 
#--length of b1 b2 b3 arrows--
((db11=b11*Lb1)); ((db12=b12*Lb1)); ((db13=b13*Lb1)); 
((db21=b21*Lb2)); ((db22=b22*Lb2)); ((db23=b23*Lb2)); 
((db31=b31*Lb3)); ((db32=b32*Lb3)); ((db33=b33*Lb3)); 
echo $hb11 $hb12 $hb13 $db11 $db12 $db13 > $recipfile
echo $hb21 $hb22 $hb23 $db21 $db22 $db23 >> $recipfile
echo $hb31 $hb32 $hb33 $db31 $db32 $db33 >> $recipfile
((tb11=hb11+b11*Tb1)); ((tb12=hb12+b12*Tb1)); ((tb13=hb13+b13*Tb1));
((tb21=hb21+b21*Tb2)); ((tb22=hb22+b22*Tb2)); ((tb23=hb23+b23*Tb2));
((tb31=hb31+b31*Tb3)); ((tb32=hb32+b32*Tb3)); ((tb33=hb33+b33*Tb3));

#---kpts---
$IFS read Nkpts srest
rm -f $kptsfile
for ((i=0;i<Nkpts;i++))
do
read kptsx[$i] kptsy[$i] kptsz[$i] klab[$i]
echo ${kptsx[$i]} ${kptsy[$i]} ${kptsz[$i]} ${klab[$i]} >> $kptsfile
done

#---kpath---
read srest
read srest
$IFS read N srest
rm -f $kpathfile
for ((i=0;i<N;i++))
do
read srest
echo $srest >> $kpathfile
done
#---irbz non kpath---
$IFS read Nnonkpath srest
rm -f $nonkpathfile
for ((i=0;i<Nnonkpath;i++))
do
read srest
echo $srest >> $nonkpathfile
done
#---front bz---
$IFS read N srest
rm -f $fbzfile
for ((i=0;i<N;i++))
do
read srest
echo $srest >> $fbzfile
done
#---back bz---
$IFS read N srest
rm -f $bbzfile
for ((i=0;i<N;i++))
do
read srest
echo $srest >> $bbzfile
done

#getting the structure name and icsd number
name=$(pwd)
echo $name | sed "s/\//\n/g" | grep '_ICSD_' > wahyutmp
exec 0< wahyutmp
$IFS read name
rm -f wahyutmp

echo "set terminal postscript eps color enhanced \"Helvetica\" 18 "> plotbz.gnu
echo "set output \"bz_$name.eps\" ">> plotbz.gnu
echo "set key off ">> plotbz.gnu
echo "unset border ">> plotbz.gnu
echo "unset xtics ">> plotbz.gnu
echo "unset ytics ">> plotbz.gnu
echo "unset ztics ">> plotbz.gnu
echo "set view equal xyz ">> plotbz.gnu
echo "set view $xrotview,$zrotview ">> plotbz.gnu
echo "set title \"$glatt  path: $gpath\" ">> plotbz.gnu
for ((i=0;i<Nkpts;i++))
do
echo "set label \"${klab[$i]}\" at ${kptsx[i]},${kptsy[i]},${kptsz[i]}" >> plotbz.gnu
done
echo "set label \"b{/Symbol _1}\" at $tb11,$tb12,$tb13 font \"Helvetica,14\" ">> plotbz.gnu
echo "set label \"b{/Symbol _2}\" at $tb21,$tb22,$tb23 font \"Helvetica,14\" ">> plotbz.gnu
echo "set label \"b{/Symbol _3}\" at $tb31,$tb32,$tb33 font \"Helvetica,14\" ">> plotbz.gnu
echo "splot '$recipfile' using 1:2:3:4:5:6 with vectors head filled lt 1 lw 1 lc rgb \"#000000\", \\" >> plotbz.gnu
echo "'$fbzfile' using 1:2:3:4:5:6 with vector nohead lt 1 lw 1 lc rgb \"#000000\", \\" >> plotbz.gnu
echo "'$bbzfile' using 1:2:3:4:5:6 with vector nohead lt 0 lw 1.5 lc rgb \"#000000\", \\" >> plotbz.gnu
echo "'$kpathfile' using 1:2:3:4:5:6 with vector nohead lt 1 lw 1 lc rgb \"#FF0000\", \\" >> plotbz.gnu
if ((Nnonkpath>0))
then
echo "'$nonkpathfile' using 1:2:3:4:5:6 with vector nohead lt 2 lw 1 lc rgb \"#FF0000\", \\" >> plotbz.gnu
fi
echo "'$kptsfile' using 1:2:3 with points pt 7 ps 0.8 lc rgb \"#FF0000\" " >> plotbz.gnu

gnuplot plotbz.gnu
convert -density 200 bz_$name.eps bz_$name.png
rm -f plotbz.gnu $datfile $recipfile $kptsfile $kpathfile $nonkpathfile $fbzfile $bbzfile
