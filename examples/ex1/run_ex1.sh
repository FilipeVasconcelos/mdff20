#!/bin/bash
# files :
#                 config/control.s  (MDFF)
#                 config/POSFF      (MDFF)
#                 config/CONFIG     (DLPOLY)
#                 config/CONTROL    (DLPOLY)
#                 config/FIELD      (DLPOLY)
#
# ====================================================
# user settings

EXEDLPOLY=/home/filipe/dev/dl_poly/execute/DLPOLY.X
EXEMDFF=/home/filipe/dev/mdff20/src/mdff20.x

mdff=true
do_dlpoly=true
# =====================================================

sep="=================================================="
echo $sep
echo "#Example 1: LJ (argon) fcc structure at low temperature"
echo "The configuration is readed in POSFF"
echo "more info in control.s and stdout"

if $mdff; then
    rm -rf mdff/
    mkdir mdff
    cd mdff
    echo 
    echo "# mdff calculation requested"
    cp ../config/control.s .
    cp ../config/POSFF .
    $EXEMDFF control.s > stdout
    #poszi.py -i OSZIFF -n 
    grep e+ OSZIFF | sed -n 'p;n' | awk '{print $2,$3}' > ene
    cd ..
fi

if $do_dlpoly; then
        echo "# dlpoly calculation requested "
        rm -rf dl_poly
        mkdir dl_poly
        cd dl_poly
        cp ../config/CONTROL .
        cp ../config/FIELD .
        cp ../config/CONFIG .
        $EXEDLPOLY
	/home/filipe/dev/mdff20/bin/read_statis > ene
        cd ..
fi

lr=$(grep "long range correction (energy)" mdff/stdout |  awk '{print $NF}')
echo "long range correction $lr"
cat > plot << eof
#!/usr/bin/gnuplot -persist 
set term x11 
set title "Total energy MD ex1."
set xlabel "time"
set ylabel "E_{tot}"
p 'dl_poly/ene' u 1:2 w l title "DLPOLY", 'mdff/ene' u 1:(\$2+${lr}) w l title "MDFF20"
eof
chmod u+x plot
./plot
