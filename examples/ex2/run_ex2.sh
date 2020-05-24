#!/bin/bash

EXEMDFF=mdff20.x

echo "EXAMPLE 2 : calculation of total energy of LJ clusters"
echo "reference structures and energies from :"
echo "The Cambridge Cluster Database"
echo "http://www-wales.ch.cam.ac.uk/CCD.html" 
echo ""

echo "N (atoms)      mdff.x       CCD          diff"
for (( cluster=3;cluster<150;cluster++)) 
do

	echo "$cluster"               > tmp.file
	echo "CLUSTER_LJ"            >> tmp.file
	echo "100.0   0.0   0.0"     >> tmp.file
	echo "  0.0 100.0   0.0"     >> tmp.file
	echo "  0.0   0.0 100.0"     >> tmp.file
	echo "1"                     >> tmp.file
	echo "A"                     >> tmp.file
	echo $cluster                >> tmp.file
	echo "Cartesian"             >> tmp.file
	awk '{print "A",$1,$2,$3}' config/$cluster >> tmp.file 
	mv tmp.file POSFF


	$EXEMDFF config/control.s > stdout
        mdffres=$(grep "Etot" stdout | awk '{print $3}'| sed -e 's/[eE]+*/\*10\^/'|bc -l)
        CCDres=$(grep " $cluster " config/REFERENCE |awk '{print $2}')
        diff=$(echo "(${mdffres})-(${CCDres})"|bc -l)
        printf "%-8d %15.6f %15.6f %15.6f\n" ${cluster} ${mdffres} ${CCDres} ${diff}

done

echo " "
echo "Lowest energy icosahedral minima at sizes with non-icosahedral global minima. " 
echo " " 
echo "N (atoms)      mdff.x       CCD          diff"
for cluster in 38 75 76 77 98 102  103  104
do
	echo "$cluster" > tmp.file
        echo "CLUSTER_LJ" >> tmp.file
	echo "100.0   0.0   0.0 " >> tmp.file
	echo "  0.0 100.0   0.0 " >> tmp.file
	echo "  0.0   0.0 100.0 " >> tmp.file
	echo "1" >> tmp.file
        echo "A" >> tmp.file
        echo $cluster >> tmp.file
	echo "Cartesian"             >> tmp.file
        awk '{print "A",$1,$2,$3}' config/$((cluster))i >> tmp.file
        mv tmp.file POSFF

        $EXEMDFF config/control.s > stdout
        mdffres=$(grep "Etot" stdout | awk '{print $3}'| sed -e 's/[eE]+*/\*10\^/'|bc -l)
        CCDres=$(grep " ${cluster}i " config/REFERENCE |awk '{print $2}')
        diff=$(echo "(${mdffres})-(${CCDres})"|bc -l)
        printf "%-8d %15.6f %15.6f %15.6f\n" ${cluster} ${mdffres} ${CCDres} ${diff}
done


#rm CONTFF  OSZIFF  POSFF  stdout  tmp.file  TRAJFF
