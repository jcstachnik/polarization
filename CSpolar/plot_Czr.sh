#!/bin/bash

ARGS=1

if [ $# -ne "$ARGS" ]
then
	echo "Usage: `basename $0` polar_file"
	echo " polar_file = output file from CSpolar.pl"
	exit 
fi

polar=$1 # output file from CSpolar.pl
seaz=100
ts=-20
posttime=600
maxbaz=360
outps=czr_curves.ps


   # Calculate correlation values, get Azimuth at max
    #awk '{print $4/($3*$5),$2}' $polar > junk.txt # Szr/Szz*Szr
    awk '{print $4/sqrt($3*$5),$2}' $polar > Czr.txt
    mm=`minmax -I.1/1 Czr.txt`
    a=(`minmax -Eh Czr.txt`) # Azimuth at max correlation value
    S=`echo ${a[1]} | awk '{printf "%d",$1}'` # azimuth 
    #resid_S=`echo $seaz $S | awk '{printf "%.1f",$2-$1}'` # residual obs-pred
    resid_S=`echo $seaz $S | awk '{printf "%d",$1-$2}'` # residual=seaz-observed
    # zero crossing from 0 az end, not really crossing, just absolute value minimum- hack
    Z=`minmax -EL Czr.txt | awk '{printf "%d", $2}'`
    # This reverses the file, zero crossing from 360 az end
    Zb=`awk '{ a[NR]=$0 } END { for(i=NR; i; --i) print a[i] } ' Czr.txt | minmax -EL | awk '{printf "%d", $2}'`
    # ---- Figure out where 90 degrees from zero cross is appropriate
    if [ "$S" -lt "$Zb" ]
    then
    if [ "$S" -gt "$Z" ]
    then
        zero_baz=`echo $Z | awk '{print $1+90}'`
    else
        if [ "$Zb" -gt "270" ]
        then
        zero_baz=`echo $Zb | awk '{print ($1+90)-360}'`
        else
        zero_baz=`echo $Zb | awk '{print $1+90}'`
        fi
    fi
    else
        zero_baz=`echo $Z | awk '{print $1+270}'`
    fi
    echo "zero_baz: $zero_baz"
    # ----
    awk '{print $4/($3),$2}' $polar > Czr_star.txt # Szr/Szz
    mm1=`minmax -I.1/1 Czr_star.txt`
    m1=`echo $mm1 | awk -F/ '{print $2}'`
    a1=(`minmax -Eh Czr_star.txt`) # Azimuth at max of Szr/Szz
    S1=`echo ${a1[1]} | awk '{printf "%d",$1}'` # azimuth 1
    #resid_S1=`echo $seaz $S1 | awk '{printf "%.1f",$2-$1}'` # residual obs-pred
    resid_S1=`echo $seaz $S1 | awk '{printf "%d",$1-$2}'` # residual=seaz-observed

	echo "C(star)zr azimuth= $S1" echo "Czr azimuth= $S"
    ###################

   # Plot the Correlation curves
    scl=1i/4i
    # Plot a zero line
    ft="--LABEL_FONT_SIZE=12p"
    #psxy -JX$scl -R-1/1/0/$maxbaz -Bf.5a1:"S@-zr@-/@~\326@~(S@-zz@-S@-rr@-)":/a30f15g15SE -m -W1to -O -X4.2 -K $ft <<EOF >> $outps 
    psxy -JX$scl -R-1/1/0/$maxbaz -Bf.5a1:"C@-zr@-":/a30f15g15SE -m -W1to -X4.2 -K $ft <<EOF > $outps
0 0
0 $maxbaz
EOF
    # Plot Szr/sqrt(Szz*Srr)
    psxy Czr.txt -JX$scl -R-1/1/0/$maxbaz -m -W5ta -O -K >> $outps
    # Plot the answer
    psxy -JX -R -m -N -W3/0/0/0ta -O -K <<EOF >> $outps
-1 $S
1 $S
EOF
    psxy -JX -R -m -N -Sc0.15 -W3 -O -K <<EOF >> $outps
0 $Z
0 $Zb
EOF
    # Plot the second answer Szr/Szz
    psxy -JX$scl -R -m -N -W3/0/0/0 -O -K <<EOF >> $outps
1 $S1
-5 $S1
EOF
    # Plot Szr/Szz, different -R
    #psxy junk1.txt -JX$scl $mm1 -B5:"S@-zr@-/S@-zz@-":/10N -m -W2/255/0/0 -O -K $ft >> $outps  
    psxy Czr_star.txt -JX$scl $mm1 -B5:"C@+*@+@-zr@-":/10n -m -W4/0/0/0 -O -K $ft >> $outps
    echo "0 380 12 0 1 CB C@+*@+@-zr@-" |\
    pstext -R -JX$scl -G0/0/0 -N -O  >> $outps

