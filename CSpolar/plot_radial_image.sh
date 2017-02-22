#!/bin/bash

# Make sure to use the -w option in CSpolar.pl
f=radial_2.txt
mm=(`minmax -C $f`)
ts=${mm[0]}
posttime=${mm[1]}
maxbaz=${mm[3]}
outps=radial_image.ps
xtick=f25a100g100

bnd=$ts/$posttime/0/$maxbaz
dx=`awk '(NR ==1){a=$1; getline; b=$1; print b-a}' $f`
xyz2grd $f -Gradial.grd -R$bnd -I$dx/1
#grdmath radial.grd DUP MUL = rad2.grd
grd=radial.grd
gg=(`grdinfo -C $grd`)
grd2cpt $grd -Cgray -I -S0/${gg[6]}/0.05 > radial.cpt
grdimage $grd -JX4i/4i -Cradial.cpt -R$bnd -B${xtick}:"Time (sec)":/f10a30:"Azimuth (\260)":WSne > $outps

