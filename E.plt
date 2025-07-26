unset multiplot
reset
set encoding utf8
set term pngcairo background "#ffffff" font "Arial,21" size 1250,800 #fontscale 0.7 #dashlength 2
set output      'trimerE2.png'

set style line 1 lt 7 lw 1 dt '-' lc rgb "red"
set style line 2 lt 7 lw 1 lc rgb "cyan"
set style line 3 lt 7 lw 2 lc rgb "blue"

#set key top right box maxrows 4
set samples 10000
set bmargin 0.
set lmargin 2.
set rmargin 0.
set tmargin 0.

set multiplot
  set origin 0.12,0.15   #set origin 0.12,0.15
  set size 0.85,0.83     #set size 0.85,0.80
    set xlabel  "ib(blok) / 1000 koraka"
    set ylabel  "E_{trimer} /mK"
    set label 10 "N_w=500,Ns=500,Nb=1000,Rij=9.1,k=0.194" boxed  at graph -0.00, graph 0.95 #textcolor ls 1
    plot [:][:]\
         'E_RoVo_R_9.100_k_0.194.dat' u ($1):($2*1.0) w l ls 2 ti '<ES> tijekom jednog bloka',\
         'E_RoVo_R_9.100_k_0.194.dat' u ($1):($3*1.0) w l ls 3 ti 'Ukupni <ES> nakon ib blokova'

unset multiplot
unset output
reset
set term GNUTERM
