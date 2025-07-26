unset multiplot
reset
set encoding utf8
set term pngcairo background "#ffffff" font "Arial,21" size 1250,800 #fontscale 0.7 #dashlength 2
set output      'minfin.png'

set style line 1 lt 7 lw 1 dt '-' lc rgb "red"
set style line 2 lt 7 lw 1 lc rgb "cyan"
set style line 3 lt 7 lw 2 lc rgb "blue"

set key top right box maxrows 4
set samples 10000
set bmargin 0.
set lmargin 2.
set rmargin 0.
set tmargin 2.
set palette model RGB ; set title "Rij"
set multiplot
  set origin 0.12,0.15   #set origin 0.12,0.15
  set size 0.85,0.83     #set size 0.85,0.80
    set xlabel  "kij / 1/A"
    set ylabel  "E_{trimer} /mK"
    set label 10 "N_w = 100, Ns=400, Nb = 15 " boxed  at graph 0.03, graph 0.95 #textcolor ls 1
    #plot 'svi_parametri.txt' u 2:3:4  t "<E> i greška" w yerr
         
    plot [0.16:0.225] 'svi_parametri2fin.txt'  every :::::10 u 2:3:1 w lp pt 3 palette t "<E>"
unset multiplot
unset output
reset
set term GNUTERM
