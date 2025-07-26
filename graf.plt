
reset
set encoding utf8
set term pngcairo background "#ffffff" font "Arial,22" size 1250,800 #fontscale 0.7 #dashlength 2
set output      'valnaa.png'

set style line 1 lt 7 lw 1 dt '-' lc rgb "red"
set style line 2 lt 7 lw 1 lc rgb "cyan"
set style line 3 lt 7 lw 2 lc rgb "blue"

set tics
show tics
set samples 10000
set bmargin 0.
set lmargin 2.
set rmargin 0.
set tmargin 0.
n(x)=(x<7.53 && x>0) 
g(x)= x<7.535 ? -0.56 : 0
Rij=9.1
kij=0.193
f(x) = (x>0 && x<Rij) ? sin(kij*x)/x : x>Rij ? (sin(kij*Rij)/x)*exp((kij*(x-Rij))/tan(kij*Rij)) : 0
  set origin 0.12,0.15   #set origin 0.12,0.1
  set size 0.85,0.83     #set size 0.85,0.80
    set xlabel  "r / A"
    set ylabel  "V_{jama} / K"
    set label 10 "1.područje    2.područje "  at graph 0.00, graph 0.95 #textcolor ls 1
    set style rect fc lt -1 fs solid 0.15 noborder
	set obj rect from 0, graph 0 to 7.53, graph 1 fc rgbcolor "red"
	set obj rect from 7.53, graph 0 to 100, graph 1 fc rgbcolor "yellow"
    plot [0:50][-0.6:0.4] f(x) ti "valna" ,g(x) ti "prav.potencijal" w l ls 3


unset output
reset
set term GNUTERM
