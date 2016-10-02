set terminal wxt 1

#unset key


highc10="component_distribution.dat"
 lowc10="component_distribution.dat"


set xrange[0:12]
set yrange[0:1]
#set key center
plot  highc10  using 1:2  w lp lt 7 ps 0.75 title "PE",\
      highc10  using 1:3  w lp lt 8 ps 0.75 title "drug",\
#      highc10  using 1:4  w lp lt 9 ps 0.75 title "PEG(SPB1PEG)",\
#	   lowc10  using 1:2  w lp lt 5 ps 0.75 title "main chain(SPB2PEG)",\
#       lowc10  using 1:3  w lp lt 2 ps 0.75 title "linker(SPB2PEG)",\
#       lowc10  using 1:4  w lp lt 3 ps 0.75 title "PEG(SPB2PEG)",\




#set print "fit_parameters.txt"
#print a1
#print b1
#print c1
#unset print

#-----------------------------------------------------------------------------------------------------
reset
#set terminal png  lw 2  font "Arial, 15"   size 600,400
#set output "m45c20.png"

set term  postscript eps color   linewidth 1 "Helvetica" 20  enhanced
set output "potential.eps"

#unset key


highc10="potential.dat"
 lowc10="potential.dat"

set xlabel "{/Arial=25 time (step) }"   rotate by  0  offset 0, 0
set ylabel "{/Arial=25 energy}"  rotate by 90 offset 0
#set xrange[0:10]
#set yrange[-0.05:1]
set key at screen 0.95,0.5  right # center
plot  highc10  using 1:3  every ::1 w lp lt -1 pt 7  ps 1.0 linecolor rgb  "red"   title "N",\
      highc10  using 1:4  every ::1 w lp lt -1 pt 7  ps 1.0 linecolor rgb  "blue" title "VDW pairwise_energy",\
      highc10  using 1:5  every ::1 w lp lt -1 pt 7  ps 1.0 linecolor rgb  "green" title "Coulombic pairwise_energy",\
      highc10  using 1:6  every ::1 w lp lt -1 pt 7  ps 1.0 linecolor rgb  "orange" title "bond_energy",\
#	   lowc10  using 1:2  w lp lt  0 pt 65 ps 1.5 linecolor rgb  "red"   title "main chain(SPB2PEG)",\
#       lowc10  using 1:3  w lp lt  0 pt 65 ps 1.5 linecolor rgb  "cyan" title "linker(SPB2PEG)",\
#       lowc10  using 1:4  w lp lt  0 pt 65 ps 1.0 linecolor rgb  "green" title "PEG(SPB2PEG)",\




#pause -1




