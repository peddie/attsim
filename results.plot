# set term wxt font "Gill Sans,9" linewidth 4 rounded
set term pdfcairo enhanced color font "Gill Sans,9" linewidth 4 rounded
set output "quat.pdf"
# set terminal svg size 1600,1200 fname "Gill Sans" fsize 9 rounded dashed
# set output "ping.svg"

set grid

# Line style for axes
set style line 80 lt rgb "#808080"

# Line style for grid
set style line 81 lt 0  # dashed
set style line 81 lt rgb "#808080"  # grey

set grid back linestyle 81
set border 3 back linestyle 80 # Remove border on top and right.  These
                               # borders are useless and make it harder
			       # to see plotted lines near the border.
			       # Also, put it in grey; no need for so 
			       # much emphasis on a border.
set xtics nomirror
set ytics nomirror

# Line styles: try to pick pleasing colors, rather
# than strictly primary colors or hard-to-see colors
# like gnuplot's default yellow.  Make the lines thick
# so they're easy to see in small plots in papers.
set style line 1 lt rgb "#A00000" lw 2 pt 1
set style line 2 lt rgb "#00A000" lw 2 pt 6
set style line 3 lt rgb "#5060D0" lw 2 pt 2
set style line 4 lt rgb "#F25900" lw 2 pt 9

# set style line 1 lt 1
# set style line 2 lt 1
# set style line 3 lt 1
# set style line 4 lt 1
# set style line 1 lt rgb "#A00000" lw 2 pt 7
# set style line 2 lt rgb "#00A000" lw 2 pt 9
# set style line 3 lt rgb "#5060D0" lw 2 pt 5
# set style line 4 lt rgb "#F25900" lw 2 pt 13

## Time-indexed plotting

set title "Quaternion states" font "Gill Sans,9"
set ylabel "q element value" font "Gill Sans,9"
set xlabel "Time (s)" font "Gill Sans,9"

plot "sim.log" u 1:2 w lines t "q0", \
     "sim.log" u 1:3 w lines t "q1", \
     "sim.log" u 1:4 w lines t "q2", \
     "sim.log" u 1:5 w lines t "q3"

set output "omega.pdf"

set title "Angular velocity" font "Gill Sans,9"
set ylabel "Angular velocity (/s)" font "Gill Sans,9"
set xlabel "Time (s)" font "Gill Sans,9"

plot "sim.log" u 1:6 w lines t "w_x", \
     "sim.log" u 1:7 w lines t "w_y", \
     "sim.log" u 1:8 w lines t "w_z"

set output "euler.pdf"

set title "Euler angles" font "Gill Sans,9"
set ylabel "Angle (rad)" font "Gill Sans,9"
set xlabel "Time (s)" font "Gill Sans,9"

plot "sim.log" u 1:9 w lines t  "roll", \
     "sim.log" u 1:10 w lines t "pitch", \
     "sim.log" u 1:11 w lines t "yaw"

set output "work.pdf"

set title "Kinetic energy vs. time" font "Gill Sans,9"
set ylabel "Energy (J)" font "Gill Sans,9"
set xlabel "Time (s)" font "Gill Sans,9"

plot "sim.log" u 1:12 w lines t "exogenous torque", \
     "sim.log" u 1:13 w lines t "total energy", \
     "sim.log" u 1:14 w lines t "Energy conservation error"
    
set output "momentum.pdf"

set title "Angular momentum vs. time" font "Gill Sans,9"
set ylabel "Angular momentum (Nms)" font "Gill Sans,9"
set xlabel "Time (s)" font "Gill Sans,9"
 
plot "sim.log" u 1:17 w lines t "L_x", \
     "sim.log" u 1:18 w lines t "L_y", \
     "sim.log" u 1:19 w lines t "L_z", \
     "sim.log" u 1:15 w lines t "Norm of angular momentum vector", \
     "sim.log" u 1:16 w lines t "Numerically integrated angular impulse"


# set multiplot
# set size 1,0.45
# set origin 0,0.5
# unset multiplot
