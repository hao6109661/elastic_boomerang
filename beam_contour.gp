reset
set terminal png

# Loop over all files 
do for [i=0:20]{

# Non-dimensional coefficient (FSI)
I = (1.0e-3/20.0)*i

# Set the output file 
# Note that the output files will not pop up, instead they will be saved in the directory automatically 
set output sprintf("beam_contour%d.png",i)

unset colorbox
set title "Region of Steady Points Appeared (I=0...0.001)"

# Set point type to solid circle
set style fill solid 1.0 noborder

set ylabel "{/symbol a}"
set xlabel "q"
set xrange [0:0.5]
set yrange [0:3.15]
unset key
set grid back

# Set the ticks of the y-axis
set ytics( "0" 0 , "0.25π" 0.25*pi , "0.5π" 0.5*pi , "0.75π" 0.75*pi , "π" pi)

# Write the value of I in the bottom right corner of the graph
set label 1 sprintf("I=%.5f", I) at screen 0.98, 0.02 right front font ",15"

# Set file handle
file = sprintf("beam_contour_theta_eq%d.dat",i)

# Plot scatter plot with smaller solid circles
plot file using 2:3:5 with points pointsize 1.2 pointtype 7

}
