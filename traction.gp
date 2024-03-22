reset
set term png

set grid back

file_ref=sprintf("beam0.dat")

do for [i=0:10]{
Q=i

set output sprintf("beam_traction%d.png",i)
# Set the plot range
    set yrange [-1:5]
    set xrange [-2:6]

# Set labels for the axes
set xlabel "x"
set ylabel "y"

# Set title for the plot
set title "Plot of Traction"


#set label 1 sprintf("~{/symbol g}{1\\.} =%.2f, Scale=10^{-4}.", (1.0e-2)*i) at screen 0.97, 0.02 right front font ",15"

# Plot the first subplot
file = sprintf("beam%d.dat",i)

set label 1 sprintf("{/symbol q} = π/6,V = 0.5, U_0 = 0.8, X_0 = 1.5, Y_0 = 2.5") at screen 0.98, 0.02 right front font ",10"
set label 2 sprintf("Q=%.1f*10^{-7}",Q) at screen 0.2, 0.02 right front font ",12"

# Calculate the endpoint coordinates of the vectors
#set datafile separator ","  # Assuming data file is comma-separated
plot file every 10::1 using 1:2:9:10 with vectors head filled size screen 0.02,15,45 lw 2 lc "green" title "traction f0", \
     file using 1:2 with lines linewidth 2.7 lc rgb "red" title "Beam", \
     file every 10::1 using 3:4:11:12 with vectors head filled size screen 0.02,15,45 lw 2 lc "orange" title "traction f", \
     file using 3:4 with lines linewidth 2.7 lc rgb "blue" title "Beam",\
     file_ref using 1:2 with lines linewidth 2.3 lc rgb "black" title "reference",\
     file_ref using 3:4 with lines linewidth 2.3 lc rgb "black" title "reference"
}

