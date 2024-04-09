reset

# Set the size of the plot
set terminal png size 1200,600

# Give the value to the q
q=0.4
 
# Loop over all the files 
do for [i=0:20] {
    
    # Set the output file 
    set output sprintf("combine_elastic_alpha_theta_q_0.4_%d.png", i)

    # Set up the multiplot layout
    set multiplot layout 1,2
    
    # Non-dimensional coefficient (FSI)
    I = (1.0e-3/20.0)*i

    # Plot the first subplot
    file = sprintf("elastic_beam_alpha_theta_q_0.4_%d.dat", i)
    
    set xlabel "{/symbol a}"
    set ylabel "{/symbol q}_{eq}"
    set xrange [0:3.15]
    set yrange [-2:2]
    set grid back
    set title 'Fixed Points (I=0...0.001)'
    
    # Set the ticks of the x-axis and y-axis 
    set xtics("0" 0, "0.25π" 0.25*pi, "0.5π" 0.5*pi, "0.75π" 0.75*pi, "π" pi)
    set ytics("-2" -2, "-1.5" -1.5, "-1" -1, "-0.5" -0.5, "0" 0, "0.5" 0.5, "1" 1, "1.5" 1.5, "2" 2)

    # Plot the first plot  
    plot file using 3:(strcol(4) eq "#" ? 1/0 : $5) with points pointsize 1.2 pointtype 7 lc rgb "red" title sprintf('Frame %d, q=%.2f, I=%.5f', i, q, I), \
         file using 3:(strcol(6) eq "#" ? 1/0 : $7) with points pointsize 1.2 pointtype 7 lc rgb "red" title ""

    # Add a line to indicate the value of q in the first subplot
    set arrow 1 from q, graph 0 to q, graph 1 nohead lc rgb "orange" lw 3.7 front

    # Add information in the right bottom corner
    set label 1 sprintf("I=%.5f, q=%.2f", I, q) at screen 0.98, 0.02 right front font ",15"


    # Plot the second subplot
    unset colorbox
    set title "Region of Steady Points Appeared (I=0...0.001)"

    # Set point type to solid circle
    set style fill solid 1.0 noborder

    set ylabel "{/symbol a}"
    set xlabel "q"
    set xrange [0:0.5]
    set yrange [0:3.15]

    # Set the ticks of the x-axis and y-axis
    set ytics( "0" 0 , "0.25π" 0.25*pi , "0.5π" 0.5*pi , "0.75π" 0.75*pi , "π" pi)
    set xtics("0" 0, "0.1" 0.1, "0.2" 0.2, "0.3" 0.3, "0.4" 0.4, "0.5" 0.5)

    # Plot the second plot 
    file2 = sprintf("beam_contour_theta_eq%d.dat",i)
    plot file2 using 2:3:5 with points pointsize 1.2 pointtype 7 title ""

    # Remove arrow after the second subplot
    unset arrow 1

    # Remove label for the next frame
    unset label 1

    unset multiplot
}

# Save the final multiplot layout to the output file
set output
