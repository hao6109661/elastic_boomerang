reset
set terminal png

# Give the value to the q
q=0.4
    
set xlabel 'I'
set ylabel "{/symbol q}_{eq}"
 set xrange [0:0.001]
 set yrange [-0.6:0.3]
set title "Fixed Points ({/symbol a}=0.02π...0.98π)"

set grid back

# Loop over all the files 
do for [i=1:49] {

    # Update the value of alpha
    alpha=(i/50.0)*pi
    
    # Set the output file 
    set output sprintf("elastic_beam_I_theta_q_0.4_%d.png", i)
    
    # Set file handle
    file = sprintf("elastic_beam_I_theta_q_0.4_%d.dat", i)

    # Plot the picture
    plot file using 1:(strcol(5) eq "#" ? 1/0 : $5) with points pointsize 1 pointtype 7 lc rgb "red" title sprintf('Frame %d, q=%.2f, {/symbol a}=%.3fπ', i, q, alpha/pi), \
         file using 1:(strcol(7) eq "#" ? 1/0 : $7) with points pointsize 1 pointtype 7 lc rgb "red" title ""
    
}
