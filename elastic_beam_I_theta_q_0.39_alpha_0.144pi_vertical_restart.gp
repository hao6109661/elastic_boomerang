reset
# Set the size of the plot
set terminal pngcairo size 1200,1200

# Give the value to the q,alpha and initial value
q=0.39
alpha=0.144*pi
initial1=-4.80

# Set file handle for the general information which contains the I, theta_eq and step number 
file1 = sprintf("RESLT_q_%.2f_alpha_%.3fpi/elastic_beam_I_theta_q_%.2f_alpha_%.3fpi_initial_%.2f.dat", q, alpha/pi, q, alpha/pi, initial1)
file2 = sprintf("RESLT_q_%.2f_alpha_%.3fpi_restart/elastic_beam_I_theta_q_%.2f_alpha_%.3fpi_initial_%.2f.dat", q, alpha/pi, q, alpha/pi, initial1)

# Got the max value of I from the general file 
max_value = system(sprintf("awk '{if(NR==1 || $1 > max) max=$1} END {print max}' \"%s\"", file1))
#print "Max value in the first column:", max_value

# Got the number of I values in the general file 
stats file1 nooutput
num_lines = STATS_records
#print "Number of lines in the file:", num_lines

# Set the number the increment for I    
step=1000

# Set the tolerace for the small interval |I_from_the_general_file-target|<epsilon
epsilon = 0.5*1.0e-5

# Loop over different I value within the interval [0, I_max], I_max is "max_value" in the previous text
do for [i=0:step] {

# This is the target value for I
target = (max_value/step)*i

# Set the output file 
set output sprintf("RESLT_q_%.2f_alpha_%.3fpi_plot/combine_elastic_beam_I_theta_q_%.2f_alpha_%.3fpi_initial_%.2f_%d.png", q, alpha/pi, q, alpha/pi, initial1, i)

# Set up the multiplot layout
set multiplot layout 2,2

# Plot the first subplot
set xlabel 'I'
set ylabel "{/symbol q}_{eq}"
set xrange [-0.00001:0.02]
set yrange [-3.45:-3.1]
set title "Fixed Points"

set grid back

# Define label size and offset for Newton iteration number
label_size = 15
label_offset_x = 1.0
label_offset_y = 0.5

# Show the value of I in the plot
set label 1 sprintf("I = %.9f", target) at screen 0.97, 0.01 right front font ",15"

# Produce a tempprary file to record the data (I_value, theta_eq, step) which satisfies the condition |I_from_the_general_file-target|<epsilon
set table 'filtered_data.dat'
plot file1 using 1:3:(abs($1 - target) <= epsilon ? $9 : "#"):(abs($1 - target) <= epsilon ? $8 : "#") with table
unset table

# Produce a tempprary file for the restart inner closed loop solution
# Change its step number to float type to identify its result directory
set table 'filtered_data2.dat'
plot file2 using 1:3:(abs($1 - target) <= epsilon ? ($9)+0.1 : "#"):(abs($1 - target) <= epsilon ? $8 : "#") with table
unset table

# Set the scanning line to the first plot  
set arrow 1 from target, graph 0 to target, graph 1 nohead lc rgb "blue" lw 1.5 front

# Plot the picture
plot file1 using 1:(strcol(3) eq "#" ? 1/0 : $3) with lines linewidth 2 lc rgb "red" title sprintf('q=%.2f, {/symbol a}=%.3fπ', q, alpha/pi), \
     file1 using 1:3 every ::0::0 with points pointsize 1.2 pointtype 7 lc rgb "black" title "I=0" ,\
     file1 every ::(num_lines-1)::(num_lines-1) using 1:3 with points pointsize 1.2 pointtype 7 lc rgb "black" title "",\
     file2 using 1:(strcol(3) eq "#" ? 1/0 : $3) with lines linewidth 2 lc rgb "red" title "", \
     "filtered_data.dat" using 1:2:3 with points pointsize 1 pointtype 7 lc rgb "blue" title "point considered",\
     "filtered_data2.dat" using 1:2:3 with points pointsize 1 pointtype 7 lc rgb "blue" title ""
     #"filtered_data.dat" using 1:(strcol(4) eq "#" ? 1/0 : $2):(strcol(4) eq "#" ? 1/0 : $4) with labels offset label_offset_x,label_offset_y font ",label_size" tc rgb "blue" title ""
         
# Remove arrow after the first subplot
unset arrow 1
         
# Plot the second subplot
#set size ratio -1

set grid back

# Set the plot range
set xrange [-1.1:0.2]
set yrange [-0.6:0.7]

# Set labels for the axes
set xlabel "x"
set ylabel "y"

# Set title for the plot
set title "Plot of Beam with two arms"

# Got the step numbers for different I which are the cases satisfied the condition |I_from_the_general_file-target|<epsilon
# Note that here this is in the type of string
suffix_list_str = system("awk '$3 ~ /^[0-9]+(\\.[0-9]+)?$/{print $3}' filtered_data.dat filtered_data2.dat")
#print suffix_list_str
#suffix_list_str2 = system("awk '$3 ~ /^[0-9]+(\\.[0-9]+)?$/{print $3}' filtered_data2.dat")
#print suffix_list_str2
#print suffix_list_str

# Got the length of this string
suffix_list_length = words(suffix_list_str)

# Declare an array
array suffix_list[suffix_list_length]

# Split the string into words and convert each word to a number
# Note that since we cannot manipulate the elements in a string, so it should be converted to an array 
do for [i=1:suffix_list_length] {
    num_str = word(suffix_list_str, i)
    suffix_list[i] = real(num_str)
}

# Now, suffix_list contains the numbers extracted from the string
#print "Suffix list:", suffix_list

# Got the values of theta-eq which are satisfied the condition |I_from_the_general_file-target|<epsilon
# Note that here this is in the type of string
theta_list_str = system("awk '$3 ~ /^[0-9]+(\\.[0-9]+)?$/{print $2}' filtered_data.dat filtered_data2.dat")
#print theta_list_str

# Got the length of this string
theta_list_length = words(theta_list_str)

# Declare an array
array theta_list[theta_list_length]


# Split the string into words and convert each word to a number
# Note that since we cannot manipulate the elements in a string, so it should be converted to an array 
do for [i=1:theta_list_length] {
    num_str = word(theta_list_str, i)
    theta_list[i] = real(num_str)
}
#print theta_list


# This is a bubble sort algorithm, arranging step numbers according to the size of theta_eq, 
# so that points on different solution branches correspond to the shape of boomerang
do for [i = 1:|theta_list|] {
    do for [j = 1:|theta_list|-i] {
        if (theta_list[j] > theta_list[j+1]) {
            tmp1 = suffix_list[j]
            suffix_list[j] = suffix_list[j+1]
            suffix_list[j+1] = tmp1
            
            tmp2 = theta_list[j]
            theta_list[j] = theta_list[j+1]
            theta_list[j+1] = tmp2
        }
    }
}


# This is some simple outputs to check the code
#print suffix_list
#print theta_list
#print suffix_list[1]
#print suffix_list[|suffix_list|]



####################################################################################################################
# This is just a simple test, we still leave it here in case we will use it for different cases in the future

# Declare an array to store the solutions in the same branch 
#array current_group[1]
#current_group[1] = suffix_list[1]

# Add the second element from suffix_list to current_group, and the second 
#array new_current_group[|current_group| + 1]
#do for [j=1:|current_group|] {
    #new_current_group[j] = current_group[j]
#}
#new_current_group[|current_group| + 1] = suffix_list[2]

# Reassign the modified array back to current_group
#array current_group[|new_current_group|]
# Copy elements from new_current_group to current_group
#do for [j=1:|new_current_group|] {
    #current_group[j] = new_current_group[j]
#}
#print current_group
####################################################################################################################



# Declare an array to store the solutions in the same branch 
# (such solutions are very closed, so we will not give the offset to the corresponding figures of boomerang)
array current_group[1]

# Record the number of branches we have for this scanning line 
counter_branch=0

# Iterate through the elements of suffix_list
do for [i=1:|suffix_list|] {
    # If it's the first element, add it directly to the current group
    if (i == 1) {
        current_group[1] = suffix_list[i]
    } else {
        # Calculate the difference between the current element and the previous one
        # Again, suffix_list stores the step numbers which are satisfied the condition
        diff = abs(theta_list[i] - theta_list[i-1])
        # If the difference is less than 1.0e-3, add the current element to the current group
        # It means that these adjacent solutions are in the same branches and are very closed
        if (diff < 1.0e-3) {
            # Add the second element from suffix_list to current_group
            # Since gnuplot does not have the push_back function to add new elements to an array, we will do it in the following lines
            array new_current_group[|current_group| + 1]
            do for [j=1:|current_group|] {
               new_current_group[j] = current_group[j]
            }
           new_current_group[|current_group| + 1] = suffix_list[i]

           # Reassign the modified array back to current_group
           array current_group[|new_current_group|]
           # Copy elements from new_current_group to current_group
           do for [j=1:|new_current_group|] {
              current_group[j] = new_current_group[j]
           }
        } else {  
            # If the difference (theta_eq) is greater than or equal to 1.0e-3, start a new group
            # It means that the adjacent solutions are not in the same branch
        
            # Set files to record the data for solutions in different branches  
            filename = sprintf("branch%d.dat", counter_branch)
            set table filename
            plot current_group using 2 with table
            unset table
            
            counter_branch=counter_branch+1
        
            # Start a new group to store the solutions in the same branch
            array current_group[1]
            current_group[1] = suffix_list[i]
        }
    }
}

# Record the final group (store the solutions which are closed and in the same branch) to a file 
filename = sprintf("branch%d.dat", counter_branch)
set table filename         
plot current_group using 2 with table 
unset table


# Plot the second subplot

# Initialize filelist elements to empty strings
array filelist[2*counter_branch+2]
do for [i=1:2 * counter_branch + 2] {
  filelist[i] = ""
}


# Loop over all the branch groups which are satisfied the condition 
do for [i=0:counter_branch] {
  filename = sprintf("branch%d.dat", i)
  
  # The source of the step number (normal or restart) is identified by its type. 
  # Integer: normal; float: restart
  integers = system(sprintf("awk '$1 ~ /^[0-9]+$/{print $1}' %s", filename))
  
  # When navigating the boomerang-shaped files in the directory, 
  # the step number should remain an integer, so it is output as an int type
  float_rounded = system(sprintf("awk '$1 ~ /^[0-9]*\\.[0-9]+$/{print int($1)}' %s", filename))

  # Assemble the file name of outer solutions
  do for [suffix in integers] {
     # Create file names using suffix and add them to filelist
     filename_first = sprintf("RESLT_q_%.2f_alpha_%.3fpi/beam_first_arm_initial_%.2f_%s.dat", q, alpha/pi, initial1, suffix)
     filelist[i+1] = filelist[i+1] . "\"" . filename_first . "\" "
     filename_second = sprintf("RESLT_q_%.2f_alpha_%.3fpi/beam_second_arm_initial_%.2f_%s.dat", q, alpha/pi, initial1, suffix)
     filelist[i+counter_branch+2] = filelist[i+counter_branch+2] . "\"" . filename_second . "\" "}
    
  # Assemble the file name of inner closed loop solutions  
  do for [suffix in float_rounded] {
     # Create file names using suffix and add them to filelist
      filename_first = sprintf("RESLT_q_%.2f_alpha_%.3fpi_restart/beam_first_arm_initial_%.2f_%s.dat", q, alpha/pi, initial1, suffix)
      filelist[i+1] = filelist[i+1] . "\"" . filename_first . "\" "
      filename_second = sprintf("RESLT_q_%.2f_alpha_%.3fpi_restart/beam_second_arm_initial_%.2f_%s.dat", q, alpha/pi, initial1, suffix)
      filelist[i+counter_branch+2] = filelist[i+counter_branch+2] . "\"" . filename_second . "\" "}
}
    
# Plot the figures of boomerang with offsets so that it can be shown clearly
# Note that to display the legend only once in plot, we initially plot filename_first with the legend, 
# and subsequently plot the real solutions without legend
plot filename_first using 3:(0.2*(counter_branch)-0.4+column(4)) with lines linewidth 2 lc rgb "blue" title "boomerang", \
     for [i=1:counter_branch+1] \
        for [file in filelist[i]] \
             file using 3:(0.2*(i-1)-0.4+column(4)) with lines linewidth 2 lc rgb "blue" title "", \
     for [i=counter_branch+2:2*counter_branch+2] \
        for [file in filelist[i]] \
             file using 3:(0.2*(i-counter_branch-2)-0.4+column(4)) with lines linewidth 2 lc rgb "blue" title ""


# Plot the third subplot   

set grid back
#set size ratio -1


# Set the plot range
set xrange [-1.1:0.2]
set yrange [-0.1:0.1]

# Set labels for the axes
set xlabel "x"
set ylabel "y"

# Set title for the plot
set title "Plot of Beam with two arms"  
     
filelist = ""
# Loop through the elements of the suffix_list array
do for [i=1:|suffix_list|] {
    suffix = suffix_list[i]
    
    # Identify the type of the step mumber 
    if(suffix-floor(suffix)==0){
    # Step number is integer. Navigate the files from the normal (outer sloutions) directory
    
    # Create file names using suffix and add them to filelist
    filename_first = sprintf("RESLT_q_%.2f_alpha_%.3fpi/beam_first_arm_initial_%.2f_%d.dat", q, alpha/pi, initial1, suffix)
    filelist = filelist . "\"" . filename_first . "\" "
    filename_second = sprintf("RESLT_q_%.2f_alpha_%.3fpi/beam_second_arm_initial_%.2f_%d.dat", q, alpha/pi, initial1, suffix)
    filelist = filelist . "\"" . filename_second . "\" "
    }else{
    # Step number is float. Navigate the files from the inner closed loop sloution directory
    
    # Create file names using suffix and add them to filelist
    filename_first = sprintf("RESLT_q_%.2f_alpha_%.3fpi_restart/beam_first_arm_initial_%.2f_%d.dat", q, alpha/pi, initial1, suffix)
    filelist = filelist . "\"" . filename_first . "\" "
    filename_second = sprintf("RESLT_q_%.2f_alpha_%.3fpi_restart/beam_second_arm_initial_%.2f_%d.dat", q, alpha/pi, initial1, suffix)
    filelist = filelist . "\"" . filename_second . "\" "
    }
}

# This the figures of bommerang without offset
plot filename_first using 3:4 with lines linewidth 2 lc rgb "blue" title "boomerang",\
     for [file in filelist] file using 3:4 with lines linewidth 2 lc rgb "blue" title ""

#set size noratio

# Plot the fourth subplot

# Set labels for the axes
set xlabel "x"
set ylabel "y"

set xrange [-1.1:0.2]
set yrange [-0.6:0.7]

# Set title for the plot
set title "Plot of Traction"

# Sicne the magnitude of the traction is too little, we scale it so that it can be shown obviously and clearly
scale=10

# Initialize filelist elements to empty strings
array filelist[2*counter_branch+2]
do for [i=1:2 * counter_branch + 2] {
  filelist[i] = ""
}

# Loop over all the branch groups which are satisfied the condition   
do for [i=0:counter_branch] {
  filename = sprintf("branch%d.dat", i)
  
  # The source of the step number (normal or restart) is identified by its type. 
  # Integer: normal; float: restart
  integers = system(sprintf("awk '$1 ~ /^[0-9]+$/{print $1}' %s", filename))
  
  # When navigating the boomerang-shaped files in the directory, 
  # the step number should remain an integer, so it is output as an int type
  float_rounded = system(sprintf("awk '$1 ~ /^[0-9]*\\.[0-9]+$/{print int($1)}' %s", filename))

  # Assemble the file name of outer solutions
  do for [suffix in integers] {
    # Create file names using suffix and add them to filelist
    filename_first = sprintf("RESLT_q_%.2f_alpha_%.3fpi/beam_first_arm_initial_%.2f_%s.dat", q, alpha/pi, initial1, suffix)
    filelist[i+1] = filelist[i+1] . "\"" . filename_first . "\" "
    filename_second = sprintf("RESLT_q_%.2f_alpha_%.3fpi/beam_second_arm_initial_%.2f_%s.dat", q, alpha/pi, initial1, suffix)
    filelist[i+counter_branch+2] = filelist[i+counter_branch+2] . "\"" . filename_second . "\" "}
  
  # Assemble the file name of inner closed loop solutions     
  do for [suffix in float_rounded] {
    # Create file names using suffix and add them to filelist
    filename_first = sprintf("RESLT_q_%.2f_alpha_%.3fpi_restart/beam_first_arm_initial_%.2f_%s.dat", q, alpha/pi, initial1, suffix)
    filelist[i+1] = filelist[i+1] . "\"" . filename_first . "\" "
    filename_second = sprintf("RESLT_q_%.2f_alpha_%.3fpi_restart/beam_second_arm_initial_%.2f_%s.dat", q, alpha/pi, initial1, suffix)
    filelist[i+counter_branch+2] = filelist[i+counter_branch+2] . "\"" . filename_second . "\" "}
}
      
# Plot the shape of the boomrang with the traction 
plot filename_first using 3:(0.2*(counter_branch)-0.4+column(4)) with lines linewidth 2 lc rgb "blue" title "boomerang", \
     filename_first every 5::1 using 3:(0.2*(counter_branch)-0.4+column(4)):(column(11)*scale):(column(12)*scale) with vectors head nofilled size screen 0.02,15,45 lw 2 lc "red" title "10*traction f", \
     for [i=1:counter_branch+1] \
         for [file in filelist[i]] \
             file using 3:(0.2*(i-1)-0.4+column(4)) with lines linewidth 2 lc rgb "blue" title "", \
     for [i=1:counter_branch+1] \
         for [file in filelist[i]] \
             file every 5::1 using 3:(0.2*(i-1)-0.4+column(4)):(column(11)*scale):(column(12)*scale) with vectors head nofilled size screen 0.02,15,45 lw 2 lc "red" title "", \
     for [i=counter_branch+2:2*counter_branch+2] \
         for [file in filelist[i]] \
             file using 3:(0.2*(i-counter_branch-2)-0.4+column(4)) with lines linewidth 2 lc rgb "blue" title "", \
     for [i=counter_branch+2:2*counter_branch+2] \
         for [file in filelist[i]] \
             file every 15::1 using 3:(0.2*(i-counter_branch-2)-0.4+column(4)):(column(11)*scale):(column(12)*scale) with vectors head nofilled size screen 0.02,15,45 lw 2 lc "red" title ""

     

# Remove label for the next frame
unset label 1

# Exit the mode of multiplot
unset multiplot

# Delte the temporary files
system("rm filtered_data.dat")
system("rm filtered_data2.dat")
do for [i=0:counter_branch] {
    system(sprintf("rm branch%d.dat", i))
}

}

# Exit the mode of multiplot
unset multiplot

# Save the final multiplot layout to the output file
set output
