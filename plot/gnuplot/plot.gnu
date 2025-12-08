# -------------------------------
# Usage:
#   gnuplot -c plot.gp file1.txt file2.txt file3.txt
# -------------------------------

# Check if at least one file is provided
if (!exists("ARG1")) {
    print "Error: No data files provided."
    exit
}

set terminal x11
set logscale y 10
set xlabel "Index"
set ylabel "Value"
set grid
set title "Plot of data files (1-column, logscale y)"

# Build the plot command as a string - note the quotes around filenames
plotcmd = sprintf("'%s' using ($0):1 every 1000 with linespoints title '%s'", ARG1, ARG1)
if (exists("ARG2")) plotcmd = plotcmd . sprintf(", '%s' using ($0):1 every 1000 with linespoints title '%s'", ARG2, ARG2)
if (exists("ARG3")) plotcmd = plotcmd . sprintf(", '%s' using ($0):1 every 1000 with linespoints title '%s'", ARG3, ARG3)
if (exists("ARG4")) plotcmd = plotcmd . sprintf(", '%s' using ($0):1 every 1000 with linespoints title '%s'", ARG4, ARG4)
if (exists("ARG5")) plotcmd = plotcmd . sprintf(", '%s' using ($0):1 every 1000 with linespoints title '%s'", ARG5, ARG5)
if (exists("ARG6")) plotcmd = plotcmd . sprintf(", '%s' using ($0):1 every 1000 with linespoints title '%s'", ARG6, ARG6)
if (exists("ARG7")) plotcmd = plotcmd . sprintf(", '%s' using ($0):1 every 1000 with linespoints title '%s'", ARG7, ARG7)
if (exists("ARG8")) plotcmd = plotcmd . sprintf(", '%s' using ($0):1 every 1000 with linespoints title '%s'", ARG8, ARG8)
if (exists("ARG9")) plotcmd = plotcmd . sprintf(", '%s' using ($0):1 every 1000 with linespoints title '%s'", ARG9, ARG9)

eval "plot " . plotcmd
pause mouse close
