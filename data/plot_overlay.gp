#!/usr/bin/gnuplot
#
# Gnuplot script - Overlay comparison of all filters
# Usage: gnuplot plot_overlay.gp
#

set terminal png size 1400,900 enhanced font 'Arial,12'
set output 'filter_overlay.png'

set title "NMR Spectrum - Filter Comparison (Overlay)" font 'Arial,16'
set xlabel "Chemical Shift (ppm)" font 'Arial,14'
set ylabel "Intensity" font 'Arial,14'
set grid
set key top right

# Plot all three on the same axes
plot 'shifted_data.txt' using 1:2 with lines lw 1 lc rgb 'blue' title 'Original (Unfiltered)', \
     'filtered_boxcar.txt' using 1:2 with lines lw 2 lc rgb 'red' title 'Boxcar (9-pt, 3 passes)', \
     'filtered_sg.txt' using 1:2 with lines lw 2 lc rgb 'green' title 'Savitzky-Golay (11-pt, 3 passes)'

print "Overlay plot saved to: filter_overlay.png"
