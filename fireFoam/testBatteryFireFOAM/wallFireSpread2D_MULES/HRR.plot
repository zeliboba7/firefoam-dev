plot './HRR/0/cellSource.dat' u 1:($3/1000) title "gas phase"
replot './patchPanel/0/faceSource.dat' u 1:($4/1000) title "panel"
replot './patchBurner/0/faceSource.dat' u 1:($4/1000) title "burner"
replot './patchOutlet/0/faceSource.dat' u 1:(-$4/1000) title "outlet"
replot "< paste ./HRR/0/cellSource.dat ./patchOutlet/0/faceSource.dat" u 1:(($3-$7)/1000) title "gas phase + outlet"
replot "< paste ./patchPanel/0/faceSource.dat ./patchBurner/0/faceSource.dat" u 1:(($4+$8)/1000) title "burner + panel"
set xlabel "time (s)"
set ylabel "HRR [kW]"
set xrange [0:*]
set title "Heat releast rate"
set key left top
replot

set terminal postscript eps enhanced color
set output '| epstopdf --filter --outfile=plot.HRR.pdf'
replot 
set output 
set terminal x11
