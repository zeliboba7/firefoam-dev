plot './patchPanelQr/0/faceSource.dat' u 1:($3/1000) title "panel"
replot './patchWallQr/0/faceSource.dat' u 1:($3/1000) title "wall"
replot './patchSideQr/0/faceSource.dat' u 1:($3/1000) title "side"
replot './patchOutletQr/0/faceSource.dat' u 1:($3/1000) title "outlet"
replot './patchBurnerQr/0/faceSource.dat' u 1:($3/1000) title "burner"
replot './patchGroundQr/0/faceSource.dat' u 1:($3/1000) title "ground"
replot "< paste ./patchPanelQr/0/faceSource.dat ./patchWallQr/0/faceSource.dat ./patchSideQr/0/faceSource.dat ./patchOutletQr/0/faceSource.dat ./patchBurnerQr/0/faceSource.dat ./patchGroundQr/0/faceSource.dat" u 1:(($3+$6+$9+$12+$15+$18)/1000) title "sum of above"
replot './HRR/0/cellSource.dat' u 1:($3/1000) title "gas phase HRR"
replot './HRR/0/cellSource.dat' u 1:($3/1000*0.6) title "HRR * RadFrac"
set xlabel "time [s]"
set ylabel "integrated heat flux [kW]"
set xrange [0:*]
set title "Integraed Net Radiative Heat Flux"
replot

