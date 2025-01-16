set style data dots
set nokey
set xrange [0: 4.02877]
set yrange [-22.91167 : 16.83335]
set arrow from  1.47463, -22.91167 to  1.47463,  16.83335 nohead
set arrow from  2.32601, -22.91167 to  2.32601,  16.83335 nohead
set xtics ("G"  0.00000,"M"  1.47463,"K"  2.32601,"G"  4.02877)
 plot "graphene_band.dat"
