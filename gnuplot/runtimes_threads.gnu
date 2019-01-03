	set terminal pngcairo  transparent enhanced font "arial,10" fontscale 1.0 size 600, 400 
	set output 'runtimes.png'

	set xtics("1" 1, "2" 2, "4" 3, "8" 4)
	set yrange [0:100]

	set style fill solid border -1
	set key invert
	set grid

	num_of_ksptypes=2
	set boxwidth 0.5/num_of_ksptypes
	dx=0.5/num_of_ksptypes
	offset=-0.12

	set xlabel "threads"
	set ylabel "seconds"

	plot 'data1.dat' using ($1+offset):($2+$3+$4+$5) title "SDO" linecolor rgb "#006400" with boxes, \
			 ''                   using ($1+offset):($3+$4+$5) title "BGM" linecolor rgb "#FFFF00" with boxes, \
			 ''                   using ($1+offset):($4+$5) title "TSQR" linecolor rgb "#FFA500 " with boxes, \
			 ''                   using ($1+offset):5 notitle linecolor rgb "#0000FF" with boxes, \
			 'data2.dat' using ($1+offset+dx):($2+$3) title "MGS" linecolor rgb "#8B008B" with boxes, \
			 ''                   using ($1+offset+dx):3 title "SpMV" linecolor rgb "#0000FF" with boxes