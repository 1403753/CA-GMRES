	set terminal pngcairo  transparent enhanced font "arial,10" fontscale 1.0 size 600, 400 
	set output 'residuals.png'

	set yrange [1e-10:1e-0]
	set xrange [0:1000]

	set style fill solid border -1
	set key bottom left
	set key box vertical width 2 height 1
	set key reverse Left
	
	set border 3
	set tics nomirror
	
	#set grid

	set xtics 100
	set xlabel "Iteration count"

	set ylabel "Residual 2-norm, log scale"	
	set logscale y
	set format y "%1.e"
	set ytics 1e-10,1e-1,1e+0
	set ytics add (-1,-1,-10)
	
	
	plot '../gnuplot/gmres_75.dat' using 1:2 title "GMRES(75)" linecolor "black" with lines, \
		'../gnuplot/gmres_mono_15_5.dat' using 1:2 title "Monomial-GMRES(15,5)" linecolor "blue" pt 21 ps 1, \
		'../gnuplot/gmres_newt_15_5.dat' using 1:2 title "Newton-GMRES(15,5)" linecolor "green" pt 8 ps 1, \
		'../gnuplot/gmres_mono_25_3.dat' using 1:2 title "Monomial-GMRES(25,3)" linecolor "magenta" pt 21 ps 1, \
		'../gnuplot/gmres_newt_25_3.dat' using 1:2 title "Newton-GMRES(25,3)" linecolor "cyan" pt 8 ps 1