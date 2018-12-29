	set terminal pngcairo  transparent enhanced font "arial,10" fontscale 1.0 size 600, 400 
	set output 'residuals.png'

	set yrange [1e-10:1e-0]

	set style fill solid border -1
	set key bottom right
	set grid

	set xtics 100
	set xlabel "Iteration count"

	set ylabel "Residual 2-norm, log scale"	
	set logscale y
	set format y "%1.e"
	set ytics 1e-10,1e-1,1e-0
	
	plot '../gnuplot/std_gmres.dat' using 1:2 title "standard GMRES" with lines, \
		'../gnuplot/gmres_monomial.dat' using 1:2 title "CA-GMRES mono" linecolor rgb "#006400" pt 6 ps 1, \
		'../gnuplot/gmres_newton.dat' using 1:2 title "CA-GMRES newt" linecolor "red" pt 3 ps 1