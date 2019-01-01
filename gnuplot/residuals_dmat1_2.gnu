	set terminal epslatex font 'phv'
	set output 'residuals_dmat1_2.tex'
 	
	set xrange [0:1000]
	set yrange [1e-10:1e-0]

	set style fill solid border -1
	set key at 720, 1e-06
	#set key top right
	set key box vertical width 2 height 1
	set key reverse Left
	
	set border 3
	set tics nomirror
	
	#set grid
	
	set xlabel '\footnotesize Iteration count'
	set ylabel '\footnotesize Residual 2-norm, log scale' offset 3
	set logscale y
	set format x '\footnotesize %g'
	set format y '\footnotesize %1.e'
	
	set xtics 100,100,900 
	#set ytics 1e-10,1e-1,1e+0
	
	set ytics ('\footnotesize -10' 1e-10, '\footnotesize -9' 1e-9, '\footnotesize -8' 1e-8, '\footnotesize -7' 1e-7, '\footnotesize -6' 1e-6, \
		'\footnotesize -5' 1e-5, '\footnotesize -4' 1e-4, '\footnotesize -3' 1e-3, '\footnotesize -2' 1e-2, '\footnotesize -1' 1e-1, '\footnotesize 0' 1e-0,)
	
	set title '\shortstack{\footnotesize 1e+04 $\times$ 1e+04 diagonal, logspace eigs, cond 1e+05 \\ \footnotesize Residual 2-norm, log scale}'
	
	plot '../gnuplot/gmres_standard.dat' using 1:2 title '\footnotesize GMRES(75)' linecolor "black" with lines, \
		'../gnuplot/gmres_mono_small_s.dat' using 1:2 title '\footnotesize Monomial-GMRES(15,5)' linecolor "blue" pt 6 ps 1, \
		'../gnuplot/gmres_newt_small_s.dat' using 1:2 title '\footnotesize Newton-GMRES(15,5)' linecolor "green" pt 8 ps 1, \
		'../gnuplot/gmres_mono_large_s.dat' using 1:2 title '\footnotesize Monomial-GMRES(25,3)' linecolor "magenta" pt 6 ps 1, \
		'../gnuplot/gmres_newt_large_s.dat' using 1:2 title '\footnotesize Newton-GMRES(25,3)' linecolor "cyan" pt 8 ps 1

	set out