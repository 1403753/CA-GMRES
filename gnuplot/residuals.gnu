	set terminal epslatex font 'phv'
	
	fname = system("head -1 " . '../gnuplot/gmres_standard.dat' . " | awk '{print $1}'")
	
	set output fname . '.tex'
 	
	set xrange [0:1000]
	set yrange [1e-16:1e-0]

	set style fill solid border -1
	set key at 480, 9e-06
	#set key top right opaque
	set key box vertical width -70 height 0.1
	set key reverse Left
	set key spacing 1.45
	set key samplen 2
	set border 3
	set tics nomirror
	
	#set grid
	
	set xlabel '\footnotesize Iteration count'
	set ylabel '\footnotesize Residual 2-norm, log scale' offset 3
	set logscale y
	set format x '\footnotesize %g'
	set format y '\footnotesize %1.e'
	
	set xtics 100,100,1000 
	#set ytics 1e-10,1e-1,1e+0
	
	set ytics ('\footnotesize -10' 1e-10, '\footnotesize -9' 1e-9, '\footnotesize -8' 1e-8, '\footnotesize -7' 1e-7, '\footnotesize -6' 1e-6, \
		'\footnotesize -5' 1e-5, '\footnotesize -4' 1e-4, '\footnotesize -3' 1e-3, '\footnotesize -2' 1e-2, '\footnotesize -1' 1e-1, '\footnotesize 0' 1e-0,)
	
  cmin1 = system("head -1 " . '../gnuplot/gmres_mono_small_s.dat' . " | awk '{print $4}'")
  cmax1 = system("head -1 " . '../gnuplot/gmres_mono_small_s.dat' . " | awk '{print $5}'")
  cmin2 = system("head -1 " . '../gnuplot/gmres_newt_small_s.dat' . " | awk '{print $4}'")
  cmax2 = system("head -1 " . '../gnuplot/gmres_newt_small_s.dat' . " | awk '{print $5}'")
  cmin3 = system("head -1 " . '../gnuplot/gmres_mono_large_s.dat' . " | awk '{print $4}'")
  cmax3 = system("head -1 " . '../gnuplot/gmres_mono_large_s.dat' . " | awk '{print $5}'")
  cmin4 = system("head -1 " . '../gnuplot/gmres_newt_large_s.dat' . " | awk '{print $4}'")
  cmax4 = system("head -1 " . '../gnuplot/gmres_newt_large_s.dat' . " | awk '{print $5}'")

	m = system("head -1 " . '../gnuplot/gmres_standard.dat' . " | awk '{print $2}'")	

	s1 = system("head -1 " . '../gnuplot/gmres_newt_small_s.dat' . " | awk '{print $2}'")
	t1 = system("head -1 " . '../gnuplot/gmres_newt_small_s.dat' . " | awk '{print $3}'")
	
	s2 = system("head -1 " . '../gnuplot/gmres_newt_large_s.dat' . " | awk '{print $2}'")
	t2 = system("head -1 " . '../gnuplot/gmres_newt_large_s.dat' . " | awk '{print $3}'")
	
	mname = system("head -1 " . '../gnuplot/gmres_standard.dat' . " | awk '{print $3}'")
	
	set title mname
	
	plot '../gnuplot/gmres_standard.dat' using 1:2 title '\scriptsize GMRES('. m .')' linecolor "black" with lines, \
		'../gnuplot/gmres_mono_small_s.dat' using 1:2 title '\begin{minipage}[l]{.95\textwidth} \scriptsize Monomial-GMRES(' . s1 . ',' . t1 . ') \newline \tiny min, max basis rcond \#: ' . cmin1 . ', ' . cmax1 . '\end{minipage}' linecolor "blue" pt 6 ps 1 , \
		'../gnuplot/gmres_newt_small_s.dat' using 1:2 title '\begin{minipage}[l]{.95\textwidth} \scriptsize Newton-GMRES(' . s1 . ',' . t1 . ') \newline \tiny min, max basis rcond \#: ' . cmin2 . ', ' . cmax2 . '\end{minipage}' linecolor "green" pt 8 ps 1, \
		'../gnuplot/gmres_mono_large_s.dat' using 1:2 title '\begin{minipage}[l]{.95\textwidth} \scriptsize Monomial-GMRES(' . s2 . ',' . t2 . ') \newline \tiny min, max basis rcond \#: ' . cmin3 . ', ' . cmax3 . '\end{minipage}' linecolor "magenta" pt 6 ps 1, \
		'../gnuplot/gmres_newt_large_s.dat' using 1:2 title '\begin{minipage}[l]{.95\textwidth} \scriptsize Newton-GMRES(' . s2 . ',' . t2 . ') \newline \tiny min, max basis rcond \#: ' . cmin4 . ', ' . cmax4 . '\end{minipage}' linecolor "cyan" pt 8 ps 1

	set out