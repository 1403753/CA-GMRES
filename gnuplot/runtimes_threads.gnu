	set terminal epslatex font 'phv' 8
	set output 'runtimes_threads.tex'
	#set terminal pngcairo
	#set output 'runtimes_threads.png'

	set yrange [0:0.5]	
	set offset 0, 0, 2, 0	

	set style fill solid border -1
	set key invert reverse Left top center opaque
	set key box vertical samplen 3
	set key height 0.1 width -5
	
	set grid

	num_of_ksptypes=2
	set boxwidth 0.5/num_of_ksptypes
	dx=0.5/num_of_ksptypes
	offset=-0.12

	set format y '\footnotesize %g'
	set ytics 0,1,100

	set xlabel '\footnotesize Threads'
	set ylabel '\footnotesize Relative runtime'
	
	mname = system("head -1 " . '../gnuplot/threads_gmres.dat' . " | awk '{print $6}'")
	
	set title mname
	
	plot '../gnuplot/threads_gmres_ca.dat' using ($1+offset):($2+$3+$4+$5+$6) title '\scriptsize Small dense operations' linecolor rgb '#006400' with boxes, \
			 ''          using ($1+offset):($3+$4+$5+$6) title '\scriptsize Block Gram-Schmidt' linecolor rgb '#FFFF00' with boxes, \
			 ''          using ($1+offset):($4+$5+$6) title '\scriptsize TSQR' linecolor rgb '#FFA500 ' with boxes, \
			 ''          using ($1+offset):($5+$6) title '\scriptsize SpMV' linecolor rgb '#0000FF' with boxes, \
			 ''          using ($1+offset):6 title '\scriptsize INIT' linecolor rgb 'red' with boxes, \
			 '../gnuplot/threads_gmres.dat' using ($1+offset+dx):($3+$4) title '\scriptsize MGS' linecolor rgb "#8B008B" with boxes, \
			 ''                   using ($1+offset+dx):4 title '\scriptsize SpMV' linecolor rgb "#0000FF" with boxes, \
			 ''                   using ($1):4:xtic(2) notitle with points pt -1, \
			 ''                   using ($1 - 0.15):($4 + 0.5):5 notitle with labels rotate by 90
	set out