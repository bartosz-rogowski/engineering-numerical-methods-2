set term png size 1600, 1000
# set size square
# set size 1.,1.
set out 'wyniki.png'
set key bottom right
set grid
set xrange [-3.2:3.2]
set xlabel "x"
set ylabel "{/Symbol D}u(x)"
set title "R_n = f(E)"
set multiplot layout 2,3 rowsfirst

	# set out 'kolokacja_1.png'
	set title 'baza 1 - m. kolokacji'
	plot for [idx=0:4] 'kolokacja.dat' i idx u 1:2 w l lw 2 t 'N='.(idx+6)

	# set out 'kwadraty_1.png'
	set key top right
	set title 'baza 1 - m. najmniejszych kwadratow'
	plot for [idx=0:4] 'kwadraty.dat' i idx u 1:2 w l lw 2 t 'N='.(idx+6)

	# set out 'galerkin_1.png'
	set key bottom right
	set title 'baza 1 - m. galerkina'
	plot for [idx=0:4] 'galerkin.dat' i idx u 1:2 w l lw 2 t 'N='.(idx+6)

	# set out 'kolokacja_2.png'
	set title 'baza 2 - m. kolokacji'
	plot for [idx=0:4] 'kolokacja.dat' i idx+5 u 1:2 w l lw 2 t 'N='.(idx+6)

	# set out 'kwadraty_2.png'
	set key bottom center
	set title 'baza 2 - m. najmniejszych kwadratow'
	plot for [idx=0:4] 'kwadraty.dat' i idx+5 u 1:2 w l lw 2 t 'N='.(idx+6)

	# set out 'galerkin_2.png'
	set title 'baza 2 - m. galerkina'
	plot for [idx=0:4] 'galerkin.dat' i idx+5 u 1:2 w l lw 2 t 'N='.(idx+6)


unset multiplot