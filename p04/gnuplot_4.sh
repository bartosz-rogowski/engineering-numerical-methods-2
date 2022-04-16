set term png size 1200, 400

set out 'wyniki.png'

set xrange [0:pi]
set yrange [0:pi]
set xtics ("0" 0, "0.25{/Symbol p}" pi/4, "0.5{/Symbol p}" pi/2, "0.75{/Symbol p}" 0.75*pi, "{/Symbol p}" pi)
set ytics ("0" 0, "0.25{/Symbol p}" pi/4, "0.5{/Symbol p}" pi/2, "0.75{/Symbol p}" 0.75*pi, "{/Symbol p}" pi)
set xlabel "x"
set ylabel "y"
set grid
set pm3d map
set cntrparam levels 20
unset clabel #all contour lines color set to black
set contour
set palette rgbformulae 22,13,10

set multiplot layout 1,3 rowsfirst
	set size square
	set title "nx=ny=3"
	splot "n3.dat" i 0 u 1:2:3 notitle

	set title "nx=ny=10"
	splot "n10.dat" i 0 u 1:2:3 notitle

	set title "Analytical solution"
	splot "ndokl.dat" i 0 u 1:2:3 notitle

unset multiplot