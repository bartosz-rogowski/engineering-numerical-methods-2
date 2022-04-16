set term png size 600, 600
set size square
# set size 1.,1.
# set key top right

set xrange [-6:6]
set yrange [-6:6]
set xlabel "x"
set ylabel "y"
set title "Grid"
set grid

set out 'nodes.png'
plot "nodes.dat" u 1:2 w l t ''

set term png size 2000, 650
set xrange [-5:5]
set yrange [-5:5]
set out 'modes.png'
set multiplot layout 2,5 rowsfirst
set pm3d map
set palette rgbformulae 22,13,10

do for [idx=0:9]{
	set title "{/Symbol l} = ".idx.""
	splot "wyniki.dat" i idx u 1:2:3 t ''
}

unset multiplot


reset
set term png size 1100, 1100
set out 'vals.png'
set multiplot layout 2,2 rowsfirst
# set xrange [0:9]
set xlabel "time"

set ylabel "y^kOc_2"
set title "y^kOc_2"
plot "vals.dat" u 1:2 w l t ''

set ylabel "y^kOc_3"
set title "y^kOc_3"
plot "vals.dat" u 1:3 w l t ''

set ylabel "y^kOy^k"
set title "y^kOy^k"
plot "vals.dat" u 1:4 w l t ''

set ylabel "y^kEy^k"
set title "y^kEy^k"
plot "vals.dat" u 1:5 w l t ''

unset multiplot

reset
set term gif size 600,500 animate delay 100
set output "maps.gif"
n=9
set view map # widok z gory
set size ratio -1
set cbr [-0.25:0.25]
set xrange [-5:5]
set yrange [-5:5]
set xlabel "x"
set ylabel "y"
set palette rgbformulae 22,13,10
# set title "u(x, y, t)"
do for [k=0:n] {
  set title sprintf("u(x, y, t); iteration: %i", 1000*(k+1))
  splot "maps.dat" index k u 1:2:3 w pm3d title sprintf("t=%i",k)
} 