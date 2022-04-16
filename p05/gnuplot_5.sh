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

set xrange [-5:5]
set yrange [-5:5]
set out 'wynik.png'
set pm3d map
set palette rgbformulae 22,13,10
set title "u(x,y)"
splot "wyniki.dat" i 0 u 1:2:3