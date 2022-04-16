set term png size 800, 600
set out "zad_1.png"
set xrange [0:150]
set xlabel "E"
set ylabel "R_n"
set title "R_n = f(E)"
set grid
plot "roznice.dat" i 0 u 1:2 w p lw 1 t "", \
	0 t ""


set xrange [0:1]
set out "zad_3.png"
set xlabel "r"
set ylabel "{/Symbol D}R(r)"
set title "Blad globalny"
plot "roznice.dat" i 1 u 1:2 w l lw 2 t "1. rozw.", \
	'' i 2 u 1:2 w l lw 2 t "2. rozw.", \
	'' i 3 u 1:2 w l lw 2 t "3. rozw.", \
	'' i 4 u 1:2 w l lw 2 t "4. rozw."

set xrange [0:1]
set out "test.png"
set xlabel "r"
set ylabel "R(r)"
set title "Test R(r) Numerov"
plot "numerov.dat" i 1 u 1:3 w l lw 2 t "1. rozw.", \
	'' i 2 u 1:3 w l lw 2 t "2. rozw.", \
	'' i 3 u 1:3 w l lw 2 t "3. rozw.", \
	'' i 4 u 1:3 w l lw 2 t "4. rozw."



set xrange [0:150]
set out "Numerov.png"
set xlabel "E"
set ylabel "R_n"
set title "R_n = f(E) (Numerov)"
plot "numerov.dat" i 0 u 1:2 w p lw 1 t "", \
	0 t ""

set xrange [0:1]
set out "Numerov_bledy.png"
set xlabel "r"
set ylabel "{/Symbol D}R(r)"
set title "Blad globalny (Numerov)"
plot "numerov.dat" i 1 u 1:2 w l lw 2 t "1. rozw.", \
	'' i 2 u 1:2 w l lw 2 t "2. rozw.", \
	'' i 3 u 1:2 w l lw 2 t "3. rozw.", \
	'' i 4 u 1:2 w l lw 2 t "4. rozw."