set term png size 1600, 1000
# set term png size 800, 600
# set size square
# set size 1.,1.
set out 'wyniki.png'
set key top right
set grid
set title "R_n = f(E)"
set multiplot layout 2,3 rowsfirst

	# set out 'M5_E.png'
	set xrange [0.4:2]
	set xlabel "{/Symbol a}"
	set ylabel "E_{/Symbol m}"
	set title 'M=5'
	plot for [idx=1:5] 'M5.dat' i 0 u 1:idx+1 w l lw 2 t '{/Symbol m} = '.(idx)

	# set out 'M10_E.png'
	set title 'M=10'
	plot for [idx=1:5] 'M10.dat' i 0 u 1:idx+1 w l lw 2 t '{/Symbol m} = '.(idx)
	
	# set out 'M30_E.png'
	set title 'M=30'
	set yrange [0:6.5]
	plot for [idx=1:5] 'M30.dat' i 0 u 1:idx+1 w l lw 2 t '{/Symbol m} = '.(idx)
	
	set yrange [*:*]
	# set out 'M5_u.png'
	set xrange [-6:6]
	set xlabel "x"
	set ylabel "u_{/Symbol m}(x)"
	set title 'M=5'
	set key bottom right
	plot for [idx=1:5] 'M5.dat' i idx u 1:2 w l lw 2 t '{/Symbol m} = '.(idx)

	# set out 'M10_u.png'
	set title 'M=10'
	set key top right
	plot for [idx=1:5] 'M10.dat' i idx u 1:2 w l lw 2 t '{/Symbol m} = '.(idx)
	
	# set out 'M30_u.png'
	set key bottom right
	set title 'M=30'
	plot for [idx=1:5] 'M30.dat' i idx u 1:2 w l lw 2 t '{/Symbol m} = '.(idx)
	

unset multiplot


# set term png size 800, 600
# set out 'M5_E.png'
# set xrange [0.4:2]
# set xlabel "{/Symbol a}"
# set ylabel "E_{/Symbol m}"
# set title 'M=5'
# set key top right
# plot for [idx=1:5] 'M5.dat' i 0 u 1:idx+1 w l lw 2 t '{/Symbol m} = '.(idx), \
# 	for [idx=0:4] idx+0.5 t ""

# set out 'M10_E.png'
# set title 'M=10'
# plot for [idx=1:5] 'M10.dat' i 0 u 1:idx+1 w l lw 2 t '{/Symbol m} = '.(idx), \
# 	for [idx=0:4] idx+0.5 t ""

# set out 'M30_E.png'
# set title 'M=30'
# set yrange [0:6.5]
# plot for [idx=1:5] 'M30.dat' i 0 u 1:idx+1 w l lw 2 t '{/Symbol m} = '.(idx), \
# 	for [idx=0:4] idx+0.5 t ""


# set yrange [*:*]
# set xrange [-6:6]
# set xlabel "x"
# set ylabel "u_{/Symbol m}(x)"

# do for [idx = 1:5] {
# 	outfile = sprintf('mu_%1.0f.png',idx)
#   	set output outfile

# 	set title '{/Symbol m} = '.(idx)
# 	plot 'M5.dat' i idx u 1:2 w l lw 2 t 'M=5', \
# 		'M10.dat' i idx u 1:2 w l lw 2 t 'M=10', \
# 		'M30.dat' i idx u 1:2 w l lw 2 t 'M=30', \
# 		'u_dokl.dat' i 0 u 1:idx+1 w l lw 2 t 'dokl.'
# }