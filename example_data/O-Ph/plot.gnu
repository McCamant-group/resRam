set terminal png
set output "data/model.png"

a = "data/absorption.dat"
b = "O1_abs.dat"
#c = "data_o1_T298/absorption.dat"
c = "data/fluorescence.dat"

set parametric
stokes = 17320
unset key
set xrange [10000:22000]
set yrange [0:3]
p a u 1:2 w l,b u 2:1 w l , c u 1:2 w l, stokes,t

do for [i=0:18] {
	j = i+2
	k = i+1
	outfile = sprintf('data/mode%03.0f.png',i)
	set output outfile
	p "data/profiles.dat" u 1:(column(j)*1e7) w l,"crossSections.dat" u 20:(column(k)*1.e7)
}

#set output "raman_560.png"
#d = "output/raman_spec.dat"

#set xrange [0:1000]
#p d u 2 w l

