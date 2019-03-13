set terminal png
set output "data/model.png"

a = "data/absorption.dat"
b = "O2_abs.dat"
c = "data/fluorescence.dat"

unset key
set parametric
stokes = 16920
set xrange [10000:22000]
set yrange [0:4]
p a u 1:2 w l,b u 1:2 w l, c u 1:2 w l, stokes, t

do for [i=1:26] {
	j = i+1
	k = i
	outfile = sprintf('data/mode%03.0f.png',i)
	set output outfile
	p "data/profiles.dat" u 1:(column(j)*1e7) w l,"crossSections.dat" u 27:(column(k)*1.e7)
}
#do for [i=0:3] {
#	j = i+2
#	outfile = sprintf('data/mode%03.0f.png',i)
#	set output outfile
#	p "data/profiles.dat" u 1:(column(j)*1e7) w l
#}

#set output "raman_560.png"
#d = "output/raman_spec.dat"

#set xrange [0:1000]
#p d u 2 w l

