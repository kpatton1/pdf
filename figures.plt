

set term png

set output 'ps_compare.png'

set logscale
set key bot right

set xlabel 'k [1/Mpc]'
set ylabel 'delta^2'

plot 'test_ps_lin.dat' using 1:2 w l title 'Linear PS', 'test_ps_nonlin.dat' using 1:2 w l title 'Nonlinear PS'

set output 'ps_var_compare.png'

plot 'test_ps_lin_var.dat' using 1:2 w l title 'Linear PS', 'test_ps_nonlin_var.dat' using 1:2 w l title 'Nonlinear PS'

set ylabel '2d P(k) = delta^2 / (2*pi*k^2)'

set output 'ps_2d.png'

plot 'test_ps_lin.dat' using 1:($2/$1/$1/(2*pi)) w l title 'Linear PS', 'test_ps_nonlin.dat' using 1:($2/$1/$1/(2*pi)) w l title 'Nonlinear PS'

set output 'ps_3d.png'

set ylabel '3d P(k) = delta^2 / (4*pi*k^3)'

plot 'test_ps_lin.dat' using 1:($2/$1/$1/$1/(4*pi)) w l title 'Linear PS', 'test_ps_nonlin.dat' using 1:($2/$1/$1/$1/(4*pi)) w l title 'Nonlinear PS'

set ylabel 'sigma'

plot 'test_ps_lin_var.dat' using 1:2 w l title 'Linear variance', 'test_ps_nonlin_var.dat' using 1:2 w l title 'Nonlinear variance'

set xlabel 'log(projected density)'
set ylabel 'PDF(log(projected density))'

set key top right

unset logscale

set output 'pdf.png'

plot 'pdf.dat' using 1:2 w l title 'Projected density PDF (1 - 10^15 Msolar halos)'

set xlabel 'log(projected density)'
set ylabel 'Integrated mass fraction'

set output 'pdf_integrated_mass.png'

plot 'pdf_integrated_mass.dat' using 1:($2*0.887552) w l title 'Integrated mass (1 - 10^15 Msolar halos)'

set xlabel 'log(projected density)'
set ylabel 'Integrated variance'

set output 'pdf_var.png'

plot 'pdf_var.dat' using 1:2 w l title 'Integrated variance (1 - 10^15 Msolar halos)'

set xlabel 'log(M/Msolar)'
set ylabel 'dN/dlogM [1/Mpc^3]'

set logscale y

set output 'mf.png'

set xrange [0:18]
set yrange [1e-20:1e0]

plot 'test_mf.dat' using 2:($3/exp($2)) w l title 'Halo mass function'

set output 'mf_mass.png'

unset logscale
set autoscale

set xlabel 'log(M/Msolar)'
set ylabel 'M dN/dlogM'

plot 'test_mf.dat' using 2:3 w l title 'Fractional mass in halos'



