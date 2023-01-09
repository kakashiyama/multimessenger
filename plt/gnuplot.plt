##############################################
reset
se ter postscript eps enhanced
se log
se xlabel "t [sec]"
se ylabel "radius [cm]"
se output "radius.eps"
plot \
"../op/hydro.dat" usi ($1):($3) w l title "r_{ej}",\
"../op/hydro.dat" usi ($1):($5) w l title "r_w"

reset
se ter postscript eps enhanced
se log
se xlabel "t [sec]"
se ylabel "ejecta velocity [cm/s]"
se output "v_ej.eps"
plot \
"../op/hydro.dat" usi ($1):($4) w l title "v_{ej}",\
"../op/hydro.dat" usi ($1):($6) w l title "v_w"

reset
se ter postscript eps enhanced
se log
se xlabel "t [sec]"
se ylabel "Thomson optical depth"
se output "tau.eps"
plot \
"../op/hydro.dat" usi ($1):($7) w l title ""

reset
se ter postscript eps enhanced
se log
se xlabel "t [sec]"
se ylabel "ejecta temperature [K]"
se output "T_ej.eps"
plot \
"../op/hydro.dat" usi ($1):($8) w l title ""

reset
se ter postscript eps enhanced
se log
se xlabel "t [sec]"
se ylabel "PWN magnetic field [G]"
se output "B_pwn.eps"
plot \
"../op/hydro.dat" usi ($1):($9) w l title ""

##############################################

reset
se ter postscript eps enhanced
se log
se yrange [1.0e30:*]
se xlabel "t [sec]"
se ylabel "bolometric luminosity [erg/s]"
se output "lc_SN.eps"
plot \
"../op/lc.dat" usi ($1):($2) w l title "SN",\
"../op/lc.dat" usi ($1):($4+$5) w l title "PSR",\
"../op/lc.dat" usi ($1):($6+$7) w l title "Radioactive nuclei",\
"../op/lc.dat" usi ($1):($3) w l title "GW"

reset
se ter postscript eps enhanced
se log
se yrange [1.0e30:*]
se xlabel "t [sec]"
se ylabel "hard-X-ray luminosity [erg/s]"
se output "lc_hardx.eps"
plot \
"../op/lc.dat" usi ($1):($8) w l title ""

reset
se ter postscript eps enhanced
se log
se yrange [1.0e30:*]
se xlabel "t [sec]"
se ylabel "gamma-ray luminosity [erg/s]"
se output "lc_gamma.eps"
plot \
"../op/lc.dat" usi ($1):($9) w l title ""

reset
se ter postscript eps enhanced
se log
se yrange [1.0e-10:2]
se xlabel "t [sec]"
se ylabel "Escape fraction"
se output "f_esc.eps"
plot \
"../op/lc.dat" usi ($1):($13) w l title "X",\
"../op/lc.dat" usi ($1):($15) w l title "gamma"

reset
se ter postscript eps enhanced
se log
se xlabel "t [sec]"
se ylabel "Spectral normalization factor"
se output "spec_norm.eps"
plot \
"../op/lc.dat" usi ($1):($14) w l title "X",\
"../op/lc.dat" usi ($1):($16) w l title "gamma"