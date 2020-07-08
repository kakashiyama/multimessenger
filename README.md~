**multi_messenters.c (ver.1.7.1)**
  
Original ver. created by Kazumi Kashiyama on 13/11/27.
Updated by Kohta Murase on 17/10/08.

In this ver. the energy depoisiton rate from gamma rays to the thermal bath is modeled in more detailed.
Note that albedo should change between magnetar and BH, epsilon_B can be 0.01 or 0.003 etc., note that one should change eps_e for BH or NS.
Be careful about the usage of this code when the spin-down time is too short.

Murase et al. 2019; Adiabatic loss of magnetic field in the PWN is implemented. Be careful how to set "timebin".
See line ~410.

---

**input**
<list of input parameters>
1. NS mass [MSUN]
2. NS radius [cm]
3. NS toloidal magnetic field [G]
4. NS poloidal magnetic field [G]
5. NS initial spin period [ms]
6. SN explosion energy [erg]
7. SN ejecta mass [MSUN]
8. Ni mass [Msun]
9. Rprocess-element mass [Msun]
10. SN explosion inital radius [cm] 
11. power-law index of SN-ejecta density profile : 9./4 in mergers, 1.0 in magnetar 
12. the albedo factor of x-rays : 0 or 0.5 for high-ionization
13. opacity in unit of [g^-1 cm^2]
14. luminosity distance to SN [pc]
