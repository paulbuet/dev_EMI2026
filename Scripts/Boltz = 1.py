Boltz      = 1.380658E-23
Avogadro   = 6.0221367E+23

XTT   = 273.16
EsTt  = 611.14
LvTt  = 2.5008E+6
Cl    = 4.218E+3
Mv    = 18.0153E-3
Md    = 28.9644E-3
XRV   = Avogadro * Boltz / Mv
XRD   = Avogadro * Boltz / Md
Cpv   = 4.* Rv

xgamw  = (Cl - Cpv) / XRV
xbetaw = (LvTt / XRV) + (xgamw * XTT)
xalpw  = np.log(EsTt) + (xbetaw / XTT) + (xgamw * np.log(XTT))

print(xalpw)


