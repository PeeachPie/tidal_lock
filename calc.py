PI = 3.14159265358979323846
G  = 6.6743e-11
MOON_EARTH_DIST        = 384400000 / 15
EARTH_MASS             = 5.972e24
EARTH_RADIUS           = 6371000
MOON_ANGULAR_VELOCITY  = 2.665e-6 * 150
MOON_VELOCITY          = 1023

def Ml(Me, R, r, w):
    ml = 1.0 / ((G / (2 * w * w * R)) * ((1.0 / (R + r)**2) + (1.0 / (R - r)**2))) - Me
    b = (ml * R) / (ml + Me)
    print(r/R, b / R, b)
    return ml

ml = Ml(EARTH_MASS, MOON_EARTH_DIST, EARTH_RADIUS, MOON_ANGULAR_VELOCITY)
print(ml / EARTH_MASS, ml)