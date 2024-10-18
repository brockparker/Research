'''Script to calculate the angular and physical separations between target stars.'''

import astropy.units as u
from astropy.coordinates import SkyCoord
# BP Required imports.

zeta_oph = SkyCoord(ra = 249.2898*u.degree, dec = -10.5670*u.degree, distance = 135*u.pc)
# BP Defining Zeta Ophiuchi's coordinates.
pk04 = SkyCoord(ra = 249.1036*u.degree, dec = -10.5442*u.degree, distance = 223*u.pc)
pk13 = SkyCoord(ra = 249.2367*u.degree, dec = -10.4871*u.degree, distance = 449*u.pc)
pk24 = SkyCoord(ra = 249.4317*u.degree, dec = -10.5519*u.degree, distance = 202*u.pc)
pk26 = SkyCoord(ra = 249.5049*u.degree, dec = -10.5797*u.degree, distance = 669*u.pc)
pk27 = SkyCoord(ra = 249.5108*u.degree, dec = -10.4688*u.degree, distance = 499*u.pc)
# BP Specifying target coordinates.

sep04 = zeta_oph.separation(pk04)
sep13 = zeta_oph.separation(pk13)
sep24 = zeta_oph.separation(pk24)
sep26 = zeta_oph.separation(pk26)
sep27 = zeta_oph.separation(pk27)
# BP Calculating angular separation between Zeta Oph and target stars.

print('The angular separation between $\zeta$ Oph and PK-04 is {:.3f}.'.format(sep04.to(u.arcmin)))
print('The angular separation between $\zeta$ Oph and PK-13 is {:.3f}.'.format(sep13.to(u.arcmin)))
print('The angular separation between $\zeta$ Oph and PK-24 is {:.3f}.'.format(sep24.to(u.arcmin)))
print('The angular separation between $\zeta$ Oph and PK-26 is {:.3f}.'.format(sep26.to(u.arcmin)))
print('The angular separation between $\zeta$ Oph and PK-27 is {:.3f}.\n'.format(sep27.to(u.arcmin)))

dist04 = sep04.to(u.radian) * zeta_oph.distance
dist13 = sep13.to(u.radian) * zeta_oph.distance
dist24 = sep24.to(u.radian) * zeta_oph.distance
dist26 = sep26.to(u.radian) * zeta_oph.distance
dist27 = sep27.to(u.radian) * zeta_oph.distance
# BP Calculating projected linear distances at the distance of Zeta Oph.

print('The projected linear separation between $\zeta$ Oph and PK-04 is {:.3f}.'.format(dist04.to(u.pc, equivalencies=u.dimensionless_angles())))
print('The projected linear separation between $\zeta$ Oph and PK-13 is {:.3f}.'.format(dist13.to(u.pc, equivalencies=u.dimensionless_angles())))
print('The projected linear separation between $\zeta$ Oph and PK-24 is {:.3f}.'.format(dist24.to(u.pc, equivalencies=u.dimensionless_angles())))
print('The projected linear separation between $\zeta$ Oph and PK-26 is {:.3f}.'.format(dist26.to(u.pc, equivalencies=u.dimensionless_angles())))
print('The projected linear separation between $\zeta$ Oph and PK-27 is {:.3f}.'.format(dist27.to(u.pc, equivalencies=u.dimensionless_angles())))