import astropy
import astroplan
import numpy as np
import astropy.units as u
import matplotlib.pyplot as plt
from astropy.time import Time
from astropy.utils import iers
from astropy.visualization import astropy_mpl_style, quantity_support
from astropy.coordinates import SkyCoord, EarthLocation, AltAz, get_moon, get_sun, Angle, ICRS
from astroplan import Observer
iers.conf.auto_download = True
plt.style.use(astropy_mpl_style)
quantity_support()

####################################################################################
# Enter Target Name

# TARGET NAME CANNOT CONTAIN SPACES

# Enter Target Coordinates (Right ascenion then Declination in DEGREES or in SEXAGESIMAL)
# Enter Transit Epoch in Julian Date
# Enter Target Period in Days
# Enter Transit Duration in Hours
# Enter Location to Save Graph to

# FILE LOCATION CANNOT CONTAIN SPACES
# FILE LOCATION MUST END IN '/'
####################################################################################

name = 'Hat-P-7b'
#tar = SkyCoord(331.117808, +81.565951, unit='deg')
tar = SkyCoord('19h28m59.35s 47d58m10.20s', unit=(u.hourangle, u.deg), frame='icrs')
#tar = SkyCoord('68.64564 +69.34811', unit=(u.hourangle, u.deg), frame='icrs')
epoch = 2459363.82778 #58765.6206 + 2400000.5   #t0 ref point for event such as a transit or a superior conjunction ; 58491.55337 + 2400000.5
period = 2.204737   # period of orbit in days
duration = 4.0398  #transit duration in hr
file_loc = '/d/www/brock/public_html/HatP7b/'
end_date =  2460610.50000 # look for transite between now and this end_date

####################################################################################

current_time = Time.now().jd

# Location of Observer - Change Latitude and Longitude (in DEGREES) to Change Observatory

wiro = EarthLocation(lat=41.0971*u.deg, lon=254.0235*u.deg, height=2943*u.m)
wiro_ob = Observer(longitude=254.0235*u.deg, latitude=41.0971*u.deg, elevation=2943*u.m, name="Wiro")

for i in range(0, 1600):

    next = epoch + i*period

    time = Time(float(next), format = 'jd')

    # Determine if Calculated Transit is After Today and Before End Date (Will give transits up to 10 days prior)

    if((next > (current_time - 10)) & (next < end_date) ):
        
        # Find Nearest Midnight from Given Transit Time

        midnight = Time(float(np.floor(next) + 0.5), format = 'jd')

        # Convert transit time to formated time

        transit = Time(float(next), format = 'jd')#, out_subfmt='longdate')
        transit_name = Time(float(next), format = 'jd')#, out_subfmt='longdate')

        # Change formatting of times

        time.format = 'isot'

        transit.format = 'iso'
        transit_name.format = 'fits'
        stransit_name=str(transit_name)[0:12]+'UTC'
        
        # Calculate Position of Target Object and Sun over the 48 Hours Around Tranist Time

        delta_midnight = np.linspace(-24, 24, 1000)*u.hour
        times_night = midnight + delta_midnight
        frame_night = AltAz(obstime=times_night, location=wiro)
        sunaltazs_night = get_sun(times_night).transform_to(frame_night)

        # Calculate Sunrise and Sunset Times

        sun_rise = wiro_ob.sun_rise_time(midnight, which = "nearest")
        sun_set = wiro_ob.sun_set_time(midnight, which = "nearest")

        # Calculate Morning and Evening Astronomical Twilight

        mor_twilight = wiro_ob.twilight_morning_astronomical(midnight, which = "nearest")
        eve_twilight = wiro_ob.twilight_evening_astronomical(midnight, which = "nearest")

        # If Sunrise Happens Before Sunset, Select Next Sunrise to Make Sure the Sunrise and Sunset Span the Same Night

        if(sun_rise - sun_set < 0):
            sun_rise = wiro_ob.sun_rise_time(midnight, which = "next")
            mor_twilight = wiro_ob.twilight_morning_astronomical(midnight, which = "next")

        # Only continue if transit happens after evening twilight and before morning

        if ((next > eve_twilight.jd) & (next < mor_twilight.jd) ) :
            
           # Convert Sunrise, Sunset, Morning Twilight, and Evening Twilight to UT Hours for Given Night

           sunrise = ((sun_rise.jd - midnight.jd) * 24)*u.hour
           sunset = ((sun_set.jd - midnight.jd) * 24)*u.hour

           mortwi = ((mor_twilight.jd - midnight.jd) * 24)*u.hour
           evetwi = ((eve_twilight.jd - midnight.jd) * 24)*u.hour

           # Calculate Position of Moon over the 48 Hours Around Tranist Time

           moon = get_moon(transit)
           moon_night = get_moon(times_night)
           moonaltazs_night = moon_night.transform_to(frame_night)

           # Calculate Illumination Percent of Moon on Given Night

           illumination = wiro_ob.moon_illumination(transit)*100

           # Calculate Distance Between Target and Moon

           separation = tar.separation(moon)
           separation = separation.degree

           # Convert Target Locations to Altitude as Azimuth from Given Location

           taraltazs_night = tar.transform_to(frame_night)

           # Convert Target Locations at Transit Midpoint to Airmass
           
           frame_midpoint = AltAz(obstime=transit, location=wiro)
           taraltazs_midpoint = tar.transform_to(frame_midpoint)
           tarairmass_midpoint = taraltazs_midpoint.secz

           # Calculate Transit Egress and Ingress

           addon = next*24-np.floor(next)*24

           if(addon + 12 > 24):
               addon = addon - 12
           else:
               addon = addon + 12

           ingress = addon - duration/2
           egress = addon + duration/2

           # Create Array with Ingress and Egress
           
           addon_array = np.array([addon-duration/2,addon+duration/2])

           # Convert Morning Twilight, Evening Twilight, Ingress, and Egress into Properly Formated Hours and Minutes with Leading zeroes before single minutes
        
           mor_hour = int(np.floor(mortwi.value))
           mor_minute = int(np.floor((mortwi.value - mor_hour)*60))
           mor_minute_str = str(mor_minute)
           mor_minute_str = mor_minute_str.zfill(2)

           eve_hour = int(np.floor(evetwi.value))
           eve_minute = int(np.floor((evetwi.value - eve_hour)*60))
           eve_minute_str = str(eve_minute)
           eve_minute_str = eve_minute_str.zfill(2)

           ingress_hour = int(np.floor(ingress))
           ingress_minute = int(np.floor((ingress - ingress_hour)*60))
           ingress_minute_str = str(ingress_minute)
           ingress_minute_str = ingress_minute_str.zfill(2)

           egress_hour = int(np.floor(egress))
           egress_minute = int(np.floor((egress - egress_hour)*60))
           egress_minute_str = str(egress_minute)
           egress_minute_str = egress_minute_str.zfill(2)

           # Create Array to Graph Ingress, Egress, Morning Twilight, and Evening Twilight

           tarr = [mortwi, evetwi, ingress*u.hour, egress*u.hour]

           # Create String to label Ingress, Egress, Morning Twilight and Evening Tiwlight
           
           tname = [str(mor_hour) + ':' + str(mor_minute_str) ,str(eve_hour) + ':' + str(eve_minute_str), str(ingress_hour) + ':' + str(ingress_minute_str) ,str(egress_hour) + ':' + str(egress_minute_str)]

           # Calculate Limits of Observability (Nautical Twilight)

           #if((next > 2459123.50000) & (next < 2459335.50000)):
           mor_observe = wiro_ob.twilight_morning_nautical(midnight, which = "next")
           eve_observe = wiro_ob.twilight_evening_nautical(midnight, which = "next")
           #else:
           #    mor_observe = wiro_ob.twilight_morning_nautical(midnight, which = "next")
           #    eve_observe = wiro_ob.twilight_evening_nautical(midnight, which = "nearest")

           # Convert Limits into UT Hours

           #if((next > 2459123.50000) & (next < 2459335.50000)):
           #    morobserve = ((mor_observe.jd - midnight.jd) * 24)
           #    eveobserve = ((eve_observe.jd - midnight.jd) * 24)
           #else:
           morobserve = ((mor_observe.jd - midnight.jd) * 24)
           eveobserve = ((eve_observe.jd - midnight.jd) * 24)

           # Create blank array

           observable = [0, 0]

           # Calculate Limits of Observability for Transit
    
           if(addon_array[0] > eveobserve):
               observable[0] = addon_array[0]
           else:
               observable[0] = eveobserve

           if(addon_array[1] < morobserve):
               observable[1] = addon_array[1]
           else:
               observable[1] = morobserve

           # Calculate Percentage of Transit Observable

           percent = (observable[1]-observable[0])/(addon_array[1]-addon_array[0])*100

           # Create string for file save location

           loc = str(file_loc) + '/' + str(name) + ' ' + str(transit) + '.pdf'
          
           # Plot Target, Moon, and Sun Altitudes Overlaid with Transit Light Curve

           fig, ax1 = plt.subplots()

           ax1.plot(delta_midnight+0*u.hour, sunaltazs_night.alt, color='r', label='Sun')
           ax1.plot(delta_midnight+0*u.hour, moonaltazs_night.alt, color=[0.5]*3, ls='--', label='Moon (Illumination {:.1f}%, Distance {:.1f} deg)'.format(illumination, separation))

           ax1.plot(delta_midnight+0*u.hour, taraltazs_night.alt, color = 'k', label='Target ({:.1f}% Observable, Midpoint Airmass {:.3f})'.format(percent, tarairmass_midpoint))

           ax1.fill_between(addon_array, 0*u.deg, 90*u.deg, color='#e8e8e8', zorder=0)

           ax1.axvline(mortwi+0*u.hour, linestyle = '--', color = 'k')
           ax1.axvline(evetwi+0*u.hour, linestyle = '--', color = 'k')

           ax1.legend(loc='upper left')
           ax1.set_xticks((np.arange(60)*1-12)*u.hour)
           ax1.axis([sunset+0*u.hour,sunrise+0*u.hour,0*u.deg,90*u.deg])
           plt.xlabel('Time (UTC)')
           plt.ylabel('Altitude [deg]')
           plt.title(name + ' - ' + str(transit))

           ax2 = ax1.twiny()
           ax2.set_xlim(ax1.get_xlim())
           ax2.set_xticks(tarr)
           ax2.set_xticklabels(tname)
           ax2.tick_params(length=0)
           plt.setp(ax2.get_xticklabels(), rotation=45, horizontalalignment='left')
           ax2.grid(False)
           
           fig.tight_layout()

####################################################################################

           plt.savefig(str(file_loc) + str(name) + '_' + stransit_name + '.pdf')      # Uncomment if you want to save the pdf graph to the provided location
          #plt.show()     # Uncomment if you want to simply view the graph rather than save it

####################################################################################

           ax1.clear()
           ax2.clear()
           plt.clf()
           plt.close()

