-- Orb.lua 
-- partial port of orb.js into lua lang
-- MIT License / Isana, Kashiwai 2021
-- update 2021-07-10 01:00 JST

-- name space

Orb = {}

Orb.Storage = {}

-- tools
function Orb.RoundNum(num,max)
  results = num % max
  if (results < 0) then
    results = results + max
  end
  return results;
end

function Orb.RoundAngle(num)
 return Orb.RoundNum(num,360)
end

function Orb.ZeroFill(num, length) 
  if (length)  then
    length = length;
  else
    length = string.len(num);
  end
  str = "0000000000" .. num
  str_length = string.len(str)
  return string.sub(str,str_length-length+1,str_length)
end

function Orb.Normalize(vec)
  local d = math.sqrt(vec.x*vec.x + vec.y*vec.y + vec.z*vec.z)
  return {
    x = vec.x/d,
    y = vec.y/d,
    z = vec.z/d
  }
end


Orb.Time = {}


function Orb.Time.TT(date)
  local delta_t = (37 + 32)
  local tt_in_sec = os.time(date) + delta_t
  return os.date("*t",tt_in_sec)
end


-- Jurian day

function Orb.Time.JD(date)
  local year = date.year;
  local month = date.month;;
  local day = date.day;
  local time_in_day = date.hour / 24 + date.min / 1440 + date.sec / 86400;
  if month <= 2 then
    year = year - 1;
    month = month + 12;
  end
  local julian_day = math.floor(365.25 * (year + 4716)) + math.floor(30.6001 * (month + 1)) + day - 1524.5;
  local transition_offset = 0;
  if julian_day < 2299160.5 then
    transition_offset = 0;
  else
    local tmp = math.floor(year / 100);
    transition_offset = 2 - tmp + math.floor(tmp / 4);
  end
  return julian_day + transition_offset + time_in_day;
end

function Orb.Time.JDToUTC(jd)
  local mjd = jd - 2400000.5
  local mjd0 = math.floor(mjd)
  local flac = mjd - mjd0;
  local n = mjd0 + 678881
  local a = 4*n + 3 + 4 * math.floor((3/4)*math.floor(((4*(n+1))/(146097))+1))
  local b = 5 * math.floor((a%1461)/4) + 2
  local year = math.floor(a/1461)
  local month = math.floor(b/153) + 3
  local day = math.floor((b%153)/5) +1

  if month == 13 then
    year = year + 1
    month = 1
  end

  if month == 14 then
    year = year + 1
    month = 2
  end  
  local h = flac*24
  local hour = math.floor(h)
  local m = (h - hour)*60
  local min = math.floor(m)
  local sec = math.floor((m - min)*60)
  
  return {
    year = year,
    month = Orb.ZeroFill(month,2),
    day = Orb.ZeroFill(day,2),
    hour = Orb.ZeroFill(hour,2),
    min = Orb.ZeroFill(min,2),
    sec = Orb.ZeroFill(sec,2)
  }
end

-- Greenwich Apparent Sidereal Time
function Orb.Time.gst(date) 
  local rad = math.pi/180;
  local time_in_sec = date.hour * 3600 + date.min * 60 + date.sec
  local time_in_day = date.hour / 24 + date.min / 1440 + date.sec / 86400;
  local jd = Orb.Time.JD(date);
  local jd0 = jd - time_in_day;

  -- Greenwich Mean Sidereal Time(GMST) at 0:00
  local t = (jd0 - 2451545.0) / 36525;
  local gmst_at_zero = (24110.5484 + 8640184.812866 * t + 0.093104 * t * t + 0.0000062 * t * t * t) / 3600;
  if gmst_at_zero > 24 then
     gmst_at_zero = gmst_at_zero % 24; 
  end
  -- GMST at target time
  local gmst = gmst_at_zero + (time_in_sec * 1.00273790925) / 3600;
  -- mean obliquity of the ecliptic
  local e = 23 + 26.0 / 60 + 21.448 / 3600 - 46.8150 / 3600 * t - 0.00059 / 3600 * t * t + 0.001813 / 3600 * t * t * t;
  -- Nutation in longitude
  local omega = 125.04452 - 1934.136261 * t + 0.0020708 * t * t + t * t * t / 450000;
  local long1 = 280.4665 + 36000.7698 * t;
  local long2 = 218.3165 + 481267.8813 * t;
  local phai = -17.20 * math.sin(omega * rad) - (-1.32 * math.sin(2 * long1 * rad)) - 0.23 * math.sin(2 * long2 * rad) + 0.21 * math.sin(2 * omega * rad);
  -- Greenwich Apparent Sidereal Time(GAST/GST)
  local gast = gmst + ((phai / 15) * (math.cos(e * rad))) / 3600
  gast = Orb.RoundNum(gast,24)
  return gast
end

Orb.Coord = {}

-- Equatorial Spherical(ra,dec) to Equatorial Rectangular(x,y,z)
Orb.Coord.RadecToEquatorial = function(radec)
  local rad = math.pi/180;
  local ra_deg = radec.ra * 15;
  local dec = radec.dec
  local distance = radec.distance
  local vec = {
    x = distance * math.cos(dec * rad) * math.cos(ra_deg * rad),
    y = distance * math.cos(dec * rad) * math.sin(ra_deg * rad),
    z = distance * math.sin(dec * rad)
  }
  return vec
end
Orb.Coord.RadecToXYZ = Orb.Coord.RadecToEquatorial


Orb.Coord.EquatorialToRadec = function (date,rect) 
  -- equatorial rectangular(x,y,z) to spherical(ra,dec)
  local rad = math.pi/180;
  local eqx = rect.x;
  local eqy = rect.y;
  local eqz = rect.z;
  local ra = math.atan2(eqy, eqx) / rad;
  if (ra < 0)  then
    ra = ra % 360 + 360
  end
  if (ra > 360) then
    ra = ra % 360
  end
  ra = ra / 15
  local dec = math.atan2(eqz, math.sqrt(eqx * eqx + eqy * eqy)) / rad;
  local distance = math.sqrt(eqx * eqx + eqy * eqy + eqz * eqz);
  return {
    ra = ra,
    dec = dec,
    distance = distance,
    date = date
  };
end

Orb.Coord.EquatorialToEcliptic = function (date,vec) 
  -- equatorial rectangular(x,y,z) to ecliptic rectangular(x,y,z)
  local rad = math.pi/180;
  local nao = Orb.Coord.NutationAndObliquity(date)
  local obliquity = nao.obliquity
  local equatorial = vec
  local ecliptic = {
    x = equatorial.x,
    y = math.cos(obliquity * rad) * equatorial.y + math.sin(obliquity * rad) * equatorial.z,
    z = -math.sin(obliquity * rad) * equatorial.y + math.cos(obliquity * rad) * equatorial.z
  }
  return {
    x = ecliptic.x,
    y = ecliptic.y,
    z = ecliptic.z,
    date = date
  }
end

Orb.Coord.EclipticToEquatorial = function (date,vec) 
  -- ecliptic rectangular(x,y,z) to equatorial rectangular(x,y,z)
  local rad = math.pi/180
  local ecliptic = vec
  local ep
  if Orb.Storage.earth then
    ep = Orb.Storage.earth
  else
    ep = Orb.Planet.Earth(date)
  end
  local gcx = ecliptic.x - ep.x;
  local gcy = ecliptic.y - ep.y;
  local gcz = ecliptic.z - ep.z;
  local nao = Orb.Coord.NutationAndObliquity(date)
  local ecl = nao.obliquity
  local equatorial = {
    x = gcx,
    y = gcy * math.cos(ecl * rad) - gcz * math.sin(ecl * rad),
    z = gcy * math.sin(ecl * rad) + gcz * math.cos(ecl * rad)
  }
  return {
    x = equatorial.x,
    y = equatorial.y,
    z = equatorial.z,
    date = date
  }
end

function Orb.Coord.NutationAndObliquity(date)
  local rad = math.pi/180;
  local tt = Orb.Time.TT(date)
  local jd = Orb.Time.JD(tt);
  local t = (jd - 2451545.0) / 36525;
  local omega = 125.04452 - 1934.136261 * t + 0.0020708 * t * t + (t * t * t / 450000);
  local L0 = 280.4665 + 36000.7698 * t;
  local L1 = 218.3165 + 481267.8813 * t;
  local nutation = (-17.20 / 3600) * math.sin(omega * rad) - (-1.32 / 3600) * math.sin(2 * L0 * rad) - (0.23 / 3600) * math.sin(2 * L1 * rad) + (0.21 / 3600) * math.sin(2 * omega * rad);
  local mean_obliquity = 23 + 26.0 / 60 + 21.448 / 3600 - (46.8150 / 3600) * t - (0.00059 / 3600) * t * t + (0.001813 / 3600) * t * t * t;
  local obliquity_delta = (9.20 / 3600) * math.cos(omega * rad) + (0.57 / 3600) * math.cos(2 * L0 * rad) + (0.10 / 3600) * math.cos(2 * L1 * rad) - (0.09 / 3600) * math.cos(2 * omega * rad);
  local obliquity = mean_obliquity + obliquity_delta;
  Orb.Storage.nutation = nutation
  Orb.Storage.obliquity = obliquity
  return {
    nutation = nutation,
    obliquity = obliquity
  }
end

-- Kepler orbit

Orb.Kepler = function(date, elements)
  local gm
  if (elements.gm) then
    gm = elements.gm;
  else
    gm = 2.9591220828559093*math.pow(10,-4);
  end

  local epoch
  if (elements.time_of_periapsis) then
    epoch = elements.time_of_periapsis;
  else
    epoch = elements.epoch;
  end

  if (elements.perihelion_distance) then
    elements.periapsis_distance = elements.perihelion_distance;
  end
  
  local rad = math.pi/180

  local EllipticalOrbit = function ()
    local eccentricity = elements.eccentricity;
    local semi_major_axis
    if (elements.semi_major_axis) then
      semi_major_axis = elements.semi_major_axis;
    elseif (elements.periapsis_distance) then
      semi_major_axis = (elements.periapsis_distance) / (1 - eccentricity)
    end
    local mean_motion = math.sqrt(gm / (semi_major_axis * semi_major_axis * semi_major_axis)) / rad;
    local jd = Orb.Time.JD(date)
    local elapsed_time = jd - epoch;
    local mean_anomaly,l
    if (elements.mean_anomaly) then
      mean_anomaly = elements.mean_anomaly;
      l = (mean_motion * elapsed_time) + mean_anomaly;
    elseif(elements.time_of_periapsis) then
      mean_anomaly = mean_motion * elapsed_time;
      l = mean_anomaly;
    end
    if (l > 360) then
      l = l % 360
    end
    l = l * rad
    local u = l
    local i = 0;
    local ut,delta_u
    repeat
      ut = u;
      delta_u = (l - u + (eccentricity * math.sin(u))) / (1 - (eccentricity * math.cos(u)));
      u = u + delta_u;
      if (i > 1000000) then 
        break;
      end
      i = i+1;
    until (math.abs(ut - u) > 0.0000001) 
    local eccentric_anomaly = u;
    local p = math.abs(semi_major_axis * (1 - eccentricity * eccentricity))
    local true_anomaly = 2 * math.atan(math.sqrt((1 + eccentricity) / (1 - eccentricity)) * math.tan(eccentric_anomaly / 2));
    local r = p / (1 + eccentricity * math.cos(true_anomaly));
    local orbital_plane = {
      r = r,
      x = r * math.cos(true_anomaly),
      y = r * math.sin(true_anomaly),
      xdot = -math.sqrt(gm / p) * math.sin(true_anomaly),
      ydot = math.sqrt(gm / p) * (eccentricity + math.cos(true_anomaly))
    };
    return orbital_plane
  end
  
  local orbital_plane = EllipticalOrbit()

  EclipticRectangular = function ()
    local lan = elements.longitude_of_ascending_node * rad;
    local ap = elements.argument_of_periapsis * rad;
    local inc = elements.inclination * rad;
    local op2xyz = function (opx, opy, lan, ap, inc)
      return {
        x = opx * (math.cos(lan) * math.cos(ap) - math.sin(lan) * math.cos(inc) * math.sin(ap)) - opy * (math.cos(lan) * math.sin(ap) + math.sin(lan) * math.cos(inc) * math.cos(ap)),
        y = opx * (math.sin(lan) * math.cos(ap) + math.cos(lan) * math.cos(inc) * math.sin(ap)) - opy * (math.sin(lan) * math.sin(ap) - math.cos(lan) * math.cos(inc) * math.cos(ap)),
        z = opx * math.sin(inc) * math.sin(ap) + opy * math.sin(inc) * math.cos(ap)
      }
    end
    local vec = op2xyz(orbital_plane.x, orbital_plane.y, lan, ap, inc)
    local dotvec = op2xyz(orbital_plane.xdot, orbital_plane.ydot, lan, ap, inc)
    return {
      x = vec.x,
      y = vec.y,
      z = vec.z,
      xdot = dotvec.x,
      ydot = dotvec.y,
      zdot = dotvec.z
    };
  end

  local xyz = EclipticRectangular()
  return xyz
end

Orb.Observe = {}

Orb.Observe.AtmosphericRefraction = function(elevation)
  local rad = math.pi/180
  local tmp = elevation+7.31/(elevation + 4.4)
  local ar = 0.0167*rad/(math.tan(tmp*rad))/rad
  return ar
end

Orb.Observe.LatLngToRect = function(date,latlng)
  local rad = math.pi/180
  local lat = latlng.latitude;
  local lng = latlng.longitude;
  local gst = Orb.Time.gst(date);
  local lst = gst*15 + lng;
  local a = 6378.135 + latlng.altitude;
  local f = 0.00335277945
  local sin_lat =math.sin(lat*rad);
  local c = 1/math.sqrt(1+f*(f-2)*sin_lat*sin_lat);
  local s = (1-f)*(1-f)*c;
  return {
    x = a*c*math.cos(lat*rad)*math.cos(lst*rad),
    y = a*c*math.cos(lat*rad)*math.sin(lst*rad),
    z = a*s*math.sin(lat*rad)
  }
end

Orb.Observe.RadecToHorizontal = function(date,radec,observer)
  local rad = math.pi/180
  local ra = radec.ra;
  local dec = radec.dec;
  local distance = radec.distance;
  local latitude = observer.latitude;
  local longitude = observer.longitude;
  local altitude = observer.altitude;
  dec = dec*rad
  local gst = Orb.Time.gst(date);
  local hour_angle = (gst*15 + longitude - (ra*15));
  local h = hour_angle*rad;
  local lat = latitude*rad;
  local azimuth = (math.atan2(-math.cos(dec)*math.sin(h),math.sin(dec)*math.cos(lat)-math.cos(dec)*math.sin(lat)*math.cos(h)))/rad;
  azimuth = Orb.RoundNum(azimuth,360)
  local elevation = (math.asin(math.sin(dec)*math.sin(lat)+math.cos(lat)*math.cos(dec)*math.cos(h)))/rad;
  local atmospheric_refraction = Orb.Observe.AtmosphericRefraction(elevation)
  if (azimuth<0) then
    azimuth = azimuth%360 + 360
  end
  return {
    azimuth = azimuth,
    elevation = elevation + atmospheric_refraction,
    distance = distance,
    atmospheric_refraction = atmospheric_refraction
   }
end

Orb.Observe.EclipticToHorizontal = function(date,ecliptic,observer)
  local equatorial = Orb.Coord.EclipticToEquatorial(date,ecliptic)
  local radec = Orb.Coord.EquatorialToRadec(date,equatorial)
  local horizontal = Orb.Observe.RadecToHorizontal(date,radec,observer)
  return horizontal
end

Orb.Observe.EquatorialToHorizontal = function(date,equatorial,observer)
  local radec = Orb.Coord.EquatorialToRadec(date,equatorial)
  local horizontal = Orb.Observe.RadecToHorizontal(date,radec,observer)
  return horizontal
end
 
-- Sun

Orb.Sun = {}

function Orb.Sun.EclipticLongitude(date)
  local rad = math.pi/180;
  local tt = Orb.Time.TT(date)
  local jd = Orb.Time.JD(tt)
  local t = (jd - 2451545.0) / 36525;
  local mean_longitude = 280.46646 + 36000.76983 * t + 0.0003032 * t * t;
  local mean_anomaly = 357.52911 + 35999.05029 * t - 0.0001537 * t * t;
  local eccentricity = 0.016708634 - 0.000042037 * t - 0.0000001267 * t * t;
  local equation = (1.914602 - 0.004817 * t - 0.000014 * t * t) * math.sin(mean_anomaly * rad);
  equation = equation + (0.019993 - 0.000101 * t) * math.sin(2 * mean_anomaly * rad);
  equation = equation + 0.000289 * math.sin(3 * mean_anomaly * rad);
  local true_longitude = mean_longitude + equation;
  local true_anomary = mean_anomaly + equation;
  local radius = (1.000001018 * (1 - eccentricity * eccentricity)) / (1 + eccentricity * math.cos(true_anomary * rad));

  local nao = Orb.Coord.NutationAndObliquity(date)
  local nutation = nao.nutation;
  local obliquity = nao.obliquity;
  local apparent_longitude = true_longitude + nutation;
  local longitude = apparent_longitude;
  local distance = radius * 149597870.691;
  return {
    longitude = longitude,
    distance = distance,
    obliquity = obliquity
  }
end

Orb.Sun.SphericalToRectangler = function (spherical)
  local rad = math.pi/180;
  local longitude = spherical.longitude;
  local distance = spherical.distance;
  local obliquity = spherical.obliquity;
  return {
    x = distance * math.cos(longitude * rad),
    y = distance * (math.sin(longitude * rad) * math.cos(obliquity * rad)),
    z = distance * (math.sin(longitude * rad) * math.sin(obliquity * rad))
  }
end

Orb.Sun.EclipticToEquatorial = function (ecliptic)
  local rad = math.pi/180;
  local longitude = ecliptic.longitude
  local distance = ecliptic.distance
  local obliquity = ecliptic.obliquity
  local ra = math.atan2(math.cos(obliquity * rad) * math.sin(longitude * rad), math.cos(longitude * rad))
  ra = Orb.RoundNum(ra/rad,360)
  ra = ra / 15
  local dec = math.asin(math.sin(obliquity * rad) * math.sin(longitude * rad));
  dec = dec / rad;
  return {
    ra = ra,
    dec = dec,
    distance = distance
  }
end

Orb.Sun.Radec = function(date)
  local ecliptic = Orb.Sun.EclipticLongitude(date)
  local radec = Orb.Sun.EclipticToEquatorial(ecliptic)
  return radec
end

Orb.Sun.Equatorial = function(date)
  local ecliptic = Orb.Sun.EclipticLongitude(date)
  local rect = Orb.Sun.SphericalToRectangler(ecliptic)
  return rect
end

Orb.Sun.Ecliptic = function(date)
  return {x=0,y=0,z=0}
end

Orb.Sun.Horizontal = function(date,observer)
  local radec = Orb.Sun.Radec(date)
  return Orb.Observe.RadecToHorizontal(date,radec,observer)
end


-- Moon

Orb.Luna = require("luna")

Orb.Moon = function(date)
  return Orb.Luna.equatorial(date);
end

-- planets

Orb.VSOP = require("vsop")

Orb.Planet = {}

Orb.Planet.Mercury = function(date)
  local terms = Orb.VSOP.Mercury();
  return Orb.VSOP.Exec(date,terms)
end


Orb.Planet.Venus = function(date)
  local terms = Orb.VSOP.Venus();
  return Orb.VSOP.Exec(date,terms)
end


Orb.Planet.Earth = function(date)
  local terms = Orb.VSOP.Earth();
  Orb.Storage.earth = Orb.VSOP.Exec(date,terms)
  return Orb.Storage.earth
end


Orb.Planet.Mars = function(date)
  local terms = Orb.VSOP.Mars();
  return Orb.VSOP.Exec(date,terms)
end

Orb.Planet.Jupiter = function(date)
  local terms = Orb.VSOP.Jupiter();
  return Orb.VSOP.Exec(date,terms)
end

Orb.Planet.Saturn = function(date)
  local terms = Orb.VSOP.Saturn();
  return Orb.VSOP.Exec(date,terms)
end

Orb.Planet.Uranus = function(date)
  local terms = Orb.VSOP.Uranus();
  return Orb.VSOP.Exec(date,terms)
end

Orb.Planet.Neptune = function(date)
  local terms = Orb.VSOP.Neptune();
  return Orb.VSOP.Exec(date,terms)
end


return Orb