-- Orb.lua 
-- partial port of orb.js into lua lang
-- MIT License / Isana, Kashiwai 2021
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

function Orb.Explode(explodeCode, str)
  
  local result = {}
  --区切りがない場合は、新たに配列を作成して返す
  if(string.find(str,explodeCode) == nil) then
      result[1] = str
      return result
  end

  local maxIndex = #str
  local index=1
  local resultID = 1

  while (index<=maxIndex) do
      local findIndex = string.find(str,explodeCode,index)
      if(findIndex~=nil) then
          result[resultID] = string.sub(str,index,findIndex-1)
          resultID = resultID + 1
          index = findIndex + 1
      else
          result[resultID] = string.sub(str,index)
          break
      end
  end

  return result
end


Orb.Time = {}

-- Terrestrial Time
function Orb.Time.DeltaT()
  return (37 + 32)
end

-- Jurian day

function Orb.Time.JD(date)
  local year = tonumber(date.year);
  local month = tonumber(date.month);;
  local day = tonumber(date.day);
  local time_in_day = tonumber(date.hour) / 24 + tonumber(date.min) / 1440 + tonumber(date.sec) / 86400;
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
  local hour = tonumber(date.hour)
  local min = tonumber(date.min)
  local sec = tonumber(date.sec)

  local time_in_sec = hour * 3600 + min * 60 + sec
  local time_in_day = hour / 24 + min / 1440 + sec / 86400;
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
  ra = Orb.RoundAngle(ra)
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
  local eqx = gcx;
  local eqy = gcy * math.cos(ecl * rad) - gcz * math.sin(ecl * rad);
  local eqz = gcy * math.sin(ecl * rad) + gcz * math.cos(ecl * rad);
  local ra = math.atan2(eqy, eqx) / rad;
  ra = Orb.RoundAngle(ra)
  ra = ra / 15
  local dec = math.atan2(eqz, math.sqrt(eqx * eqx + eqy * eqy)) / rad;
  local distance = math.sqrt(eqx * eqx + eqy * eqy + eqz * eqz);
  return {
    x = eqx,
    y = eqy,
    z = eqz,
    ra = ra,
    dec = dec,
    distance = distance,
    date = date
  }
end

function Orb.Coord.NutationAndObliquity(date)
  local rad = math.pi/180;
  local dt = Orb.Time.DeltaT()
  local jd = Orb.Time.JD(date) + dt/86400;
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
    obliquity = obliquity,
    date = date
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
      zdot = dotvec.z,
      date = date
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
    atmospheric_refraction = atmospheric_refraction,
    date = date
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

Orb.Sun = function(date)
  local rad = math.pi/180;
  local dt = Orb.Time.DeltaT()
  local jd = Orb.Time.JD(date) + dt/86400;
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
  local distance = radius;
  local x = distance * math.cos(longitude * rad);
  local y = distance * (math.sin(longitude * rad) * math.cos(obliquity * rad));
  local z = distance * (math.sin(longitude * rad) * math.sin(obliquity * rad));
  local ra = math.atan2(math.cos(obliquity * rad) * math.sin(longitude * rad), math.cos(longitude * rad))
  ra = Orb.RoundNum(ra/rad,360)
  ra = ra / 15
  local dec = math.asin(math.sin(obliquity * rad) * math.sin(longitude * rad));
  dec = dec / rad;
  return {
    ra = ra,
    dec = dec,
    x = x,
    y = y,
    z = z,
    longitude = longitude,
    distance = distance,
    obliquity = obliquity
  }
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


--SGP4.satellite
Orb.SGP4 = {}

Orb.SGP4.Exec = function (time, name, tle)
  local rad =  math.pi / 180
  local omm = Orb.SGP4.TLE2OMM(time,name, tle.first_line, tle.second_line) 
  local sgp4 = Orb.SGP4.SetSGP4(omm)
  local epoch_array =  Orb.Explode("T", omm.EPOCH)
  local epoch_date = Orb.Explode("-", epoch_array[1])
  local epoch_time = Orb.Explode(":", epoch_array[2])
  local now_sec = {year = time.year, month = time.month , day = time.day, hour = time.hour, min = time.min, sec = time.sec}
  local epoch_sec = {year = tonumber(epoch_date[1]), month = tonumber(epoch_date[2]) , day = tonumber(epoch_date[3]), hour = tonumber(epoch_time[1]), min = tonumber(epoch_time[2]), sec = tonumber(epoch_time[3])}
  local tsince = (os.time(now_sec) - os.time(epoch_sec)) / 60

  local xmo = sgp4.xmo
  local xmdot = sgp4.xmdot
  local omegao = sgp4.omegao
  local omgdot = sgp4.omgdot
  local xnodeo = sgp4.xnodeo
  local xnodot = sgp4.xnodot
  local xnodcf = sgp4.xnodcf
  local bstar = sgp4.bstar
  local t2cof = sgp4.t2cof
  local omgcof = sgp4.omgcof
  local isimp = sgp4.isimp
  local xmcof = sgp4.xmcof
  local eta = sgp4.eta
  local delmo = sgp4.delmo
  local c1 = sgp4.c1
  local c4 = sgp4.c4
  local c5 = sgp4.c5
  local d2 = sgp4.d2
  local d3 = sgp4.d3
  local d4 = sgp4.d4
  local sinmo = sgp4.sinmo
  local t3cof = sgp4.t3cof
  local t4cof = sgp4.t4cof
  local t5cof = sgp4.t5cof
  local aodp = sgp4.aodp
  local eo = sgp4.eo
  local xnodp = sgp4.xnodp
  local xke = sgp4.xke
  local xlcof = sgp4.xlcof
  local aycof = sgp4.aycof
  local x3thm1 = sgp4.x3thm1
  local x1mth2 = sgp4.x1mth2
  local xincl = sgp4.xincl
  local cosio = sgp4.cosio
  local sinio = sgp4.sinio
  local e6a = sgp4.e6a
  local ck2 = sgp4.ck2
  local x7thm1 = sgp4.x7thm1
  local xkmper = sgp4.xkmper
  local epoch_year = sgp4.epoch_year
  local epoch = sgp4.epoch
  local xmdf = xmo + xmdot * tsince
  local omgadf = omegao + omgdot * tsince
  local xnoddf = xnodeo + xnodot * tsince
  local omega = omgadf
  local xmp = xmdf
  local tsq = tsince * tsince
  local xnode = xnoddf + xnodcf * tsq
  local tempa = 1.0 - c1 * tsince
  local tempe = bstar * c4 * tsince
  local templ = t2cof * tsq

  if not(isimp == 1) then
    local delomg = omgcof * tsince
    local delm = xmcof * (math.pow((1.0 + eta * math.cos(xmdf)), 3) - delmo)
    local temp = delomg + delm
    local xmp = xmdf + temp
    omega = omgadf - temp
    local tcube = tsq * tsince
    local tfour = tsince * tcube
    tempa = tempa - d2 * tsq - d3 * tcube - d4 * tfour
    tempe = tempe + bstar * c5 * (math.sin(xmp) - sinmo)
    templ = templ + t3cof * tcube + tfour * (t4cof + tsince * t5cof)
  end

  local a = aodp * tempa * tempa
  local e = eo - tempe
  local xl = xmp + omega + xnode + xnodp * templ
  local beta = math.sqrt(1.0 - e * e)
  local xn = xke / math.pow(a, 1.5)

  -- long period periodics
  local axn = e * math.cos(omega)
  local temp = 1.0 / (a * beta * beta)
  local xll = temp * xlcof * axn
  local aynl = temp * aycof
  local xlt = xl + xll
  local ayn = e * math.sin(omega) + aynl

  -- solve keplers equation
  local capu = (xlt - xnode) % (2.0 * math.pi)
  local temp2 = capu;
  local sinepw
  local cosepw
  local temp3
  local temp4
  local temp5
  local temp6
  local epw

  for i = 1, 10 do
    sinepw = math.sin(temp2)
    cosepw = math.cos(temp2)
    temp3 = axn * sinepw
    temp4 = ayn * cosepw
    temp5 = axn * cosepw
    temp6 = ayn * sinepw
    epw = (capu - temp4 + temp3 - temp2) / (1.0 - temp5 - temp6) + temp2

    if math.abs(epw - temp2) <= e6a then
      break
    end

    temp2 = epw
    
  end

  -- short period preliminary quantities
  local ecose = temp5 + temp6
  local esine = temp3 - temp4
  local elsq = axn * axn + ayn * ayn
  local temp = 1.0 - elsq
  local pl = a * temp
  local r = a * (1.0 - ecose)
  local temp1 = 1.0 / r
  local rdot = xke * math.sqrt(a) * esine * temp1
  local rfdot = xke * math.sqrt(pl) * temp1
  local temp2 = a * temp1
  local betal = math.sqrt(temp)
  local temp3 = 1.0 / (1.0 + betal)
  local cosu = temp2 * (cosepw - axn + ayn * esine * temp3)
  local sinu = temp2 * (sinepw - ayn - axn * esine * temp3)
  local u = math.atan2(sinu, cosu)

  if u < 0 then
    u = u + 2 * math.pi
  end

  local sin2u = 2.0 * sinu * cosu
  local cos2u = 2.0 * cosu * cosu - 1.
  local temp = 1.0 / pl
  local temp1 = ck2 * temp
  local temp2 = temp1 * temp
  -- update for short periodics
  local rk = r * (1.0 - 1.5 * temp2 * betal * x3thm1) + 0.5 * temp1 * x1mth2 * cos2u
  local uk = u - 0.25 * temp2 * x7thm1 * sin2u
  local xnodek = xnode + 1.5 * temp2 * cosio * sin2u
  local xinck = xincl + 1.5 * temp2 * cosio * sinio * cos2u
  local rdotk = rdot - xn * temp1 * x1mth2 * sin2u
  local rfdotk = rfdot + xn * temp1 * (x1mth2 * cos2u + 1.5 * x3thm1)
  -- orientation vectors
  local sinuk = math.sin(uk)
  local cosuk = math.cos(uk)
  local sinik = math.sin(xinck)
  local cosik = math.cos(xinck)
  local sinnok = math.sin(xnodek)
  local cosnok = math.cos(xnodek)
  local xmx = -sinnok * cosik
  local xmy = cosnok * cosik
  local ux = xmx * sinuk + cosnok * cosuk
  local uy = xmy * sinuk + sinnok * cosuk
  local uz = sinik * sinuk
  local vx = xmx * cosuk - cosnok * sinuk
  local vy = xmy * cosuk - sinnok * sinuk
  local vz = sinik * cosuk
  local x = rk * ux
  local y = rk * uy
  local z = rk * uz
  local xdot = rdotk * ux + rfdotk * vx
  local ydot = rdotk * uy + rfdotk * vy
  local zdot = rdotk * uz + rfdotk * vz
  local xkm = (x * xkmper)
  local ykm = (y * xkmper)
  local zkm = (z * xkmper)
  local xdotkmps = (xdot * xkmper / 60)
  local ydotkmps = (ydot * xkmper / 60)
  local zdotkmps = (zdot * xkmper / 60)
  return {
    x = xkm,
    y = ykm,
    z = zkm,
    xdot = xdotkmps,
    ydot = ydotkmps,
    zdot = zdotkmps,
  }
end


Orb.SGP4.TLE2OMM = function (date,name, first_line, second_line) 

  local line1 = first_line;
  local line2 = second_line;
  local creation_date = date.year .. "-" .. Orb.ZeroFill(date.month , 2) .. "-" .. Orb.ZeroFill(date.day, 2) .. " " .. Orb.ZeroFill(date.hour, 2) .. ":" .. Orb.ZeroFill(date.min, 2) .. ":" .. Orb.ZeroFill(date.sec, 2)
  local id = line1:sub(10, 18)
  local epystr = ""

  if tonumber(id:sub(1, 2)) < 58 then
    epystr = "20"
  else 
    epystr = "19"
  end

  local international_designator = epystr .. id:sub(1, 2) .. "-" .. id:sub(3, 7)
  local epy = tonumber(line1:sub(19, 20))
  local epoch_year

  if epy < 57 then
    epoch_year = epy + 2000
  else 
    epoch_year = epy + 1900
  end

  local doy = tonumber(line1:sub(21, 32))
  local year2 = epoch_year
  local epoch_data = {year = year2, month = 1, day = 1, hour = 0, min = 0, sec = 0}

  local epoch_jd = Orb.Time.JD(epoch_data) + (doy) -1
  local epoch = Orb.Time.JDToUTC(epoch_jd)
  --print(epoch.year)
  --print(epoch.month)
  --print(epoch.day)
  --print(epoch.hour)
  --print(epoch.min)
  --print(epoch.sec)

  local epoch_str = epoch.year .. "-" .. Orb.ZeroFill(epoch.month, 2) .. "-" .. Orb.ZeroFill(epoch.day, 2) .. 
  "T" .. Orb.ZeroFill(epoch.hour, 2) .. ":" .. Orb.ZeroFill(epoch.min, 2) .. ":" .. Orb.ZeroFill(epoch.sec, 2)

  local bstar_mantissa = tonumber(line1:sub(54, 59)) * 1e-5;
  local bstar_exponent = tonumber("1e" .. tonumber(line1:sub(60, 61)));
  local bstar = bstar_mantissa * bstar_exponent
  local mm_ddot = Orb.Explode("-", line1:sub(46, 52))
 
  if mm_ddot[2] == nil then
    mm_ddot = Orb.Explode("+", line1:sub(46, 52))
  end
  
  local mean_motion_ddot = tonumber(mm_ddot[1]) * 10 ^ (0 - tonumber(mm_ddot[2]))

  local omm = {
    CCSDS_OMM_VERS = "2.0",
    COMMENT = "GENERATED VIA ORB.JS",
    CREATION_DATE = creation_date,
    ORIGINATOR = "",
    OBJECT_NAME = name,
    OBJECT_ID = international_designator,
    CENTER_NAME = "EARTH",
    REF_FRAME = "TEME",
    TIME_SYSTEM = "UTC",
    MEAN_ELEMENT_THEORY = "SGP4",
    EPOCH = epoch_str,
    MEAN_MOTION = tonumber(line2:sub(53, 63)),
    ECCENTRICITY = tonumber(line2:sub(27, 33)) * 1e-7,
    INCLINATION = tonumber(line2:sub(9, 16)),
    RA_OF_ASC_NODE = tonumber(line2:sub(18, 25)),
    ARG_OF_PERICENTER = tonumber(line2:sub(35, 42)),
    MEAN_ANOMALY = tonumber(line2:sub(44, 51)),
    EPHEMERIS_TYPE = tonumber(line1:sub(63, 63)),
    CLASSIFICATION_TYPE = line1:sub(8, 8),
    NORAD_CAT_ID = tonumber(line1:sub(3, 7)),
    ELEMENT_SET_NO = tonumber(line1:sub(64, 68)),
    REV_AT_EPOCH = tonumber(line2:sub(64, 68)),
    BSTAR = bstar,
    MEAN_MOTION_DOT = tonumber(line1:sub(35, 43)),
    MEAN_MOTION_DDOT = mean_motion_ddot,
    USER_DEFINED_TLE_LINE0 = "0 " .. name,
    USER_DEFINED_TLE_LINE1 = line1,
    USER_DEFINED_TLE_LINE2 = line2
  }

  return omm

end


Orb.SGP4.SetSGP4 = function (omm)

  local torad = math.pi / 180
  local ck2 = 5.413080e-4
  local ck4 = 0.62098875e-6
  local e6a = 1.0e-6
  local qoms2t = 1.88027916e-9
  local s = 1.01222928 -- 1.0+78.0/xkmper
  local tothrd = 0.66666667
  local xj3 = -0.253881e-5
  local xke = 0.743669161e-1
  local xkmper = 6378.135
  local xmnpda = 1440.0 -- min_par_day
  local ae = 1.0
  local pi = 3.14159265
  local pio2 = 1.57079633
  local twopi = 6.2831853
  local x3pio2 = 4.71238898
  local bstar = omm.BSTAR
  local xincl = omm.INCLINATION * torad
  local xnodeo = omm.RA_OF_ASC_NODE * torad
  local eo = omm.ECCENTRICITY
  local omegao = omm.ARG_OF_PERICENTER * torad
  local xmo = omm.MEAN_ANOMALY * torad
  local xno = omm.MEAN_MOTION * 2.0 * math.pi / 1440.0
  local a1 = math.pow(xke / xno, tothrd)
  local cosio = math.cos(xincl)
  local theta2 = cosio * cosio
  local x3thm1 = 3 * theta2 - 1.0
  local eosq = eo * eo
  local betao2 = 1 - eosq
  local betao = math.sqrt(betao2)
  local del1 = 1.5 * ck2 * x3thm1 / (a1 * a1 * betao * betao2)
  local ao = a1 * (1 - del1 * ((1.0 / 3.0) + del1 * (1.0 + (134.0 / 81.0) * del1)))
  local delo = 1.5 * ck2 * x3thm1 / (ao * ao * betao * betao2)
  local xnodp = xno / (1.0 + delo) --original_mean_motion
  local aodp = ao / (1.0 - delo) --semi_major_axis
  local orbital_period = 1440.0 / omm.MEAN_MOTION
  local isimp = 0

  if ((aodp * (1.0 - eo) / ae) < (220.0 / xkmper + ae)) then
    isimp = 1;
  end
  
  local s4 = s
  local qoms24 = qoms2t
  local perigee = (aodp * (1.0 - eo) - ae) * xkmper
  local apogee = (aodp * (1.0 + eo) - ae) * xkmper

  if perigee < 156.0 then
    s4 = perigee - 78.0
    if perigee <= 98.0 then
      s4 = 20.0
    else
      qoms24 = math.pow(((120.0 - s4) * ae / xkmper), 4)
      s4 = s4 / xkmper + ae;
    end
  end
  
  local pinvsq = 1.0 / (aodp * aodp * betao2 * betao2)
  local tsi = 1.0 / (aodp - s4)
  local eta = aodp * eo * tsi
  local etasq = eta * eta
  local eeta = eo * eta
  local psisq = math.abs(1.0 - etasq)
  local coef = qoms24 * math.pow(tsi, 4)
  local coef1 = coef / math.pow(psisq, 3.5)
  local c2 = coef1 * xnodp * (aodp * (1.0 + 1.5 * etasq + eeta * (4.0 + etasq)) + 0.75 * ck2 * tsi / psisq * x3thm1 * (8.0 + 3.0 * etasq * (8.0 + etasq)))
  local c1 = bstar * c2
  local sinio = math.sin(xincl)
  local a3ovk2 = -xj3 / ck2 * math.pow(ae, 3)
  local c3 = coef * tsi * a3ovk2 * xnodp * ae * sinio / eo
  local x1mth2 = 1.0 - theta2
  local c4 = 2.0 * xnodp * coef1 * aodp * betao2 * (eta * (2.0 + 0.5 * etasq) + eo * (0.5 + 2.0 * etasq) - 2.0 * ck2 * tsi / (aodp * psisq) * 
  (-3.0 * x3thm1 * (1.0 - 2.0 * eeta + etasq * (1.5 - 0.5 * eeta)) + 0.75 * x1mth2 * (2.0 * etasq - eeta * (1.0 + etasq)) * math.cos((2.0 * omegao))))
  local c5 = 2.0 * coef1 * aodp * betao2 * (1.0 + 2.75 * (etasq + eeta) + eeta * etasq)
  local theta4 = theta2 * theta2
  local temp1 = 3.0 * ck2 * pinvsq * xnodp
  local temp2 = temp1 * ck2 * pinvsq
  local temp3 = 1.25 * ck4 * pinvsq * pinvsq * xnodp
  local xmdot = xnodp + 0.5 * temp1 * betao * x3thm1 + 0.0625 * temp2 * betao * (13.0 - 78.0 * theta2 + 137.0 * theta4)
  local x1m5th = 1.0 - 5.0 * theta2
  local omgdot = -0.5 * temp1 * x1m5th + 0.0625 * temp2 * (7.0 - 114.0 * theta2 + 395.0 * theta4) + temp3 * (3.0 - 36.0 * theta2 + 49.0 * theta4)
  local xhdot1 = -temp1 * cosio
  local xnodot = xhdot1 + (0.5 * temp2 * (4.0 - 19.0 * theta2) + 2.0 * temp3 * (3.0 - 7.0 * theta2)) * cosio
  local omgcof = bstar * c3 * math.cos(omegao)
  local xmcof = -tothrd * coef * bstar * ae / eeta
  local xnodcf = 3.5 * betao2 * xhdot1 * c1
  local t2cof = 1.5 * c1
  local xlcof = 0.125 * a3ovk2 * sinio * (3.0 + 5.0 * cosio) / (1.0 + cosio)
  local aycof = 0.25 * a3ovk2 * sinio;
  local delmo = math.pow((1.0 + eta * math.cos(xmo)), 3)
  local sinmo = math.sin(xmo);
  local x7thm1 = 7.0 * theta2 - 1.0
  local c1sq
  local d2
  local temp
  local d3
  local d4
  local t3cof
  local t4cof
  local t5cof
  if not(isimp == 1) then
    c1sq = c1 * c1
    d2 = 4.0 * aodp * tsi * c1sq
    temp = d2 * tsi * c1 / 3.0
    d3 = (17.0 * aodp + s4) * temp
    d4 = 0.5 * temp * aodp * tsi * (221.0 * aodp + 31.0 * s4) * c1
    t3cof = d2 + 2.0 * c1sq
    t4cof = 0.25 * (3.0 * d3 + c1 * (12.0 * d2 + 10.0 * c1sq))
    t5cof = 0.2 * (3.0 * d4 + 12.0 * c1 * d3 + 6.0 * d2 * d2 + 15.0 * c1sq * (2.0 * d2 + c1sq))
  end
  --set accesser
  return {
    omm = omm,
    apogee = apogee,
    perigee = perigee,
    orbital_period = orbital_period,
    xmo = xmo,
    xmdot = xmdot,
    omegao = omegao,
    omgdot = omgdot,
    xnodeo = xnodeo,
    xnodot = xnodot,
    xnodcf = xnodcf,
    bstar = bstar,
    t2cof = t2cof,
    omgcof = omgcof,
    isimp = isimp,
    xmcof = xmcof,
    eta = eta,
    delmo = delmo,
    c1 = c1,
    c4 = c4,
    c5 = c5,
    d2 = d2,
    d3 = d3,
    d4 = d4,
    sinmo = sinmo,
    t3cof = t3cof,
    t4cof = t4cof,
    t5cof = t5cof,
    aodp = aodp,
    eo = eo,
    xnodp = xnodp,
    xke = xke,
    xlcof = xlcof,
    aycof = aycof,
    x3thm1 = x3thm1,
    x1mth2 = x1mth2,
    xincl = xincl,
    cosio = cosio,
    sinio = sinio,
    e6a = e6a,
    ck2 = ck2,
    x7thm1 = x7thm1,
    xkmper = xkmper
  }

end


Orb.SGP4.RectangularToGeographic = function (time, rect)

  local time = time
  local xkm = rect.x
  local ykm = rect.y
  local zkm = rect.z
  local xdotkmps = rect.xdot
  local ydotkmps = rect.ydot
  local zdotkmps = rect.zdot
  local rad = math.pi / 180
  local gmst = Orb.Time.gst(time)
  local lst = gmst * 15;
  local f = 0.00335277945 --Earth's flattening term in WGS-72 (= 1/298.26)
  local a = 6378.135  --Earth's equational radius in WGS-72 (km)
  local r = math.sqrt(xkm * xkm + ykm * ykm)
  local lng = math.atan2(ykm, xkm) / rad - lst

  if lng > 360 then 
    lng = lng % 360
  end

  if lng < 0 then
    lng = lng % 360 + 360
  end

  if lng > 180 then
    lng = lng - 360
  end

  local lat = math.atan2(zkm, r)
  local e2 = f * (2 - f)
  local tmp_lat = 0

  tmp_lat = lat;
  local sin_lat = math.sin(tmp_lat)
  local c = 1 / math.sqrt(1 - e2 * sin_lat * sin_lat);
  lat = math.atan2(zkm + a * c * e2 * (math.sin(tmp_lat)), r)
  
  while math.abs(lat - tmp_lat) > 0.0001 do

    tmp_lat = lat;
    sin_lat = math.sin(tmp_lat)
    c = 1 / math.sqrt(1 - e2 * sin_lat * sin_lat);
    lat = math.atan2(zkm + a * c * e2 * (math.sin(tmp_lat)), r)

  end

  local alt = r / math.cos(lat) - a * c;
  local v = math.sqrt(xdotkmps * xdotkmps + ydotkmps * ydotkmps + zdotkmps * zdotkmps)
  return {
    longitude = lng,
    latitude = lat / rad,
    altitude = alt,
    velocity = v
  }
end


Orb.SGP4.LongitudeLatitude = function(LongitudeLatitude)

  local LongitudeLatitude = math.abs(LongitudeLatitude)
  local h_LongitudeLatitude = math.floor(LongitudeLatitude)
  local m_LongitudeLatitude = math.floor((LongitudeLatitude  - h_LongitudeLatitude) * 60)
  local s_LongitudeLatitude = ((LongitudeLatitude  - h_LongitudeLatitude) * 60 * 60) - (m_LongitudeLatitude * 60)

  if h_LongitudeLatitude > 360 then 
    h_LongitudeLatitude = h_LongitudeLatitude % 360
  end

  if h_LongitudeLatitude < 0 then
    h_LongitudeLatitude = h_LongitudeLatitude % 360 + 360
  end

  if h_LongitudeLatitude > 180 then
    h_LongitudeLatitude = h_LongitudeLatitude - 360
  end
  h_LongitudeLatitude = math.abs(h_LongitudeLatitude)

  return{
    h = h_LongitudeLatitude,
    m = m_LongitudeLatitude,
    s = s_LongitudeLatitude
  }

end


Orb.SGP4.satellite = function(date, name, tle)

  local rect = Orb.SGP4.Exec(date, name, tle)
  local geo = Orb.SGP4.RectangularToGeographic(date, rect)
  local long = Orb.SGP4.LongitudeLatitude(geo.longitude)
  local lat = Orb.SGP4.LongitudeLatitude(geo.latitude)

  local longitude_EW = "W"
  local latitude_NS = "N"
  
  if geo.longitude <= 0 then
    latitude_NS = "S"
  end
  
  if geo.latitude <= 0 then
    longitude_EW = "E"
  end

  local NS = Orb.ZeroFill(lat.h, 3) .. "˚" .. Orb.ZeroFill(lat.m, 2) .. "'" .. Orb.ZeroFill(lat.s, 2)
  local EW = Orb.ZeroFill(long.h, 3) .. "˚" .. Orb.ZeroFill(long.m, 2) .. "'" .. Orb.ZeroFill(long.s, 2)

  return {
    x = rect.x,
    y = rect.y,
    z = rect.z,
    xdot = rect.xdot,
    ydot = rect.ydot,
    zdot = rect.zdot,
    longitude = geo.longitude,
    latitude = geo.latitude,
    altitude = geo.altitude,
    velocity = geo.velocity,
    NS = latitude_NS,
    EW = longitude_EW,
    lat = NS,
    lng = EW
  }

end

function Orb.info()

  local Luna = Orb.Luna.info()
  local VSOP = Orb.VSOP.info()
  local Orb = {name = "AkinoMizuki", update = "2023-04-02 17:00 JST"}

  return{
    Lunainfo = Luna,
    VSOPinfo = VSOP,
    Orbinfo = Orb
  }

end

return Orb