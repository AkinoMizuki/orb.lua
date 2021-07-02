local Orb = require("orb")

-- time conversion
local t = {year=2021,month=1,day=1,hour=0,min=0,sec=0}

local utc_str = os.date('%Y-%m-%dT%H:%M:%S', os.time(t))
print("UTC: " .. utc_str)

local jd = Orb.Time.jd(t)
print("JD: " .. jd)

local gst = Orb.Time.gst(t)
print("GST: ".. gst)


-- Equatorial Spherical(ra,dec) to Equatorial Rectangular(x,y,z)
-- Note: Right Ascension(ra) must be hours(not degree)

local sirius = {
  ra=6.75257,
  dec = -16.7131,
  distance = 543300
}

local vec = Orb.RadecToXYZ(sirius.ra,sirius.dec,1)
print("Sirius(equatorial)\n x:" .. vec.x .. ", y: " .. vec.y .. ", z: " .. vec.z)

local observer = {
  latitude = 35.658,
  longitude = 139.741,
  altitude = 0
}

local obs = Orb.Observe.RadecToHorizontal(t,sirius,observer)
print("Sirius(horizontal)\n azimuth:" .. obs.azimuth .. ", elevation: " .. obs.elevation)

-- Planet Position in Ecliptic Rectangular Coordinate (au) via VSOP Theory
local earth = Orb.Planet.Earth(t)
print("Earth(ecliptic)\n x:" .. earth.x .. ", y: " .. earth.y .. ", z: " .. earth.z)

local mars = Orb.Planet.Mars(t)
print("Mars(ecliptic)\n x:" .. mars.x .. ", y: " .. mars.y .. ", z: " .. mars.z)

-- Ecliptic Coordinate(Earth Centered) Mars Position
local ecm = Orb.EclipticToEquatorial(t,mars)
print("Mars(equatorial)\n x:" .. ecm.x .. ", y: " .. ecm.y .. ", z: " .. ecm.z)

-- Ecliptic Coordinate(Normalized)
local normalized = Orb.Normalize(ecm)
print("Mars(ecliptic/normalized)\n x:" .. normalized.x .. ", y: " .. normalized.y .. ", z: " .. normalized.z)

-- Horizontal Coordinate
local obs = Orb.Observe.RectToHorizontal(t,ecm,observer)
print("Mars(horizontal)\n azimuth:" .. obs.azimuth .. ", elevation: " .. obs.elevation)

