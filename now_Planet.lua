-- Test code for orb.lua
--https://github.com/AkinoMizuki/MoonSharpInspectionでの動作確認用
-- update 2024-0-18JST

local Orb = require("orb")

local auToKm = 149597870.700

-- Pluto elements from JPL Horizons
-- Unit: au & au/d
-- Reference frame : Ecliptic/J2000.0
-- GM:2.9591220828411951E-04 (Sun)
-- 2459396.500000000 = A.D. 2021-Jul-01 00:00:00.0000 TDB 
--  EC= 2.515592767106864E-01 QR= 2.967588460857290E+01 IN= 1.729107375790253E+01
--  OM= 1.103546798589673E+02 W = 1.141843508986021E+02 Tp=  2447879.317561375909
--  N = 3.947613990481508E-03 MA= 4.546539052564052E+01 TA= 7.086411718013353E+01
--  A = 3.965028049001756E+01 AD= 4.962467637146222E+01 PR= 9.119432671685543E+04

local pluto_elements = {
  gm = 2.9591220828411951E-04,
  eccentricity = 0.2519446,
  inclination = 17.09860,
  longitude_of_ascending_node = 110.29702,
  argument_of_periapsis = 115.37952,
  time_of_periapsis = 2448031.24959,
  semi_major_axis = 39.8362800
}

--距離計算
function CalculateDistance(a, b)

  local dx = b.x - a.x;
  local dy = b.y - a.y;
  local dz = b.z - a.z; 

  return math.sqrt(dx * dx + dy * dy + dz * dz)

end

while true do

  -- now
  local t = os.date("!*t")
  local jd = Orb.Time.JD(t)

  --惑星の位置計算
  local Mercury_Pos = Orb.Planet.Mercury(t)
  local Venus_Pos = Orb.Planet.Venus(t)
  local Earth_Pos = Orb.Planet.Earth(t)
  local Mars_Pos = Orb.Planet.Mars(t)
  local Jupiter_Pos = Orb.Planet.Jupiter(t)
  local Saturn_Pos = Orb.Planet.Saturn(t)
  local Uranus_Pos = Orb.Planet.Uranus(t)
  local Neptune_Pos = Orb.Planet.Neptune(t)
  local pluto_ecliptic = Orb.Kepler(t,pluto_elements)

  --月
  Moon_Pos = Orb.Moon(t)
  local moon_ecliptic = Orb.Coord.EquatorialToEcliptic(t,Moon_Pos)

  local EarthKm = {
    x = 0,
    y = 0,
    z = 0
  }

  local moonToEarthKm = CalculateDistance(EarthKm, moon_ecliptic)

  local Result = "" ..
  "JD :".. jd .. "\r\n" ..
  "Mercury  x:" .. Mercury_Pos.x .. ", y: " .. Mercury_Pos.z .. ", z: " .. Mercury_Pos.y .. "\r\n" ..
  "Venus  x:" .. Venus_Pos.x .. ", y: " .. Venus_Pos.z .. ", z: " .. Venus_Pos.y .. "\r\n" ..
  "Earth  x:" .. Earth_Pos.x .. ", y: " .. Earth_Pos.z .. ", z: " .. Earth_Pos.y .. "\r\n" ..
  "Mars  x:" .. Mars_Pos.x .. ", y: " .. Mars_Pos.z .. ", z: " .. Mars_Pos.y .. "\r\n" ..
  "Jupiter  x:" .. Jupiter_Pos.x .. ", y: " .. Jupiter_Pos.z .. ", z: " .. Jupiter_Pos.y .. "\r\n" ..
  "Saturn  x:" .. Saturn_Pos.x .. ", y: " .. Saturn_Pos.z .. ", z: " .. Saturn_Pos.y .. "\r\n" ..
  "Uranus  x:" .. Uranus_Pos.x .. ", y: " .. Uranus_Pos.z .. ", z: " .. Uranus_Pos.y .. "\r\n" ..
  "Neptune  x:" .. Neptune_Pos.x .. ", y: " .. Neptune_Pos.z .. ", z: " .. Neptune_Pos.y .. "\r\n" ..
  "Pluto x:" .. pluto_ecliptic.x .. ", y:" .. pluto_ecliptic.y .. ", z:" .. pluto_ecliptic.z .. "\r\n" ..
  "Moon x:" .. string.format("%3f",moon_ecliptic.x) .. ", y:" .. string.format("%3f",moon_ecliptic.y) .. ", z:" .. string.format("%3f",moon_ecliptic.z) .. "\r\n" ..
  "Moon Distance:" .. moonToEarthKm .. "\r\n" ..
  "Moon Distance:" .. moonToEarthKm / auToKm .. "\r\n" ..
  "Moon x:" .. string.format("%3f",moon_ecliptic.x / auToKm) .. ", y:" .. string.format("%3f",moon_ecliptic.y / auToKm) .. ", z:" .. string.format("%3f",moon_ecliptic.z / auToKm) .. "\r\n" ..
  ""
  WinForm.LED(false)
  printf(Result)

end
