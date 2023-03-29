
local Luna = {}

function Luna.info()
--更新日

  return{
    name = "Isana, Kashiwai",
    update = "2021-07-11 23:00 JST",
  }

end--END_更新日

Luna.latlng =  function (date)
  local rad = math.pi/180;
  local deg = 180 / math.pi;
  local dt = Orb.Time.DeltaT()
  local jd = Orb.Time.JD(date) + dt/86400;

  --ephemeris days from the epch J2000.0
  local t = (jd - 2451545.0) / 36525;
  local t2 = t * t;
  local t3 = t * t * t;
  local t4 = t * t * t * t;
  local e = 1 - 0.002516 * t - 0.0000074 * t2;
  local L1 = (218.3164477 + 481267.88123421 * t - 0.0015786 * t2 + t3 / 538841 - t4 / 65194000);
  L1 = Orb.RoundAngle(L1) * rad;
  local D0 = (297.8501921 + 445267.1114034 * t - 0.0018819 * t2 + t3 / 545868 - t4 / 113065000);
  D0 = Orb.RoundAngle(D0) * rad;
  local M0 = (357.5291092 + 35999.0502909 * t - 0.0001536 * t2 + t3 / 24490000);
  M0 = Orb.RoundAngle(M0) * rad;
  local M1 = (134.9633964 + 477198.8675055 * t + 0.0087414 * t2 + t3 / 69699 - t4 / 14712000);
  M1 = Orb.RoundAngle(M1) * rad;
  local F0 = (93.2720950 + 483202.0175233 * t - 0.0036539 * t2 - t3 / 3526000 + t4 / 863310000);
  F0 = Orb.RoundAngle(F0) * rad;
  local A1 = (119.75 + 131.849 * t);
  A1 = Orb.RoundAngle(A1) * rad;
  local A2 = (53.09 + 479264.290 * t);
  A2 = Orb.RoundAngle(A2) * rad;
  local A3 = (313.45 + 481266.484 * t);
  A3 = Orb.RoundAngle(A3) * rad;

  local SigmaL = function ()
    local result = 0;
    local terms = Luna.Terms.LR();
    local terms_length = #terms;
    for i = 1, terms_length do
      local coef = terms[i][5];
      local multi = {terms[i][1], terms[i][2], terms[i][3], terms[i][4]}
      local e_coef
      if math.abs(multi[2]) == 1 then
        e_coef = e;
      elseif math.abs(multi[2]) == 2 then
        e_coef = e * e;
      else
        e_coef = 1;
      end
      local asin = multi[1] * D0 + multi[2] * M0 + multi[3] * M1 + multi[4] * F0;
      result = result + coef * math.sin(asin) * e_coef;
    end

    result = result + 3958 * math.sin(A1)
    result = result +  1962 * math.sin(L1 - F0)
    result = result +  318 * math.sin(A2)
    return result;
  end

  local SigmaR = function ()
    local result = 0;
    local terms = Luna.Terms.LR();
    local terms_length = #terms;
    for i = 1, terms_length do
      local coef = terms[i][6];
      local e_coef
      local multi = {terms[i][1], terms[i][2], terms[i][3], terms[i][4]}
      if math.abs(multi[2]) == 1 then
        e_coef = e;
     elseif math.abs(multi[2]) == 2 then
        e_coef = e * e;
     else
        e_coef = 1;
     end
      local acos = multi[1] * D0 + multi[2] * M0 + multi[3] * M1 + multi[4] * F0
      result = result + coef * math.cos(acos) * e_coef;
    end
    return result;
  end

  local SigmaB = function ()
    local result = 0;
    local terms = Luna.Terms.B();
    local terms_length = #terms;
    for i = 1, terms_length do
      local coef = terms[i][5];
      local e_coef
      local multi = {terms[i][1], terms[i][2], terms[i][3], terms[i][4]}
      if math.abs(multi[2]) == 1 then
        e_coef = e;
      elseif math.abs(multi[2]) == 2 then
        e_coef = e * e;
      else
        e_coef = 1;
      end
      local asin = multi[1] * D0 + multi[2] * M0 + multi[3] * M1 + multi[4] * F0
      result = result + coef * math.sin(asin) * e_coef;
    end
    result = result +  -2235 * math.sin(L1)
    result = result +  382 * math.sin(A3)
    result = result +  175 * math.sin(A1 - F0)
    result = result +  175 * math.sin(A1 + F0)
    result = result +  127 * math.sin(L1 - M1)
    result = result +  -115 * math.sin(L1 + M1)
    return result;
  end

  local sigma_l = SigmaL();
  local sigma_r = SigmaR();
  local sigma_b = SigmaB();
  local true_longitude = (L1 / rad) % 360 + (sigma_l) / 1000000
  local latitude = (sigma_b) / 1000000
  local distance = 385000.56 + sigma_r / 1000
  local nao = Orb.Coord.NutationAndObliquity(date)
  local nutation = nao.nutation;
  local obliquity = nao.obliquity;
  local apparent_longitude = true_longitude + nutation;
  local longitude = apparent_longitude;
  return {
    latitude = latitude,
    longitude = longitude,
    distance = distance,
    obliquity = obliquity,
    date = date
  }
end

Luna.equatorial = function (date)
    local latlng = Luna.latlng(date);
    local rad = math.pi/180;
    local latitude = latlng.latitude
    local longitude = latlng.longitude
    local distance = latlng.distance
    local obliquity = latlng.obliquity
    local ra_deg = math.atan2(math.sin(longitude * rad) * math.cos(obliquity * rad) - math.tan(latitude * rad) * math.sin(obliquity * rad), math.cos(longitude * rad)) / rad;
    local dec = math.asin(math.sin(latitude * rad) * math.cos(obliquity * rad) + math.cos(latitude * rad) * math.sin(obliquity * rad) * math.sin(longitude * rad)) / rad;
    local vec = {
      x = distance * math.cos(dec * rad) * math.cos(ra_deg * rad),
      y = distance * math.cos(dec * rad) * math.sin(ra_deg * rad),
      z = distance * math.sin(dec * rad)
    }
    ra = Orb.RoundAngle(ra_deg) / 15;  
    return {
      ra = ra,
      dec = dec,
      x = vec.x,
      y = vec.y,
      z = vec.z,
      distance = distance,
      obliquity = obliquity,
      date = date
    }
  end

  Luna.parallax = function (date)
    local latlng = Luna.latlng(date);
    local rad = math.pi/180;
    return math.asin(6378.14 / latlng.distance) / rad
  end

  Luna.phase = function (date)
    local rad = math.pi/180;
    local jd = Orb.Time.JD(date)
    local date_first = new Date(date.year, 0, 1, 0, 0, 0);
    local date_last = new Date(date.year, 11, 31, 11, 59, 59, 999);
    local since_new_year = (now - date_first) / (date_last - date_first);
    local y = date.year + since_new_year;

    local k = math.floor((y - 2000) * 12.3685);
    local t = k / 1236.85;
    local t2 = t * t;
    local t3 = t * t * t;
    local t4 = t * t * t * t;
    local jde0 = 2451550.09766 + 29.530588861 * k + 0.00015437 * t2 - 0.000000150 * t3 + 0.00000000073 * t4;
    local e = 1 - 0.002516 * t - 0.0000074 * t2;
    e = Orb.RoundAngle(e);
    --Sun's mean anomary at the time;
    local m0 = 2.5534 + 29.10535670 * k - 0.0000014 * t2 - 0.00000011 * t3;
    m0 = Orb.RoundAngle(m0);
    --Moon's mean anomary at the time;
    local m1 = 201.5643 + 385.81693528 * k + 0.0107582 * t2 + 0.00001238 * t3 - 0.000000011 * t4;
    m1 = Orb.RoundAngle(m1);
    --Moon's argument of latitude
    local f = 160.7108 + 390.67050284 * k - 0.0016118 * t2 - 0.00000227 * t3 + 0.000000011 * t4;
    f = Orb.RoundAngle(f);
    --Longitude of the ascending node of lunar orbit
    local omega = 124.7746 - 1.56375588 * k + 0.0020672 * t2 + 0.00000215 * t3;
    omega = Orb.RoundAngle(omega);
    local c1 = 0;
    c1 = c1 - 0.40720 * math.sin(m1 * rad);
    c1 = c1 + 0.17241 * e * math.sin(m0 * rad);
    c1 = c1 + 0.01608 * math.sin(2 * m1 * rad);
    c1 = c1 + 0.01039 * math.sin(2 * f * rad);
    c1 = c1 + 0.00739 * e * math.sin((m1 - m0) * rad);
    c1 = c1 - 0.00514 * e * math.sin((m1 + m0) * rad);
    c1 = c1 + 0.00208 * e * e * math.sin(2 * m0 * rad);
    c1 = c1 - 0.00111 * math.sin((m1 - 2 * f) * rad)
    c1 = c1 - 0.00057 * math.sin((m1 + 2 * f) * rad)
    c1 = c1 + 0.00056 * e * math.sin((2 * m1 + m0) * rad);
    c1 = c1 - 0.00042 * math.sin(3 * m1 * rad);
    c1 = c1 + 0.00042 * e * math.sin((m0 + 2 * f) * rad)
    c1 = c1 + 0.00038 * e * math.sin((m0 - 2 * f) * rad)
    c1 = c1 - 0.00024 * e * math.sin((2 * m1 - m0) * rad);
    c1 = c1 - 0.00017 * math.sin(omega * rad);
    c1 = c1 - 0.00007 * math.sin((m1 + 2 * m0) * rad);
    c1 = c1 + 0.00004 * math.sin((2 * m1 - 2 * f) * rad);
    c1 = c1 + 0.00004 * math.sin(3 * m0 * rad);
    c1 = c1 + 0.00003 * math.sin((m1 + m0 - 2 * f) * rad);
    c1 = c1 + 0.00003 * math.sin((2 * m1 + 2 * f) * rad);
    c1 = c1 - 0.00003 * math.sin((m1 + m0 + 2 * f) * rad);
    c1 = c1 + 0.00003 * math.sin((m1 - m0 + 2 * f) * rad);
    c1 = c1 - 0.00002 * math.sin((m1 - m0 - 2 * f) * rad);
    c1 = c1 - 0.00002 * math.sin((3 * m1 + m0) * rad);
    c1 = c1 + 0.00002 * math.sin(4 * m1 * rad);
    local a1 = 299.77 + 0.107408 * k - 0.009173 * t2;
    local a2 = 251.88 + 0.016321 * k;
    local a3 = 251.83 + 26.651886 * k;
    local a4 = 349.42 + 36.412478 * k;
    local a5 = 84.66 + 18.206239 * k;
    local a6 = 141.74 + 53.303771 * k;
    local a7 = 207.14 + 2.453732 * k;
    local a8 = 154.84 + 7.306860 * k;
    local a9 = 34.52 + 27.261239 * k;
    local a10 = 207.19 + 0.121824 * k;
    local a11 = 291.34 + 1.844379 * k;
    local a12 = 161.72 + 24.198154 * k;
    local a13 = 239.56 + 25.513099 * k;
    local a14 = 331.55 + 3.592518 * k;
    local c2 = 0;
    c2 = c2 + 0.000325 * math.sin(a1 * rad);
    c2 = c2 + 0.000165 * math.sin(a2 * rad);
    c2 = c2 + 0.000164 * math.sin(a3 * rad);
    c2 = c2 + 0.000126 * math.sin(a4 * rad);
    c2 = c2 + 0.000110 * math.sin(a5 * rad);
    c2 = c2 + 0.000062 * math.sin(a6 * rad);
    c2 = c2 + 0.000060 * math.sin(a7 * rad);
    c2 = c2 + 0.000056 * math.sin(a8 * rad);
    c2 = c2 + 0.000047 * math.sin(a9 * rad);
    c2 = c2 + 0.000042 * math.sin(a10 * rad);
    c2 = c2 + 0.000040 * math.sin(a11 * rad);
    c2 = c2 + 0.000037 * math.sin(a12 * rad);
    c2 = c2 + 0.000035 * math.sin(a13 * rad);
    c2 = c2 + 0.000023 * math.sin(a14 * rad);
    local jde = jde0 + c1 + c2;
    local phase_of_the_moon = jd - jde;
    return phase_of_the_moon;
  end


Luna.Terms = {}
Luna.Terms.LR = function()
  return {
    {0, 0, 1, 0, 6288774, -20905335},
    {2, 0, -1, 0, 1274027, -3699111},
    {2, 0, 0, 0, 658314, -2955968},
    {0, 0, 2, 0, 213618, -569925},
    {0, 1, 0, 0, -185116, 48888},
    {0, 0, 0, 2, -114332, -3149},
    {2, 0, -2, 0, 58793, 246158},
    {2, -1, -1, 0, 57066, -152138},
    {2, 0, 1, 0, 53322, -170733},
    {2, -1, 0, 0, 45758, -204586},
    {0, 1, -1, 0, -40923, -129620},
    {1, 0, 0, 0, -34720, 108743},
    {0, 1, 1, 0, -30383, 104755},
    {2, 0, 0, -2, 15327, 10321},
    {0, 0, 1, 2, -12528, 0},
    {0, 0, 1, -2, 10980, 79661},
    {4, 0, -1, 0, 10675, -34782},
    {0, 0, 3, 0, 10034, -23210},
    {4, 0, -2, 0, 8548, -21636},
    {2, 1, -1, 0, -7888, 24208},
    {2, 1, 0, 0, -6766, 30824},
    {1, 0, -1, 0, -5163, -8379},
    {1, 1, 0, 0, 4987, -16675},
    {2, -1, 1, 0, 4036, -12831},
    {2, 0, 2, 0, 3994, -10445},
    {4, 0, 0, 0, 3861, -11650},
    {2, 0, -3, 0, 3665, 14403},
    {0, 1, -2, 0, -2689, -7003},
    {2, 0, -1, 2, -2602, 0},
    {2, -1, -2, 0, 2390, 10056},
    {1, 0, 1, 0, -2348, 6322},
    {2, -2, 0, 0, 2236, -9884},
    {0, 1, 2, 0, -2120, 5751},
    {0, 2, 0, 0, -2069, 0},
    {2, -2, -1, 0, 2048, -4950},
    {2, 0, 1, -2, -1773, 4130},
    {2, 0, 0, 2, -1595, 0},
    {4, -1, -1, 0, 1215, -3958},
    {0, 0, 2, 2, -1110, 0},
    {3, 0, -1, 0, -892, 3258},
    {2, 1, 1, 0, -810, 2616},
    {4, -1, -2, 0, 759, -1897},
    {0, 2, -1, 0, -713, -2117},
    {2, 2, -1, 0, -700, 2354},
    {2, 1, -2, 0, 691, 0},
    {2, -1, 0, -2, 596, 0},
    {4, 0, 1, 0, 549, -1423},
    {0, 0, 4, 0, 537, -1117},
    {4, -1, 0, 0, 520, -1571},
    {1, 0, -2, 0, -487, -1739},
    {2, 1, 0, -2, -399, 0},
    {0, 0, 2, -2, -381, -4421},
    {1, 1, 1, 0, 351, 0},
    {3, 0, -2, 0, -340, 0},
    {4, 0, -3, 0, 330, 0},
    {2, -1, 2, 0, 327, 0},
    {0, 2, 1, 0, -323, 1165},
    {1, 1, -1, 0, 299, 0},
    {2, 0, 3, 0, 294, 0},
    {2, 0, -1, -2, 0, 8752}
  }
end

  Luna.Terms.B = function()
    return{
    {0, 0, 0, 1, 5128122},
    {0, 0, 1, 1, 280602},
    {0, 0, 1, -1, 277693},
    {2, 0, 0, -1, 173237},
    {2, 0, -1, 1, 55413},
    {2, 0, -1, -1, 46271},
    {2, 0, 0, 1, 32573},
    {0, 0, 2, 1, 17198},
    {2, 0, 1, -1, 9266},
    {0, 0, 2, -1, 8822},
    {2, -1, 0, -1, 8216},
    {2, 0, -2, -1, 4324},
    {2, 0, 1, 1, 4200},
    {2, 1, 0, -1, -3359},
    {2, -1, -1, 1, 2463},
    {2, -1, 0, 1, 2211},
    {2, -1, -1, -1, 2065},
    {0, 1, -1, -1, -1870},
    {4, 0, -1, -1, 1828},
    {0, 1, 0, 1, -1794},
    {0, 0, 0, 3, -1749},
    {0, 1, -1, 1, -1565},
    {1, 0, 0, 1, -1491},
    {0, 1, 1, 1, -1475},
    {0, 1, 1, -1, -1410},
    {0, 1, 0, -1, -1344},
    {1, 0, 0, -1, -1335},
    {0, 0, 3, 1, 1107},
    {4, 0, 0, -1, 1021},
    {4, 0, -1, 1, 833},
    {0, 0, 1, -3, 777},
    {4, 0, -2, 1, 671},
    {2, 0, 0, -3, 607},
    {2, 0, 2, -1, 596},
    {2, -1, 1, -1, 491},
    {2, 0, -2, 1, -451},
    {0, 0, 3, -1, 439},
    {2, 0, 2, 1, 422},
    {2, 0, -3, -1, 421},
    {2, 1, -1, 1, -366},
    {2, 1, 0, 1, -351},
    {4, 0, 0, 1, 331},
    {2, -1, 1, 1, 315},
    {2, -2, 0, -1, 302},
    {0, 0, 1, 3, -283},
    {2, 1, 1, -1, -229},
    {1, 1, 0, -1, 223},
    {1, 1, 0, 1, 223},
    {0, 1, -2, -1, -220},
    {2, 1, -1, -1, -220},
    {1, 0, 1, 1, -185},
    {2, -1, -2, -1, 181},
    {0, 1, 2, 1, -177},
    {4, 0, -2, -1, 176},
    {4, -1, -1, -1, 166},
    {1, 0, 1, -1, -164},
    {4, 0, 1, -1, 132},
    {1, 0, -1, -1, -119},
    {4, -1, 0, -1, 115},
    {2, -2, 0, 1, 107}
    }
    end

return Luna