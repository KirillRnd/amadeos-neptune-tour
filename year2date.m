function [year,month,day] = year2date(year)
%CALC_TIME  По числу YEAR определяет день, месяц и год.

if year < 1000, year = year + 2000; end

isleap = Leap_year(fix(year));

partyear = year - fix(year);

for i = 0:(365+isleap)
    if (partyear >= i/(365+isleap))  && (partyear < (i+1)/(365+isleap))
        numday = i + 1;
        break
    end
end

if (numday >= 1) && (numday <= 31)
    day   = numday;
    month = 1;
elseif (numday > 31) && (numday <= 59+isleap)
    day   = numday - 31;
    month = 2;
elseif (numday > 59+isleap) && (numday <= 90+isleap)
    day   = numday - 59 - isleap;
    month = 3;
elseif (numday > 90+isleap) && (numday <= 120+isleap)
    day   = numday - 90 - isleap;
    month = 4;
elseif (numday > 120+isleap) && (numday <= 151+isleap)
    day   = numday - 120 - isleap;
    month = 5;
elseif (numday > 151+isleap) && (numday <= 181+isleap)
    day   = numday - 151 - isleap;
    month = 6;
elseif (numday > 181+isleap) && (numday <= 212+isleap)
    day   = numday - 181 - isleap;
    month = 7;
elseif (numday > 212+isleap) && (numday <= 243+isleap)
    day   = numday - 212 - isleap;
    month = 8;
elseif (numday > 243+isleap) && (numday <= 273+isleap)
    day   = numday - 243 - isleap;
    month = 9;
elseif (numday > 273+isleap) && (numday <= 304+isleap)
    day   = numday - 273 - isleap;
    month = 10;
elseif (numday > 304+isleap) && (numday <= 334+isleap)
    day   = numday - 304 - isleap;
    month = 11;
elseif (numday > 334+isleap) && (numday <= 365+isleap)
    day   = numday - 334 - isleap;
    month = 12;
end

year = fix(year);