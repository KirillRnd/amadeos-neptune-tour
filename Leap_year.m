function out = Leap_year(YEAR)
%LEAP_YEAR  Определяет по году YEAR, является ли год високосным или нет

if (fix(YEAR/400) == YEAR/400)
	%Год високосный
	out = 1;
elseif (fix(YEAR/4) == YEAR/4) && (fix(YEAR/100) ~= YEAR/100)
	%Год високосный
    out = 1;
else
	%Год невисокосный
	out = 0;
end
