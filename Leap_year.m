function out = Leap_year(YEAR)
%LEAP_YEAR  ���������� �� ���� YEAR, �������� �� ��� ���������� ��� ���

if (fix(YEAR/400) == YEAR/400)
	%��� ����������
	out = 1;
elseif (fix(YEAR/4) == YEAR/4) && (fix(YEAR/100) ~= YEAR/100)
	%��� ����������
    out = 1;
else
	%��� ������������
	out = 0;
end
