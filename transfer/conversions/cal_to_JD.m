function JD = cal_to_JD(YY, MM, DD, hh, mm, ss)
    JD = datenum(YY, MM, DD, hh, mm, ss) + 1721058.5;
end