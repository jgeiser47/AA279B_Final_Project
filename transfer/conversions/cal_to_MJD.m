function MJD = cal_to_MJD(YY, MM, DD, hh, mm, ss)
    JD = datenum(YY, MM, DD, hh, mm, ss) + 1721058.5;
    MJD = JD - 2400000.5;
end