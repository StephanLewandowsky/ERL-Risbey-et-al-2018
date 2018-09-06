wget ftp://ftp.ncdc.noaa.gov/pub/data/anomalies/usingGHCNMv2/monthly.land_ocean.90S.90N.df_1901-2000mean.dat

foreach name (*.txt.?? *.txt.??? *.dat)
set date = `awk '/^[12]/{n1 = $1; n2 = $2}END{printf "%04d%02d",n1,n2}' $name`
awk '/^[12]/{if ($3>-99) print $1+$2/12.-1/24.,$3}' $name > NOAA.$date.temp
end

foreach name (from_noaa/*.asc)
set date = `awk '/^[12]/{n1 = $1; n2 = $2}END{printf "%04d%02d",n1,n2}' $name`
awk '/^[12]/{if ($3>-99) print $1+$2/12.-1/24.,$3}' $name > NOAA.$date.temp
end

@ y = 2015
while ($y < 2017)
curl http://www.ncdc.noaa.gov/cag/time-series/global/globe/land_ocean/p12/12/1880-$y.csv > 1880-$y.csv
@ y = $y + 1
end

foreach name (*.csv)
set date = `awk -F , '/^[12]/{n1 = $1}END{print n1}' $name`
awk -F , '/^[12]/{if ($2>-99) print substr($1,1,4)+substr($1,5,2)/12.-1/24.,$2}' $name > NOAA.$date.temp
end
