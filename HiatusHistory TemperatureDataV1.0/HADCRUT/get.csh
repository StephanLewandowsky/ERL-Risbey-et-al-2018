wget http://www.metoffice.gov.uk/hadobs/hadcrut3/data/HadCRUT3.gz
wget http://www.metoffice.gov.uk/hadobs/hadcrut3/data_dec_2009/HadCRUT3_archive_dec_2009.gz
wget http://www.metoffice.gov.uk/hadobs/hadcrut4/data/current/time_series/HadCRUT.4.5.0.0.monthly_ns_avg.txt
wget http://www.metoffice.gov.uk/hadobs/hadcrut4/data/4.4.0.0/time_series/HadCRUT.4.4.0.0.monthly_ns_avg.txt
wget http://www.metoffice.gov.uk/hadobs/hadcrut4/data/4.3.0.0/time_series/HadCRUT.4.3.0.0.monthly_ns_avg.txt
wget http://www.metoffice.gov.uk/hadobs/hadcrut4/data/4.2.0.0/time_series/HadCRUT.4.2.0.0.monthly_ns_avg.txt
wget http://www.metoffice.gov.uk/hadobs/hadcrut4/data/4.1.1.0/time_series/HadCRUT.4.1.1.0.monthly_ns_avg.txt
wget http://www.metoffice.gov.uk/hadobs/hadcrut4/data/4.0.0.0/time_series/hadcrut4_monthly_ns_avg.txt

gunzip HadCRUT3.gz
gunzip HadCRUT3_archive_dec_2009.gz
python tempns.py HadCRUT3.txt | awk '{printf "%8.3f %8.4f\n",$1,$2}' > HadCRUT3_201405.temp
python tempns.py HadCRUT3_archive_dec_2009.txt | awk '{printf "%8.3f %8.4f\n",$1,$2}' > HadCRUT3_200912.temp

foreach name (HadCRUT.4.*.txt)
awk '{print substr($0,1,4)+substr($0,6,2)/12.0-1/24.0, $2}' $name > $name.temp
end

