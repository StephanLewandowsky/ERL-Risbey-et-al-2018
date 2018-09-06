\rm Land_and_Ocean_complete.*

wget http://web.archive.org/web/20140503051102/http://berkeleyearth.lbl.gov/auto/Global/Land_and_Ocean_complete.txt
wget http://web.archive.org/web/20150217123234/http://berkeleyearth.lbl.gov/auto/Global/Land_and_Ocean_complete.txt
wget http://web.archive.org/web/20150325094341/http://berkeleyearth.lbl.gov/auto/Global/Land_and_Ocean_complete.txt
wget http://web.archive.org/web/20150406185813/http://berkeleyearth.lbl.gov/auto/Global/Land_and_Ocean_complete.txt
wget http://web.archive.org/web/20150511073225/http://berkeleyearth.lbl.gov/auto/Global/Land_and_Ocean_complete.txt
wget http://web.archive.org/web/20150525014933/http://berkeleyearth.lbl.gov/auto/Global/Land_and_Ocean_complete.txt
wget http://web.archive.org/web/20150717222850/http://berkeleyearth.lbl.gov/auto/Global/Land_and_Ocean_complete.txt
wget http://web.archive.org/web/20150906073402/http://berkeleyearth.lbl.gov/auto/Global/Land_and_Ocean_complete.txt
wget http://web.archive.org/web/20160304030546/http://berkeleyearth.lbl.gov/auto/Global/Land_and_Ocean_complete.txt
wget http://web.archive.org/web/20160328191324/http://berkeleyearth.lbl.gov/auto/Global/Land_and_Ocean_complete.txt
wget http://web.archive.org/web/20160406211459/http://berkeleyearth.lbl.gov/auto/Global/Land_and_Ocean_complete.txt
wget http://web.archive.org/web/20160620073552/http://berkeleyearth.lbl.gov/auto/Global/Land_and_Ocean_complete.txt
wget http://web.archive.org/web/20160725185749/http://berkeleyearth.lbl.gov/auto/Global/Land_and_Ocean_complete.txt
wget http://web.archive.org/web/20161110213315/http://berkeleyearth.lbl.gov/auto/Global/Land_and_Ocean_complete.txt/
wget http://web.archive.org/web/20161229074317/http://berkeleyearth.lbl.gov/auto/Global/Land_and_Ocean_complete.txt
wget http://web.archive.org/web/20170201143221/http://berkeleyearth.lbl.gov/auto/Global/Land_and_Ocean_complete.txt


foreach name (Land_and_Ocean_complete.*)
set date = `awk '/^ *[12]/{n1 = $1; n2 = $2}END{printf "%04d%02d",n1,n2}' $name`
awk '/^ *[12]/{if ($3>-99) print $1+$2/12.-1/24.,$3}/from Water Temperatures/{exit}' $name > BERKELEY.$date.temp
end
