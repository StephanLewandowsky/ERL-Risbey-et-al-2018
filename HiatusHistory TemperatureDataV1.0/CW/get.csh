wget http://www-users.york.ac.uk/~kdc3/papers/coverage2013/had4_krig_long_v1_0_0.txt
wget http://www-users.york.ac.uk/~kdc3/papers/coverage2013/had4_krig_v2_0_0.txt
foreach name (*.txt)
set date = `awk '/^[12]/{n=$1;n1=int(n);n2=int(12*(n-n1)+1)}END{printf "%04d%02d",n1,n2}' $name`
awk '{print $1,$2}' $name > $name:r.$date.temp
end
