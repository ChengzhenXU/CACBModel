#!/bin/bash
Here=`pwd`

FILENAME=$1
ionic=$2

if [ ! -e ${Here}/${FILENAME}/TCI ]; then 
mkdir ${Here}/${FILENAME}/TCI
fi
mkdir ${Here}/${FILENAME}/TCI/${ionic}
mkdir ${Here}/${FILENAME}/TCI/${ionic}/Pphi
#dihn=$3
for dihn in {0..50}
do
#DIH=${Here}/${FILENAME}/dih_WHAM/WHAM/dih_${dihn}
DIH=${Here}/${FILENAME}/re_F/${ionic}/WHAM/dih_${dihn}
mkdir ${Here}/${FILENAME}/TCI/${ionic}/Pphi/dih_${dihn}
cat > ${Here}/${FILENAME}/TCI/${ionic}/Pphi/dih_${dihn}/Sigmoidf.m <<EOF
%%%loading data
dataFilePath = '${DIH}/A_Q.csv';
y = dlmread(dataFilePath, '', 0, 0);
y = y(:); % 确保 y 是列向量

min_x = 20 * 8.314 / 1000;
max_x = 320 * 8.314 / 1000;

n = length(y);

x = linspace(min_x, max_x, n).'; % 确保 x 是列向量

x_f = x(1:10);
y_f = y(1:10);
p_f = polyfit(x_f, y_f, 1);
f_f_fit = p_f(1)*x + p_f(2);
%%%Unfolding Curves
x_u = x(end-9:end); % Assuming 1501 is the total number (row) of the data
y_u = y(end-9:end);
p_u = polyfit(x_u, y_u, 1);
f_u_fit = p_u(1)*x + p_u(2);

%%%Fitting function
sigfunc = @(A, x)((p_f(2) + p_f(1)*(x-0)) .* 1./(1 + exp((-A(1) + x.*A(1)./A(2)) ./x)) + (p_u(2) + p_u(1).*x) - (p_u(2) + p_u(1).*(x-0)) .* 1./(1 + exp((-A(1) + x.*A(1)./A(2)) ./x)));

A0 = [2 0.8];
A_fit = nlinfit(x, y, sigfunc, A0);
F = sigfunc(A_fit, x);

F_p = 1.00 ./(1 + exp((-A_fit(1) + x*A_fit(1)/A_fit(2)) ./x));
f_f_fit = p_f(1)*x + p_f(2);
f_u_fit = p_u(1)*x + p_u(2);

Raw = [x, y];
save('${Here}/${FILENAME}/TCI/${ionic}/Pphi/dih_${dihn}/Raw.dat','Raw','-ascii');
Baseline_u = [x, f_u_fit];
save('${Here}/${FILENAME}/TCI/${ionic}/Pphi/dih_${dihn}/Baseline_u.dat','Baseline_u','-ascii');
Baseline_f = [x, f_f_fit];
save('${Here}/${FILENAME}/TCI/${ionic}/Pphi/dih_${dihn}/Baseline_f.dat','Baseline_f','-ascii');
F_fit = [x, F];
save('${Here}/${FILENAME}/TCI/${ionic}/Pphi/dih_${dihn}/F_fit.dat','F_fit','-ascii');
F_p = [x, F_p];
save('${Here}/${FILENAME}/TCI/${ionic}/Pphi/dih_${dihn}/F_p.dat','F_p','-ascii');
save('${Here}/${FILENAME}/TCI/${ionic}/Pphi/dih_${dihn}/Fitting.dat','A_fit','-ascii');

EOF


done

cat > ${Here}/${FILENAME}/TCI/${ionic}/TCI.py <<EOF
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import  MultipleLocator
from pandas import DataFrame
import math
from cycler import cycler
import seaborn as sns
from sklearn import datasets 
from matplotlib.ticker import  MultipleLocator
from pandas import DataFrame

P=np.zeros(51)
EOF

for i in {0..50}
do
cat >> ${Here}/${FILENAME}/TCI/${ionic}/TCI.py <<EOF
F_${i}=[]
T_${i}=[]
num=0
sum=0
with open("${Here}/${FILENAME}/TCI/${ionic}/Pphi/dih_${i}/F_p.dat") as f:
    for line in f:
        cols = line.split()
        T_${i}.append(float(cols[0]))
        F_${i}.append(float(cols[1]))
        num=num+1
EOF
done


cat >> ${Here}/${FILENAME}/TCI/${ionic}/TCI.py <<EOF
c=np.zeros([51,51,num])
EOF


for((i=0;i<51;i++))
do
for((j=0;j<51;j++))
do
cat >> ${Here}/${FILENAME}/TCI/${ionic}/TCI.py <<EOF
for k in range(num):
    c[${i}][${j}] = F_${i}[k]-F_${j}[k]
EOF
done
done


cat >> ${Here}/${FILENAME}/TCI/${ionic}/TCI.py <<EOF
m=0
mi=0

TCI=np.zeros([51,51])
for i in range(51):
    for j in range(51):
        if i!=j :
            I=np.dot(c[i][j][:],c[i][j][:])
            if I==0:
                I=1e-100
                print(i,j)
            TCI[i][j]=-math.log(math.sqrt(I))
        if i==j:
            TCI[i][j]==0

        if TCI[i][j]<mi:
            mi=TCI[i][j] 
        if TCI[i][j]>m:
            m=TCI[i][j]
print("${FILENAME},${ionic}:",m,mi)
df1=DataFrame(TCI)
df1.to_csv('${Here}/${FILENAME}/TCI/${ionic}/TCI.csv',index=False)

EOF



