#!/bin/bash
Here=`pwd`
FILENAME=$1
POT=${Here}/${FILENAME}/pot_all
DIH=${Here}/${FILENAME}/dih
if [ ! -e ${Here}/${FILENAME}/dih_WHAM ]; then
mkdir ${Here}/${FILENAME}/dih_WHAM
fi

cat > ${POT}/all_calc.py <<EOF

import sys
import math

file_path = sys.argv[1]

with open(file_path, 'r') as file:
    numbers = [float(line.strip()) for line in file]

max_num = max(numbers)
min_num = min(numbers)

min_rounded = math.floor(min_num) 
max_rounded = math.ceil(max_num)

print(f"{min_rounded} {max_rounded}")
EOF

find ${POT} -type f -name '*.xvg' -print0 | xargs -0 cat > ${POT}/all_numbers.txt
module load anaconda3 

read pot_min pot_max <<< $(python3 ${POT}/all_calc.py ${POT}/all_numbers.txt)

# 在Shell中使用捕获的变量
echo "Min rounded down to nearest ten: $pot_min"
echo "Max rounded up to nearest ten: $pot_max"

if [ -z "$pot_max" ]; then
    echo "Error: Failed to read max rounded number"
    exit 1
fi

declare -a Temp_seq
N_t=0
for line in `cat $Here/T.csv`
do
Temp_seq[$N_t]=$line
N_t=$[$N_t+1]
done
mkdir ${Here}/${FILENAME}/dih_WHAM/count
mkdir ${Here}/${FILENAME}/dih_WHAM/WHAM
for dihn in {0..50}
do


mkdir ${Here}/${FILENAME}/dih_WHAM/count/dih_${dihn}
for(( Tn=0;Tn<${#Temp_seq[@]};Tn++))
do
Temp=${Temp_seq[$Tn]}
cat > ${Here}/${FILENAME}/dih_WHAM/count/dih_${dihn}/count_${Tn}.m <<EOF
Temp=${Temp}*8.314/1000;
Q_ini = 0;
Q_bins = 100;
Q_step = 0.01;

pot_init = ${pot_min};
pot_bins = ${pot_max}-(${pot_min});
pot_step = 1;

pot = dlmread('${POT}/${Tn}.xvg');
Q = dlmread('${DIH}/T_${Tn}/dih_all/dih_${dihn}.xvg');
p = length(pot);
count = zeros(Q_bins, pot_bins);

for i = 1:p
    m = 1;
    Qnum = floor(Q(i) / Q_step) + 1;
    Potnum = floor((pot(i) - pot_init) / pot_step) + 1;
    if Potnum > pot_bins
        Potnum = pot_bins;
    end
    if Qnum > Q_bins
        Qnum = Q_bins;
    end
    count(Qnum, Potnum) = count(Qnum, Potnum) + m;
end


writematrix(count, '${Here}/${FILENAME}/dih_WHAM/count/dih_${dihn}/count_${Tn}.csv');

EOF
#module load matlab2022 
#matlab -nosplash -nodesktop -r "run  ${Here}/${FILENAME}/dih_WHAM/count/dih_${dihn}/count_${Tn}.m;clear scriptname;quit;"
done


mkdir ${Here}/${FILENAME}/dih_WHAM/WHAM/dih_${dihn}

cat > ${Here}/${FILENAME}/dih_WHAM/WHAM/dih_${dihn}/wham.m <<EOF
count='${Here}/${FILENAME}/dih_WHAM/count/dih_${dihn}/';
reweighting(count);

function reweighting(count)
    temp_filename = fullfile('${Here}/T.csv');
    if ~isfile(temp_filename)
        error('温度文件不存在：%s', temp_filename);
    end
    Temp_seq = readmatrix(temp_filename);
    if isempty(Temp_seq)
        error('温度文件读取结果为空：%s', temp_filename);
    end
    if any(isnan(Temp_seq))
        error('温度文件包含NaN：%s', temp_filename);
    end

    filename = fullfile(count, 'count_0.csv');

    if ~isfile(filename)
        error('第一个直方图文件不存在：%s', filename);
    end

    fid = fopen(filename, 'rt');
    if fid == -1
        error('无法打开文件：%s', filename);
    else
        fclose(fid); % 如果可以打开，记得关闭文件
    end

    try
        data = readmatrix(filename, 'Range', 'B1');
    catch ME
        error('读取文件时出现错误：%s\n错误信息：%s', filename, ME.message);
    end


    if isempty(data)
        error('文件内容为空：%s', filename);
    elseif any(isnan(data), 'all')
        error('文件内容包含NaN：%s，已输出数据到%s', filename, output_filename);
    end

    dis = zeros(size(data, 1)-1, size(data, 2) , length(Temp_seq)); % 初始化 dis
    dims = size(dis);
    for kk = 1:length(Temp_seq)
        filename = fullfile(count, sprintf('count_%d.csv', kk-1));
        if ~isfile(filename)
            error('直方图文件不存在：%s', filename);
        end
        data = readmatrix(filename, 'Range', 'B2');
        if isempty(data) || any(isnan(data), 'all')
            error('直方图文件读取失败或包含NaN：%s', filename);
        end
        dis(:,:,kk) = data(:, :); % 假设第一列是不需要的索引列
        
    end


    slice = dis(:,:,1);

    filename = fullfile('${Here}/${FILENAME}/dih_WHAM/WHAM/dih_${dihn}','dis_0.list');

    dlmwrite(filename, slice, 'delimiter', ',', 'precision', 9);

    disp('Array dimensions:');
    disp(['Dimension 1 (rows): ', num2str(size(dis, 1))]);
    disp(['Dimension 2 (columns): ', num2str(size(dis, 2))]);
    disp(['Dimension 3 (slices): ', num2str(size(dis, 3))]);

    total_sum = sum(sum(dis(:,:,1)));



    Q_num = size(dis, 1); 
    pot_num = size(dis, 2);
    temperatures = Temp_seq * 8.314 / 1000;
    
    [n, New_p, A_Q, cv, F] = nnn('${Here}/${FILENAME}/dih_WHAM/WHAM/dih_${dihn}',temperatures, dis, Q_num, pot_num);
    
    writematrix(A_Q, fullfile('${Here}/${FILENAME}/dih_WHAM/WHAM/dih_${dihn}', 'A_Q.csv'));
    writematrix(n, fullfile('${Here}/${FILENAME}/dih_WHAM/WHAM/dih_${dihn}', 'n.csv'));
    writematrix(cv, fullfile('${Here}/${FILENAME}/dih_WHAM/WHAM/dih_${dihn}', 'cv.csv'));
    writematrix(F, fullfile('${Here}/${FILENAME}/dih_WHAM/WHAM/dih_${dihn}', 'F.csv'));
    disp('Reweighting completed.');
end



function [n, New_p, A_Q , cv, F] = nnn(Here,temperatures, dis, Q_num, pot_num)
    max_iter = 10000;
    toln = 1e-10;
    tolf = 1e-10;
    T0 = 1; % 参考温度，根据您的系统设置

    f = ones(size(temperatures)); % 初始权重f设置为1
    n = zeros(Q_num, pot_num); % 初始化态密度分布n

    N1 = squeeze(sum(dis, 3)); % 按照第三维（温度）进行求和
    N2 = squeeze(sum(sum(dis, 1), 2)); % 按照第一维和第二维（Q和pot）进行求和
    for iter = 1:max_iter
        old_n = n;
        old_f = f;

        for Qn = 1:Q_num
            for En = 1:pot_num
                NN_n = 0;
                for Tn = 1:length(temperatures)
                    NN_n = NN_n + N2(Tn)  * exp(f(Tn)-((1/temperatures(Tn)-1/T0)) * (En - 1 + 0.5));
                end
                n(Qn, En) = N1(Qn, En) / NN_n;
            end
        end

        if any(isnan(n), 'all') || any(isinf(n), 'all')
            error('在迭代过程中n包含NaN或Inf。');
            output_filename = fullfile(Here, 'n.csv');
            writematrix(n, output_filename);
            error('在迭代过程中n包含NaN或Inf。已输出数据到%s',  output_filename);



        end

        n_sum = sum(n, 'all');
        n = n / n_sum;





        for Tn = 1:length(temperatures)
            NN = 0;
            for Qn = 1:Q_num
                for En = 1:pot_num
                    NN = NN + n(Qn, En) * exp(-((1/temperatures(Tn))-1/T0) * ( En - 1 + 0.5));
                    if isnan(NN) 
                        disp(['NN: ', num2str(NN)]);
                        disp(n_value);
                        disp(exponential_term);
                    end

                end
            end
            
            f(Tn) = -log(NN);
        end



        if any(isinf(f))
            error('在迭代过程中f包含Inf');
        end
        if any(isnan(f))
            error('在迭代过程中f包含NaN');
       
        end

        condition1 = max(abs(old_n - n), [], 'all') < toln;
        condition2 = max(abs(old_f - f), [], 'all') < tolf;

        if condition1 && condition2
            nuu=nuu+1
            delta_n = max(abs(old_n - n), [], 'all');
            delta_f = max(abs(old_f - f), [], 'all');
            fprintf('Iteration %d: condition1=%d, condition2=%d, delta_n: %f delta_f: %f\n', iter, condition1, condition2, delta_n, delta_f);
            if nuu==3
                fprintf('Converged after %d iterations.\n', iter);
                break;
            end
        else
            delta_n = max(abs(old_n - n), [], 'all');
            delta_f = max(abs(old_f - f), [], 'all');
            nuu=0;
            fprintf('Iteration %d: condition1=%d, condition2=%d, delta_n: %f delta_f: %f\n', iter, condition1, condition2, delta_n, delta_f);
        end
    end

    T_bins = 40000; 
    T_init = min(temperatures);
    T_max = max(temperatures);
    New_T = linspace(T_init, T_max, T_bins);
    New_p = zeros(length(New_T), Q_num, pot_num);
    for Tn = 1:length(New_T)
        for Qn = 1:Q_num
            for En = 1:pot_num
                New_p(Tn, Qn, En) = n(Qn, En).*exp(-(1./New_T(Tn)-1/T0).*(En - 1 + 0.5));
            end
        end
        New_p(Tn, :, :) = New_p(Tn, :, :) ./ sum(New_p(Tn, :, :), 'all');
    end
    if any(isinf(New_p))
        error('New_p包含Inf');
    end
    if any(isnan(New_p))
        error('New_p包含NaN');

    end
    A_Q = zeros(length(New_T), 1);
    Q = (0:(Q_num-1)) * 0.01 + 0.005;

    for Tn = 1:length(New_T)
        for Qn = 1:Q_num
            A_Q(Tn) = A_Q(Tn) + Q(Qn) * sum(New_p(Tn, Qn, :));
        end
    end

    dis_Q = zeros(length(New_T), Q_num);
    dis_E = zeros(length(New_T), pot_num);
    for Tn = 1:length(New_T)
        for Qn = 1:Q_num
            dis_Q(Tn,Qn) = sum(New_p(Tn,Qn,:));
            if dis_Q(Tn,Qn)==0
                dis_Q(Tn,Qn)=1e-50;
            end
        end
    end
    disp(dis_Q(1,:));
    for Tn = 1:length(New_T)
        for En = 1:pot_num
            dis_E(Tn,En) = sum(New_p(Tn,:,En));
        end
    dis_E(Tn, :) = dis_E(Tn, :) ./ sum(dis_E(Tn,:), 'all');
    end
    disp(dis_E(1,:));
    

    cv=zeros(length(New_T),1);
    for Tn = 1:length(New_T)
        E2=0;
        Emean=0;
        for En = 1:pot_num
            E2=E2+dis_E(Tn,En)*( En - 1 + 0.5+(${pot_min}))*( En - 1 + 0.5 +(${pot_min}));
            Emean=Emean+dis_E(Tn,En)*( En - 1 + 0.5 +(${pot_min}));
        end
        cv(Tn)=(E2-Emean*Emean)/(New_T(Tn)*New_T(Tn));
    end

    F=zeros(length(New_T),Q_num);
    for Tn = 1:length(New_T)
        F(Tn,1)=0;
        minss=100;
        for Qn = 2:Q_num
            F(Tn,Qn)=F(Tn,1)-New_T(Tn)*log(dis_Q(Tn,Qn)/dis_Q(Tn,1));
            if F(Tn,Qn)<minss
                minss=F(Tn,Qn);
            end
        end
        for Qn = 1:Q_num
            F(Tn,Qn)=F(Tn,Qn)-minss;
        end
    end




end




EOF
done
#module load matlab2022
#matlab -nosplash -nodesktop -r "run ${FILENAME}/dih_WHAM/WHAM/dih_${dihn}/wham.m;clear scriptname;quit;"

