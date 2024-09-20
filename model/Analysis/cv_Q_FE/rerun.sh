#!/bin/bash
Here=`pwd`
filename=$1
FILEDIR=${filename}
ionic=$2
#####origin:0.10
GROMACSDIR="/hpc/users/CONNECT/cxu807/opt/gromacs-4.5.7/bin"
Tablefile=$Here/../tablefiles/Table_file_DH_C_0.${ionic}.xvg
LD_LIBRARY_PATH=/hpc/users/CONNECT/cxu807/opt/gromacs-4.5.7/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH
mkdir ${Here}/${filename}/rerun/
mkdir ${Here}/${filename}/rerun/${ionic}
declare -a Temp_seq
N=0
for line in `cat $Here/T.csv`
do
Temp_seq[$N]=$line
N=$[$N+1]
done

i=1
while true; do
if [[ -e "${Here}/../${i}_${FILEDIR}" ]]
then
echo "Found ${Here}/../${i}_${FILEDIR}"
RUNROOT=${Here}/../${i}_${FILEDIR}
mkdir ${Here}/${filename}/rerun/${ionic}/${i}
for(( Tn=0;Tn<${#Temp_seq[@]};Tn++))
do
mkdir ${Here}/${filename}/rerun/${ionic}/${i}/T_${Tn} 
cat > ${Here}/${filename}/rerun/${ionic}/${i}/T_${Tn}/Jobs_rerun.sh <<EOF
#!/bin/bash             
#BSUB -J rerun_${i}_${ionic}_${Tn}
#BSUB -n 1                      
#BSUB -W 72:00
#BSUB -q long_cpu                                         
#BSUB -o out.%J                
#BSUB -e err.%J    
#BSUB -cwd ${Here}/${filename}/rerun/${ionic}/${i}/T_${Tn} 
######################集群调取cpu信息####################
source /hpc/jhinno/unischeduler/exec/unisched
#       \$NCPU:  Number of CPU cores                                  #
#       \$HOST_FILE:  List of computer hostfiles                      #
########################################################
LD_LIBRARY_PATH=/hpc/users/CONNECT/cxu807/opt/gromacs-4.5.7/lib:\$LD_LIBRARY_PATH
export LD_LIBRARY_PATH
${GROMACSDIR}/mdrun_mpi -rerun ${RUNROOT}/traj${Tn}.xtc -s ${RUNROOT}/md${Tn}.tpr  -table $Tablefile -tablep $Tablefile 
echo 6 8 |$GROMACSDIR/g_energy_mpi -f  ${Here}/${filename}/rerun/${ionic}/${i}/T_${Tn}/ener.edr     -o ${Here}/${filename}/rerun/${ionic}/pot_${Tn}_${i}.xvg -sum
done
EOF
bsub < ${Here}/${filename}/rerun/${ionic}/${i}/T_${Tn}/Jobs_rerun.sh
done

else
echo "Not found ${Here}/../x_${i}, exiting loop."
break # 如果文件不存在，退出循环
fi
((i++)) # 增加数字，准备检查下一个
done

echo "Done checking."


