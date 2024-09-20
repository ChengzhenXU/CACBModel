#!/bin/bash
Here=`pwd`
FILEDIR=$1
TOPFILE=${Here}/../files/1ENH.top
PDBfile=${Here}/../files/1ENH.pdb
if [ ! -f ${FILEDIR} ]; then
mkdir ${FILEDIR}
fi
echo "######reading Temperature######"
declare -a Temp_seq

N_t=0
for line in `cat ${Here}/T.csv`
do
Temp_seq[$N_t]=$line
N_t=$[$N_t+1]
done


echo "######getting pairs.txt###### "

cat > get_pairs.sh <<EOF

outputfile="pairs.txt"

copying=false

while IFS= read -r line
do
    if [[ "\$line" == "[ pairs ]"* ]]; then
        copying=true
        continue
    fi

    if \$copying && [[ "\$line" == "["* ]]; then
        break
    fi

    if \$copying && [[ "\$line" != ";"* ]]; then
        echo "\$line" >> "\$outputfile"
    fi
done < "$TOPFILE"

echo "Pairs extracted to \$outputfile"
EOF

if [ ! -f "$Here/pairs.txt" ]; then
    sh get_pairs.sh
else
    echo "File pairs.txt exists."
fi

rm  get_pairs.sh


echo "#######Calculating Q3######"
if [ ! -f ${Here}/${FILEDIR}/Q ]; then
mkdir ${Here}/${FILEDIR}/Q
fi
if [ ! -f ${Here}/${FILEDIR}/pot ]; then
mkdir ${Here}/${FILEDIR}/pot
fi
if [ ! -f ${Here}/${FILEDIR}/col_pot ]; then
mkdir ${Here}/${FILEDIR}/col_pot
fi
GROMACSDIR=/hpc/users/CONNECT/cxu807/opt/gromacs-4.5.7/bin
RUNROOT="${Here}/.."

i=1
while true; do
if [[ -e "${RUNROOT}/${i}_${FILEDIR}" ]]
then
echo "Found ${RUNROOT}/${i}_${FILEDIR}"
for(( Tn=0;Tn<${#Temp_seq[@]};Tn++))
do
Temp=${Temp_seq[$Tn]}
if [ ! -f ${Here}/${FILEDIR}/Q/${Tn}_${i}.csv ]; then
cat > ${Here}/${FILEDIR}/Q/Q_${Tn}_${i}.py <<EOF
import MDAnalysis as mda
import numpy as np
import re
from pandas import DataFrame
import math
pairfile_path = '${Here}/pairs.txt'
beta=5.0
contact = []
current_contact = None
pdb_path= '${PDBfile}'
ATOMNUM=0
with open(pdb_path,'r') as pdb:
        X0,Y0,Z0,R0=[],[],[],[]
        for line in pdb:
                cols = line.split()
                if cols[0]=="ATOM":
                        X0.append((float(cols[6])))
                        Y0.append((float(cols[7])))
                        Z0.append((float(cols[8])))
                        ATOMNUM=ATOMNUM+1


CONNUM=0
with open(pairfile_path, 'r') as file:
    for line in file:
        current_contact = contact
        if line.startswith(';') or not line.strip():
            continue
        parts = line.split()
        if len(parts) >= 5:
            n1=int(float(parts[0])-1)
            n2=int(float(parts[1])-1)
            dis=math.sqrt(( ( float ( X0 [ int(n1)] )-float(X0[int(n2)]))**2.0 +(float(Y0[int(n1)])-float(Y0[int(n2)]))**2.0 + (float(Z0[int(n1)])-float(Z0[int(n2)]))** 2.0 ))
            data = [int(parts[0]), int(parts[1]),dis]
            current_contact.append(data)
            CONNUM=CONNUM+1
print("number of contact pairs:",CONNUM,len(contact[0]))
contact = list(map(list, zip(*contact)))
Q=[]


u = mda.Universe(pdb_path, '${Here}/../${i}_${FILEDIR}/traj${Tn}.xtc')
for ts in u.trajectory:
    Qi=0
    for i in range(CONNUM):
        a1=int(contact[0][i]-1)
        a2=int(contact[1][i]-1)
        r0=contact[2][i]
        atom1 = u.atoms[a1]  
        atom2 = u.atoms[a2]
        distance = np.linalg.norm(atom1.position - atom2.position)
        EXP_THRESHOLD = 700
        exp_argument = beta * (distance - 1.2 * float(r0))
        if exp_argument > EXP_THRESHOLD:
            qi = 0
        else:
            qi = 1.0 / (1.0 + math.exp(exp_argument))
        Qi=Qi+qi
    Q.append(Qi)
Q_fin=[]
for i in range(len(Q)):
    Q_fin.append(float(Q[i]/CONNUM))

df=DataFrame(Q_fin)
df.to_csv('${Here}/${FILEDIR}/Q/${Tn}_${i}.csv',index=False)
print('\r')
EOF
cat > ${Here}/${FILEDIR}/Q/job_${Tn}_${i}.sh <<EOF
#!/bin/bash            
#BSUB -J Q_${FILEDIR}_${i}_${Tn}
#BSUB -n 1                      
#BSUB -q long_cpu                                        
#BSUB -o out.%J                
#BSUB -e err.%J    
#BSUB -cwd ${Here}/${FILEDIR}/Q
######################集群调取cpu信息####################
source /hpc/jhinno/unischeduler/exec/unisched
#       \$NCPU: Number of CPU cores                                  #
#       \$HOST_FILE:  List of computer hostfiles                     #
########################################################
LD_LIBRARY_PATH=/hpc/users/CONNECT/cxu807/opt/gromacs-4.5.7/lib:\$LD_LIBRARY_PATH
export LD_LIBRARY_PATH
module load anaconda3
##########软件运行命令 ###################################
python3 Q_${Tn}_${i}.py 
EOF
bsub < ${Here}/${FILEDIR}/Q/job_${Tn}_${i}.sh
fi
done
echo "#######Q done######"

echo "#######Potential#####"

if [ ! -f ${Here}/${FILEDIR}/pot/pot_0_${i}.xvg ]; then

cat >  ${Here}/${FILEDIR}/pot/potential_${FILEDIR}_${i}.sh <<BBB
#!/bin/bash
#BSUB -J pot_${FILEDIR}_${i}   
#BSUB -n 1                      
#BSUB -q long_cpu                                         
#BSUB -o out.%J                
#BSUB -e err.%J    
#BSUB -cwd  ${Here}/${FILEDIR}/pot

######################集群调取cpu信息####################
source /hpc/jhinno/unischeduler/exec/unisched
#       \$NCPU:  Number of CPU cores                                  #
#       \$HOST_FILE:  List of computer hostfiles                      #
########################################################
LD_LIBRARY_PATH=/hpc/users/CONNECT/cxu807/opt/gromacs-4.5.7/lib:\$LD_LIBRARY_PATH
export LD_LIBRARY_PATH
##########软件运行命令 ###################################

for(( Num=0;Num<${#Temp_seq[@]};Num++))
do
echo 3 4 5 6 7 8 |$GROMACSDIR/g_energy_mpi -f  ${Here}/../${i}_${FILEDIR}/ener\${Num}.edr     -o ${Here}/${FILEDIR}/pot/pot_\${Num}_${i}.xvg -sum
done
BBB
bsub < ${Here}/${FILEDIR}/pot/potential_${FILEDIR}_${i}.sh
fi
if [ ! -f ${Here}/${FILEDIR}/col_pot/col_pot_0_${i}.xvg ]; then
cat >  ${Here}/${FILEDIR}/col_pot/potential_${FILEDIR}_${i}.sh <<BBB
#!/bin/bash
#BSUB -J col_${FILEDIR}_${i}   
#BSUB -n 1                      
#BSUB -q long_cpu                                         
#BSUB -o out.%J                
#BSUB -e err.%J    
#BSUB -cwd  ${Here}/${FILEDIR}/col_pot

######################集群调取cpu信息####################
source /hpc/jhinno/unischeduler/exec/unisched
#       \$NCPU:  Number of CPU cores                                  #
#       \$HOST_FILE:  List of computer hostfiles                      #
########################################################
LD_LIBRARY_PATH=/hpc/users/CONNECT/cxu807/opt/gromacs-4.5.7/lib:\$LD_LIBRARY_PATH
export LD_LIBRARY_PATH
##########软件运行命令 ###################################

for(( Num=0;Num<${#Temp_seq[@]};Num++))
do
echo 6 8 |$GROMACSDIR/g_energy_mpi -f  ${Here}/../${i}_${FILEDIR}/ener\${Num}.edr     -o ${Here}/${FILEDIR}/col_pot/col_pot_\${Num}_${i}.xvg -sum
done
BBB
bsub < ${Here}/${FILEDIR}/col_pot/potential_${FILEDIR}_${i}.sh
fi
echo "######potential done####"
else
echo "Not found ${RUNROOT}/x_${i}, exiting loop."
break 
fi
((i++)) 
done

echo "Done checking."


