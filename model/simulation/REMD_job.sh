#!/bin/bash

ionic=0.10


RUNROOT=`pwd`
GROMACSDIR="/hpc/users/CONNECT/cxu807/opt/gromacs-4.5.7/bin"
Filename=1ENH
TOPfile=$RUNROOT/../files/R01.top
GROfile=$RUNROOT/../files/${Filename}.gro

Tablefile=$RUNROOT/../tablefiles/Table_file_DH_C_${ionic}.xvg
mpirun="/usr/mpi/gcc/openmpi-4.1.2a1/bin/mpirun"

Steps=2000000000 #1000ns 

mkdir mdpfiles

LD_LIBRARY_PATH=/hpc/users/CONNECT/cxu807/opt/gromacs-4.5.7/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH

sed "s/XSteps/${Steps}/g" $RUNROOT/../CACB.mdp > mdpfiles/tmp.mdp


declare -a Temp_seq
N=0
for line in `cat $RUNROOT/T.csv`
do
Temp_seq[$N]=$line
N=$[$N+1]
done


for(( i=0;i<${#Temp_seq[@]};i++))
do
    Temp=${Temp_seq[i]}
	Num=$i
	sed "s/XTemp/$Temp/g" $RUNROOT/mdpfiles/tmp.mdp > ${RUNROOT}/mdpfiles/md${Num}.mdp
	MDPfile="$RUNROOT/mdpfiles/md${Num}.mdp"
	grompp_mpi -f $MDPfile -c $GROfile -p $TOPfile -o $RUNROOT/md${Num}.tpr
done

cat > Jobs_REMD.sh <<EOF
#!/bin/bash             
#BSUB -J ${ionic}_R01
#BSUB -n $N                      
#BSUB -W 72:00
#BSUB -q long_cpu                                         
#BSUB -o out.%J                
#BSUB -e err.%J    
#BSUB -cwd $RUNROOT

######################集群调取cpu信息####################
source /hpc/jhinno/unischeduler/exec/unisched
#       \$NCPU:  Number of CPU cores                                  #
#       \$HOST_FILE:  List of computer hostfiles                      #
########################################################
LD_LIBRARY_PATH=/hpc/users/CONNECT/cxu807/opt/gromacs-4.5.7/lib:\$LD_LIBRARY_PATH
export LD_LIBRARY_PATH

$mpirun -np $N $GROMACSDIR/mdrun_mpi -s md.tpr  -replex 2000 -multi $N  -table $Tablefile -tablep $Tablefile  -cpi state.cpt –append 
EOF

bsub < Jobs_REMD.sh
