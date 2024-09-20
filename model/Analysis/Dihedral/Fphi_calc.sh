#!/bin/bash
Here=`pwd`

FILEDIR=$1
TOPFILE=${Here}/../files/1ENH.top
PDBfile=${Here}/../files/1ENH.pdb
GROfile=${Here}/../files/1ENH.gro
GROMACSDIR=/hpc/users/CONNECT/cxu807/opt/gromacs-4.5.7/bin
RUNROOT="${Here}/.."



declare -a Temp_seq
N_t=0
for line in `cat ${Here}/T.csv`
do
Temp_seq[$N_t]=$line
N_t=$[$N_t+1]
done


mkdir ${Here}/${FILEDIR}/dih

if [ ! -f "$PDBfile" ]; then
    echo "File not found!"
    exit 1
fi
ca_atom_count=$(grep -E '^ATOM|^HETATM' "$PDBfile" | awk '{if ($3 == "CA") print $0}' | wc -l)
dih_num=`expr ${ca_atom_count} - 3`
echo "二面角一共有 $dih_num 个"


i=1
while true; do
if [[ -e "${RUNROOT}/${i}_${FILEDIR}" ]]
then
echo "Found ${RUNROOT}/${i}_${FILEDIR}"
for(( Tn=0;Tn<${#Temp_seq[@]};Tn++))
do
mkdir ${Here}/${FILEDIR}/dih/T_${Tn}
Temp=${Temp_seq[$Tn]}
cat > ${Here}/${FILEDIR}/dih/T_${Tn}/dih_${i}.py <<EOF
import MDAnalysis as mda
import numpy as np
import math
import os
u = mda.Universe('${GROfile}')
ca_atoms = u.select_atoms('name CA')
n_atoms = len(ca_atoms)
X0 = np.zeros(n_atoms)
Y0 = np.zeros(n_atoms)
Z0 = np.zeros(n_atoms)
for atom_index, atom in enumerate(ca_atoms):
    X0[atom_index] = atom.position[0]
    Y0[atom_index] = atom.position[1]
    Z0[atom_index] = atom.position[2]
ATOMNUM=n_atoms
a_0=np.zeros(ATOMNUM-1)
b_0=np.zeros(ATOMNUM-1)
c_0=np.zeros(ATOMNUM-1)
r_0=np.zeros(ATOMNUM-1)
phi0=np.zeros(ATOMNUM-3)
F_0=np.zeros([ATOMNUM-2,3])
f_0=np.zeros([ATOMNUM-2,3])
xl_0=np.zeros([ATOMNUM-1,3])
rf_0=np.zeros(ATOMNUM-2)
for i in range(ATOMNUM-1):
    a_0[i]=(X0[i+1]-X0[i])
    b_0[i]=(Y0[i+1]-Y0[i])
    c_0[i]=(Z0[i+1]-Z0[i])
    r_0[i]=(a_0[i]**2.0+b_0[i]**2.0+c_0[i]**2.0)**(0.5)
    xl_0[i,0]=a_0[i]
    xl_0[i,1]=b_0[i]
    xl_0[i,2]=c_0[i]
for i in range(ATOMNUM-2):
    F_0[i]=(np.cross(xl_0[i,:],xl_0[i+1,:]))
    rf_0[i]=(F_0[i,0]**2+F_0[i,1]**2+F_0[i,2]**2)**(0.5)
    f_0[i,0]=F_0[i,0]/rf_0[i]
    f_0[i,1]=F_0[i,1]/rf_0[i]
    f_0[i,2]=F_0[i,2]/rf_0[i] 
m=np.zeros(3)
for i in range(ATOMNUM-3):
    m=np.cross(f_0[i,:],xl_0[i+1,:])/r_0[i+1]
    y=np.dot(m,f_0[i+1,:])
    x=np.dot(f_0[i,:],f_0[i+1,:])
    phi0_t=math.atan2(y,x)*360/(2*np.pi)
    if(phi0_t>0):
        phi0[i]=180-phi0_t+360
    elif(phi0_t<0):
        phi0[i]=180-phi0_t
print("get phi0")

u = mda.Universe('${GROfile}', '${RUNROOT}/${i}_${FILEDIR}/traj${Tn}.xtc')
n_frames = len(u.trajectory)
ca_atoms = u.select_atoms('name CA')
n_atoms = len(ca_atoms)

X0 = np.zeros(n_atoms)
Y0 = np.zeros(n_atoms)
Z0 = np.zeros(n_atoms)

X = np.zeros((n_atoms, n_frames))
Y = np.zeros((n_atoms, n_frames))
Z = np.zeros((n_atoms, n_frames))

for frame_index, ts in enumerate(u.trajectory):
    for atom_index, atom in enumerate(ca_atoms):
        X[atom_index, frame_index] = atom.position[0]
        Y[atom_index, frame_index] = atom.position[1]
        Z[atom_index, frame_index] = atom.position[2]

ATOMNUM=n_atoms
phi=np.zeros([n_frames,ATOMNUM-3])
for nt in range(n_frames):
    a=np.zeros(ATOMNUM-1)
    b=np.zeros(ATOMNUM-1)
    c=np.zeros(ATOMNUM-1)
    r=np.zeros(ATOMNUM-1)
    F=np.zeros([ATOMNUM-2,3])
    f=np.zeros([ATOMNUM-2,3])
    xl=np.zeros([ATOMNUM-1,3])
    rf=np.zeros(ATOMNUM-2)
    for i in range(ATOMNUM-1):
        a[i]=(X[i+1,nt]-X[i,nt])
        b[i]=(Y[i+1,nt]-Y[i,nt])
        c[i]=(Z[i+1,nt]-Z[i,nt])
        r[i]=(a[i]**2.0+b[i]**2.0+c[i]**2.0)**(0.5)
        xl[i,0]=a[i]
        xl[i,1]=b[i]
        xl[i,2]=c[i]
    for i in range(ATOMNUM-2):
        F[i]=(np.cross(xl[i,:],xl[i+1,:]))
        rf[i]=(F[i,0]**2+F[i,1]**2+F[i,2]**2)**(0.5)
        f[i,0]=F[i,0]/rf[i]
        f[i,1]=F[i,1]/rf[i]
        f[i,2]=F[i,2]/rf[i] 
    m=np.zeros(3)
    for i in range(ATOMNUM-3):
        m=np.cross(f[i,:],xl[i+1,:])/r[i+1]
        y=np.dot(m,f[i+1,:])
        x=np.dot(f[i,:],f[i+1,:])
        phi_t=math.atan2(y,x)*360/(2*np.pi)
        if(phi_t>0):
            phi[nt][i]=180-phi_t+360
        elif(phi_t<0):
            phi[nt][i]=180-phi_t


sigma=60.0
F=np.zeros([n_frames,ATOMNUM-3])
for i in range(n_frames):
    for k in range(ATOMNUM-3):
        delta0=(phi[i][k]-phi0[k])
        if delta0 >=-180 and delta0 < 180:
            delta=delta0
        elif delta0 < -180 :
            delta=delta0 + 360.0
        elif delta0 >= 180 :
            delta=delta0 - 360.0 
        F[i][k]=math.exp(-(delta*delta)/(sigma*sigma))
print("complete F")
for i in range(ATOMNUM-3):
    file="${Here}/${FILEDIR}/dih/T_${Tn}/"+"dih_"+'%s'%i+"_${i}.xvg"
    if os.path.exists(file):
        os.remove(file)

    try:
        os.mknod(file)
    except OSError as e:
        print(f"An error occurred while creating the file: {e}")
    Fphi = open(file,"a")
    for nt in range(n_frames):
        Fphi.write(str(F[nt][i]))
        Fphi.write('\n')
    Fphi.close() 
EOF

cat > ${Here}/${FILEDIR}/dih/T_${Tn}/dih_${i}.sh <<EOF

#!/bin/bash
#BSUB -J ${FILEDIR}_${Tn}_${i}   
#BSUB -n 1                      
#BSUB -q long_cpu                                         
#BSUB -o out.%J                
#BSUB -e err.%J    
#BSUB -cwd  ${Here}/${FILEDIR}/dih/T_${Tn}

######################集群调取cpu信息####################
source /hpc/jhinno/unischeduler/exec/unisched
#       \$NCPU:  Number of CPU cores                                  #
#       \$HOST_FILE:  List of computer hostfiles                      #
########################################################
module load anaconda3
##########软件运行命令 #####################################
python3  ${Here}/${FILEDIR}/dih/T_${Tn}/dih_${i}.py

EOF
bsub <  ${Here}/${FILEDIR}/dih/T_${Tn}/dih_${i}.sh 
done

else
echo "Not found ${RUNROOT}/x_${i}, exiting loop."
break 
fi
((i++)) 
done