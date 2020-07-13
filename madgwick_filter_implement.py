# -*- coding: utf-8 -*-
"""
Created on Mon Jul 13 15:13:32 2020
@author: kanno
Madgwick filetr
"""
# import
import pandas as pd
import numpy as np
import math
from numpy.linalg import norm
from scipy.spatial.transform import Rotation as R
import matplotlib.pyplot as plt

# input
# input file path
file_path_i = "input/腕振り　右手系.csv"
# output file path
file_path_o = "output/腕振り　右手系_euler_test1.csv"
# save image path
img_path_o = "image/腕振り　右手系.png"

csv_input1 = pd.read_csv(filepath_or_buffer= file_path_i, encoding="ms932", sep=",")
array = csv_input1.values

time = array[:,0]
acc = array[:,1:4]
gyro_deg = array[:,4:7]
gyro = gyro_deg * (np.pi/180) 

# setting
Ts_imu = 0.05
N_imu = len(gyro)
#t = range(0,(N_imu-1)*Ts_imu,Ts_imu)
q = [1,0,0,0]

# Madgwick filter
q_mad = []
q_mad.append(q)
beta = math.sqrt(3/4)*math.pi*(5/180)
for i in range(1,N_imu)[:200]:
    if norm(acc[i,:]) != 0:
        acc_norm = acc[i,:]/norm(acc[i,:])
    else:
        acc_norm = acc[i,:]
    f_g = np.array([
        2*(q[1]*q[3] - q[0]*q[2]) - acc_norm[0],
        2*(q[0]*q[1] + q[2]*q[3]) - acc_norm[1],
        2*(0.5 - q[1]**2 - q[2]**2) - acc_norm[2]
    ])
    
    j_g = np.array([
        [-2*q[2], 2*q[3], -2*q[0], 2*q[1]],
        [2*q[1], 2*q[0], 2*q[3], 2*q[2]],
        [0, -4*q[1], -4*q[2], 0]
    ])
    
    nabla_f = np.dot(j_g.T,f_g)
    
    #44
    q_ep = nabla_f / norm(nabla_f) 
    
    #43
    imu_quatlized = np.array([
        [0,gyro[i,0],gyro[i,1],gyro[i,2]],
        [-gyro[i,0],0,-gyro[i,2],gyro[i,1]],
        [-gyro[i,1],-gyro[i,2],0,-gyro[i,0]],
        [-gyro[i,2],gyro[i,1],gyro[i,0],0]
        ])
    q_dot_t_omega = 0.5*np.dot(q_mad[i-1],imu_quatlized)
    q_dot_est = q_dot_t_omega - beta * q_ep.T 
    
    #42
    #q_new = q_mad[i-1] + q_dot_est*Ts_imu
    q_new = (q_mad[i-1] + q_dot_est*Ts_imu) / norm(q_mad[i-1] + q_dot_est*Ts_imu) #norm
    q_mad.append(q_new)
    
# to euler
def quaternion_to_euler_zyx(q):
    r = R.from_quat([q[0], q[1], q[2], q[3]])
    return r.as_euler('zyx', degrees=True)

euler = []
euler_z = []
euler0 = quaternion_to_euler_zyx(q)
euler.append(euler0)

for j in range(i):
    euler.append(quaternion_to_euler_zyx(q_mad[j]))
    
    # 1周に変換
    if euler[j][2] > 0:
        euler_z.append(euler[j][2]-180) 
    else:
        euler_z.append(euler[j][2]+180)

# plot
plt.plot(euler_z)
#plt.show()
plt.savefig(img_path_o)

# output to csv


"""
#q_new[np.isneginf(q_new)] = 100
#q_new = np.nan_to_num(q_new)
"""