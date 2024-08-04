import matplotlib.pyplot as plt
import numpy as np
from glob import glob
import os
from datetime import date
import shutil


def graph_raw_video(plot_ped_sig=True):
    files = glob("./*.csv")
    current_date = date.today()
    
    base_directory_name  = current_date.strftime("%Y-%m-%d")


    if not os.path.exists(base_directory_name):
        os.mkdir(base_directory_name)
    
    else:
        suffix = 2
        while True:
            base_directory_name = f"{base_directory_name}_{suffix}"
            if not os.path.exists(base_directory_name):
                os.mkdir(base_directory_name)
                break
        suffix += 1
    

    for i,f in enumerate(files):
        ped,pix_val = np.loadtxt(f,skiprows=0,usecols=(2,3),delimiter=',',unpack=True)
        pix_val = -1*pix_val
        range=np.linspace(0,len(pix_val),len(pix_val))
        fig, ax = plt.subplots(figsize=(8,8))
        plt.plot(range,pix_val,linewidth=1)
        plt.grid()
        plt.title("HDU={}".format(i), fontsize=20)
        plt.tick_params(labelsize=16)
        plt.ylabel("Pixel Value", fontsize=20)
        plt.xlabel("Range",fontsize=20)
    
    plt.show()

    if plot_ped_sig:
        norm_factor=1e5
        for i,f in enumerate(files):
            ped,pix_val = np.loadtxt(f,skiprows=0,usecols=(2,3),delimiter=',',unpack=True)
            pix_val=-1*pix_val
            range=np.linspace(0,len(pix_val),len(pix_val))
            fig, ax = plt.subplots(figsize=(8,8))
            plt.plot(range,pix_val/norm_factor,linewidth=1)
            plt.plot(range,ped,linewidth=1, color='red')
            plt.grid()
            plt.title("HDU={}".format(i), fontsize=20)
            plt.tick_params(labelsize=16)
            plt.ylabel("Pixel Value", fontsize=20)
            plt.xlabel("Range",fontsize=20)
    
    plt.show()

    
    files_to_move = os.listdir()
    for file in files_to_move:
        if os.path.isfile(file) and not file.endswith((".sh",".py")):
            shutil.move(file, os.path.join(base_directory_name,file))

if __name__ == "__main__":
    graph_raw_video()

