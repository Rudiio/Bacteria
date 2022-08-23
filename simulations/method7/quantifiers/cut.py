import os
import shutil
import glob

f = "simu12"

all_files = glob.glob("./energies/" + f + "/*.npy")


for file in all_files:
    print(file[-6:-4])
    if(file[-6:-4]=="El"):
        shutil.move(file,"./energies/" + f + "/El")
    if(file[-6:-4]=="Et"):
        shutil.move(file,"./energies/" + f + "/Et")