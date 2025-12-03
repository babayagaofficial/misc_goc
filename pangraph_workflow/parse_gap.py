import glob
import os
import pandas as pd
import math

def parse_gap():
    status = {"solved":{},"objective":{}, "bound":{}, "gap":{}}
    for file in glob.glob("/home/daria/Documents/projects/ABC/pangraph_workflow/spp_dcj_sol/sol/*.log"):
        name = os.path.basename(file).replace(".log","").replace("_","\_").replace("community","c").replace("sub","s")
        with open(file, "r") as f:
            lines = f.readlines()[-2:]
            exit_line = lines[0].split(" ")
            if exit_line[0]=="Optimal":
                status["solved"][name]=True
            else:
                status["solved"][name]=False
            values = lines[1].split(" ")
            status["objective"][name]=int(float(values[2].replace(",","")))
            status["bound"][name]=int(math.ceil(float(values[5].replace(",",""))))
            try:
                gap = (status["objective"][name]-status["bound"][name])/status["objective"][name]*100
            except:
                gap = 0
            status["gap"][name]="{:.4f}".format(gap)+"\%"


    status_df = pd.DataFrame(status)  
    status_df.sort_values(by="gap", inplace=True)
    status_df.to_latex("status.tex")
