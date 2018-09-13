from commons.utils import file2stringlist

def dump_template(ORIG, outfile,SAT_OBS,ARGO):
    LINES=[]
    for line in ORIG:
        newline=line
        if (line.find("@@SAT_OBS@@") != -1):  newline=line.replace("@@SAT_OBS@@",SAT_OBS)
        if (line.find("@@ARGO@@")   != -1):   newline=line.replace("@@ARGO@@"   ,  ARGO)
        LINES.append(newline + "\n")    

    fid=open(filename,"w")
    fid.writelines(LINES)
    fid.close()





ORIG=file2stringlist('satfloat.template')
SAT_OBS="1"
ARGO   ="1"
filename="satfloat.20150102.nml"

dump_template(ORIG, filename, SAT_OBS, ARGO)





