import sys

exitlist=[]
hideOrder =[]
with open(sys.argv[1],'r') as gtdbF, open(sys.argv[2],'r') as myF,open(sys.argv[3],'r') as orderF, open(sys.argv[4],'w') as outF:
    for a in myF:
        exitlist.append(a.strip().lstrip("p__"))
    n = 0
    for o in orderF:
        hideOrder.append(o.strip().lstrip("o__"))


    for b in gtdbF:
        if n ==0 :
            outF.write(b)
        else:
            bl=b.strip().split("\t")[0].split("|")
            if bl[1] in exitlist:
               if  bl[3] not in hideOrder:
                    outF.write(b)
        n+=1

