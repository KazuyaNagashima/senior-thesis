import csv
from itertools import islice
import pandas as pd
import mygene

#dataframeの定義
cols = ["cluster"]
cv = pd.DataFrame(index=[], columns=cols)

#Rで取り出した各clusterの遺伝子リストは別々のファイルなのでfor文でmygeneを回す
for i in range(1,11):
    i = str(i)
    msig_file = "trimmed_genename_cluster"+i+".txt"
    sl = list()

    with open(msig_file, 'r') as fl:
        reader = csv.reader(fl,delimiter='\n')
        for row in islice(reader, 0, None):
            sl.extend(row)
    
    mg = mygene.MyGeneInfo()
    ref = mg.querymany(sl, 'name,symbol',as_dataframe=True,species='pig')  
    
    df = pd.DataFrame(ref)
    df = df[df['_id'].str.match('^\d*$', na=False)]
    df = df[~df['symbol'].str.match('^LOC', na=False)]
    df = df[df.index == df['symbol']]
    df["cluster"] = i

    df.to_csv("cluster"+i+"_list.csv") 

    cv = cv.append(df)
    

cv.to_csv("allcluster_list.csv")
