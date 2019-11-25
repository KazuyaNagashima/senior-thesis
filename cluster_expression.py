import csv
import pandas as pd

#RNAseqで得られたデータをDataframeにする。stringtieでアセンブリしたため、Non-cordingRNA領域なども含まれるのでトリミングを行う。
df = pd.read_csv('ERR3153937.abund.tab', sep = '\t')
df = df[~df['Gene Name'].str.match('^ENSSSCG', na=False)]
df = df[~df['Gene Name'].str.match('^RF\d*$', na=False)]
df = df[~df['Gene Name'].str.match('^ssc-.*', na=False)]
df = df[~df['Gene Name'].str.match('-', na=False)]

#アノテーションされたcluster listをDataframeにする
cl = pd.read_csv('allcluster_list.csv', sep = ',')
cl = cl.rename(columns={'symbol':'Gene Name'})

#Gene Nameをkeyに２つのDataframeをmergeする
al = pd.merge(cl,df,on='Gene Name')

#dataframeの定義
cv = pd.DataFrame(index=['expression'], columns=[])
ab = pd.DataFrame(index=[])

#clusterごとのTPM,FPKMの平均値を求め、Dataframeにする
for i in range(1,11):
    dk = al[al['cluster'] == i]
    cv['TPM'] = dk['TPM'].mean()
    cv['FPKM'] = dk['FPKM'].mean()
    i = str(i)
    cv['cluster'] = i
    cv = cv.set_index('cluster')
    
    ab = ab.append(cv)

#画像データの出力
import matplotlib.pyplot as plt

label = ab.index
x = range(1, 11)
plt.bar(x,ab['TPM'], tick_label = label, width=0.5)
plt.rcParams["font.size"] = 8
plt.title("expression")
plt.ylabel("average TPM")
plt.xlabel("cluster")
plt.savefig("cluster_TPM.png", bbox_inches='tight')

label = ab.index
x = range(1, 11)
plt.bar(x,ab['FPKM'], tick_label = label, width=0.5)
plt.rcParams["font.size"] = 8
plt.title("expression")
plt.ylabel("average FPKM")
plt.xlabel("cluster")
plt.savefig("cluster_FPKM.png", bbox_inches='tight')
