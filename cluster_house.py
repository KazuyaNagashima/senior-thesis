import pandas as pd

#ヒトhousekeeping genesのデータをDataframeにする
#\tの前にスペースがあるので注意
df1 = pd.read_csv('human_housekeeping.csv', sep = ' \t',header = None)
df1 = df1.rename(columns={0:'Gene Name'})
df1 = df1.rename(columns={1:'ENSID'})

#clusterのアノテーションデータをDataframeにする
df2 = pd.read_csv('allcluster_list.csv', sep = ',')
df2 = df2.rename(columns={'Unnamed: 0':'Gene Name'})
df2 = df2.set_index('Gene Name')

#Dataframeの定義
df4 = pd.DataFrame(index=[0], columns=['Gene Name'])

#ハウスキーピング遺伝子のデータを一行ずつ抽出しkeyとする。keyをもとにclusterのアノテーションデータを検索し、Dataframeにする
for index, row in df1.iterrows():
    key = row[0]
    key = str(key)
    df3 = df2[df2.index == key]
    df4 = df4.append(df3,sort=False)
 
#pd.options.display.float_format = '{:,.0f}'.format

#clusterごとのハウスキーピング遺伝子の数のリストをDataframeにする
cl1 = df4['cluster'].value_counts()
#今回はcluster5のデータが0なので追加
cl1.loc[5] = 0
cl1.sort_index(inplace=True)
cl1 = pd.DataFrame(cl1)
cl1 = cl1.rename(columns={'cluster':'house_genes'})

#clusterごとの遺伝子の数のリストをDataframeにする
cl2 = df2['cluster'].value_counts()
cl2.sort_index(inplace=True)
cl2 = pd.DataFrame(cl2)
cl2 = cl2.rename(columns={'cluster':'cluster_genes'})

#clusterごとの遺伝子数、ハウスキーピング遺伝子のDataframeをmergeする
cl3 = pd.merge(cl1,cl2,left_index=True, right_index=True)

#ハウスキーピング遺伝子の数/遺伝子数　で比率をだす
cl3['house_rate'] = cl3['house_genes'] / cl3['cluster_genes']

#画像データを出力する
import matplotlib.pyplot as plt

x = range(1, 11)
plt.bar(x,cl3['house_rate'], tick_label = x, width=0.5)
plt.rcParams["font.size"] = 8
plt.title("housekeeping genes")
plt.ylabel("rate of housekeeping genes")
plt.xlabel("cluster")
plt.savefig("cluster_house.png", bbox_inches='tight')

