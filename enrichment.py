import pandas as pd

import matplotlib.pyplot as plt
matplotlib.use('Agg')
%matplotlib inline

#データのセット
df = pd.read_csv('fibro3min-190515_peaks.annotated.txt', sep = '\t')
ip = pd.read_csv('ERR3154128_peaks.annotated.txt', sep="\t")
# pandasはデフォルトでは１行目をヘッダーにする。今回のようにそれが不要の場合はheader = Noneにする

#無駄な文字を消す(sample)
df['Annotation'].replace('\(.*?\)', '', regex=True, inplace=True)
df['Annotation'].replace('\.\d*', '', regex=True, inplace=True)
df['Annotation'].replace(' $', '' ,regex = True, inplace=True)
df['Annotation'].replace('promoter-', '' ,regex = True, inplace=True)
df['Annotation'].unique()

#無駄な文字を消す(input)
ip['Annotation'].replace('\(.*?\)', '', regex=True, inplace=True)
ip['Annotation'].replace('\.\d*', '', regex=True, inplace=True)
ip['Annotation'].replace(' $', '' ,regex = True, inplace=True)
ip['Annotation'].replace('promoter-', '' ,regex = True, inplace=True)
ip['Annotation'].unique()

#sample,inputのアノテーション数のDataframeを作成
dc = pd.DataFrame({'sample_peaks':df['Annotation'].value_counts(), 'input_peaks':ip['Annotation'].value_counts()})

#enrichment = sampleの個数の比率/inputの個数の比率
dc['sample_rate'] = dc['sample_peaks'] / dc['sample_peaks'].sum()
dc['input_rate'] = dc['input_peaks'] / dc['input_peaks'].sum()
dc['enrichment'] = dc['sample_rate'] / dc['input_rate']

#画像データの出力
label = dc.index
x = range(0, 8)
plt.bar(x,dc['enrichment'], tick_label = label, width=0.5)
plt.rcParams["font.size"] = 8
plt.title("fibro3min190515")
plt.ylabel("enrichment")
plt.xlabel("region name")
plt.savefig("fibro190515.enrichment.png", bbox_inches='tight')
