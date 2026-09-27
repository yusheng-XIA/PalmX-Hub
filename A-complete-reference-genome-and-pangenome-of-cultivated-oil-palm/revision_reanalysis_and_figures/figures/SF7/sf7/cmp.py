import pandas as pd, numpy as np
X='../../deliver/Source_Data_split/Source_Data_Supplementary_Figs.xlsx'
raw=pd.read_excel(X,sheet_name='SF7_raw_dCt')
n=pd.read_csv('../qpcr/work/negdct_workbook_all.tsv',sep='\t')
n['assay']=n.primer.str[1].astype(int)
m=raw.merge(n,left_on=['Stage','Material','Primer assay'],right_on=['stage','material','assay'],how='outer',indicator=True)
print(m._merge.value_counts())
m['d']=m['log₂(target/ACTIN) = −ΔCt']-m.negdct
m['dn']=m['Target-gene wells (n)']-m.n
m['dmean']=m['Mean target-gene Ct']-m['mean']
bad=m[(m.d.abs()>1e-9)|(m.dn!=0)|(m.dmean.abs()>1e-9)]
print(bad[['Stage','Material','Primer assay','Target-gene Ct values','Mean target-gene Ct','mean','n','log₂(target/ACTIN) = −ΔCt','negdct','d']].to_string())
print('max abs diff others', m.drop(bad.index).d.abs().max())
