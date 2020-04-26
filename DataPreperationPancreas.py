import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split

data=pd.read_excel('GTEx_pancreas_liver_images_liverfat_pancreasfat.xlsx')
NNdata=data.drop(columns=['Subject.ID','Tissue.Sample.ID_pancreas', 'Pathology.Notes_liver','Pathology.Notes_pancreas', 'Tissue.Sample.ID_liver'])


NNdata = NNdata.rename(columns={'Fat.Percentage_liver':'Fat_Percentage_liver', 'Fat.Percentage_pancreas': 'Fat_Percentage_pancreas'})

NNdata = NNdata.rename(columns={'Age.Bracket':'Age_Bracket'})

NNdata = NNdata.rename(columns={'Hardy.Scale':'Hardy_Scale'})

NNdata = NNdata.rename(columns={'Pathology.Categories_pancreas':'Pathology_Categories_pancreas', 'Pathology.Categories_liver': 'Pathology_Categories_liver'})

NNdataPan=NNdata.copy()
NNdataPan=NNdataPan.drop(columns=['Fat_Percentage_liver', 'Pathology_Categories_liver'])
short_data2=NNdataPan.copy()

from sklearn.preprocessing import MultiLabelBinarizer
mlb = MultiLabelBinarizer()
mask = short_data2['Pathology_Categories_pancreas'].notnull()
arr = mlb.fit_transform(short_data2.loc[mask, 'Pathology_Categories_pancreas'].dropna().str.split(', '))
df = (pd.DataFrame(arr, index=short_data2.index[mask], columns=mlb.classes_).reindex(short_data2.index, fill_value=0))

# # check all values are 0 or 1
print(pd.unique(df.values.ravel('K')))

Final_Data = NNdataPan.join(df)

Final_Data = Final_Data.join(pd.get_dummies(Final_Data['Sex'],prefix='Sex'))
Final_Data = Final_Data.join(pd.get_dummies(Final_Data['Hardy_Scale'],prefix='Hardy'))
Final_Data=Final_Data.join(pd.get_dummies(Final_Data['Age_Bracket'],prefix='Age'))
Final_Data=Final_Data.drop(columns=['Age_Bracket','Hardy_Scale','Sex', 'Pathology_Categories_pancreas'])

# Final_Data.to_csv(r'insert path here', index=False)
