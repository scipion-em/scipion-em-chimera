from google.colab import files
import pandas as pd

df = pd.DataFrame(data={'col1': [1, 2], 'col2': [3, 4]})
df.to_csv('dataframe.csv')
!zip -q results.zip dataframe.cvs 
files.download('result.zip ')
