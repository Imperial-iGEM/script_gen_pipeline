# Testing exception warning
import pandas as pd 

f = {'prefixes': [1,1,1,1,1], 'parts': [2,2,2,2,2], 'suffixes': [3,3,3,3,3]}
a = []
for i in range(len(f)):
    a.append(pd.DataFrame.from_dict(f))

c = pd.concat(a)
print(c)