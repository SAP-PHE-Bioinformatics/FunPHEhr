#!/usr/bin/env python3
import pandas, pathlib, subprocess
kfiles = f"{input}".split()
print(kfiles)
id_table = pandas.DataFrame()
print(id_table)
for k in kfiles:
    kraken = pathlib.Path(k)
    df = pandas.read_csv(kraken, sep = "\t", header =None, names = ['percentage', 'frag1', 'frag2','code','taxon','name'])
    df['percentage'] = df['percentage'].apply(lambda x:float(x.strip('%')) if isinstance(x, str) == True else float(x))
    df = df.sort_values(by = ['percentage'], ascending = False)
    df = df[df['code'].isin(['U','S'])]
    df = df.reset_index(drop = True)
    print(df)
    tempdf = pandas.DataFrame()
    d = {}
    t = df['name'].count()
    print(t)
    for i in range((t if t < 3 else 3)):
        d.update({'Isolate': f"{kraken.parts[1]}",
            f"#{i+1} Match": df.loc[i,'name'].strip(), f"%{i+1}": df.loc[i,'percentage']
            })
        print(d)
    #   d = {'Isolate': f"{kraken.parts[1]}",
    #           '#1 Match': df.loc[0,'name'].strip(), '%1': df.loc[0,'percentage'],
    #           '#2 Match': df.loc[1,'name'].strip(), '%2': df.loc[1,'percentage'],
    #           '#3 Match': df.loc[2,'name'].strip(), '%3': df.loc[2,'percentage']
    #        }

    tempdf = pandas.DataFrame(data = d, index= [0])
    if id_table.empty:
            id_table = tempdf
    else:
            id_table = id_table.append(tempdf, sort = True)
cols_list = ['Isolate', '#1 Match', '%1', '#2 Match', '%2', '#3 Match', '%3']
id_table = id_table.reindex(cols_list, axis = 'columns')
id_table.to_csv(f"{output}", sep = "\t", index = False)
subprocess.run("sed -i 's/%[0-9]/%/g' {output}", shell=True)
