import pandas as pd

# retrieve data from the txt files of compositions
df_03atm = {
    '20C': pd.read_table('output/compositions/xy_20C_03atm.txt', header=0),
    '25C': pd.read_table('output/compositions/xy_25C_03atm.txt', header=0),
    '30C': pd.read_table('output/compositions/xy_30C_03atm.txt', header=0),
    '35C': pd.read_table('output/compositions/xy_35C_03atm.txt', header=0),
    '40C': pd.read_table('output/compositions/xy_40C_03atm.txt', header=0),
    '45C': pd.read_table('output/compositions/xy_45C_03atm.txt', header=0),
    '50C': pd.read_table('output/compositions/xy_50C_03atm.txt', header=0),
    '55C': pd.read_table('output/compositions/xy_55C_03atm.txt', header=0),
}
df_04atm = {
    '20C': pd.read_table('output/compositions/xy_20C_04atm.txt', header=0),
    '25C': pd.read_table('output/compositions/xy_25C_04atm.txt', header=0),
    '30C': pd.read_table('output/compositions/xy_30C_04atm.txt', header=0),
    '35C': pd.read_table('output/compositions/xy_35C_04atm.txt', header=0),
    '40C': pd.read_table('output/compositions/xy_40C_04atm.txt', header=0),
    '45C': pd.read_table('output/compositions/xy_45C_04atm.txt', header=0),
    '50C': pd.read_table('output/compositions/xy_50C_04atm.txt', header=0),
    '55C': pd.read_table('output/compositions/xy_55C_04atm.txt', header=0),
}
df_05atm = {
    '20C': pd.read_table('output/compositions/xy_20C_05atm.txt', header=0),
    '25C': pd.read_table('output/compositions/xy_25C_05atm.txt', header=0),
    '30C': pd.read_table('output/compositions/xy_30C_05atm.txt', header=0),
    '35C': pd.read_table('output/compositions/xy_35C_05atm.txt', header=0),
    '40C': pd.read_table('output/compositions/xy_40C_05atm.txt', header=0),
    '45C': pd.read_table('output/compositions/xy_45C_05atm.txt', header=0),
    '50C': pd.read_table('output/compositions/xy_50C_05atm.txt', header=0),
    '55C': pd.read_table('output/compositions/xy_55C_05atm.txt', header=0),
}
df_07atm = {
    '20C': pd.read_table('output/compositions/xy_20C_07atm.txt', header=0),
    '25C': pd.read_table('output/compositions/xy_25C_07atm.txt', header=0),
    '30C': pd.read_table('output/compositions/xy_30C_07atm.txt', header=0),
    '35C': pd.read_table('output/compositions/xy_35C_07atm.txt', header=0),
    '40C': pd.read_table('output/compositions/xy_40C_07atm.txt', header=0),
    '45C': pd.read_table('output/compositions/xy_45C_07atm.txt', header=0),
    '50C': pd.read_table('output/compositions/xy_50C_07atm.txt', header=0),
    '55C': pd.read_table('output/compositions/xy_55C_07atm.txt', header=0),
}
df_10atm = {
    '20C': pd.read_table('output/compositions/xy_20C_10atm.txt', header=0),
    '25C': pd.read_table('output/compositions/xy_25C_10atm.txt', header=0),
    '30C': pd.read_table('output/compositions/xy_30C_10atm.txt', header=0),
    '35C': pd.read_table('output/compositions/xy_35C_10atm.txt', header=0),
    '40C': pd.read_table('output/compositions/xy_40C_10atm.txt', header=0),
    '45C': pd.read_table('output/compositions/xy_45C_10atm.txt', header=0),
    '50C': pd.read_table('output/compositions/xy_50C_10atm.txt', header=0),
    '55C': pd.read_table('output/compositions/xy_55C_10atm.txt', header=0),
}
df_12atm = {
    '20C': pd.read_table('output/compositions/xy_20C_12atm.txt', header=0),
    '25C': pd.read_table('output/compositions/xy_25C_12atm.txt', header=0),
    '30C': pd.read_table('output/compositions/xy_30C_12atm.txt', header=0),
    '35C': pd.read_table('output/compositions/xy_35C_12atm.txt', header=0),
    '40C': pd.read_table('output/compositions/xy_40C_12atm.txt', header=0),
    '45C': pd.read_table('output/compositions/xy_45C_12atm.txt', header=0),
    '50C': pd.read_table('output/compositions/xy_50C_12atm.txt', header=0),
    '55C': pd.read_table('output/compositions/xy_55C_12atm.txt', header=0),
}
df_15atm = {
    '20C': pd.read_table('output/compositions/xy_20C_15atm.txt', header=0),
    '25C': pd.read_table('output/compositions/xy_25C_15atm.txt', header=0),
    '30C': pd.read_table('output/compositions/xy_30C_15atm.txt', header=0),
    '35C': pd.read_table('output/compositions/xy_35C_15atm.txt', header=0),
    '40C': pd.read_table('output/compositions/xy_40C_15atm.txt', header=0),
    '45C': pd.read_table('output/compositions/xy_45C_15atm.txt', header=0),
    '50C': pd.read_table('output/compositions/xy_50C_15atm.txt', header=0),
    '55C': pd.read_table('output/compositions/xy_55C_15atm.txt', header=0),
}
df_16atm = {
    '20C': pd.read_table('output/compositions/xy_20C_16atm.txt', header=0),
    '25C': pd.read_table('output/compositions/xy_25C_16atm.txt', header=0),
    '30C': pd.read_table('output/compositions/xy_30C_16atm.txt', header=0),
    '35C': pd.read_table('output/compositions/xy_35C_16atm.txt', header=0),
    '40C': pd.read_table('output/compositions/xy_40C_16atm.txt', header=0),
    '45C': pd.read_table('output/compositions/xy_45C_16atm.txt', header=0),
    '50C': pd.read_table('output/compositions/xy_50C_16atm.txt', header=0),
    '55C': pd.read_table('output/compositions/xy_55C_16atm.txt', header=0),
}

# create a list for temperature indexes
T = [float(key[:-1]) for key in df_16atm.keys()]

# concat all the compositions at different temperature
df_03_atm = pd.concat(
    [df_03atm['20C'], df_03atm['25C'], df_03atm['30C'],
     df_03atm['35C'], df_03atm['40C'], df_03atm['45C'],
     df_03atm['50C'], df_03atm['55C']], ignore_index=True)
df_04_atm = pd.concat(
    [df_04atm['20C'], df_04atm['25C'], df_04atm['30C'],
     df_04atm['35C'], df_04atm['40C'], df_04atm['45C'],
     df_04atm['50C'], df_04atm['55C']], ignore_index=True)
df_05_atm = pd.concat(
    [df_05atm['20C'], df_05atm['25C'], df_05atm['30C'],
     df_05atm['35C'], df_05atm['40C'], df_05atm['45C'],
     df_05atm['50C'], df_05atm['55C']], ignore_index=True)
df_07_atm = pd.concat(
    [df_07atm['20C'], df_07atm['25C'], df_07atm['30C'],
     df_07atm['35C'], df_07atm['40C'], df_07atm['45C'],
     df_07atm['50C'], df_07atm['55C']], ignore_index=True)
df_10_atm = pd.concat(
    [df_10atm['20C'], df_10atm['25C'], df_10atm['30C'],
     df_10atm['35C'], df_10atm['40C'], df_10atm['45C'],
     df_10atm['50C'], df_10atm['55C']], ignore_index=True)
df_12_atm = pd.concat(
    [df_12atm['20C'], df_12atm['25C'], df_12atm['30C'],
     df_12atm['35C'], df_12atm['40C'], df_12atm['45C'],
     df_12atm['50C'], df_12atm['55C']], ignore_index=True)
df_15_atm = pd.concat(
    [df_15atm['20C'], df_15atm['25C'], df_15atm['30C'],
     df_15atm['35C'], df_15atm['40C'], df_15atm['45C'],
     df_15atm['50C'], df_15atm['55C']], ignore_index=True)
df_16_atm = pd.concat(
    [df_16atm['20C'], df_16atm['25C'], df_16atm['30C'],
     df_16atm['35C'], df_16atm['40C'], df_16atm['45C'],
     df_16atm['50C'], df_16atm['55C']], ignore_index=True)

# reset the indexes as temperature
df_03_atm = df_03_atm.assign(Temp = T); df_03_atm.set_index('Temp', inplace=True)
df_04_atm = df_04_atm.assign(Temp = T); df_04_atm.set_index('Temp', inplace=True)
df_05_atm = df_05_atm.assign(Temp = T); df_05_atm.set_index('Temp', inplace=True)
df_07_atm = df_07_atm.assign(Temp = T); df_07_atm.set_index('Temp', inplace=True)
df_10_atm = df_10_atm.assign(Temp = T); df_10_atm.set_index('Temp', inplace=True)
df_12_atm = df_12_atm.assign(Temp = T); df_12_atm.set_index('Temp', inplace=True)
df_15_atm = df_15_atm.assign(Temp = T); df_15_atm.set_index('Temp', inplace=True)
df_16_atm = df_16_atm.assign(Temp = T); df_16_atm.set_index('Temp', inplace=True)

# write on csv files the dataframes compositions
df_03_atm.to_csv('output/compositions/csv/xy_03atm.csv', index=True)
df_04_atm.to_csv('output/compositions/csv/xy_04atm.csv', index=True)
df_05_atm.to_csv('output/compositions/csv/xy_05atm.csv', index=True)
df_07_atm.to_csv('output/compositions/csv/xy_07atm.csv', index=True)
df_10_atm.to_csv('output/compositions/csv/xy_10atm.csv', index=True)
df_12_atm.to_csv('output/compositions/csv/xy_12atm.csv', index=True)
df_15_atm.to_csv('output/compositions/csv/xy_15atm.csv', index=True)
df_16_atm.to_csv('output/compositions/csv/xy_16atm.csv', index=True)