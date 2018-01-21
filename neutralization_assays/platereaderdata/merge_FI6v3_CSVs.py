import glob
import pandas

dfs = []
for f in glob.glob('FI6v3_*.csv'):
    df = pandas.read_csv(f)
    sample = f.split('_', 1)[1].split('.')[0].replace(
            '_HA2', '-HA2').replace('-8T', '(-8T)')
    df = df.rename(columns=dict(
            [(str(i), '{0}-{1}'.format(sample, i)) for i in [1, 2, 3]]))
    dfs.append(df.set_index('concentration'))
pandas.concat(dfs, axis=1).to_csv('FI6v3_neutcurves.csv')
