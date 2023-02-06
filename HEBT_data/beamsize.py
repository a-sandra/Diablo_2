import pandas as pd
from scipy.interpolate import InterpolatedUnivariateSpline

d = {'Ekinetic': [50, 60, 70, 90 ,110, 130, 150, 170, 190, 210, 230],
     'SIGMAF': [0.001, 0.001,  0.001,  0.001,  0.001,  0.001, 0.00122, 0.00164,0.00194,0.00219,0.00238]}
df = pd.DataFrame(data=d)

order =1
s = InterpolatedUnivariateSpline(df["Ekinetic"], df["SIGMAF"], k=order)

    def compute_beam_size(f, value):
        return f(value)


