import pandas as pd
from decimal import Decimal

num_files = 9 # Number of output files to read
for i in range(num_files):
    df = pd.read_csv(f"outputCleaned{i+1}.txt", sep=" ", header=None)
    df.columns = ["Time", "GrowthRate", "Cov(CO)", "Cov(O2)"]

    # Plot the coverage of CO and O2
    res1 = df.plot(x="Time", y=["Cov(CO)", "Cov(O2)"], title=f"Coverage of CO and O2 for stoichiometric coefficient {round(Decimal(0.1 + i * 0.1), 2)}").get_figure()
    res1.savefig(f"coverage{i+1}.png")
    # Plot the growth rate
    res2 = df.plot(x="Time", y="GrowthRate", title=f"Growth rate for stoichiometric coefficient {round(Decimal(0.1 + i * 0.1), 2)}")
    res2.get_figure().savefig(f"growthRate{i+1}.png")