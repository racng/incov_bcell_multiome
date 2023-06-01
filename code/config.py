from pathlib import Path

import seaborn as sns
sns.set_context('paper')

from matplotlib import rcParams

rcParams["savefig.transparent"] = True
rcParams['figure.facecolor'] = (1,1,1,0)
rcParams['axes.facecolor'] = (1,1,1,0)
rcParams["figure.dpi"] = 300
rcParams["savefig.dpi"] = 300
rcParams["savefig.format"] = 'png'
rcParams["savefig.bbox"] = 'tight'
rcParams["axes.grid"] = False
rcParams['figure.figsize'] = (3,3)

datadir = Path("/home/rng/data/incov_bcell_multiome/")
outputdir = Path("../output/")
outputdir.mkdir(exist_ok=True)