"""updates the encyclopedia section in the mkdocs.yml
 should be run whenever a file is removed or added into the directory"""

import os
import yaml

try:
    from yaml import CLoader as Loader, CDumper as Dumper
except ImportError:
    from yaml import Loader, Dumper

ABSFILEPATH = os.path.dirname(os.path.realpath(__file__))
FILEARRAY = os.listdir(ABSFILEPATH + "/docs/Encyclopedia/pages")
FILEARRAY = sorted(FILEARRAY, key=str.lower)

with open(ABSFILEPATH + "/mkdocs.yml", "r") as f:
    doc = yaml.load(f, Loader=Loader)

encycdict = next(
    d for (index, d) in enumerate(doc["nav"]) if d.get("Encyclopedia", False) != False
)

newlist = []

for x in range(len(FILEARRAY)):
    if FILEARRAY[x][-3:] == ".md":
        tempdict = {
            FILEARRAY[x][:-3].replace("_", " "): "".join(
                ["Encyclopedia/pages/", FILEARRAY[x][:-3], ".md"]
            )
        }
        newlist.append(tempdict)

encycdict["Encyclopedia"] = newlist

with open(ABSFILEPATH + "/mkdocs.yml", "w+") as f:
    f.write(
        yaml.dump(
            doc,
            default_flow_style=False,
            explicit_start=True,
            indent=2,
            width=80,
            Dumper=Dumper,
        )
    )
