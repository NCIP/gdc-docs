"""updates the encyclopedia section in the mkdocs.yml
 should be run whenever a file is removed or added into the directory"""

import os
import yaml

ABSFILEPATH = os.path.dirname(os.path.realpath(__file__))
FILEARRAY = os.listdir(ABSFILEPATH + '/../../../docs/Encyclopedia/pages')

with open(ABSFILEPATH + '/../../../mkdocs.yml', 'r') as f:
    doc = yaml.load(f)

encycdict = next(d for (index, d) in enumerate(doc['pages']) \
             if d.get('EncyclopediaEntries', False) != False)

newlist = []

for x in range(len(FILEARRAY)):
    if FILEARRAY[x][-3:] == ".md":
        tempdict = {FILEARRAY[x][:-3]:"".join(['Encyclopedia/pages/', FILEARRAY[x][:-3], '.md'])}
        newlist.append(tempdict)

encycdict['EncyclopediaEntries'] = newlist

with open(ABSFILEPATH + '/../../../mkdocs.yml', 'w+') as f:
    f.write(yaml.dump(doc, default_flow_style=False))
