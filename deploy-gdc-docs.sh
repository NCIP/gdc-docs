#!/bin/bash
cd gdc-docs
echo "$(date +'%d %B %Y - %k:%M'): Pulling latest version from git"
git pull
echo "$(date +'%d %B %Y - %k:%M'): Data_Portal UG: Building pandoc document"
mkdocs2pandoc -f Data_Portal_UG.yml -o docs/Data_Portal/PDF/Data_portal_UG.pd
echo "$(date +'%d %B %Y - %k:%M'): Data_Portal UG: Replacing strings in pandoc document "
sed -i -e 's/# / /g' docs/Data_Portal/PDF/Data_portal_UG.pd
sed -i -e 's/### /## /g' docs/Data_Portal/PDF/Data_portal_UG.pd
sed -i -e 's/\/site\//\/docs\//g' docs/Data_Portal/PDF/Data_portal_UG.pd
echo "$(date +'%d %B %Y - %k:%M'): Data_Portal UG: Building PDF from pandoc document "
pandoc --toc -V documentclass=report -V geometry:"top=2cm, bottom=1.5cm, left=1cm, right=1cm" -f markdown+grid_tables+table_captions -o docs/Data_Portal/PDF/Data_portal_UG.pdf docs/Data_Portal/PDF/Data_portal_UG.pd

echo "$(date +'%d %B %Y - %k:%M'): Data_Submission_Portal UG: Building pandoc document"
mkdocs2pandoc -f Data_Portal_UG.yml -o docs/Data_Submission_Portal/PDF/Data_Submission_Portal_UG.pd
echo "$(date +'%d %B %Y - %k:%M'): Data_Submission_Portal UG: Replacing strings in pandoc document "
sed -i -e 's/# / /g' docs/Data_Submission_Portal/PDF/Data_Submission_Portal_UG.pd
sed -i -e 's/### /## /g' docs/Data_Submission_Portal/PDF/Data_Submission_Portal_UG.pd
sed -i -e 's/\/site\//\/docs\//g' docs/Data_Submission_Portal/PDF/Data_Submission_Portal_UG.pd
echo "$(date +'%d %B %Y - %k:%M'): Data_Submission_Portal UG: Building PDF from pandoc document "
pandoc --toc -V documentclass=report -V geometry:"top=2cm, bottom=1.5cm, left=1cm, right=1cm" -f markdown+grid_tables+table_captions -o docs/Data_Portal/PDF/Data_Submission_Portal_UG.pdf docs/Data_Portal/PDF/Data_Submission_Portal_UG.pd


#Add other guides and mkdocs deploy
