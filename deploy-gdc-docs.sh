#!/bin/bash
cd ~/gdc-docs
echo "$(date +'%d %B %Y - %k:%M'): Pulling latest version from git"
/usr/bin/git pull
/usr/bin/git checkout develop

if [ ! -d "docs/Data_Portal/PDF/" ]; then
   mkdir docs/Data_Portal/PDF/
fi
echo "$(date +'%d %B %Y - %k:%M'): Data_Portal UG: Building pandoc document"
/usr/local/bin/mkdocs2pandoc -f Data_Portal_UG.yml -o docs/Data_Portal/PDF/Data_portal_UG.pd
echo "$(date +'%d %B %Y - %k:%M'): Data_Portal UG: Replacing strings in pandoc document "
/bin/sed -i -e 's/# / /g' docs/Data_Portal/PDF/Data_portal_UG.pd
/bin/sed -i -e 's/### /## /g' docs/Data_Portal/PDF/Data_portal_UG.pd
/bin/sed -i -e 's/\/site\//\/docs\//g' docs/Data_Portal/PDF/Data_portal_UG.pd
echo "$(date +'%d %B %Y - %k:%M'): Data_Portal UG: Building PDF from pandoc document "
/usr/bin/pandoc --toc -V documentclass=report -V geometry:"top=2cm, bottom=1.5cm, left=1cm, right=1cm" -f markdown+grid_tables+table_captions docs/Data_Portal/PDF/Data_Portal_Title.txt -o docs/Data_Portal/PDF/Data_portal_UG.pdf docs/Data_Portal/PDF/Data_portal_UG.pd

if [ ! -d "docs/Data_Portal/PDF/" ]; then
   mkdir docs/Data_Submission_Portal_UG/PDF/
fi
echo "$(date +'%d %B %Y - %k:%M'): Data_Submission_Portal UG: Building pandoc document"
/usr/local/bin/mkdocs2pandoc -f Data_Submission_Portal_UG.yml -o docs/Data_Submission_Portal/PDF/Data_Submission_Portal_UG.pd
echo "$(date +'%d %B %Y - %k:%M'): Data_Submission_Portal UG: Replacing strings in pandoc document "
/bin/sed -i -e 's/# / /g' docs/Data_Submission_Portal/PDF/Data_Submission_Portal_UG.pd
/bin/sed -i -e 's/### /## /g' docs/Data_Submission_Portal/PDF/Data_Submission_Portal_UG.pd
/bin/sed -i -e 's/\/site\//\/docs\//g' docs/Data_Submission_Portal/PDF/Data_Submission_Portal_UG.pd
echo "$(date +'%d %B %Y - %k:%M'): Data_Submission_Portal UG: Building PDF from pandoc document "
/usr/bin/pandoc --toc -V documentclass=report -V geometry:"top=2cm, bottom=1.5cm, left=1cm, right=1cm" -f markdown+grid_tables+table_captions docs/Data_Submission_Portal/PDF/Data_Submission_Portal_Title.txt -o docs/Data_Submission_Portal/PDF/Data_Submission_Portal_UG.pdf docs/Data_Submission_Portal/PDF/Data_Submission_Portal_UG.pd

if [ ! -d "docs/API/PDF/" ]; then
   mkdir docs/API/PDF/
fi
echo "$(date +'%d %B %Y - %k:%M'): API UG: Building pandoc document"
/usr/local/bin/mkdocs2pandoc -f API_UG.yml -o docs/API/PDF/API_UG.pd
echo "$(date +'%d %B %Y - %k:%M'): API UG: Replacing strings in pandoc document "
/bin/sed -i -e 's/# / /g' docs/API/PDF/API_UG.pd
/bin/sed -i -e 's/### /## /g' docs/API/PDF/API_UG.pd
/bin/sed -i -e 's/\/site\//\/docs\//g' docs/API/PDF/API_UG.pd
echo "$(date +'%d %B %Y - %k:%M'): API UG: Building PDF from pandoc document "
/usr/bin/pandoc --toc -V documentclass=report -V geometry:"top=2cm, bottom=1.5cm, left=1cm, right=1cm" -f markdown+grid_tables+table_captions docs/API/PDF/API_Title.txt -o docs/API/PDF/API_UG.pdf docs/API/PDF/API_UG.pd

if [ ! -d "docs/Data_Transfer_Tool/PDF/" ]; then
   mkdir docs/Data_Transfer_Tool/PDF/
fi
echo "$(date +'%d %B %Y - %k:%M'): Data_Transfer_Tool UG: Building pandoc document"
/usr/local/bin/mkdocs2pandoc -f Data_Transfer_Tool_UG.yml -o docs/Data_Transfer_Tool/PDF/Data_Transfer_Tool_UG.pd
echo "$(date +'%d %B %Y - %k:%M'): Data_Transfer_Tool UG: Replacing strings in pandoc document "
/bin/sed -i -e 's/# / /g' docs/Data_Transfer_Tool/PDF/Data_Transfer_Tool_UG.pd
/bin/sed -i -e 's/### /## /g' docs/Data_Transfer_Tool/PDF/Data_Transfer_Tool_UG.pd
/bin/sed -i -e 's/\/site\//\/docs\//g' docs/Data_Transfer_Tool/PDF/Data_Transfer_Tool_UG.pd
echo "$(date +'%d %B %Y - %k:%M'): Data_Transfer_Tool UG: Building PDF from pandoc document "
/usr/bin/pandoc --toc -V documentclass=report -V geometry:"top=2cm, bottom=1.5cm, left=1cm, right=1cm" -f markdown+grid_tables+table_captions docs/Data_Transfer_Tool/PDF/Data_Transfer_Tool_Title.txt -o docs/Data_Transfer_Tool/PDF/Data_Transfer_Tool_UG.pdf docs/Data_Transfer_Tool/PDF/Data_Transfer_Tool_UG.pd

/usr/local/bin/mkdocs build -v --site-dir /var/www/gdc-docs-dev.nci.nih.gov/

#Temporary fix to address a link that would be broken otherwise in the submission portal due to a change of the dictionary name
rm /var/www/gdc-docs-dev.nci.nih.gov/Dictionary/ -R
ln -sfn /var/www/gdc-docs-dev.nci.nih.gov/Data_Dictionary /var/www/gdc-docs-dev.nci.nih.gov/Dictionary
