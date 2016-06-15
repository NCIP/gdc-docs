#!/bin/bash
# This script pull the latest version from git, build PDF and mkdocs.
# Takes the following arguments:
# - dev
# - qa
while [[ $# > 1 ]]
do
key="$1"

case $key in
    -e|--environment)
    ENVIRONMENT="$2"
    shift
    ;;
    --default)
    ENVIRONMENT='dev'
    ;;
    *)

    ;;
esac
shift
done
if [ "$ENVIRONMENT" != "dev" ] && [ "$ENVIRONMENT" != "qa" ]; then
   echo "$(date +'%d %B %Y - %k:%M'): ${ENVIRONMENT}: Incorrect environment, needs to be dev or qa"
   exit;
fi

echo "$(date +'%d %B %Y - %k:%M'): ${ENVIRONMENT}: Building script for ${ENVIRONMENT}"

if [ -d "~/gdc-docs-${ENVIRONMENT}/" ]; then
   echo "$(date +'%d %B %Y - %k:%M'): ${ENVIRONMENT}: Directory exists, removing"
   sudo rm ~/gdc-docs-${ENVIRONMENT}/ -R
fi
echo "$(date +'%d %B %Y - %k:%M'): ${ENVIRONMENT}: Creating directory and cloning git repo"
mkdir ~/gdc-docs-${ENVIRONMENT}/
git clone git@github.com:NCI-GDC/gdc-docs.git ~/gdc-docs-${ENVIRONMENT}/
cd ~/gdc-docs-${ENVIRONMENT}/

echo "$(date +'%d %B %Y - %k:%M'): ${ENVIRONMENT}: Switching to correct branch"
if [ "$ENVIRONMENT" = "dev" ] ; then
   /usr/bin/git checkout develop
elif [ "$ENVIRONMENT" = "qa" ] ; then
   /usr/bin/git checkout master
else
   exit;
fi

echo "$(date +'%d %B %Y - %k:%M'): ${ENVIRONMENT}: Veryfing if all MARKDOWN files are UTF-8 encoded"
countWrongFiles=$(for f in `find docs/ | egrep -v Eliminate`; do echo "$f" ' -- ' `file -bi "$f"` ; done | grep ".md" | grep -v "utf-8" | wc -l)
echo "$(date +'%d %B %Y - %k:%M'): ${ENVIRONMENT}: Number of incorrectly encoded files: ${countWrongFiles}"

if [ "$countWrongFiles" -gt 0 ] ; then
   echo "$(date +'%d %B %Y - %k:%M'): ${ENVIRONMENT}: ERROR the following files are not encoded in UTF-8"
   for f in `find docs/ | egrep -v Eliminate`; do echo "$f" ' -- ' `file -bi "$f"` ; done | grep ".md" | grep -v "utf-8"
   exit
fi

echo "$(date +'%d %B %Y - %k:%M'): ${ENVIRONMENT}: Looking for User Guides"
userGuides=()
for i in $( ls *_UG.yml ); do
   userGuides+=(${i::${#i}-7})
done
echo "$(date +'%d %B %Y - %k:%M'): ${ENVIRONMENT}: Number of User Guides found: ${#userGuides[@]}"

for userGuide in "${userGuides[@]}"; do
   echo "$(date +'%d %B %Y - %k:%M'): ${ENVIRONMENT}: ${userGuide}: Starting creation"
   if [ ! -d "docs/Data_Portal/PDF/" ]; then
      echo "$(date +'%d %B %Y - %k:%M'): ${ENVIRONMENT}: ${userGuide}: PDF Directory does not exists, creating ..."
      mkdir docs/${userGuide}/PDF/
   fi
   echo "$(date +'%d %B %Y - %k:%M'): ${ENVIRONMENT}: ${userGuide}: Building pandoc document"
   /usr/local/bin/mkdocs2pandoc -f ${userGuide}_UG.yml -o docs/${userGuide}/PDF/${userGuide}_UG.pd
   echo "$(date +'%d %B %Y - %k:%M'): ${ENVIRONMENT}: ${userGuide}: Replacing strings in pandoc document "
   /bin/sed -i -e 's/# / /g' docs/${userGuide}/PDF/${userGuide}_UG.pd
   /bin/sed -i -e 's/### /## /g' docs/${userGuide}/PDF/${userGuide}_UG.pd
   /bin/sed -i -e 's/\/site\//\/docs\//g' docs/${userGuide}/PDF/${userGuide}_UG.pd
   /bin/sed -i -e "s/(images/(https:\/\/gdc-docs.nci.nih.gov\/"$userGuide"\/Users_Guide\/images/g" docs/${userGuide}/PDF/${userGuide}_UG.pd #To make images clickable in the PDF
   echo "$(date +'%d %B %Y - %k:%M'): ${ENVIRONMENT}: ${userGuide}: Building PDF from pandoc document "
   /usr/bin/pandoc --toc -V documentclass=report -V geometry:"top=2cm, bottom=1.5cm, left=1cm, right=1cm" -f markdown+grid_tables+table_captions docs/${userGuide}/PDF/${userGuide}_Title.txt -o docs/${userGuide}/PDF/${userGuide}_UG.pdf docs/${userGuide}/PDF/${userGuide}_UG.pd
done

echo "$(date +'%d %B %Y - %k:%M'): ${ENVIRONMENT}: Cleaning previous website directory (rm)"
sudo rm /var/www/gdc-docs-${ENVIRONMENT}.nci.nih.gov/* -R

echo "$(date +'%d %B %Y - %k:%M'): ${ENVIRONMENT}: Deploying new version to /var/www/gdc-docs-${ENVIRONMENT}.nci.nih.gov/"
/usr/local/bin/mkdocs build -v --site-dir /var/www/gdc-docs-${ENVIRONMENT}.nci.nih.gov/

if [ -f /tmp/buildlog.txt ]; then
   echo "$(date +'%d %B %Y - %k:%M'): ${ENVIRONMENT}: Copying log file"
   cp /tmp/${ENVIRONMENT}-buildlog.txt /var/www/gdc-docs-${ENVIRONMENT}.nci.nih.gov/buildlog.txt
fi

#echo "$(date +'%d %B %Y - %k:%M'): ${ENVIRONMENT}: Temporarily creating symlink"
#Temporary fix to address a link that would be broken otherwise in the submission portal due to a change of the dictionary name
#rm /var/www/gdc-docs-${ENVIRONMENT}.nci.nih.gov/Dictionary/ -R
#ln -sfn /var/www/gdc-docs-${ENVIRONMENT}.nci.nih.gov/Data_Dictionary /var/www/gdc-docs-${ENVIRONMENT}.nci.nih.gov/Dictionary
