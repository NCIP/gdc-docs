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

# Identify the various user guides from the yml index file
# Strip the last part of the script to keep only the user guide name
# Data_Transfer_Tool_UG.yml => Data_Transfer_Tool
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
   echo "$(date +'%d %B %Y - %k:%M'): ${ENVIRONMENT}: ${userGuide}: Building PDF from pandoc document "
   /usr/bin/pandoc --toc -V documentclass=report -V geometry:"top=2cm, bottom=1.5cm, left=1cm, right=1cm" -f markdown+grid_tables+table_captions docs/${userGuide}/PDF/${userGuide}_Title.txt -o docs/${userGuide}/PDF/${userGuide}_UG.pdf docs/${userGuide}/PDF/${userGuide}_UG.pd
done
