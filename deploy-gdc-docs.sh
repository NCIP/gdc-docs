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

if [ -d "/home/ubuntu/gdc-docs-${ENVIRONMENT}/" ]; then
   echo "$(date +'%d %B %Y - %k:%M'): ${ENVIRONMENT}: Directory exists, removing"
   sudo rm /home/ubuntu/gdc-docs-${ENVIRONMENT}/ -R
	# Added by ray 2018-07-11 since it was erroring trying to find a directory ?
   sleep 10s
   echo "$(date +'%d %B %Y - %k:%M'): ${ENVIRONMENT}: Directory exists, removing"
else
   echo "$(date +'%d %B %Y - %k:%M'): ${ENVIRONMENT}: Directory ~/gdc-docs-${ENVIRONMENT}/ does not exist"
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

# Building virtualenv and installing dependencies
echo "$(date +'%d %B %Y - %k:%M'): ${ENVIRONMENT}: Building virtualenv for ${ENVIRONMENT}"
python3 -m virtualenv -p python3.8 venv
source venv/bin/activate
pip install pip-tools
pip-sync requirements.txt

#iconv --verbose -f ascii -t utf-8 -o /tmp/test docs/Data_Portal/PDF/Data_Portal_UG.pd
hasEncodingError=false
echo "$(date +'%d %B %Y - %k:%M'): ${ENVIRONMENT}: Veryfing if all MARKDOWN have no encoding issue"
for scanFile in $( find docs/ | grep md | egrep -v Eliminate ); do
   iconv -f ascii -t utf-8 -o /tmp/test ${scanFile} >> /tmp/${ENVIRONMENT}-buildlog.txt
   if [ "$?" -gt 0 ] ; then
      echo "$(date +'%d %B %Y - %k:%M'): ${ENVIRONMENT}: NOK: ${scanFile} has an encoding error, to address open in vim and use :goto POSITION"
      hasEncodingError=true
   else
      echo "$(date +'%d %B %Y - %k:%M'): ${ENVIRONMENT}: OK: ${scanFile} "
   fi
done

if $hasEncodingError  ; then
   echo "$(date +'%d %B %Y - %k:%M'): ${ENVIRONMENT}: ERROR: Some of the files have encoding errors, not building the site"
   if [ -f /tmp/${ENVIRONMENT}-buildlog.txt ]; then
      echo "$(date +'%d %B %Y - %k:%M'): ${ENVIRONMENT}: Copying log file"
      cp /tmp/${ENVIRONMENT}-buildlog.txt /var/www/gdc-docs-${ENVIRONMENT}.nci.nih.gov/buildlog.txt
   fi
   exit
fi

# Generating PDFs from *_UG.yml files.
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
   echo "$(date +'%d %B %Y - %k:%M'): ${ENVIRONMENT}: ${userGuide}: Building PDF documents"
   ENABLE_PDF_EXPORT=1 venv/bin/mkdocs build -f ${userGuide}_UG.yml
done

echo "$(date +'%d %B %Y - %k:%M'): ${ENVIRONMENT}: Cleaning previous website directory (rm)"
sudo rm /var/www/gdc-docs-${ENVIRONMENT}.nci.nih.gov/* -R
   # Added by ray 2018-07-11 since it was erroring trying to find a directory ?
   sleep 10s

echo "$(date +'%d %B %Y - %k:%M'): ${ENVIRONMENT}: Build Encyclopedia"
# python buildencyclopedia.py

echo "$(date +'%d %B %Y - %k:%M'): ${ENVIRONMENT}: Deploying new version to /var/www/gdc-docs-${ENVIRONMENT}.nci.nih.gov/"
venv/bin/mkdocs build -v --site-dir /var/www/gdc-docs-${ENVIRONMENT}.nci.nih.gov/

if [ -f /tmp/${ENVIRONMENT}-buildlog.txt ]; then
   echo "$(date +'%d %B %Y - %k:%M'): ${ENVIRONMENT}: Copying log file"
   cp /tmp/${ENVIRONMENT}-buildlog.txt /var/www/gdc-docs-${ENVIRONMENT}.nci.nih.gov/buildlog.txt
fi

# After build files will have ubuntu:ubuntu permissions giving error on requests. Fixing permissions here.
sudo chown -R ubuntu:www-data /var/www/gdc-docs-${ENVIRONMENT}.nci.nih.gov

