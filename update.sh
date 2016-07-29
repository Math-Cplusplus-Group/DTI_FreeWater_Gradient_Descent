#!/bin/bash

clear
echo "stagig ..."
echo git add --all

echo "checking the status:"
echo git status

echo "Do you want to proceed pushing the changes?"
read answer

if [ answer = "y" ]; then
	echo "commit ..."
	echo git commit -m \"update\"
	echo "pushing to master"
	echo git push origin master
fi
