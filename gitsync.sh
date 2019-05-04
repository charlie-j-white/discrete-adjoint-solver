
# prepare date-time string

string=$(date +"%a")
day=${string^^}
string=$(date +"%d-%b-%Y")
date=${string^^}
string=$(date +"%T")
time=${string^^}

# list of files to add, combine with .gitignore

git add ./

# rest of the sync commands

git commit -m 'auto sync:'\ $day\ $date\ $time
git push -u origin master





