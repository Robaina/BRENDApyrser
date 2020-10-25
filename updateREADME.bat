rem This is a batch file to update README.md on Github based on README.ipynb
rem Will be called by a git hook upon commit of README.ipynb

"jupyter nbconvert --to markdown README.ipynb"
"git add daily_prices.json"
"git commit -m 'Update data'"
"git push origin master"
