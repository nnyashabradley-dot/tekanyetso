@echo off
cd /d "C:\Users\Admin\Documents\Developing stuff\Pula peg"
echo =============== %date% %time% >> run.log
py fetch.py >> run.log 2>&1
py estimate.py --window 90 --log >> run.log 2>&1
py scoreboard.py >> run.log 2>&1

git add data predictions.csv fetch.log README.md
git diff --cached --quiet
if errorlevel 1 (
  git commit -q -m "daily run"
  git push -q
)
