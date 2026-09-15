@echo off
cd /d "C:\Users\Admin\Documents\Developing stuff\Pula peg"
echo ================ %date% %time% >> run.log
py fetch.py >> run.log 2>&1
py estimate.py --window 90 --log >> run.log 2>&1
