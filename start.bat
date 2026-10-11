@echo off
chcp 65001 >nul
cd /d "%~dp0"
if not exist ".venv\Scripts\python.exe" (
    py -3.12 -m venv .venv
    if errorlevel 1 goto failed
)
.venv\Scripts\python.exe -c "import rdkit, PIL, requests" >nul 2>&1
if errorlevel 1 (
    .venv\Scripts\python.exe -m pip install -r requirements.txt
    if errorlevel 1 goto failed
)
.venv\Scripts\python.exe chem_drawer.py
if errorlevel 1 goto failed
exit /b 0
:failed
echo 启动失败。请安装 Python 3.12 64位，并检查网络和上方错误。
pause
exit /b 1
