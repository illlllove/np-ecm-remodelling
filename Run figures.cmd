@echo off
setlocal
if exist "%~dp0Environment\Python\python.exe" (
    "%~dp0Environment\Python\python.exe" -I -B -X utf8 "%~dp0Run figures.py" %*
) else (
    python -I -B -X utf8 "%~dp0Run figures.py" %*
)
exit /b %errorlevel%
