@echo off
cd /d %~dp0
python ncbi_downloader.py
if errorlevel 1 (
  echo.
  echo NCBI Downloader closed with an error.
  pause
)
