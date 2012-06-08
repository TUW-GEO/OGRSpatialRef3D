@echo off
IF "%1"=="" GoTo usage
IF "%2"=="" GoTo usage

rem if desitinatiion file doesnt' exists, copy it
IF NOT EXIST %1 (
	echo source file doesn't exist
	goto end
)
IF NOT EXIST %2 GOTO copy_file

rem if file sizes differ, copy it
if %~z1 NEQ %~z2 goto copy_file

rem both files exist and have identical size so we must analyse in detail
set TEMPFILE=~compare.tmp

FC %1 %2 > %TEMPFILE%
rem check if files are perfectly identical
findstr "FC:" %TEMPFILE% >nul
IF %ERRORLEVEL% EQU 0 GOTO nothing

rem lib files differ a little bit accept this differences as identical (there must be a date information in there)
echo analyse differences...
SET count=0
FOR /F "usebackq delims=" %%a IN (%TEMPFILE%) DO SET /A count+=1
echo %count%
IF %count% LSS 20 GOTO nothing

:copy_file
ECHO Copy file
copy /y %1 %2 >nul
goto end

:nothing
echo Files are identical
goto end

:usage
echo %0 copy the source file to the destination file if files differ or the destination file doesn't exist
echo usage: %0 <source-file> <dest-file>

:end
