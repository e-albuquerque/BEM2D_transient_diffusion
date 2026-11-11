@echo off
setlocal
set "PATH=C:\msys64\ucrt64\bin;%PATH%"

echo Compilando DCC Moderno (Fortran 90 + OpenMP + LAPACK)...
gfortran -fopenmp -O3 -fdefault-real-8 -Wall -Wextra -static ^
    data_mod.f90 solver_mod.f90 assembly_mod.f90 io_mod.f90 main.f90 ^
    -o dcc_modern.exe -llapack -lblas

if %ERRORLEVEL% NEQ 0 (
    echo [ERRO] Falha na compilacao!
    exit /b %ERRORLEVEL%
)

del /Q /F *.mod 2>nul
echo [OK] dcc_modern.exe compilado com sucesso.
endlocal
