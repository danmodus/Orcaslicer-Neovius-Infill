@REM OrcaSlicer build script for Windows (VS 2026) - Version 2.2 (Canary)
@echo off
set WP=%CD%
echo Current Directory: %WP%
echo Script Version: 2.2 (Canary)

set debug=OFF
set debuginfo=OFF
if "%1"=="debug" set debug=ON
if "%2"=="debug" set debug=ON
if "%1"=="debuginfo" set debuginfo=ON
if "%2"=="debuginfo" set debuginfo=ON
if "%debug%"=="ON" (
    set build_type=Debug
    set build_dir=build-dbg
) else (
    if "%debuginfo%"=="ON" (
        set build_type=RelWithDebInfo
        set build_dir=build-dbginfo
    ) else (
        set build_type=Release
        set build_dir=build
    )
)
echo Build type: %build_type%
echo Build dir: %build_dir%

@REM Try to detect the correct generator for VS 2026
if not defined GENERATOR set GENERATOR=Visual Studio 18 2026
echo Using generator: %GENERATOR%

setlocal DISABLEDELAYEDEXPANSION 

if "%1"=="slicer" (
    GOTO :slicer
)

echo ##########################################
echo # Building Dependencies...
echo ##########################################
echo.

if not exist "%WP%\deps\CMakeLists.txt" (
    echo ERROR: "%WP%\deps\CMakeLists.txt" NOT FOUND!
    echo Directory listing of %WP%\deps:
    dir "%WP%\deps"
    exit /b 1
)

cd deps
if not exist %build_dir% mkdir %build_dir%
cd %build_dir%
set "SIG_FLAG="
if defined ORCA_UPDATER_SIG_KEY set "SIG_FLAG=-DORCA_UPDATER_SIG_KEY=%ORCA_UPDATER_SIG_KEY%"

echo on
set CMAKE_POLICY_VERSION_MINIMUM=3.5
cmake "%WP%\deps" -G "%GENERATOR%" -A x64 -DCMAKE_BUILD_TYPE=%build_type%
if %errorlevel% neq 0 (
    echo "Dependency configuration failed!"
    exit /b %errorlevel%
)

cmake --build . --config %build_type% --target deps -- -m
if %errorlevel% neq 0 (
    echo "Dependency build failed! Please check the output above."
    exit /b %errorlevel%
)
@echo off

if "%1"=="deps" exit /b 0

:slicer
echo.
echo ##########################################
echo # Building OrcaSlicer...
echo ##########################################
echo.

cd %WP%
if not exist %build_dir% mkdir %build_dir%
cd %build_dir%

@REM Set the prefix path using forward slashes for CMake compatibility
set "PREFIX_PATH=%WP%\deps\%build_dir%\OrcaSlicer_dep\usr\local"
set "PREFIX_PATH=%PREFIX_PATH:\=/%"
set "WXWIN=%PREFIX_PATH%"

echo Checking for dependencies in: %PREFIX_PATH%
set "BOOST_CHECK_DIR=%WP%\deps\%build_dir%\OrcaSlicer_dep\usr\local\include\boost"
if not exist "%BOOST_CHECK_DIR%" (
    echo.
    echo WARNING: Boost not found in dependencies folder!
    echo Expected: %BOOST_CHECK_DIR%
    echo Please make sure the 'deps' build finished successfully.
    echo.
) else (
    echo Found Boost headers. Searching for BoostConfig.cmake...
    for /d %%d in ("%WP%\deps\%build_dir%\OrcaSlicer_dep\usr\local\lib\cmake\Boost-*") do (
        set "BOOST_DIR_SET=%%d"
    )
)

if defined BOOST_DIR_SET (
    set "BOOST_DIR_SET=%BOOST_DIR_SET:\=/%"
    echo Found Boost Config at: %BOOST_DIR_SET%
    set "BOOST_FLAGS=-DBoost_DIR=%BOOST_DIR_SET%"
) else (
    set "BOOST_FLAGS="
)

echo on
set CMAKE_POLICY_VERSION_MINIMUM=3.5
cmake "%WP%" -G "%GENERATOR%" -A x64 -DORCA_TOOLS=ON %SIG_FLAG% -DCMAKE_BUILD_TYPE=%build_type% -DCMAKE_PREFIX_PATH="%PREFIX_PATH%" -DWXWIN="%WXWIN%" %BOOST_FLAGS% -DBoost_DEBUG=ON
if %errorlevel% neq 0 (
    echo "OrcaSlicer configuration failed!"
    exit /b %errorlevel%
)

cmake --build . --config %build_type% --target ALL_BUILD -- -m
if %errorlevel% neq 0 (
    echo "OrcaSlicer build failed!"
    exit /b %errorlevel%
)
@echo off
cd ..
if exist scripts\run_gettext.bat call scripts\run_gettext.bat
cd %build_dir%
cmake --build . --target install --config %build_type%

