# How to Build on Windows (64-bit)

Building OrcaSlicer from source on Windows requires several tools and can take 20–40 minutes depending on your hardware.

## 1. Prerequisites
Install these tools using a terminal (like PowerShell) or download them manually:

*   **Visual Studio**: [VS 2019 or 2022](https://visualstudio.microsoft.com/downloads/) (Community Edition is free).
    *   *Command*: `winget install --id=Microsoft.VisualStudio.Community -e`
*   **CMake**: [CMake 3.x or 4.x](https://cmake.org/download/)
    *   *Command*: `winget install --id=Kitware.CMake -e`
*   **Git & Git LFS**: Essential for cloning and downloading large tools.
    *   *Command*: `winget install --id=Git.Git -e`
    *   *Command*: `winget install --id=GitHub.GitLFS -e`
*   **Strawberry Perl**: Required for some dependencies.
    *   *Command*: `winget install --id=StrawberryPerl.StrawberryPerl -e`

## 2. Prepare the Source
1.  **Clone the Repository**:
    ```bash
    git clone https://github.com/danmodus/Orcaslicer-Neovius-Infill.git
    cd Orcaslicer-Neovius-Infill
    ```
2.  **Download Large Files (Git LFS)**:
    This is a critical step for Windows users to get the necessary build tools:
    ```bash
    git lfs pull
    ```

## 3. Build the Project
1.  Open the **"x64 Native Tools Command Prompt for VS 2022"** (search for it in the Start Menu).
2.  Navigate to your cloned folder:
    ```cmd
    cd C:\path\to\Orcaslicer-Neovius-Infill
    ```
3.  Run the build script:
    ```cmd
    build_release_vs2022.bat
    ```
    *(This script will automatically build all dependencies and then OrcaSlicer itself).*

## 4. Run the Application
Once the build finished successfully, the executable will be located here:
`build\src\Release\OrcaSlicer.exe`

---
> [!TIP]
> **Troubleshooting**: If the build fails, try deleting the `build` and `deps/build` folders and running the `.bat` script again. For more details, see the [official OrcaSlicer Wiki](https://github.com/SoftFever/OrcaSlicer/wiki/How-to-build).
