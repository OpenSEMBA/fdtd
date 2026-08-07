# 🐞 Issue: Build Problem (remember to assign tag: 'Compilation Issue')

## 📌 Environment Information

- **Intel oneAPI Version (if applicable):**  
  > e.g., 2024.1.0

- **CMake Version:**  
  > e.g., 3.28.2

- **Operating System:**  
  > e.g., Ubuntu 22.04 / Windows 11 / WSL2

- **Fortran Compiler:**  
  > e.g., ifx 2024.1 / ifort / gfortran 13

- **CMake Generator:**  
  > e.g., Ninja / Unix Makefiles / Visual Studio 17 2022

- **IDE:**
  > e.g., VSCode

---

## ⚙️ Build Configuration

- **Execution Type:**  
  - [ ] Preset  
  - [ ] Direct Command

---

## 📄 Preset or Command Used

### ▶️ Preset (if applicable)

```json
{
  "configurePresets": [
    {
      "name": "",
      "generator": "",
      "binaryDir": "",
      "cacheVariables": {}
    }
  ]
}

```
### ▶️ Preset (if applicable)
```cmd
cmake -S . -B build \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_Fortran_COMPILER=ifx
cmake --build build
```

### Error log or explanation
> e.g., Missing variable is: CMAKE_Fortran_PREPROCESS_SOURCE
