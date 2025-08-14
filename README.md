# Métodos numéricos para EDP 2025-2

**Andrés Contreras, Melany Daza, Amelia Hoyos**

Este repositorio es el hogar de los tres miniproyectos realizados durante el curso de *Métodos Numéricos* en el semestre *2025-2*.

Para finales del año 2025 contará con cuatro carpetas:
- **Miniproyecto 1**
- **Miniproyecto 2**
- **Miniproyecto 3**
- **Numerical Methods**

Cada carpeta tendrá un documento explicando cada uno de los puntos y el código correspondiente con el que se solucionaron, excepto la carpeta de **Numerical Methods** que contendrá la implementación de los métodos numéricos vistos en clase que se usarán en los miniproyectos. 

La estructura de este proyecto se está gestionando con `uv` y contiene los siguientes archivos estructurales:

- `.gitignore` — Indica a Git qué archivos no rastrear.
- `.python-version` — Documenta con qué versión de Python se desarrolló el proyecto.
- `README.md` — Este archivo.
- `pyproject.toml` — Define metadatos del proyecto, dependencias y configuración de herramientas en un solo lugar.
- `uv.lock` — Archivo de bloqueo de dependencias.

---

## Instalación de la biblioteca editable `numericalMethods`

### 1️⃣ Crear y activar el entorno virtual
Desde la raíz del proyecto, crea un entorno virtual con `uv`:

```bash
uv venv
````

Activa el entorno virtual según tu sistema operativo:

* **Windows (PowerShell):**

```powershell
.venv\Scripts\Activate.ps1
```

* **Windows (CMD):**

```cmd
.venv\Scripts\activate.bat
```

* **macOS/Linux (bash/zsh):**

```bash
source .venv/bin/activate
```

---

### 2️⃣ Instalar la biblioteca en modo editable

Con el entorno virtual activado, instala el paquete `numericalMethods`:

```bash
uv pip install -e .
```

Esto permite que cualquier cambio en el código fuente se refleje inmediatamente sin necesidad de reinstalar.

---