# 🧬 CRISPR Toolkit - Project Organization Rules

This document defines the organizational rules and standards for the CRISPR Toolkit project to maintain a clean, professional structure.

## 📁 **Root Directory Rules**

### ✅ **ALLOWED in Root:**
- **README.md** - Main project documentation
- **LICENSE** - Project license file
- **run.sh** - Main automation script (ONLY automation script allowed)
- **setup.py** / **pyproject.toml** - Package configuration
- **requirements*.txt** - Dependency files
- **docker-compose.yml** - Docker orchestration (if not in /docker/)
- **.gitignore** / **.editorconfig** - Version control and editor configuration
- **Core directories only**: src/, tests/, docs/, docker/, data/, assets/, notebooks/, etc.

### ❌ **PROHIBITED in Root:**
- **.env files** - Use Docker environment variables instead
- **Test files (.py)** - ALL tests must go in /tests/ subdirectories
- **Script files (.py)** - ALL scripts must go in /tests/scripts/ or appropriate subdirectories
- **Documentation (.md)** - ALL markdown files except README.md must go in /docs/
- **Verification scripts** - Must go in /tests/scripts/
- **Individual Docker files** - Must go in /docker/ directory

## 🗂️ **File Organization Standards**

### **1. Test Files (.py)**
```
tests/
├── clinical/           # Clinical-related tests
├── demos/             # Demonstration tests
├── integration/       # Integration tests
├── scripts/           # All test scripts and verification scripts
├── unit/              # Unit tests
└── validation/        # Validation tests
```

**Rule**: ANY .py file that is a test, verification script, or utility script MUST be in the appropriate /tests/ subdirectory.

### **2. Documentation Files (.md)**
```
docs/
├── REORGANIZATION_SUCCESS.md
├── CHANGELOG.md
├── COMPLETION_REPORTS.md
├── API_DOCS.md
└── user_guide.md
```

**Rule**: ALL markdown files except README.md MUST be in /docs/ directory.

### **3. Docker Configuration**
```
docker/
├── Dockerfile
├── docker-compose.yml
├── .dockerignore
└── docker-compose.*.yml
```

**Rule**: ALL Docker-related files MUST be in /docker/ directory.

### **4. Source Code**
```
src/
└── crispr_toolkit/
    ├── cli/           # Command-line interfaces
    ├── analysis/      # Analysis modules
    ├── models/        # ML models
    ├── data/          # Data handling
    ├── pipelines/     # Analysis pipelines
    └── utils/         # Utility functions
```

**Rule**: ALL source code MUST be in /src/ directory following proper Python package structure.

## 🐳 **Docker Environment Rules**

### **No .env Files Policy**
- **Prohibition**: NO .env files in the project
- **Alternative**: Use Docker environment variables with default values
- **Pattern**: `${VARIABLE_NAME:-default_value}` in docker-compose.yml

### **Environment Variable Examples**
```yaml
environment:
  - DATABASE_URL=${DATABASE_URL:-postgresql://crispr:password@db:5432/crispr_toolkit}
  - REDIS_URL=${REDIS_URL:-redis://redis:6379/0}
  - CRISPR_ENV=${CRISPR_ENV:-development}
  - LOG_LEVEL=${LOG_LEVEL:-INFO}
```

## 🧹 **Cleanup Enforcement**

### **Automated Cleanup**
The `run.sh` script includes cleanup functions that:
1. Remove __pycache__ directories
2. Clean temporary files
3. Enforce file organization rules

### **Pre-commit Checks**
- Verify no .py test files in root
- Verify no .md files in root (except README.md)
- Verify no .env files exist
- Verify Docker files are in /docker/

## ✅ **Compliance Checklist**

Before any commit, ensure:

- [ ] No .py test/script files in root directory
- [ ] No .md files in root (except README.md)
- [ ] No .env files anywhere in project
- [ ] All Docker files in /docker/ directory
- [ ] All test files in appropriate /tests/ subdirectories
- [ ] All documentation in /docs/ directory
- [ ] Root directory contains only essential project files (<25 files)

## 🔧 **Maintenance Commands**

### **Check Compliance**
```bash
./run.sh clean          # Remove cache and temporary files
./run.sh test           # Run tests to verify organization
```

### **Manual Cleanup**
```bash
# Remove any .env files (prohibited)
find . -name ".env*" -type f -delete

# Move any misplaced test files
find . -maxdepth 1 -name "test_*.py" -exec mv {} tests/ \;

# Move any misplaced .md files (except README.md)
find . -maxdepth 1 -name "*.md" ! -name "README.md" -exec mv {} docs/ \;
```

## 🎯 **Benefits of This Organization**

1. **Professional Appearance**: Clean root directory suitable for deployment
2. **Clear Separation**: Easy to find files based on their purpose
3. **Docker-First**: Environment configuration through Docker
4. **Maintainable**: Consistent structure across all development phases
5. **Scalable**: Easy to add new features without cluttering root

## 📚 **References**

- [Python Packaging User Guide](https://packaging.python.org/)
- [Docker Best Practices](https://docs.docker.com/develop/best-practices/)
- [Git Repository Organization](https://git-scm.com/book/en/v2)

---

**Note**: These rules are enforced automatically through automation scripts and should be followed by all contributors to maintain project quality and professionalism.
