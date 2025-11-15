#!/bin/bash
# Helper script to properly update Python bindings
#
# Issue: maturin develop doesn't update the .so file in site-packages
# Solution: Manually copy the built library after maturin develop
#
# Usage: ./update_python_bindings.sh [ClassName]
#        If ClassName is provided, verifies it exists in the module

set -e

echo "============================================================"
echo "Updating Python Bindings"
echo "============================================================"

# Step 1: Build with maturin
echo ""
echo "1. Building with maturin develop..."
maturin develop --release --features python

# Step 2: Get site-packages path (try user site-packages first, then system)
SITE_PACKAGES=$(python3 -c "
import site
import os

# Try user site-packages first (--user installs)
user_site = site.getusersitepackages()
if os.path.exists(os.path.join(user_site, 'biometal')):
    print(user_site)
else:
    # Fall back to system site-packages
    for sp in site.getsitepackages():
        if os.path.exists(os.path.join(sp, 'biometal')):
            print(sp)
            break
    else:
        print(user_site)  # Default to user site even if not found
")

TARGET_SO="$SITE_PACKAGES/biometal/biometal.cpython-314-darwin.so"
SOURCE_LIB="target/release/libbiometal.dylib"

echo ""
echo "2. Copying .so file..."
echo "   Source: $SOURCE_LIB"
echo "   Target: $TARGET_SO"

# Check source exists
if [ ! -f "$SOURCE_LIB" ]; then
    echo "   ❌ ERROR: Source library not found!"
    echo "   Run 'cargo build --release --features python' first"
    exit 1
fi

# Copy the file
cp "$SOURCE_LIB" "$TARGET_SO"
echo "   ✓ Copied successfully"

# Step 3: Verify timestamps
echo ""
echo "3. Verifying timestamps..."
SOURCE_TIME=$(stat -f "%Sm" "$SOURCE_LIB")
TARGET_TIME=$(stat -f "%Sm" "$TARGET_SO")
echo "   Source: $SOURCE_TIME"
echo "   Target: $TARGET_TIME"

# Step 4: Test import (optional class verification)
echo ""
echo "4. Testing Python import..."
python3 << EOF
import sys
# Force reimport
for mod in list(sys.modules.keys()):
    if 'biometal' in mod:
        del sys.modules[mod]

import biometal
print("   ✓ biometal imported successfully")

# If a class name was provided, verify it exists
import sys
if len(sys.argv) > 1:
    class_name = sys.argv[1]
    if hasattr(biometal, class_name):
        print(f"   ✓ {class_name} found in module")
    else:
        print(f"   ❌ {class_name} NOT found in module")
        print(f"   Available classes: {[c for c in dir(biometal) if not c.startswith('_')][:10]}...")
        sys.exit(1)
EOF

if [ $# -eq 1 ]; then
    # Verify the specific class exists
    python3 -c "import biometal; assert hasattr(biometal, '$1'), '$1 not found'; print('   ✓ $1 verified')"
fi

echo ""
echo "============================================================"
echo "✓ Python bindings updated successfully!"
echo "============================================================"
