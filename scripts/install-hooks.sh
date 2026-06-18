#!/usr/bin/env bash
set -euo pipefail
root="$(cd "$(dirname "$0")/.." && pwd)"
ln -sf ../../scripts/hooks/pre-commit "$root/.git/hooks/pre-commit"
chmod +x "$root/scripts/hooks/pre-commit"
echo "Installed pre-commit hook."
