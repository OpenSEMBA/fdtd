# Output file
OUTFILE=${1:-git_info.txt}

# Check if git command exists and we're in a git repository
if command -v git >/dev/null 2>&1 && [ -e .git ]; then
    commit=$(git rev-parse --short HEAD 2>/dev/null || echo "UNKNOWN")
else
    commit="UNKNOWN"
fi

# Write file
cat > "$OUTFILE" <<EOF
GIT_COMMIT="$commit"
EOF

echo "Wrote commit info to $OUTFILE (commit: $commit)"