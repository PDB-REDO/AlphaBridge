#!/bin/bash
#set -euo pipefail

# --- Defaults ---
PLOT=False
PARALLEL_JOBS=1   # default to sequential execution
SKIP_EXISTING=false

# --- Require positional arguments ---
if [ $# -lt 2 ]; then
    echo "Usage: $0 <alphabridge_folder> <input_directory> [--plot True|False] [--parallel N] [--skip-existing]"
    echo "  <alphabridge_folder> : folder where AlphaBridge is installed (define_interfaces.py is inside)"
    echo "  <input_directory>    : folder containing subdirectories to process"
    echo "  --plot True|False    : optional, default is False"
    echo "  --parallel N         : run N jobs in parallel (optional, default is 1)"
    echo "  --skip-existing      : skip subdirs that already contain AlphaBridge/alphabridge_data.json"
    exit 1
fi

ALPHABRIDGE_DIR="$1"
INPUT_DIR="$2"
shift 2

# --- Path to define_interfaces.py ---
SCRIPT_PATH="$ALPHABRIDGE_DIR/define_interfaces.py"
if [ ! -f "$SCRIPT_PATH" ]; then
    echo "Error: $SCRIPT_PATH not found!"
    exit 1
fi

# --- Validate INPUT_DIR ---
if [ ! -d "$INPUT_DIR" ]; then
    echo "Error: $INPUT_DIR is not a directory"
    exit 1
fi

# --- Parse optional flags ---
while [[ $# -gt 0 ]]; do
    case "$1" in
        --parallel)
            PARALLEL_JOBS="$2"
            if ! [[ "$PARALLEL_JOBS" =~ ^[0-9]+$ ]] || [ "$PARALLEL_JOBS" -lt 1 ]; then
                echo "Error: --parallel requires a positive integer >= 1"
                exit 1
            fi
            shift 2
            ;;
        --skip-existing)
            SKIP_EXISTING=false
            shift
            ;;
        --plot)
            if [[ $# -lt 2 ]]; then
                echo "Error: --plot requires an argument: True or False"
                exit 1
            fi
            PLOT="$2"
            if [[ "$PLOT" != "True" && "$PLOT" != "False" ]]; then
                echo "Error: --plot must be 'True' or 'False'"
                exit 1
            fi
            shift 2
            ;;
        *)
            echo "Unknown option: $1"
            exit 1
            ;;
    esac
done

# --- Find all subdirectories and store in an array ---
mapfile -d '' SUBDIRS < <(find "$INPUT_DIR" -mindepth 1 -maxdepth 1 -type d -print0)

# --- Function to process a subdirectory ---
run_job() {
    local subdir="$1"

    if [ "$SKIP_EXISTING" = "true" ] && [ -f "$subdir/AlphaBridge/alphabridge_data.json" ]; then
        echo "Skipping (already processed): $subdir"
        return 0
    fi

    echo "Processing: $subdir"
    python3 "$SCRIPT_PATH" -i "$subdir" -p $PLOT
}

# --- Unified processing loop ---
jobs=0
echo "Processing ${#SUBDIRS[@]} subdirectories with $PARALLEL_JOBS parallel jobs:"
for subdir in "${SUBDIRS[@]}"; do
    run_job "$subdir" &

    ((jobs++))
    while (( jobs >= PARALLEL_JOBS )); do
        wait -n
        ((jobs--))
    done
done
wait
