#!/bin/bash
#
# Disk Usage Monitor for BugBuster Pipeline
# Monitors work directory and output directory sizes during pipeline execution
#
# Usage: ./monitor_disk_usage.sh [work_dir] [output_dir] [interval_seconds]
#

set -euo pipefail

# Default values
WORK_DIR="${1:-./work}"
OUTPUT_DIR="${2:-./results}"
INTERVAL="${3:-60}"
LOG_FILE="disk_usage_monitor.log"
CSV_FILE="disk_usage_monitor.csv"

# Colors for terminal output
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
NC='\033[0m' # No Color

# Function to get directory size in bytes
get_dir_size() {
    local dir="$1"
    if [ -d "$dir" ]; then
        du -sb "$dir" 2>/dev/null | cut -f1 || echo "0"
    else
        echo "0"
    fi
}

# Function to convert bytes to human readable
bytes_to_human() {
    local bytes=$1
    if [ "$bytes" -ge 1073741824 ]; then
        echo "scale=2; $bytes / 1073741824" | bc | awk '{printf "%.2f GB", $1}'
    elif [ "$bytes" -ge 1048576 ]; then
        echo "scale=2; $bytes / 1048576" | bc | awk '{printf "%.2f MB", $1}'
    elif [ "$bytes" -ge 1024 ]; then
        echo "scale=2; $bytes / 1024" | bc | awk '{printf "%.2f KB", $1}'
    else
        echo "${bytes} B"
    fi
}

# Function to get percentage change
get_change() {
    local prev=$1
    local curr=$2
    if [ "$prev" -eq 0 ]; then
        echo "N/A"
    else
        local change=$(echo "scale=2; (($curr - $prev) * 100) / $prev" | bc)
        if (( $(echo "$change > 0" | bc -l) )); then
            echo "+${change}%"
        else
            echo "${change}%"
        fi
    fi
}

# Initialize log files
echo "=== BugBuster Pipeline Disk Usage Monitor ===" > "$LOG_FILE"
echo "Started: $(date)" >> "$LOG_FILE"
echo "Work Directory: $WORK_DIR" >> "$LOG_FILE"
echo "Output Directory: $OUTPUT_DIR" >> "$LOG_FILE"
echo "Monitoring Interval: ${INTERVAL}s" >> "$LOG_FILE"
echo "" >> "$LOG_FILE"

# CSV header
echo "Timestamp,WorkDir_Bytes,WorkDir_GB,OutputDir_Bytes,OutputDir_GB,Total_GB,WorkDir_Change,OutputDir_Change" > "$CSV_FILE"

# Initial values
PREV_WORK_SIZE=0
PREV_OUTPUT_SIZE=0
ITERATION=0

echo -e "${GREEN}=== BugBuster Disk Usage Monitor ===${NC}"
echo "Work Directory: $WORK_DIR"
echo "Output Directory: $OUTPUT_DIR"
echo "Monitoring interval: ${INTERVAL}s"
echo "Logs: $LOG_FILE, $CSV_FILE"
echo ""
echo "Press Ctrl+C to stop monitoring"
echo ""

# Trap Ctrl+C
trap 'echo -e "\n${YELLOW}Monitoring stopped${NC}"; exit 0' INT

# Main monitoring loop
while true; do
    ITERATION=$((ITERATION + 1))
    TIMESTAMP=$(date +"%Y-%m-%d %H:%M:%S")
    
    # Get current sizes
    WORK_SIZE=$(get_dir_size "$WORK_DIR")
    OUTPUT_SIZE=$(get_dir_size "$OUTPUT_DIR")
    TOTAL_SIZE=$((WORK_SIZE + OUTPUT_SIZE))
    
    # Convert to GB for display
    WORK_GB=$(echo "scale=2; $WORK_SIZE / 1073741824" | bc)
    OUTPUT_GB=$(echo "scale=2; $OUTPUT_SIZE / 1073741824" | bc)
    TOTAL_GB=$(echo "scale=2; $TOTAL_SIZE / 1073741824" | bc)
    
    # Calculate changes
    if [ $ITERATION -gt 1 ]; then
        WORK_CHANGE=$(get_change $PREV_WORK_SIZE $WORK_SIZE)
        OUTPUT_CHANGE=$(get_change $PREV_OUTPUT_SIZE $OUTPUT_SIZE)
    else
        WORK_CHANGE="N/A"
        OUTPUT_CHANGE="N/A"
    fi
    
    # Write to CSV
    echo "$TIMESTAMP,$WORK_SIZE,$WORK_GB,$OUTPUT_SIZE,$OUTPUT_GB,$TOTAL_GB,$WORK_CHANGE,$OUTPUT_CHANGE" >> "$CSV_FILE"
    
    # Write to log
    {
        echo "[$TIMESTAMP]"
        echo "  Work Dir:   $(bytes_to_human $WORK_SIZE) ($WORK_CHANGE)"
        echo "  Output Dir: $(bytes_to_human $OUTPUT_SIZE) ($OUTPUT_CHANGE)"
        echo "  Total:      $(bytes_to_human $TOTAL_SIZE)"
        echo ""
    } >> "$LOG_FILE"
    
    # Terminal output with color coding
    echo -e "${YELLOW}[$TIMESTAMP]${NC}"
    
    # Color code work dir based on change
    if [ "$WORK_CHANGE" != "N/A" ] && [[ "$WORK_CHANGE" == -* ]]; then
        echo -e "  Work Dir:   ${GREEN}$(bytes_to_human $WORK_SIZE)${NC} (${GREEN}$WORK_CHANGE${NC})"
    else
        echo -e "  Work Dir:   $(bytes_to_human $WORK_SIZE) ($WORK_CHANGE)"
    fi
    
    # Color code output dir
    if [ "$OUTPUT_CHANGE" != "N/A" ] && [[ "$OUTPUT_CHANGE" == +* ]]; then
        echo -e "  Output Dir: ${GREEN}$(bytes_to_human $OUTPUT_SIZE)${NC} (${GREEN}$OUTPUT_CHANGE${NC})"
    else
        echo -e "  Output Dir: $(bytes_to_human $OUTPUT_SIZE) ($OUTPUT_CHANGE)"
    fi
    
    echo -e "  Total:      $(bytes_to_human $TOTAL_SIZE)"
    echo ""
    
    # Store previous values
    PREV_WORK_SIZE=$WORK_SIZE
    PREV_OUTPUT_SIZE=$OUTPUT_SIZE
    
    # Wait for next iteration
    sleep "$INTERVAL"
done
