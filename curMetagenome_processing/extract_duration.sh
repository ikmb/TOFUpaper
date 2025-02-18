grep -r -l "Workflow execution completed successfully" ./*batch*/pipeline_info/TOFU-MAaPO*.html | while read -r file; do
        duration=$(grep -oP '(duration: <strong>).*?(?=\))' "$file" | sed -E 's/duration: <strong>//g; s/<\/strong>//g')
    if [[ -n $duration ]]; then
        echo "$file,$duration"
    fi
done

