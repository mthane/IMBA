# Run the Python MVP tracker on example_video.mp4 (experiment output for R Shiny / scripts).
# Writes stdout.log, stderr.log, metadata.txt, <timestamp>-data/detections.csv.
# Place example_video.mp4 in legacy/data/ or pass a path as the first argument.

$ErrorActionPreference = "Stop"
$here = Split-Path -Parent $MyInvocation.MyCommand.Path
$repoRoot = Split-Path -Parent $here
$video = if ($args.Count -ge 1) { $args[0] } else { Join-Path $repoRoot "legacy\data\example_video.mp4" }

if (-not (Test-Path $video)) {
    Write-Host "Video not found: $video"
    Write-Host "Copy example_video.mp4 into legacy\data\ or pass the full path as the first argument."
    exit 1
}

Set-Location $repoRoot
python -m imba_tracker.cli --experiment-output --output-dir (Split-Path -Parent $video) $video
exit $LASTEXITCODE
