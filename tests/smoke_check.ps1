$root = Split-Path -Parent $PSScriptRoot

$required = @(
    (Join-Path $root "README.md"),
    (Join-Path $root "LICENSE"),
    (Join-Path $root "CITATION.cff"),
    (Join-Path $root "CHANGELOG.md"),
    (Join-Path $root "VERSION"),
    (Join-Path $root "bin\ume"),
    (Join-Path $root "config\ume.env.example"),
    (Join-Path $root "docs\USAGE.md"),
    (Join-Path $root "examples\toy-workflow\README.md"),
    (Join-Path $root "code\UMCAL\src\UME_RCALL_V2.1.sh")
)

$missing = @()
foreach ($path in $required) {
    if (-not (Test-Path $path)) {
        $missing += $path
    }
}

if ($missing.Count -gt 0) {
    Write-Error ("Missing required files:`n" + ($missing -join "`n"))
    exit 1
}

Write-Output "UME repository smoke check passed."
