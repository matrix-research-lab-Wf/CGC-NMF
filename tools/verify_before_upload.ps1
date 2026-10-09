[CmdletBinding()]
param()

$ErrorActionPreference = 'Stop'
$repoRoot = Split-Path -Parent $PSScriptRoot
$failures = New-Object System.Collections.Generic.List[string]

$required = @(
    'README.md', 'LICENSE', 'LICENSE_SCOPE.md', 'CITATION.cff',
    'DATA_PROVENANCE.md', 'FINAL_MANUSCRIPT_RESULT_AUDIT.md',
    'manuscript\Evidence-guided conservative graph calibration.tex',
    'manuscript\Evidence-guided conservative graph calibration.pdf',
    'run_smoke_test.m', 'run_public_smoke_test.m', 'THIRD_PARTY_NOTICES.md',
    'data\README.md', 'data\main_six\Optdigits_Han.mat',
    'results\cgc_gnmfld_transfer\CGC_GNMFLD_exact_wilcoxon_holm.csv',
    'results\main_comparison\CGC_minus_GOCNMF_bootstrap_10000.csv'
)

foreach ($item in $required) {
    if (-not (Test-Path -LiteralPath (Join-Path $repoRoot $item))) {
        $failures.Add("Missing required file: $item")
    }
}

$prohibited = @(
    'data\main_six\PIE.mat','data\main_six\YaleB.mat',
    'data\main_six\COIL20_Obj.mat','data\main_six\COIL100_Obj.mat',
    'data\main_six\MNIST_Han.mat','data\neu_cls\NEU_CLS_32x32.mat',
    'data\neu_cls\NEU_CLS_LBP59_4x4.mat','code\core\hungarian.m'
)
foreach ($item in $prohibited) {
    if (Test-Path -LiteralPath (Join-Path $repoRoot $item)) {
        $failures.Add("License-blocked file is present: $item")
    }
}

$provenancePath = Join-Path $repoRoot 'DATA_PROVENANCE.md'
if (Test-Path -LiteralPath $provenancePath) {
    $unresolved = Select-String -LiteralPath $provenancePath -Pattern '\bunresolved\b'
    if ($unresolved) {
        $failures.Add("DATA_PROVENANCE.md still contains $($unresolved.Count) unresolved line(s).")
    }
}

$oversized = Get-ChildItem -LiteralPath $repoRoot -Recurse -File -Force |
    Where-Object { $_.FullName -notlike "*\.git\*" -and $_.Length -gt 95MB }
foreach ($file in $oversized) {
    $failures.Add("File exceeds the 95 MB safety limit: $($file.FullName)")
}

$texPath = Join-Path $repoRoot 'manuscript\Evidence-guided conservative graph calibration.tex'
if (Test-Path -LiteralPath $texPath) {
    $falseSentence = Select-String -LiteralPath $texPath -SimpleMatch `
        'SNMFWLP reached the stopping threshold in only a subset of runs'
    if ($falseSentence) {
        $failures.Add('The audited false SNMFWLP stopping sentence is still present.')
    }
}

if ($failures.Count -gt 0) {
    Write-Host 'PRE-UPLOAD VERIFICATION: FAIL' -ForegroundColor Red
    $failures | ForEach-Object { Write-Host " - $_" }
    exit 1
}

Write-Host 'PRE-UPLOAD VERIFICATION: PASS' -ForegroundColor Green
exit 0
