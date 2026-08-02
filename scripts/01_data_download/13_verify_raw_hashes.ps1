########################################
# Script: 13_verify_raw_hashes.ps1
# Purpose: Verify final raw-file SHA-256 values against the locked baseline.
# Input: Raw-data directory and baseline manifest.
# Output: Integrity report and final hash manifest.
# Software: PowerShell
# Version: 7.x or Windows PowerShell 5.1
# Random seed: Not applicable (deterministic)
# Author: Study authors
########################################
$ErrorActionPreference = 'Stop'
$projectRoot = (Get-Location).Path
$outputRoot = Join-Path $projectRoot '13_methodological_reanalysis_v9'
$baselinePath = Join-Path $outputRoot 'raw_data_sha256_baseline.csv'
$finalPath = Join-Path $outputRoot 'raw_data_sha256_final.csv'
$reportPath = Join-Path $outputRoot 'raw_data_integrity_check.md'

$baseline = Import-Csv -LiteralPath $baselinePath
$rows = foreach ($row in $baseline) {
    $absolutePath = Join-Path $projectRoot $row.relative_path
    if (Test-Path -LiteralPath $absolutePath -PathType Leaf) {
        try {
            $item = Get-Item -LiteralPath $absolutePath
            $hash = (Get-FileHash -LiteralPath $absolutePath -Algorithm SHA256).Hash.ToLowerInvariant()
            [pscustomobject]@{
                relative_path = $row.relative_path
                baseline_sha256 = $row.sha256
                final_sha256 = $hash
                baseline_size_bytes = $row.size_bytes
                final_size_bytes = $item.Length
                hash_status = if ($hash -eq $row.sha256 -and [int64]$item.Length -eq [int64]$row.size_bytes) { 'unchanged' } else { 'changed' }
                error = ''
            }
        } catch {
            [pscustomobject]@{
                relative_path = $row.relative_path
                baseline_sha256 = $row.sha256
                final_sha256 = ''
                baseline_size_bytes = $row.size_bytes
                final_size_bytes = ''
                hash_status = 'hash_error'
                error = $_.Exception.Message
            }
        }
    } else {
        [pscustomobject]@{
            relative_path = $row.relative_path
            baseline_sha256 = $row.sha256
            final_sha256 = ''
            baseline_size_bytes = $row.size_bytes
            final_size_bytes = ''
            hash_status = 'missing'
            error = 'File not found during final verification'
        }
    }
}
$rows | Export-Csv -LiteralPath $finalPath -NoTypeInformation -Encoding UTF8
$changed = @($rows | Where-Object { $_.hash_status -eq 'changed' })
$missing = @($rows | Where-Object { $_.hash_status -eq 'missing' })
$errors = @($rows | Where-Object { $_.hash_status -eq 'hash_error' })
$unchanged = @($rows | Where-Object { $_.hash_status -eq 'unchanged' })

$lines = @(
    '# Final raw-data integrity check',
    '',
    ('Generated: ' + (Get-Date -Format 'yyyy-MM-dd HH:mm:ss zzz')),
    '',
    ('- Baseline raw files: ' + $baseline.Count),
    ('- Unchanged by SHA-256 and byte size: ' + $unchanged.Count),
    ('- Changed: ' + $changed.Count),
    ('- Missing: ' + $missing.Count),
    ('- Hash errors: ' + $errors.Count),
    ('- raw_data_modified = ' + $(if ($changed.Count -eq 0 -and $missing.Count -eq 0 -and $errors.Count -eq 0) { 'FALSE' } else { 'UNRESOLVED' })),
    ''
)
if ($changed.Count -gt 0) { $lines += 'Changed files:'; $lines += ($changed.relative_path | ForEach-Object { '- ' + $_ }) }
if ($missing.Count -gt 0) { $lines += 'Missing files:'; $lines += ($missing.relative_path | ForEach-Object { '- ' + $_ }) }
if ($errors.Count -gt 0) { $lines += 'Hash errors:'; $lines += ($errors.relative_path | ForEach-Object { '- ' + $_ }) }
$lines | Set-Content -LiteralPath $reportPath -Encoding UTF8
