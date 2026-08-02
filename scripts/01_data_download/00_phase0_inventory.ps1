########################################
# Script: 00_phase0_inventory.ps1
# Purpose: Inventory the locked project and calculate pre-analysis file metadata.
# Input: Local project tree staged according to config.yaml.
# Output: Inventory and manifest records.
# Software: PowerShell
# Version: 7.x or Windows PowerShell 5.1
# Random seed: Not applicable (deterministic)
# Author: Study authors
########################################
param(
  [string]$ProjectRoot = (Resolve-Path (Join-Path $PSScriptRoot '..\..')).Path,
  [string]$OutputDir = (Resolve-Path (Join-Path $PSScriptRoot '..')).Path
)

$ErrorActionPreference = 'Stop'
$analysisDirName = Split-Path -Leaf $OutputDir

function Get-Category([string]$relativePath) {
  $ext = [IO.Path]::GetExtension($relativePath).ToLowerInvariant()
  if ($relativePath -like '01_raw_data\*') { return 'raw_data' }
  if ($relativePath -like '04_intermediate\expression_matrices\*') { return 'processed_expression' }
  if ($relativePath -like '04_intermediate\normalized_data\*') { return 'normalized_expression' }
  if ($relativePath -like '02_metadata\*') { return 'metadata_or_annotation' }
  if ($ext -in @('.r','.py','.ps1','.sh','.ipynb')) { return 'code' }
  if ($ext -in @('.docx','.md','.txt')) { return 'document' }
  if ($ext -in @('.csv','.tsv','.xlsx','.xls')) { return 'table' }
  if ($ext -in @('.pdf','.png','.svg','.tiff','.jpg','.jpeg')) { return 'figure_or_render' }
  return 'other'
}

$files = Get-ChildItem -LiteralPath $ProjectRoot -Recurse -File -Force |
  Where-Object { $_.FullName -notlike (Join-Path $OutputDir '*') } |
  Sort-Object FullName
$rootPrefix = $ProjectRoot.TrimEnd([IO.Path]::DirectorySeparatorChar, [IO.Path]::AltDirectorySeparatorChar) + [IO.Path]::DirectorySeparatorChar

$manifest = foreach ($file in $files) {
  $relative = if ($file.FullName.StartsWith($rootPrefix, [StringComparison]::OrdinalIgnoreCase)) {
    $file.FullName.Substring($rootPrefix.Length)
  } else {
    $file.FullName
  }
  $hash = $null
  $hashStatus = 'ok'
  $hashError = ''
  try {
    $hash = (Get-FileHash -LiteralPath $file.FullName -Algorithm SHA256 -ErrorAction Stop).Hash.ToLowerInvariant()
  } catch {
    $hashStatus = 'unreadable_or_locked'
    $hashError = $_.Exception.Message
  }
  [pscustomobject]@{
    relative_path = $relative
    category = Get-Category $relative
    extension = $file.Extension.ToLowerInvariant()
    size_bytes = $file.Length
    last_write_time = $file.LastWriteTime.ToString('yyyy-MM-dd HH:mm:ss zzz')
    sha256 = $hash
    hash_status = $hashStatus
    hash_error = $hashError
  }
}

$manifestPath = Join-Path $OutputDir '02_file_manifest_sha256.csv'
$manifest | Export-Csv -LiteralPath $manifestPath -NoTypeInformation -Encoding UTF8
$manifest | Where-Object { $_.category -eq 'raw_data' } |
  Export-Csv -LiteralPath (Join-Path $OutputDir 'raw_data_sha256_baseline.csv') -NoTypeInformation -Encoding UTF8

$groups = $manifest | Group-Object category | Sort-Object Name
$extGroups = $manifest | Group-Object extension | Sort-Object Count -Descending
$largest = $manifest | Sort-Object size_bytes -Descending | Select-Object -First 30

$lines = [System.Collections.Generic.List[string]]::new()
$lines.Add('# Project inventory')
$lines.Add('')
$lines.Add(('Generated: {0}' -f (Get-Date).ToString('yyyy-MM-dd HH:mm:ss zzz')))
$lines.Add('')
$lines.Add(('Project root: `{0}`' -f $ProjectRoot))
$lines.Add('')
$lines.Add(('Baseline scope: all pre-existing files outside `{0}`.' -f $analysisDirName))
$lines.Add('')
$lines.Add(('Total files: {0}' -f $manifest.Count))
$lines.Add('')
$lines.Add('## File classes')
$lines.Add('')
$lines.Add('| Category | Files | Bytes |')
$lines.Add('|---|---:|---:|')
foreach ($g in $groups) {
  $bytes = ($g.Group | Measure-Object size_bytes -Sum).Sum
  $lines.Add(('| {0} | {1} | {2} |' -f $g.Name, $g.Count, $bytes))
}
$lines.Add('')
$lines.Add('## Extensions')
$lines.Add('')
$lines.Add('| Extension | Files |')
$lines.Add('|---|---:|')
foreach ($g in $extGroups) {
  $name = if ([string]::IsNullOrWhiteSpace($g.Name)) {'[none]'} else {$g.Name}
  $lines.Add(('| {0} | {1} |' -f $name, $g.Count))
}
$lines.Add('')
$lines.Add('## Largest files')
$lines.Add('')
$lines.Add('| Relative path | Size (GiB) | Category |')
$lines.Add('|---|---:|---|')
foreach ($f in $largest) {
  $lines.Add(('| `{0}` | {1:N3} | {2} |' -f $f.relative_path, ($f.size_bytes/1GB), $f.category))
}
$lines.Add('')
$lines.Add('## Complete path-level manifest')
$lines.Add('')
$lines.Add('The complete path, size, timestamp, category and SHA-256 inventory is in `02_file_manifest_sha256.csv`.')
$lines | Set-Content -LiteralPath (Join-Path $OutputDir '01_project_inventory.md') -Encoding UTF8

Write-Output ("WROTE {0} records to {1}" -f $manifest.Count, $manifestPath)
