"""Static repository checks only; never import/execute study analysis code."""
import json
import re
import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
tracked = subprocess.check_output(['git', 'ls-files', '-z'], cwd=root).decode().split('\0')
errors = []
for name in filter(None, tracked):
    path = root / name
    if path.suffix.lower() not in {'.r', '.md', '.json', '.csv', '.yml', '.yaml', '.py'} and path.name not in {'.gitignore', '.gitattributes'}:
        errors.append((name, 'unexpected file type'))
        continue
    text = path.read_text(encoding='utf-8-sig')
    if '\0' in text:
        errors.append((name, 'NUL bytes'))
    if not name.startswith('legacy/'):
        for n, line in enumerate(text.splitlines(), 1):
            if re.search(r'''(?i)(?:api_key|password|access_token|secret_key)\s*(?:=|<-)\s*["'][A-Za-z0-9_/-]{12,}''', line):
                errors.append((name, f'possible credential assignment at line {n}'))
    if path.suffix == '.md':
        for link in re.findall(r'\]\(([^)]+)\)', text):
            target = link.split('#')[0]
            if target and not re.match(r'https?://|mailto:', target) and not (path.parent / target).exists():
                errors.append((name, 'broken relative link: ' + target))
rows = json.loads((root / 'docs/script_inventory.json').read_text(encoding='utf-8'))
for row in rows:
    if row['reference'] and not (root / row['reference']).exists():
        errors.append((row['source'], 'selected reference absent'))
print(f"Inventory: {len(rows)} scripts; {sum(r['status']=='damaged' for r in rows)} damaged; "
      f"{sum(r['status']=='selected_reference' for r in rows)} selected references")
for item in errors:
    print(item)
if errors:
    raise SystemExit(1)
print('PASS: tracked file types, NUL/credential patterns, documentation links and references. No analyses executed.')
