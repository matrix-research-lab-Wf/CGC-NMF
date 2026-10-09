import json
import re
import sys


session_path = sys.argv[1]
target_name = sys.argv[2]
maximum_line = int(sys.argv[3])
patches = []

with open(session_path, "r", encoding="utf-8") as stream:
    for line_number, line in enumerate(stream, 1):
        if line_number > maximum_line:
            break
        record = json.loads(line)
        payload = record.get("payload", {})
        source = payload.get("input", "") if isinstance(payload, dict) else ""
        if "const patch =" not in source or target_name not in source:
            continue
        match = re.search(r'const patch = ("(?:\\.|[^"])*")', source, re.S)
        if not match:
            continue
        patches.append({"line": line_number, "patch": json.loads(match.group(1))})

print(json.dumps(patches, ensure_ascii=True))
