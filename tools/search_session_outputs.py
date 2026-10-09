import json
import re
import sys

session_path = sys.argv[1]
pattern = re.compile(sys.argv[2], re.IGNORECASE)

def walk(value):
    if isinstance(value, dict):
        for item in value.values():
            yield from walk(item)
    elif isinstance(value, list):
        for item in value:
            yield from walk(item)
    elif isinstance(value, str):
        yield value

with open(session_path, "r", encoding="utf-8") as handle:
    for line_number, line in enumerate(handle, 1):
        try:
            record = json.loads(line)
        except json.JSONDecodeError:
            continue
        for text in walk(record):
            if pattern.search(text):
                matches = [row for row in text.splitlines() if pattern.search(row)]
                for row in matches:
                    print(f"{line_number}: {row[:2000]}")
