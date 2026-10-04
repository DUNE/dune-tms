#!/usr/bin/env python3
"""Build the stage-by-stage detector-simulation page: page_template.html + page_data.json (stage_page_data.py) +
notes.json (the prose that depends on the numbers; keys INTRO, META, N1-N6, N4X, OPEN, REGEN).

    python3 build_stage_page.py <page_data.json> <notes.json> <output.html>
"""
import json
import os
import sys

here = os.path.dirname(os.path.abspath(__file__))
data = json.load(open(sys.argv[1]))
data["notes"] = json.load(open(sys.argv[2]))
html = open(os.path.join(here, "page_template.html")).read()
html = html.replace("/*__DATA__*/null", json.dumps(data, separators=(",", ":")).replace("</", "<\\/"))
open(sys.argv[3], "w").write(html)
print("wrote", sys.argv[3], len(html), "bytes")
