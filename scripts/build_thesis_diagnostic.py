"""Transcribe Smith thesis Tables 5.1/5.2 into a separate diagnostic config.

The accessible PDF text is OCR: LV elastance reads 100 x 10^6 N/m^5.
Visual verification was unavailable. This is an explicitly provisional candidate,
not a correction to the frozen baseline or a confirmed paper parameter set.
"""
import json
from pathlib import Path
from scipy.constants import mmHg

ROOT = Path(__file__).resolve().parents[1]
c = json.loads((ROOT/"configs/baseline.yaml").read_text())
p = c["source_parameters"]
# 1 mL = 1e-6 m^3; divide Pa by Pa/mmHg for all pressure factors.
for component, suffix, elastance, pressure_scale in (
    ("lvf_calculator", "lvf", 100, 10),
    ("rvf_calculator", "rvf", 54, 10),
    ("septum", "spt", 6500, 148)):
    p[component]["E_es_"+suffix] = elastance/mmHg
    p[component]["P_0_"+suffix] = pressure_scale/mmHg
p["pericardium"]["P_0_pcd"] = 66.7/mmHg
for component,suffix,e in (("vena_cava","vc",1.3),("pulmonary_artery","pa",72),
                           ("pulmonary_vein","pu",1.9),("aorta","ao",98)):
    p[component]["E_es_"+suffix] = e/mmHg
for suffix,r in {"mt":6.1,"av":2.75,"tc":1.,"pv":1.,"pul":9.4,"sys":170.}.items():
    p["heart_parameters"]["R_"+suffix] = r/mmHg
for suffix,l in {"mt":1.3e4,"av":5e4,"tc":1.3e4,"pv":2e4}.items():
    p["heart_parameters"]["L_"+suffix] = l*1e-6/mmHg
c["status"] = "Provisional OCR transcription of thesis Tables 5.1/5.2; not frozen or visually verified"
c["source"] = "https://ir.canterbury.ac.nz/bitstreams/fdab391f-f441-4288-a7c0-ca22d8624287/download"
c["driver"]["C_s"] = .27
c["transcription_notes"] = ["PDF pages 89/90 (printed 74/75); LV E reads 100, unverified visually",
    "A=1 B=80 C=.27 from thesis Eq 2.4; period .75s from journal's explicit 80 bpm",
    "Initial allocation inherited only for warm-up; stressed inventory remains 1500 mL",
    "All three thoracic offsets retained: thesis Eqs 2.19/2.20 explicitly include Pth",
    "SI conversion uses scipy.constants.mmHg; no parameter fitted to the 9 percent result"]
(ROOT/"configs/m3-thesis-candidate.json").write_text(json.dumps(c,indent=2)+"\n")
protocol={"version":"2026-10-04-thesis-ocr-v1", "base_config":"configs/m3-thesis-candidate.json",
    "warmup_phase_intervals":100,"consecutive_periodic_comparisons":3,"measurement_cycles":10,"sample_step_s":.001,
    "cases":[{"name":"thesis_baseline"},{"name":"thesis_thoracic_zero","knobs":{"thoracic_pressure_mmHg":0}}]}
(ROOT/"configs/m3-thesis-diagnostic.json").write_text(json.dumps(protocol,indent=2)+"\n")
