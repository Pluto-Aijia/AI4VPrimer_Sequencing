"""
app.py
Flask web application for AI4VPrimer Amplicon-Sanger Suite.
Provides a modern no-code browser interface for running Sanger primer designs.
"""

from flask import Flask, render_template, request, jsonify
import webbrowser
import threading
import os
from sanger_designer.pipeline import SangerAmpliconPipeline

app = Flask(__name__, template_folder="templates")


import re

UPLOAD_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "uploads")
os.makedirs(UPLOAD_DIR, exist_ok=True)


def sanitize_input_path(raw_path: str) -> str:
    """Normalize file paths from Windows, macOS, Linux, or file:// URLs."""
    if not raw_path:
        return ""
    # Strip quotes and whitespace
    p = raw_path.strip().strip('"').strip("'").strip()
    # Strip file:// URI scheme
    if p.startswith("file://"):
        p = p[7:]
        # On Windows, file:///C:/path -> C:/path
        if len(p) > 2 and p[0] == '/' and p[2] == ':':
            p = p[1:]
    # Expand user home dir ~
    p = os.path.expanduser(p)
    # Check direct or normalized path
    if os.path.exists(p):
        return os.path.abspath(p)
    norm = os.path.normpath(p)
    if os.path.exists(norm):
        return os.path.abspath(norm)
    return p


@app.route("/")
def index():
    return render_template("index.html")


@app.route("/api/upload", methods=["POST"])
def upload_file_api():
    """Direct file upload endpoint supporting Windows, Mac, and Linux clients."""
    try:
        if "file" not in request.files:
            return jsonify({"status": "ERROR", "message": "No file included in upload request"}), 400
        file = request.files["file"]
        if not file or file.filename == "":
            return jsonify({"status": "ERROR", "message": "No file selected"}), 400
        
        # Sanitize filename (handling Windows and Unix client paths)
        raw_name = file.filename.replace('\\', '/').split('/')[-1]
        clean_name = re.sub(r'[^a-zA-Z0-9_.-]', '_', raw_name)
        if not clean_name:
            clean_name = "uploaded_sequence.fasta"
        
        dest_path = os.path.join(UPLOAD_DIR, clean_name)
        file.save(dest_path)
        
        return jsonify({
            "status": "SUCCESS",
            "filename": clean_name,
            "filepath": dest_path
        })
    except Exception as e:
        return jsonify({"status": "ERROR", "message": f"Upload failed: {str(e)}"}), 500


@app.route("/api/run", methods=["POST"])
def run_pipeline_api():
    try:
        data = request.get_json() or {}
        raw_fasta = data.get("fasta_path", "")
        fasta_path = sanitize_input_path(raw_fasta)

        if not fasta_path or not os.path.exists(fasta_path):
            return jsonify({
                "status": "ERROR",
                "message": (
                    f"Input FASTA file not found: '{raw_fasta}'.\n"
                    "Tip: If you are connecting from a Windows PC or client machine, "
                    "please use the 'Upload File' tab to upload your FASTA file directly."
                )
            }), 400

        fwd_pcr = data.get("fwd_pcr")
        rev_pcr = data.get("rev_pcr")
        target_start = data.get("target_start")
        target_end = data.get("target_end")
        min_cov = float(data.get("min_cov", 80.0))
        min_tm = float(data.get("min_tm", 55.0))
        max_tm = float(data.get("max_tm", 60.0))
        out_path = sanitize_input_path(data.get("out_path", "")) or None

        pipeline = SangerAmpliconPipeline(
            fasta_path=fasta_path,
            fwd_pcr_primer=fwd_pcr,
            rev_pcr_primer=rev_pcr,
            target_subregion_start=target_start,
            target_subregion_end=target_end,
            min_coverage_pct=min_cov,
            min_tm=min_tm,
            max_tm=max_tm,
            output_report_path=out_path
        )

        res = pipeline.run()

        return jsonify({
            "status": "SUCCESS",
            "report": res["report_content"],
            "pcr_info": res["pcr_info"],
            "tiling_info": res["tiling_info"]
        })

    except Exception as e:
        return jsonify({
            "status": "ERROR",
            "message": str(e)
        }), 500


def open_browser():
    webbrowser.open_new("http://127.0.0.1:5001")


if __name__ == "__main__":
    port = int(os.environ.get("PORT", 5001))
    print(f"Starting Amplicon-Sanger Suite on http://127.0.0.1:{port}")
    # Run Flask
    app.run(host="127.0.0.1", port=port, debug=False)
